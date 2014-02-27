/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jit/mips/MoveEmitter-mips.h"

using namespace js;
using namespace js::jit;

MoveEmitterMIPS::MoveEmitterMIPS(MacroAssemblerMIPSCompat &masm)
  : inCycle_(false),
    masm(masm),
    pushedAtCycle_(-1),
    pushedAtSpill_(-1),
    spilledReg_(InvalidReg),
    spilledFloatReg_(InvalidFloatReg)
{
    pushedAtStart_ = masm.framePushed();
}

void
MoveEmitterMIPS::emit(const MoveResolver &moves)
{
    if (moves.hasCycles()) {
        // Reserve stack for cycle resolution
        masm.reserveStack(sizeof(double));
        pushedAtCycle_ = masm.framePushed();
    }

    for (size_t i = 0; i < moves.numMoves(); i++)
        emit(moves.getMove(i));
}

MoveEmitterMIPS::~MoveEmitterMIPS()
{
    assertDone();
}

Address
MoveEmitterMIPS::cycleSlot() const
{
    int offset = masm.framePushed() - pushedAtCycle_;
    MOZ_ASSERT(Imm16::isInSignedRange(offset));
    return Address(StackPointer, offset);
}

int32_t
MoveEmitterMIPS::getAdjustedOffset(const MoveOperand &operand)
{
    MOZ_ASSERT(operand.isMemoryOrEffectiveAddress());
    if (operand.base() != StackPointer)
        return operand.disp();

    // Adjust offset if stack pointer has been moved.
    return operand.disp() + masm.framePushed() - pushedAtStart_;
}

Address
MoveEmitterMIPS::getAdjustedAddress(const MoveOperand &operand)
{
    return Address(operand.base(), getAdjustedOffset(operand));
}


Register
MoveEmitterMIPS::tempReg()
{
    // We will not use spill reg. Ve will use Macroassembler's second scratch
    // and be carefull not to use macro instructions that clobber it.
    spilledReg_ = t8;
    return t8;
}

void
MoveEmitterMIPS::breakCycle(const MoveOperand &from, const MoveOperand &to, MoveOp::Type type)
{
    // There is some pattern:
    //   (A -> B)
    //   (B -> A)
    //
    // This case handles (A -> B), which we reach first. We save B, then allow
    // the original move to continue.
    switch (type) {
      case MoveOp::FLOAT32:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ls(temp, getAdjustedAddress(to));
            masm.ma_ss(temp, cycleSlot());
        } else {
            masm.ma_ss(to.floatReg(), cycleSlot());
        }
        break;
      case MoveOp::DOUBLE:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ld(temp, getAdjustedAddress(to));
            masm.ma_sd(temp, cycleSlot());
        } else {
            masm.ma_sd(to.floatReg(), cycleSlot());
        }
        break;
      case MoveOp::INT32:
      case MoveOp::GENERAL:
        // an non-vfp value
        if (to.isMemory()) {
            Register temp = tempReg();
            masm.ma_lw(temp, getAdjustedAddress(to));
            masm.ma_sw(temp, cycleSlot());
        } else {
            // Register t8 should not be moved by MoveEmitter.
            MOZ_ASSERT(to.reg() != spilledReg_);
            masm.ma_sw(to.reg(), cycleSlot());
        }
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("Unexpected move type");
    }
}

void
MoveEmitterMIPS::completeCycle(const MoveOperand &from, const MoveOperand &to, MoveOp::Type type)
{
    // There is some pattern:
    //   (A -> B)
    //   (B -> A)
    //
    // This case handles (B -> A), which we reach last. We emit a move from the
    // saved value of B, to A.
    switch (type) {
      case MoveOp::FLOAT32:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ls(temp, cycleSlot());
            masm.ma_ss(temp, getAdjustedAddress(to));
        } else {
            masm.ma_ls(to.floatReg(), cycleSlot());
        }
        break;
      case MoveOp::DOUBLE:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ld(temp, cycleSlot());
            masm.ma_sd(temp, getAdjustedAddress(to));
        } else {
            masm.ma_ld(to.floatReg(), cycleSlot());
        }
        break;
      case MoveOp::INT32:
      case MoveOp::GENERAL:
        if (to.isMemory()) {
            Register temp = tempReg();
            masm.ma_lw(temp, cycleSlot());
            masm.ma_sw(temp, getAdjustedAddress(to));
        } else {
            // Register t8 should not be moved by MoveEmitter.
            MOZ_ASSERT(to.reg() != spilledReg_);
            masm.ma_lw(to.reg(), cycleSlot());
        }
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("Unexpected move type");
    }
}

void
MoveEmitterMIPS::emitMove(const MoveOperand &from, const MoveOperand &to)
{
    if (from.isGeneralReg()) {
        // Register t8 should not be moved by MoveEmitter.
        MOZ_ASSERT(from.reg() != spilledReg_);

        if (to.isGeneralReg())
            masm.ma_move(to.reg(), from.reg());
        else if (to.isMemory())
            masm.ma_sw(from.reg(), getAdjustedAddress(to));
        else
            MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
    } else if (from.isMemory()) {
        if (to.isGeneralReg()) {
            masm.ma_lw(to.reg(),  getAdjustedAddress(from));
        } else if (to.isMemory()) {
            masm.ma_lw(tempReg(),  getAdjustedAddress(from));
            masm.ma_sw(tempReg(),  getAdjustedAddress(to));
        } else {
            MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
        }
    } else if (from.isEffectiveAddress()) {
        if (to.isGeneralReg()) {
            masm.ma_addu(to.reg(),  from.base(), Imm32(getAdjustedOffset(from)));
        } else if (to.isMemory()) {
            masm.ma_addu(tempReg(),  from.base(), Imm32(getAdjustedOffset(from)));
            masm.ma_sw(tempReg(),  getAdjustedAddress(to));
        } else {
            MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
        }
    } else {
        MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
    }
}

void
MoveEmitterMIPS::emitFloat32Move(const MoveOperand &from, const MoveOperand &to)
{
    // Ensure that we can use ScratchFloatReg in memory move.
    MOZ_ASSERT_IF(from.isFloatReg(), from.floatReg() != ScratchFloatReg);
    MOZ_ASSERT_IF(to.isFloatReg(), to.floatReg() != ScratchFloatReg);

    if (from.isFloatReg()) {
        if (to.isFloatReg()) {
            masm.as_movs(to.floatReg(), from.floatReg());
        } else if (to.isGeneralReg()) {
            // This should only be used when passing float parameter in a1,a2,a3
            MOZ_ASSERT(to.reg() == a1 || to.reg() == a2 || to.reg() == a3);
            masm.as_mfc1(to.reg(), from.floatReg());
        } else {
            MOZ_ASSERT(to.isMemory());
            masm.ma_ss(from.floatReg(), getAdjustedAddress(to));
        }
    } else if (to.isFloatReg()) {
        MOZ_ASSERT(from.isMemory());
        masm.ma_ls(to.floatReg(), getAdjustedAddress(from));
    } else if (to.isGeneralReg()) {
        MOZ_ASSERT(from.isMemory());
        // This should only be used when passing float parameter in a1,a2,a3
        MOZ_ASSERT(to.reg() == a1 || to.reg() == a2 || to.reg() == a3);
        masm.ma_lw(to.reg(), getAdjustedAddress(from));
    } else {
        MOZ_ASSERT(from.isMemory());
        MOZ_ASSERT(to.isMemory());
        masm.ma_ls(ScratchFloatReg, getAdjustedAddress(from));
        masm.ma_ss(ScratchFloatReg, getAdjustedAddress(to));
    }
}

void
MoveEmitterMIPS::emitDoubleMove(const MoveOperand &from, const MoveOperand &to)
{
    // Ensure that we can use ScratchFloatReg in memory move.
    MOZ_ASSERT_IF(from.isFloatReg(), from.floatReg() != ScratchFloatReg);
    MOZ_ASSERT_IF(to.isFloatReg(), to.floatReg() != ScratchFloatReg);

    if (from.isFloatReg()) {
        if (to.isFloatReg()) {
            masm.as_movd(to.floatReg(), from.floatReg());
        } else if (to.isGeneralReg()) {
            // This should only be used when passing double parameter in a2,a3
            if(to.reg() == a2)
                masm.as_mfc1(a2, from.floatReg());
            else if(to.reg() == a3)
                masm.as_mfc1_Odd(a3, from.floatReg());
            else
                MOZ_ASSUME_UNREACHABLE("Invalid emitDoubleMove arguments.");
        } else {
            MOZ_ASSERT(to.isMemory());
            masm.ma_sd(from.floatReg(), getAdjustedAddress(to));
        }
    } else if (to.isFloatReg()) {
        MOZ_ASSERT(from.isMemory());
        masm.ma_ld(to.floatReg(), getAdjustedAddress(from));
    } else if (to.isGeneralReg()) {
        MOZ_ASSERT(from.isMemory());
        // This should only be used when passing double parameter in a2,a3
        if(to.reg() == a2)
            masm.ma_lw(a2, getAdjustedAddress(from));
        else if(to.reg() == a3)
            masm.ma_lw(a3, Address (from.base(), getAdjustedOffset(from) + sizeof(uint32_t)));
        else
            MOZ_ASSUME_UNREACHABLE("Invalid emitDoubleMove arguments.");
    } else {
        MOZ_ASSERT(from.isMemory());
        MOZ_ASSERT(to.isMemory());
        masm.ma_ld(ScratchFloatReg, getAdjustedAddress(from));
        masm.ma_sd(ScratchFloatReg, getAdjustedAddress(to));
    }
}

void
MoveEmitterMIPS::emit(const MoveOp &move)
{
    const MoveOperand &from = move.from();
    const MoveOperand &to = move.to();

    if (move.isCycleEnd()) {
        MOZ_ASSERT(inCycle_);
        completeCycle(from, to, move.type());
        inCycle_ = false;
        return;
    }

    if (move.isCycleBegin()) {
        MOZ_ASSERT(!inCycle_);
        breakCycle(from, to, move.endCycleType());
        inCycle_ = true;
    }

    switch (move.type()) {
      case MoveOp::FLOAT32:
        emitFloat32Move(from, to);
        break;
      case MoveOp::DOUBLE:
        emitDoubleMove(from, to);
        break;
      case MoveOp::INT32:
      case MoveOp::GENERAL:
        emitMove(from, to);
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("Unexpected move type");
    }
}

void
MoveEmitterMIPS::assertDone()
{
    MOZ_ASSERT(!inCycle_);
}

void
MoveEmitterMIPS::finish()
{
    assertDone();

    masm.freeStack(masm.framePushed() - pushedAtStart_);
}
