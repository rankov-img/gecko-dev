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
    JS_ASSERT(Imm16::isInSignedRange(offset));
    return Address(StackPointer, offset);
}

int32_t MoveEmitterMIPS::getAdjustedOffset(const MoveOperand &operand)
{
    JS_ASSERT(operand.isMemoryOrEffectiveAddress());
    if (operand.base() != StackPointer)
        return operand.disp();

    // Adjust offset if stack pointer has been moved.
    return operand.disp() + masm.framePushed() - pushedAtStart_;
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
            masm.ma_ls(temp, to.base(), getAdjustedOffset(to));
            masm.ma_ss(temp, cycleSlot().base, cycleSlot().offset);
        } else {
            masm.ma_ss(to.floatReg(), cycleSlot().base, cycleSlot().offset);
        }
        break;
      case MoveOp::DOUBLE:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ld(temp, to.base(), getAdjustedOffset(to));
            masm.ma_sd(temp, cycleSlot().base, cycleSlot().offset);
        } else {
            masm.ma_sd(to.floatReg(), cycleSlot().base, cycleSlot().offset);
        }
        break;
      case MoveOp::INT32:
      case MoveOp::GENERAL:
        // an non-vfp value
        if (to.isMemory()) {
            Register temp = tempReg();
            masm.ma_lw(temp, to.base(), getAdjustedOffset(to));
            masm.ma_sw(temp, cycleSlot().base, cycleSlot().offset);
        } else {
            if (to.reg() == spilledReg_) {
                MOZ_ASSUME_UNREACHABLE("Register t8 should not be moved by MoveEmitter.");
            }
            masm.ma_sw(to.reg(), cycleSlot().base, cycleSlot().offset);
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
            masm.ma_ls(temp, cycleSlot().base, cycleSlot().offset);
            masm.ma_ss(temp, to.base(), getAdjustedOffset(to));
        } else {
            masm.ma_ls(to.floatReg(), cycleSlot().base, cycleSlot().offset);
        }
        break;
      case MoveOp::DOUBLE:
        if (to.isMemory()) {
            FloatRegister temp = ScratchFloatReg;
            masm.ma_ld(temp, cycleSlot().base, cycleSlot().offset);
            masm.ma_sd(temp, to.base(), getAdjustedOffset(to));
        } else {
            masm.ma_ld(to.floatReg(), cycleSlot().base, cycleSlot().offset);
        }
        break;
      case MoveOp::INT32:
      case MoveOp::GENERAL:
        if (to.isMemory()) {
            Register temp = tempReg();
            masm.ma_lw(temp, cycleSlot().base, cycleSlot().offset);
            masm.ma_sw(temp, to.base(), getAdjustedOffset(to));
        } else {
            if (to.reg() == spilledReg_) {
                MOZ_ASSUME_UNREACHABLE("Register t8 should not be moved by MoveEmitter.");
            }
            masm.ma_lw(to.reg(), cycleSlot().base, cycleSlot().offset);
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
        if (from.reg() == spilledReg_) {
            MOZ_ASSUME_UNREACHABLE("Register t8 should not be moved by MoveEmitter.");
        }

        if (to.isGeneralReg())
            masm.ma_move(to.reg(), from.reg());
        else if (to.isMemory())
            masm.ma_sw(from.reg(), to.base(), getAdjustedOffset(to));
        else
            MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
    } else if (from.isMemory()) {
        if (to.isGeneralReg()) {
            masm.ma_lw(to.reg(), from.base(), getAdjustedOffset(from));
        } else if (to.isMemory()) {
            masm.ma_lw(tempReg(), from.base(), getAdjustedOffset(from));
            masm.ma_sw(tempReg(), to.base(), getAdjustedOffset(to));
        } else {
            MOZ_ASSUME_UNREACHABLE("Invalid emitMove arguments.");
        }
    } else if (from.isEffectiveAddress()) {
        if (to.isGeneralReg()) {
            masm.ma_addu(to.reg(), from.base(), Imm32(getAdjustedOffset(from)));
        } else if (to.isMemory()) {
            masm.ma_addu(tempReg(), from.base(), Imm32(getAdjustedOffset(from)));
            masm.ma_sw(tempReg(), to.base(), getAdjustedOffset(to));
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
    JS_ASSERT_IF(from.isFloatReg(), from.floatReg() != ScratchFloatReg);
    JS_ASSERT_IF(to.isFloatReg(), to.floatReg() != ScratchFloatReg);

    if (from.isFloatReg()) {
        if (to.isFloatReg()) {
            masm.as_movs(to.floatReg(), from.floatReg());
        } else if (to.isGeneralReg()) {
            if(to.reg() == a1)
                masm.as_mfc1(a1, from.floatReg());
            else if(to.reg() == a2)
                masm.as_mfc1(a2, from.floatReg());
            else if(to.reg() == a3)
                masm.as_mfc1(a3, from.floatReg());
            else
                MOZ_ASSUME_UNREACHABLE("Invalid emitFloat32Move arguments.");
        } else {
            JS_ASSERT(to.isMemory());
            masm.ma_ss(from.floatReg(), to.base(), getAdjustedOffset(to));
        }
    } else if (to.isFloatReg()) {
        JS_ASSERT(from.isMemory());
        masm.ma_ls(to.floatReg(), from.base(), getAdjustedOffset(from));
    } else if (to.isGeneralReg()) {
        JS_ASSERT(from.isMemory());
        if(to.reg() == a1)
            masm.ma_lw(a1, from.base(), getAdjustedOffset(from));
        else if(to.reg() == a2)
            masm.ma_lw(a2, from.base(), getAdjustedOffset(from));
        else if(to.reg() == a3)
            masm.ma_lw(a3, from.base(), getAdjustedOffset(from));
        else
            MOZ_ASSUME_UNREACHABLE("Invalid emitFloat32Move arguments.");
    } else {
        JS_ASSERT(from.isMemory());
        JS_ASSERT(to.isMemory());
        masm.ma_ls(ScratchFloatReg, from.base(), getAdjustedOffset(from));
        masm.ma_ss(ScratchFloatReg, to.base(), getAdjustedOffset(to));
    }
}

void
MoveEmitterMIPS::emitDoubleMove(const MoveOperand &from, const MoveOperand &to)
{
    // Ensure that we can use ScratchFloatReg in memory move.
    JS_ASSERT_IF(from.isFloatReg(), from.floatReg() != ScratchFloatReg);
    JS_ASSERT_IF(to.isFloatReg(), to.floatReg() != ScratchFloatReg);

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
            JS_ASSERT(to.isMemory());
            masm.ma_sd(from.floatReg(), to.base(), getAdjustedOffset(to));
        }
    } else if (to.isFloatReg()) {
        JS_ASSERT(from.isMemory());
        masm.ma_ld(to.floatReg(), from.base(), getAdjustedOffset(from));
    } else if (to.isGeneralReg()) {
        JS_ASSERT(from.isMemory());
        // This should only be used when passing double parameter in a2,a3
        if(to.reg() == a2)
            masm.ma_lw(a2, from.base(), getAdjustedOffset(from));
        else if(to.reg() == a3)
            masm.ma_lw(a3, from.base(), getAdjustedOffset(from) + sizeof(intptr_t));
        else
            MOZ_ASSUME_UNREACHABLE("Invalid emitDoubleMove arguments.");
    } else {
        JS_ASSERT(from.isMemory());
        JS_ASSERT(to.isMemory());
        masm.ma_ld(ScratchFloatReg, from.base(), getAdjustedOffset(from));
        masm.ma_sd(ScratchFloatReg, to.base(), getAdjustedOffset(to));
    }
}

void
MoveEmitterMIPS::emit(const MoveOp &move)
{
    const MoveOperand &from = move.from();
    const MoveOperand &to = move.to();

    if (move.isCycleEnd()) {
        JS_ASSERT(inCycle_);
        completeCycle(from, to, move.type());
        inCycle_ = false;
        return;
    }

    if (move.isCycleBegin()) {
        JS_ASSERT(!inCycle_);
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
    JS_ASSERT(!inCycle_);
}

void
MoveEmitterMIPS::finish()
{
    assertDone();

    masm.freeStack(masm.framePushed() - pushedAtStart_);
}
