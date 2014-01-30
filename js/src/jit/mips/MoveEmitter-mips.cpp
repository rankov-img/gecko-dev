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
    MOZ_ASSUME_UNREACHABLE("NYI");
}

MoveEmitterMIPS::~MoveEmitterMIPS()
{
    assertDone();
}

Address
MoveEmitterMIPS::cycleSlot() const
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return Address(StackPointer, 0);
}

int32_t MoveEmitterMIPS::getAdjustedOffset(const MoveOperand &operand)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return 0;
}

Register
MoveEmitterMIPS::tempReg()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return t8;
}

void
MoveEmitterMIPS::breakCycle(const MoveOperand &from, const MoveOperand &to, MoveOp::Type type)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::completeCycle(const MoveOperand &from, const MoveOperand &to, MoveOp::Type type)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::emitMove(const MoveOperand &from, const MoveOperand &to)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::emitFloat32Move(const MoveOperand &from, const MoveOperand &to)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::emitDoubleMove(const MoveOperand &from, const MoveOperand &to)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::emit(const MoveOp &move)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
MoveEmitterMIPS::assertDone()
{
    JS_ASSERT(!inCycle_);
}

void
MoveEmitterMIPS::finish()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}
