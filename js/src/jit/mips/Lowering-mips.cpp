/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "mozilla/MathAlgorithms.h"

#include "jit/mips/Assembler-mips.h"
#include "jit/Lowering.h"
#include "jit/MIR.h"

#include "jit/shared/Lowering-shared-inl.h"

using namespace js;
using namespace js::jit;

using mozilla::FloorLog2;

bool
LIRGeneratorMIPS::useBox(LInstruction *lir, size_t n, MDefinition *mir,
                         LUse::Policy policy, bool useAtStart)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::useBoxFixed(LInstruction *lir, size_t n, MDefinition *mir, Register reg1,
                              Register reg2)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

LAllocation
LIRGeneratorMIPS::useByteOpRegister(MDefinition *mir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return useRegister(mir);
}

LAllocation
LIRGeneratorMIPS::useByteOpRegisterOrNonDoubleConstant(MDefinition *mir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return useRegisterOrNonDoubleConstant(mir);
}

bool
LIRGeneratorMIPS::lowerConstantDouble(double d, MInstruction *mir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerConstantFloat32(float d, MInstruction *mir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return define(new(alloc()) LFloat32(d), mir);
}

bool
LIRGeneratorMIPS::visitConstant(MConstant *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitBox(MBox *box)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitUnbox(MUnbox *unbox)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitReturn(MReturn *ret)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

// x = !y
bool
LIRGeneratorMIPS::lowerForALU(LInstructionHelper<1, 1, 0> *ins,
                              MDefinition *mir, MDefinition *input)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

// z = x+y
bool
LIRGeneratorMIPS::lowerForALU(LInstructionHelper<1, 2, 0> *ins, MDefinition *mir,
                              MDefinition *lhs, MDefinition *rhs)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerForFPU(LInstructionHelper<1, 1, 0> *ins, MDefinition *mir,
                              MDefinition *input)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerForFPU(LInstructionHelper<1, 2, 0> *ins, MDefinition *mir,
                              MDefinition *lhs, MDefinition *rhs)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerForBitAndAndBranch(LBitAndAndBranch *baab, MInstruction *mir,
                                          MDefinition *lhs, MDefinition *rhs)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::defineUntypedPhi(MPhi *phi, size_t lirIndex)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

void
LIRGeneratorMIPS::lowerUntypedPhiInput(MPhi *phi, uint32_t inputPosition,
                                       LBlock *block, size_t lirIndex)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
LIRGeneratorMIPS::lowerForShift(LInstructionHelper<1, 2, 0> *ins, MDefinition *mir,
                                MDefinition *lhs, MDefinition *rhs)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerDivI(MDiv *div)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerMulI(MMul *mul, MDefinition *lhs, MDefinition *rhs)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerModI(MMod *mod)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitPowHalf(MPowHalf *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

LTableSwitch *
LIRGeneratorMIPS::newLTableSwitch(const LAllocation &in, const LDefinition &inputCopy,
                                  MTableSwitch *tableswitch)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return new(alloc()) LTableSwitch(in, inputCopy, temp(), tableswitch);
}

LTableSwitchV *
LIRGeneratorMIPS::newLTableSwitchV(MTableSwitch *tableswitch)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return new(alloc()) LTableSwitchV(temp(), tempFloat32(), temp(), tableswitch);
}

bool
LIRGeneratorMIPS::visitGuardShape(MGuardShape *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitGuardObjectType(MGuardObjectType *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerUrshD(MUrsh *mir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSNeg(MAsmJSNeg *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerUDiv(MDiv *div)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerUMod(MMod *mod)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSUnsignedToDouble(MAsmJSUnsignedToDouble *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSUnsignedToFloat32(MAsmJSUnsignedToFloat32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSLoadHeap(MAsmJSLoadHeap *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSStoreHeap(MAsmJSStoreHeap *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::visitAsmJSLoadFuncPtr(MAsmJSLoadFuncPtr *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerTruncateDToInt32(MTruncateToInt32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
LIRGeneratorMIPS::lowerTruncateFToInt32(MTruncateToInt32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


bool
LIRGeneratorMIPS::visitStoreTypedArrayElementStatic(MStoreTypedArrayElementStatic *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

