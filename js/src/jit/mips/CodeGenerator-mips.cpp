/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jit/mips/CodeGenerator-mips.h"

#include "mozilla/MathAlgorithms.h"

#include "jscntxt.h"
#include "jscompartment.h"
#include "jsnum.h"

#include "jit/CodeGenerator.h"
#include "jit/IonFrames.h"
#include "jit/JitCompartment.h"
#include "jit/MIR.h"
#include "jit/MIRGraph.h"
#include "vm/Shape.h"

#include "jsscriptinlines.h"

#include "jit/shared/CodeGenerator-shared-inl.h"

using namespace js;
using namespace js::jit;

using mozilla::FloorLog2;
using mozilla::NegativeInfinity;
using JS::GenericNaN;

// shared
CodeGeneratorMIPS::CodeGeneratorMIPS(MIRGenerator *gen, LIRGraph *graph, MacroAssembler *masm)
  : CodeGeneratorShared(gen, graph, masm)
{
}

bool
CodeGeneratorMIPS::generatePrologue()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::generateEpilogue()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

template void
CodeGeneratorMIPS::emitBranch(Register lhs, Register rhs, Assembler::Condition cond,
                              MBasicBlock *mirTrue, MBasicBlock *mirFalse);

template void
CodeGeneratorMIPS::emitBranch(Register lhs, Imm32 rhs, Assembler::Condition cond,
                              MBasicBlock *mirTrue, MBasicBlock *mirFalse);

template <typename T> void
CodeGeneratorMIPS::emitBranch(Register lhs, T rhs, Assembler::Condition cond,
                              MBasicBlock *mirTrue, MBasicBlock *mirFalse)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

template <typename T> void
CodeGeneratorMIPS::branchToBlock(Register lhs, T rhs, MBasicBlock *mir, Assembler::Condition cond)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

void
CodeGeneratorMIPS::branchToBlock(FloatRegister lhs, FloatRegister rhs, MBasicBlock *mir,
                                 Assembler::DoubleCondition cond, bool isDouble)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
OutOfLineBailout::accept(CodeGeneratorMIPS *codegen)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return codegen->visitOutOfLineBailout(this);
}

bool
CodeGeneratorMIPS::visitTestIAndBranch(LTestIAndBranch *test)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompare(LCompare *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareAndBranch(LCompareAndBranch *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::generateOutOfLineCode()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

template bool
CodeGeneratorMIPS::bailoutIf(Register lhs, Register rhs,
                             Assembler::Condition c, LSnapshot *snapshot);

template bool
CodeGeneratorMIPS::bailoutIf(Register lhs, Imm32 rhs,
                             Assembler::Condition c, LSnapshot *snapshot);

template <typename T1, typename T2> bool
CodeGeneratorMIPS::bailoutIf(T1 lhs, T2 rhs, Assembler::Condition c, LSnapshot *snapshot)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::bailoutFrom(Label *label, LSnapshot *snapshot)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::bailout(LSnapshot *snapshot)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitOutOfLineBailout(OutOfLineBailout *ool)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitMinMaxD(LMinMaxD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAbsD(LAbsD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAbsF(LAbsF *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitSqrtD(LSqrtD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitSqrtF(LSqrtF *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAddI(LAddI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitSubI(LSubI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitMulI(LMulI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitDivI(LDivI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitDivPowTwoI(LDivPowTwoI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitModI(LModI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitModPowTwoI(LModPowTwoI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitModMaskI(LModMaskI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}
bool
CodeGeneratorMIPS::visitBitNotI(LBitNotI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitBitOpI(LBitOpI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitShiftI(LShiftI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitUrshD(LUrshD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitPowHalfD(LPowHalfD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

MoveOperand
CodeGeneratorMIPS::toMoveOperand(const LAllocation *a) const
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return MoveOperand();
}


class js::jit::OutOfLineTableSwitch : public OutOfLineCodeBase<CodeGeneratorMIPS>
{
    MTableSwitch *mir_;
    CodeLabel jumpLabel_;

    bool accept(CodeGeneratorMIPS *codegen) {
        MOZ_ASSUME_UNREACHABLE("NYI");
        return codegen->visitOutOfLineTableSwitch(this);
    }

  public:
    OutOfLineTableSwitch(MTableSwitch *mir)
      : mir_(mir)
    {}

    MTableSwitch *mir() const {
        MOZ_ASSUME_UNREACHABLE("NYI");
        return mir_;
    }

    CodeLabel *jumpLabel() {
        MOZ_ASSUME_UNREACHABLE("NYI");
        return &jumpLabel_;
    }
};

bool
CodeGeneratorMIPS::visitOutOfLineTableSwitch(OutOfLineTableSwitch *ool)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::emitTableSwitchDispatch(MTableSwitch *mir, const Register &index,
                                           const Register &address)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitMathD(LMathD *math)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitMathF(LMathF *math)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitFloor(LFloor *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitFloorF(LFloorF *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitRound(LRound *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitTruncateDToInt32(LTruncateDToInt32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return emitTruncateDouble(ToFloatRegister(ins->input()), ToRegister(ins->output()));
}

bool
CodeGeneratorMIPS::visitTruncateFToInt32(LTruncateFToInt32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return emitTruncateFloat32(ToFloatRegister(ins->input()), ToRegister(ins->output()));
}

static const uint32_t FrameSizes[] = { 128, 256, 512, 1024 };

FrameSizeClass
FrameSizeClass::FromDepth(uint32_t frameDepth)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return FrameSizeClass::None();
}

FrameSizeClass
FrameSizeClass::ClassLimit()
{
    return FrameSizeClass(JS_ARRAY_LENGTH(FrameSizes));
}

uint32_t
FrameSizeClass::frameSize() const
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return FrameSizes[class_];
}

ValueOperand
CodeGeneratorMIPS::ToValue(LInstruction *ins, size_t pos)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return ValueOperand();
}

ValueOperand
CodeGeneratorMIPS::ToOutValue(LInstruction *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return ValueOperand();
}

ValueOperand
CodeGeneratorMIPS::ToTempValue(LInstruction *ins, size_t pos)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return ValueOperand();
}

bool
CodeGeneratorMIPS::visitValue(LValue *value)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}



bool
CodeGeneratorMIPS::visitBox(LBox *box)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitBoxFloatingPoint(LBoxFloatingPoint *box)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitUnbox(LUnbox *unbox)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitDouble(LDouble *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitFloat32(LFloat32 *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return(true);
}

Register
CodeGeneratorMIPS::splitTagForTest(const ValueOperand &value)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return value.typeReg();
}

bool
CodeGeneratorMIPS::visitTestDAndBranch(LTestDAndBranch *test)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitTestFAndBranch(LTestFAndBranch *test)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareD(LCompareD *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareF(LCompareF *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


bool
CodeGeneratorMIPS::visitCompareDAndBranch(LCompareDAndBranch *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareFAndBranch(LCompareFAndBranch *comp)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareB(LCompareB *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareBAndBranch(LCompareBAndBranch *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareV(LCompareV *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitCompareVAndBranch(LCompareVAndBranch *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitBitAndAndBranch(LBitAndAndBranch *baab)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSUInt32ToDouble(LAsmJSUInt32ToDouble *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSUInt32ToFloat32(LAsmJSUInt32ToFloat32 *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitNotI(LNotI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitNotD(LNotD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitNotF(LNotF *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitLoadSlotV(LLoadSlotV *load)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitLoadSlotT(LLoadSlotT *load)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitStoreSlotT(LStoreSlotT *store)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitLoadElementT(LLoadElementT *load)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

void
CodeGeneratorMIPS::storeElementTyped(const LAllocation *value, MIRType valueType,
                                     MIRType elementType, const Register &elements,
                                     const LAllocation *index)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
CodeGeneratorMIPS::visitGuardShape(LGuardShape *guard)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitGuardObjectType(LGuardObjectType *guard)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;

}

bool
CodeGeneratorMIPS::visitGuardClass(LGuardClass *guard)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitImplicitThis(LImplicitThis *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitInterruptCheck(LInterruptCheck *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::generateInvalidateEpilogue()
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

void
DispatchIonCache::initializeAddCacheState(LInstruction *ins, AddCacheState *addState)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
CodeGeneratorMIPS::visitLoadTypedArrayElementStatic(LLoadTypedArrayElementStatic *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitStoreTypedArrayElementStatic(LStoreTypedArrayElementStatic *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadHeap(LAsmJSLoadHeap *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSStoreHeap(LAsmJSStoreHeap *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSPassStackArg(LAsmJSPassStackArg *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitUDiv(LUDiv *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitUMod(LUMod *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitEffectiveAddress(LEffectiveAddress *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadGlobalVar(LAsmJSLoadGlobalVar *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSStoreGlobalVar(LAsmJSStoreGlobalVar *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadFuncPtr(LAsmJSLoadFuncPtr *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadFFIFunc(LAsmJSLoadFFIFunc *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitNegI(LNegI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGeneratorMIPS::visitNegD(LNegD *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


// Methods that have been copied from CodeGenerator.cpp because they need
// special implementation for MIPS.

bool
CodeGenerator::visitAbsI(LAbsI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


bool
CodeGenerator::visitBoundsCheck(LBoundsCheck *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGenerator::visitBoundsCheckRange(LBoundsCheckRange *lir)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGenerator::visitBoundsCheckLower(LBoundsCheckLower *lir)
{
     MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGenerator::emitStoreHoleCheck(Register elements, const LAllocation *index, LSnapshot *snapshot)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
CodeGenerator::visitMinMaxI(LMinMaxI *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}
