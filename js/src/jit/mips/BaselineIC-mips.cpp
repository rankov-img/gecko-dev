/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jsiter.h"

#include "jit/BaselineCompiler.h"
#include "jit/BaselineHelpers.h"
#include "jit/BaselineIC.h"
#include "jit/BaselineJIT.h"
#include "jit/IonLinker.h"

#include "jsboolinlines.h"

using namespace js;
using namespace js::jit;

namespace js {
namespace jit {

// ICCompare_Int32

bool
ICCompare_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    // Guard that R0 is an integer and R1 is an integer.
    Label failure;
    Label conditionTrue;
    masm.branchTestInt32(Assembler::NotEqual, R0, &failure);
    masm.branchTestInt32(Assembler::NotEqual, R1, &failure);

    // Compare payload regs of R0 and R1.
    Assembler::Condition cond = JSOpToCondition(op, /* signed = */true);
    masm.ma_cmp_set(R0.payloadReg(), R0.payloadReg(), R1.payloadReg(), cond);

    masm.tagValue(JSVAL_TYPE_BOOLEAN, R0.payloadReg(), R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);

    return true;
}

bool
ICCompare_Double::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure, isNaN;
    masm.ensureDouble(R0, FloatReg0, &failure);
    masm.ensureDouble(R1, FloatReg1, &failure);

    Register dest = R0.scratchReg();

    Assembler::DoubleCondition doubleCond = JSOpToDoubleCondition(op);

    masm.ma_cmp_set_double(dest, FloatReg0, FloatReg1, doubleCond);

    masm.tagValue(JSVAL_TYPE_BOOLEAN, dest, R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

// ICBinaryArith_Int32

bool
ICBinaryArith_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    // Guard that R0 is an integer and R1 is an integer.
    Label failure;
    masm.branchTestInt32(Assembler::NotEqual, R0, &failure);
    masm.branchTestInt32(Assembler::NotEqual, R1, &failure);

    // Add R0 and R1.  Don't need to explicitly unbox, just use R2's payloadReg.
    Register scratchReg = R2.payloadReg();

    // DIV and MOD need an extra non-volatile ValueOperand to hold R0.
    GeneralRegisterSet savedRegs = availableGeneralRegs(2);
    savedRegs = GeneralRegisterSet::Intersect(GeneralRegisterSet::NonVolatile(), savedRegs);
    ValueOperand savedValue = savedRegs.takeAnyValue();

    Label goodMul, divTest1, divTest2;
    switch(op_) {
      case JSOP_ADD:
        // We know R0.typeReg() already contains the integer tag. No boxing
        // required.
        masm.ma_addTestOverflow(R0.payloadReg(), R0.payloadReg(), R1.payloadReg(), &failure);
        break;
      case JSOP_SUB:
        masm.ma_subTestOverflow(R0.payloadReg(), R0.payloadReg(), R1.payloadReg(), &failure);
        break;
      case JSOP_MUL: {
        masm.ma_mul_branch_overflow(scratchReg, R0.payloadReg(), R1.payloadReg(), &failure);

        masm.ma_b(scratchReg, Imm32(0), &goodMul, Assembler::NotEqual, true);

        // Result is -0 if operands have differnet signs.
        masm.as_xor(t8, R0.payloadReg(), R1.payloadReg());
        masm.ma_b(t8, Imm32(0), &failure, Assembler::LessThan, true);

        masm.bind(&goodMul);
        masm.ma_move(R0.payloadReg(), scratchReg);
        break;
      }
      case JSOP_DIV:
      case JSOP_MOD: {
        // Check for INT_MIN / -1, it results in a double.
        masm.ma_b(R0.payloadReg(), Imm32(INT_MIN), &divTest1, Assembler::NotEqual, true);
        masm.ma_b(R1.payloadReg(), Imm32(-1), &failure, Assembler::Equal, true);
        masm.bind(&divTest1);

        // Check for division by zero
        masm.ma_b(R1.payloadReg(), Imm32(0), &failure, Assembler::Equal, true);

        // Check for 0 / X with X < 0 (results in -0).
        masm.ma_b(R0.payloadReg(), Imm32(0), &divTest2, Assembler::NotEqual, true);
        masm.ma_b(R1.payloadReg(), Imm32(0), &failure, Assembler::LessThan, true);
        masm.bind(&divTest2);

        masm.as_div(R0.payloadReg(), R1.payloadReg());

        if (op_ == JSOP_DIV) {
            // Result is a double if the remainder != 0.
            masm.as_mfhi(scratchReg);
            masm.branch32(Assembler::NotEqual, scratchReg, Imm32(0), &failure);
            masm.as_mflo(scratchReg);
            masm.tagValue(JSVAL_TYPE_INT32, scratchReg, R0);
        } else {
            Label done;
            // If X % Y == 0 and X < 0, the result is -0.
            masm.as_mfhi(scratchReg);
            masm.branch32(Assembler::NotEqual, scratchReg, Imm32(0), &done);
            masm.branch32(Assembler::LessThan, R0.payloadReg(), Imm32(0), &failure);
            masm.bind(&done);
            masm.tagValue(JSVAL_TYPE_INT32, scratchReg, R0);
        }
        break;
      }
      case JSOP_BITOR:
        masm.ma_or(R0.payloadReg() , R0.payloadReg(), R1.payloadReg());
        break;
      case JSOP_BITXOR:
        masm.ma_xor(R0.payloadReg() , R0.payloadReg(), R1.payloadReg());
        break;
      case JSOP_BITAND:
        masm.ma_and(R0.payloadReg() , R0.payloadReg(), R1.payloadReg());
        break;
      case JSOP_LSH:
        // MIPS will happily try to shift by more than 0x1f.
        masm.ma_sll(R0.payloadReg(), R0.payloadReg(), R1.payloadReg());
        break;
      case JSOP_RSH:
        masm.ma_sra(R0.payloadReg(), R0.payloadReg(), R1.payloadReg());
        break;
      case JSOP_URSH:
        masm.ma_srl(scratchReg, R0.payloadReg(), R1.payloadReg());
        if (allowDouble_) {
            Label toUint;
            masm.ma_b(scratchReg, Imm32(0), &toUint, Assembler::LessThan, true);

            // Move result and box for return.
            masm.ma_move(R0.payloadReg(), scratchReg);
            EmitReturnFromIC(masm);

            masm.bind(&toUint);
            // NOTE: I'm not sure if we can use FloatReg0, but we can't use
            // ScratchFloatReg like ARM and X86 do.
            masm.convertUInt32ToDouble(scratchReg, FloatReg1);
            masm.boxDouble(FloatReg1, R0);
        } else {
            masm.ma_b(scratchReg, Imm32(0), &failure, Assembler::LessThan, true);
            // Move result for return.
            masm.ma_move(R0.payloadReg(), scratchReg);
        }
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("Unhandled op for BinaryArith_Int32.");
    }

    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);

    return true;
}

bool
ICUnaryArith_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    masm.branchTestInt32(Assembler::NotEqual, R0, &failure);

    switch (op) {
      case JSOP_BITNOT:
        masm.ma_not(R0.payloadReg(), R0.payloadReg());
        break;
      case JSOP_NEG:
        // Guard against 0 and MIN_INT, both result in a double.
        masm.branchTest32(Assembler::Zero, R0.payloadReg(), Imm32(0x7fffffff), &failure);

        // Compile -x as 0 - x.
        masm.ma_negu(R0.payloadReg(), R0.payloadReg());
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("Unexpected op");
        return false;
    }

    EmitReturnFromIC(masm);

    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}


// Mehtods that have been copied from BaselineIC.cpp because they need special
// implementation for MIPS.

bool
ICBinaryArith_BooleanWithInt32::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    if (lhsIsBool_)
        masm.branchTestBoolean(Assembler::NotEqual, R0, &failure);
    else
        masm.branchTestInt32(Assembler::NotEqual, R0, &failure);

    if (rhsIsBool_)
        masm.branchTestBoolean(Assembler::NotEqual, R1, &failure);
    else
        masm.branchTestInt32(Assembler::NotEqual, R1, &failure);

    Register lhsReg = lhsIsBool_ ? masm.extractBoolean(R0, ExtractTemp0)
                      : masm.extractInt32(R0, ExtractTemp0);
    Register rhsReg = rhsIsBool_ ? masm.extractBoolean(R1, ExtractTemp1)
                      : masm.extractInt32(R1, ExtractTemp1);

    JS_ASSERT(op_ == JSOP_ADD || op_ == JSOP_SUB ||
              op_ == JSOP_BITOR || op_ == JSOP_BITXOR || op_ == JSOP_BITAND);

    switch(op_) {
      case JSOP_ADD: {
        Label fixOverflow;

        masm.ma_addTestOverflow(lhsReg, lhsReg, rhsReg, &failure);
        masm.tagValue(JSVAL_TYPE_INT32, lhsReg, R0);
        EmitReturnFromIC(masm);

        break;
      }
      case JSOP_SUB: {
        Label fixOverflow;

        masm.ma_subTestOverflow(lhsReg, lhsReg, rhsReg, &failure);
        masm.tagValue(JSVAL_TYPE_INT32, lhsReg, R0);
        EmitReturnFromIC(masm);

        break;
      }
      case JSOP_BITOR: {
        masm.orPtr(rhsReg, lhsReg);
        masm.tagValue(JSVAL_TYPE_INT32, lhsReg, R0);
        EmitReturnFromIC(masm);
        break;
      }
      case JSOP_BITXOR: {
        masm.xorPtr(rhsReg, lhsReg);
        masm.tagValue(JSVAL_TYPE_INT32, lhsReg, R0);
        EmitReturnFromIC(masm);
        break;
      }
      case JSOP_BITAND: {
        masm.andPtr(rhsReg, lhsReg);
        masm.tagValue(JSVAL_TYPE_INT32, lhsReg, R0);
        EmitReturnFromIC(masm);
        break;
      }
      default:
       MOZ_ASSUME_UNREACHABLE("Unhandled op for BinaryArith_BooleanWithInt32.");
       return false;
    }

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// Compare_Boolean
//

bool
ICCompare_Boolean::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    masm.branchTestBoolean(Assembler::NotEqual, R0, &failure);
    masm.branchTestBoolean(Assembler::NotEqual, R1, &failure);

    Register left = masm.extractInt32(R0, ExtractTemp0);
    Register right = masm.extractInt32(R1, ExtractTemp1);

    // Compare payload regs of R0 and R1.
    Assembler::Condition cond = JSOpToCondition(op, /* signed = */true);
    // MIPS specific implementation
    masm.ma_cmp_set(left, left, right, cond);

    // Box the result and return
    masm.tagValue(JSVAL_TYPE_BOOLEAN, left, R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// Compare_Int32WithBoolean
//

bool
ICCompare_Int32WithBoolean::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    ValueOperand int32Val;
    ValueOperand boolVal;
    if (lhsIsInt32_) {
        int32Val = R0;
        boolVal = R1;
    } else {
        boolVal = R0;
        int32Val = R1;
    }
    masm.branchTestBoolean(Assembler::NotEqual, boolVal, &failure);
    masm.branchTestInt32(Assembler::NotEqual, int32Val, &failure);

    if (op_ == JSOP_STRICTEQ || op_ == JSOP_STRICTNE) {
        // Ints and booleans are never strictly equal, always strictly not
        // equal.
        masm.moveValue(BooleanValue(op_ == JSOP_STRICTNE), R0);
        EmitReturnFromIC(masm);
    } else {
        Register boolReg = masm.extractBoolean(boolVal, ExtractTemp0);
        Register int32Reg = masm.extractInt32(int32Val, ExtractTemp1);

        // Compare payload regs of R0 and R1.
        Assembler::Condition cond = JSOpToCondition(op_, /* signed = */true);
        // MIPS specific implementation
        masm.ma_cmp_set(R0.scratchReg(), lhsIsInt32_ ? int32Reg : boolReg,
                        lhsIsInt32_ ? boolReg : int32Reg, cond);

        // Box the result and return
        masm.tagValue(JSVAL_TYPE_BOOLEAN, R0.scratchReg(), R0);
        EmitReturnFromIC(masm);
    }

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// ToBool_Int32
//

bool
ICToBool_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    masm.branchTestInt32(Assembler::NotEqual, R0, &failure);

    Label ifFalse;
    // MIPS specific implementation
    masm.branchTestInt32Truthy(false, R0, &ifFalse);

    masm.moveValue(BooleanValue(true), R0);
    EmitReturnFromIC(masm);

    masm.bind(&ifFalse);
    masm.moveValue(BooleanValue(false), R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// ToBool_String
//

bool
ICToBool_String::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;
    masm.branchTestString(Assembler::NotEqual, R0, &failure);

    Label ifFalse;
    // MIPS specific implementation
    masm.branchTestStringTruthy(false, R0, &ifFalse);

    masm.moveValue(BooleanValue(true), R0);
    EmitReturnFromIC(masm);

    masm.bind(&ifFalse);
    masm.moveValue(BooleanValue(false), R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// ToBool_Double
//

bool
ICToBool_Double::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure, ifTrue;
    masm.branchTestDouble(Assembler::NotEqual, R0, &failure);
    masm.unboxDouble(R0, FloatReg0);

    // MIPS specific implementation
    masm.branchTestDoubleTruthy(true, FloatReg0, &ifTrue);

    masm.moveValue(BooleanValue(false), R0);
    EmitReturnFromIC(masm);

    masm.bind(&ifTrue);
    masm.moveValue(BooleanValue(true), R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// ToBool_Object
//

bool
ICToBool_Object::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure, ifFalse, slowPath;
    masm.branchTestObject(Assembler::NotEqual, R0, &failure);

    Register objReg = masm.extractObject(R0, ExtractTemp0);
    Register scratch = R1.scratchReg();
    // MIPS specific implementation
    masm.branchTestObjectTruthy(false, objReg, scratch, &slowPath, &ifFalse);

    // If object doesn't emulate undefined, it evaulates to true.
    masm.moveValue(BooleanValue(true), R0);
    EmitReturnFromIC(masm);

    masm.bind(&ifFalse);
    masm.moveValue(BooleanValue(false), R0);
    EmitReturnFromIC(masm);

    masm.bind(&slowPath);
    masm.setupUnalignedABICall(1, scratch);
    masm.passABIArg(objReg);
    masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, js::EmulatesUndefined));
    masm.convertBoolToInt32(ReturnReg, ReturnReg);
    masm.xor32(Imm32(1), ReturnReg);
    masm.tagValue(JSVAL_TYPE_BOOLEAN, ReturnReg, R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}

//
// IteratorMore_Native
//

bool
ICIteratorMore_Native::Compiler::generateStubCode(MacroAssembler &masm)
{
    Label failure;

    Register obj = masm.extractObject(R0, ExtractTemp0);

    GeneralRegisterSet regs(availableGeneralRegs(1));
    Register nativeIterator = regs.takeAny();
    Register scratch = regs.takeAny();

    masm.branchTestObjClass(Assembler::NotEqual, obj, scratch,
                            &PropertyIteratorObject::class_, &failure);
    masm.loadObjPrivate(obj, JSObject::ITER_CLASS_NFIXED_SLOTS, nativeIterator);

    masm.branchTest32(Assembler::NonZero, Address(nativeIterator, offsetof(NativeIterator, flags)),
                      Imm32(JSITER_FOREACH), &failure);

    // Set output to true if props_cursor < props_end.
    masm.loadPtr(Address(nativeIterator, offsetof(NativeIterator, props_end)), scratch);
    // MIPS specific implementation
    Address addr = Address(nativeIterator, offsetof(NativeIterator, props_cursor));
    masm.ma_cmp_set(scratch, addr, scratch, Assembler::LessThan);

    masm.tagValue(JSVAL_TYPE_BOOLEAN, scratch, R0);
    EmitReturnFromIC(masm);

    // Failure case - jump to next stub
    masm.bind(&failure);
    EmitStubGuardFailure(masm);
    return true;
}


} // namespace jit
} // namespace js
