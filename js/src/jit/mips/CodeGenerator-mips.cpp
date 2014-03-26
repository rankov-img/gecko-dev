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
    if (gen->compilingAsmJS()) {
        masm.Push(ra);
        // Note that this automatically sets MacroAssembler::framePushed().
        masm.reserveStack(frameDepth_);
    } else {
        // Note that this automatically sets MacroAssembler::framePushed().
        masm.reserveStack(frameSize());
        masm.checkStackAlignment();
    }

    return true;
}

bool
CodeGeneratorMIPS::generateEpilogue()
{
    masm.bind(&returnLabel_);
#if JS_TRACE_LOGGING
    masm.tracelogStop();
#endif
    if (gen->compilingAsmJS()) {
        // Pop the stack we allocated at the start of the function.
        masm.freeStack(frameDepth_);
        masm.Pop(ra);
        masm.as_jr(ra);
        masm.as_nop();
        MOZ_ASSERT(masm.framePushed() == 0);
    } else {
        // Pop the stack we allocated at the start of the function.
        masm.freeStack(frameSize());
        MOZ_ASSERT(masm.framePushed() == 0);
        masm.ret();
    }
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
    if (isNextBlock(mirFalse->lir())) {
        branchToBlock(lhs, rhs, mirTrue, cond);
    } else {
        branchToBlock(lhs, rhs, mirFalse, Assembler::InvertCondition(cond));
        jumpToBlock(mirTrue);
    }
}

template <typename T> void
CodeGeneratorMIPS::branchToBlock(Register lhs, T rhs, MBasicBlock *mir, Assembler::Condition cond)
{
    if (Label *oolEntry = labelForBackedgeWithImplicitCheck(mir)) {
        // Note: the backedge is initially a jump to the next instruction.
        // It will be patched to the target block's label during link().
        RepatchLabel rejoin;

        CodeOffsetJump backedge;
        Label skip;
        masm.ma_b(lhs, rhs, &skip, Assembler::InvertCondition(cond), ShortJump);

        backedge = masm.jumpWithPatch(&rejoin);
        masm.bind(&rejoin);
        masm.bind(&skip);

        if (!patchableBackedges_.append(PatchableBackedgeInfo(backedge, mir->lir()->label(),
                                        oolEntry)))
            MOZ_CRASH();
    } else {
        masm.ma_b(lhs, rhs, mir->lir()->label(), cond);
    }
}

void
CodeGeneratorMIPS::branchToBlock(FloatRegister lhs, FloatRegister rhs, MBasicBlock *mir,
                                 Assembler::DoubleCondition cond, bool isDouble)
{
    if (Label *oolEntry = labelForBackedgeWithImplicitCheck(mir)) {
        // Note: the backedge is initially a jump to the next instruction.
        // It will be patched to the target block's label during link().
        RepatchLabel rejoin;

        CodeOffsetJump backedge;
        Label skip;
        if (isDouble)
            masm.ma_bc1d(lhs, rhs, &skip, Assembler::InvertCondition(cond), ShortJump);
        else
            masm.ma_bc1s(lhs, rhs, &skip, Assembler::InvertCondition(cond), ShortJump);

        backedge = masm.jumpWithPatch(&rejoin);
        masm.bind(&rejoin);
        masm.bind(&skip);

        if (!patchableBackedges_.append(PatchableBackedgeInfo(backedge, mir->lir()->label(),
                                        oolEntry)))
            MOZ_CRASH();
    } else {
        if (isDouble)
            masm.ma_bc1d(lhs, rhs, mir->lir()->label(), cond);
        else
            masm.ma_bc1s(lhs, rhs, mir->lir()->label(), cond);
    }
}

bool
OutOfLineBailout::accept(CodeGeneratorMIPS *codegen)
{
    return codegen->visitOutOfLineBailout(this);
}

bool
CodeGeneratorMIPS::visitTestIAndBranch(LTestIAndBranch *test)
{
    const LAllocation *opd = test->getOperand(0);
    MBasicBlock *ifTrue = test->ifTrue();
    MBasicBlock *ifFalse = test->ifFalse();

    emitBranch(ToRegister(opd), Imm32(0), Assembler::NonZero, ifTrue, ifFalse);
    return true;
}

bool
CodeGeneratorMIPS::visitCompare(LCompare *comp)
{
    Assembler::Condition cond = JSOpToCondition(comp->mir()->compareType(), comp->jsop());
    const LAllocation *left = comp->getOperand(0);
    const LAllocation *right = comp->getOperand(1);
    const LDefinition *def = comp->getDef(0);

    if (right->isConstant())
        masm.ma_cmp_set(ToRegister(def), ToRegister(left), Imm32(ToInt32(right)), cond);
    else if (right->isGeneralReg())
        masm.ma_cmp_set(ToRegister(def), ToRegister(left), ToRegister(right), cond);
    else
        masm.ma_cmp_set(ToRegister(def), ToRegister(left), ToAddress(right), cond);

    return true;
}

bool
CodeGeneratorMIPS::visitCompareAndBranch(LCompareAndBranch *comp)
{
    Assembler::Condition cond = JSOpToCondition(comp->cmpMir()->compareType(), comp->jsop());
    if (comp->right()->isConstant()) {
        emitBranch(ToRegister(comp->left()), Imm32(ToInt32(comp->right())), cond,
                   comp->ifTrue(), comp->ifFalse());
    } else if (comp->right()->isGeneralReg()) {
        emitBranch(ToRegister(comp->left()), ToRegister(comp->right()), cond,
                   comp->ifTrue(), comp->ifFalse());
    } else {
        emitBranch(ToRegister(comp->left()), ToAddress(comp->right()), cond,
                   comp->ifTrue(), comp->ifFalse());
    }

    return true;
}

bool
CodeGeneratorMIPS::generateOutOfLineCode()
{
    if (!CodeGeneratorShared::generateOutOfLineCode())
        return false;

    if (deoptLabel_.used()) {
        // All non-table-based bailouts will go here.
        masm.bind(&deoptLabel_);

        // Push the frame size, so the handler can recover the IonScript.
        // Frame size is stored in 'ra' and pushed by GenerateBailoutThunk
        // We have to use 'ra' because generateBailoutTable will implicitly do
        // the same.
        masm.ma_li(ra, Imm32(frameSize()));

        JitCode *handler = gen->jitRuntime()->getGenericBailoutHandler();

        masm.branch(handler);
    }

    return true;
}

bool
CodeGeneratorMIPS::bailoutFrom(Label *label, LSnapshot *snapshot)
{
    if (masm.bailed())
        return false;
    MOZ_ASSERT(label->used());
    MOZ_ASSERT(!label->bound());

    CompileInfo &info = snapshot->mir()->block()->info();
    switch (info.executionMode()) {
      case ParallelExecution: {
        // in parallel mode, make no attempt to recover, just signal an error.
        OutOfLineAbortPar *ool = oolAbortPar(ParallelBailoutUnsupported,
                                             snapshot->mir()->block(),
                                             snapshot->mir()->pc());
        masm.retarget(label, ool->entry());
        return true;
      }
      case SequentialExecution:
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("No such execution mode");
    }

    if (!encode(snapshot))
        return false;

    // Though the assembler doesn't track all frame pushes, at least make sure
    // the known value makes sense. We can't use bailout tables if the stack
    // isn't properly aligned to the static frame size.
    MOZ_ASSERT_IF(frameClass_ != FrameSizeClass::None(),
                 frameClass_.frameSize() == masm.framePushed());

    // We don't use table bailouts because retargeting is easier this way.
    OutOfLineBailout *ool = new(alloc()) OutOfLineBailout(snapshot, masm.framePushed());
    if (!addOutOfLineCode(ool)) {
        return false;
    }

    masm.retarget(label, ool->entry());

    return true;
}

bool
CodeGeneratorMIPS::bailout(LSnapshot *snapshot)
{
    Label label;
    masm.ma_b(&label);
    return bailoutFrom(&label, snapshot);
}

bool
CodeGeneratorMIPS::visitOutOfLineBailout(OutOfLineBailout *ool)
{
    // Push snapshotOffset and make sure stack is aligned.
    masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(void *)));
    masm.ma_sw(Imm32(ool->snapshot()->snapshotOffset()), Address(StackPointer, 0));

    masm.ma_b(&deoptLabel_);
    return true;
}

bool
CodeGeneratorMIPS::visitMinMaxD(LMinMaxD *ins)
{
    FloatRegister first = ToFloatRegister(ins->first());
    FloatRegister second = ToFloatRegister(ins->second());
    FloatRegister output = ToFloatRegister(ins->output());

    MOZ_ASSERT(first == output);

    Assembler::DoubleCondition cond = ins->mir()->isMax()
                                      ? Assembler::DoubleLessThanOrEqual
                                      : Assembler::DoubleGreaterThanOrEqual;
    Label nan, equal, returnSecond, done;

    // First or second is NaN, result is NaN.
    masm.ma_bc1d(first, second, &nan, Assembler::DoubleUnordered, ShortJump);
    // Make sure we handle -0 and 0 right.
    masm.ma_bc1d(first, second, &equal, Assembler::DoubleEqual, ShortJump);
    masm.ma_bc1d(first, second, &returnSecond, cond, ShortJump);
    masm.ma_b(&done, ShortJump);

    // Check for zero.
    masm.bind(&equal);
    masm.ma_lid(ScratchFloatReg, 0.0);
    // First wasn't 0 or -0, so just return it.
    masm.ma_bc1d(first, ScratchFloatReg, &done, Assembler::DoubleNotEqualOrUnordered, ShortJump);

    // So now both operands are either -0 or 0.
    if (ins->mir()->isMax()) {
        // -0 + -0 = -0 and -0 + 0 = 0.
        masm.addDouble(second, first);
    } else {
        masm.negateDouble(first);
        masm.subDouble(second, first);
        masm.negateDouble(first);
    }
    masm.ma_b(&done, ShortJump);

    masm.bind(&nan);
    masm.loadConstantDouble(GenericNaN(), output);
    masm.ma_b(&done, ShortJump);

    masm.bind(&returnSecond);
    masm.as_movd(output, second);

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitAbsD(LAbsD *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    MOZ_ASSERT(input == ToFloatRegister(ins->output()));
    masm.as_absd(input, input);
    return true;
}

bool
CodeGeneratorMIPS::visitAbsF(LAbsF *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    MOZ_ASSERT(input == ToFloatRegister(ins->output()));
    masm.as_abss(input, input);
    return true;
}

bool
CodeGeneratorMIPS::visitSqrtD(LSqrtD *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    FloatRegister output = ToFloatRegister(ins->output());
    masm.as_sqrtd(output, input);
    return true;
}

bool
CodeGeneratorMIPS::visitSqrtF(LSqrtF *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    FloatRegister output = ToFloatRegister(ins->output());
    masm.as_sqrts(output, input);
    return true;
}

bool
CodeGeneratorMIPS::visitAddI(LAddI *ins)
{
    const LAllocation *lhs = ins->getOperand(0);
    const LAllocation *rhs = ins->getOperand(1);
    const LDefinition *dest = ins->getDef(0);

    MOZ_ASSERT(rhs->isConstant() || rhs->isGeneralReg());

    // If there is no snapshot, we don't need to check for overflow
    if (!ins->snapshot()) {
        if (rhs->isConstant())
            masm.ma_addu(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)));
        else
            masm.as_addu(ToRegister(dest), ToRegister(lhs), ToRegister(rhs));
        return true;
    }

    Label overflow;
    if (rhs->isConstant())
        masm.ma_addTestOverflow(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)), &overflow);
    else
        masm.ma_addTestOverflow(ToRegister(dest), ToRegister(lhs), ToRegister(rhs), &overflow);

    if (!bailoutFrom(&overflow, ins->snapshot()))
        return false;

    return true;
}

bool
CodeGeneratorMIPS::visitSubI(LSubI *ins)
{
    const LAllocation *lhs = ins->getOperand(0);
    const LAllocation *rhs = ins->getOperand(1);
    const LDefinition *dest = ins->getDef(0);

    MOZ_ASSERT(rhs->isConstant() || rhs->isGeneralReg());

    // If there is no snapshot, we don't need to check for overflow
    if (!ins->snapshot()) {
        if (rhs->isConstant())
            masm.ma_subu(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)));
        else
            masm.as_subu(ToRegister(dest), ToRegister(lhs), ToRegister(rhs));
        return true;
    }

    Label overflow;
    if (rhs->isConstant())
        masm.ma_subTestOverflow(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)), &overflow);
    else
        masm.ma_subTestOverflow(ToRegister(dest), ToRegister(lhs), ToRegister(rhs), &overflow);

    if (!bailoutFrom(&overflow, ins->snapshot()))
        return false;

    return true;
}

bool
CodeGeneratorMIPS::visitMulI(LMulI *ins)
{
    const LAllocation *lhs = ins->lhs();
    const LAllocation *rhs = ins->rhs();
    Register dest = ToRegister(ins->output());
    MMul *mul = ins->mir();

    MOZ_ASSERT_IF(mul->mode() == MMul::Integer, !mul->canBeNegativeZero() && !mul->canOverflow());

    if (rhs->isConstant()) {
        int32_t constant = ToInt32(rhs);
        Register src = ToRegister(lhs);

        // Bailout on -0.0
        if (mul->canBeNegativeZero() && constant <= 0) {
            Label bailNegativeZero;
            masm.ma_b(src, Imm32(0), &bailNegativeZero,
                      (constant == 0) ? Assembler::LessThan : Assembler::Equal);
            if (!bailoutFrom(&bailNegativeZero, ins->snapshot()))
                return false;
        }

        switch (constant) {
          case -1:
            if (mul->canOverflow()) {
                if (!bailoutCmpPtr(Assembler::Equal, src, Imm32(INT32_MIN), ins->snapshot()))
                    return false;
            }
            masm.ma_negu(dest, src);
            break;
          case 0:
            masm.ma_li(dest, Imm32(0));
            break;
          case 1:
            masm.ma_move(dest, src);
            break;
          case 2:
            if (mul->canOverflow()) {
                Label mulTwoOverflow;
                masm.ma_addTestOverflow(dest, src, src, &mulTwoOverflow);

                if (!bailoutFrom(&mulTwoOverflow, ins->snapshot()))
                    return false;
            } else {
                masm.as_addu(dest, src, src);
            }
            break;
          default:
            uint32_t shift = FloorLog2(constant);

            if (!mul->canOverflow() && (constant > 0)) {
                // If it cannot overflow, we can do lots of optimizations.
                uint32_t rest = constant - (1 << shift);

                // See if the constant has one bit set, meaning it can be
                // encoded as a bitshift.
                if ((1 << shift) == constant) {
                    masm.ma_sll(dest, src, Imm32(shift));
                    return true;
                }

                // If the constant cannot be encoded as (1<<C1), see if it can
                // be encoded as (1<<C1) | (1<<C2), which can be computed
                // using an add and a shift.
                uint32_t shift_rest = FloorLog2(rest);
                if (src != dest && (1u << shift_rest) == rest) {
                    masm.ma_sll(dest, src, Imm32(shift - shift_rest));
                    masm.as_addu(dest, dest, src);
                    if (shift_rest != 0)
                        masm.ma_sll(dest, dest, Imm32(shift_rest));
                    return true;
                }
            }

            if (mul->canOverflow() && (constant > 0) && (src != dest)) {
                // To stay on the safe side, only optimize things that are a
                // power of 2.

                if ((1 << shift) == constant) {
                    // dest = lhs * pow(2, shift)
                    masm.ma_sll(dest, src, Imm32(shift));
                    // At runtime, check (lhs == dest >> shift), if this does
                    // not hold, some bits were lost due to overflow, and the
                    // computation should be resumed as a double.
                    masm.ma_sra(ScratchRegister, dest, Imm32(shift));
                    if (!bailoutCmpPtr(Assembler::NotEqual, src, ScratchRegister, ins->snapshot()))
                        return false;
                    return true;
                }
            }

            if (mul->canOverflow()) {
                Label mulConstOverflow;
                masm.ma_mul_branch_overflow(dest, ToRegister(lhs), Imm32(ToInt32(rhs)),
                                            &mulConstOverflow);

                if (!bailoutFrom(&mulConstOverflow, ins->snapshot()))
                    return false;
            } else {
                masm.ma_mult(src, Imm32(ToInt32(rhs)));
                masm.as_mflo(dest);
            }
            break;
        }
    } else {
        Label multRegOverflow;

        if (mul->canOverflow()) {
            masm.ma_mul_branch_overflow(dest, ToRegister(lhs), ToRegister(rhs), &multRegOverflow);
            if (!bailoutFrom(&multRegOverflow, ins->snapshot()))
                return false;
        } else {
            masm.as_mult(ToRegister(lhs), ToRegister(rhs));
            masm.as_mflo(dest);
        }

        if (mul->canBeNegativeZero()) {
            Label done;
            masm.ma_b(dest, dest, &done, Assembler::NonZero, ShortJump);

            // Result is -0 if lhs or rhs is negative.
            // In that case result must be double value so bailout
            Label bailNegativeZero;
            masm.ma_b(ToRegister(lhs), ToRegister(lhs), &bailNegativeZero, Assembler::Signed);
            if (!bailoutFrom(&bailNegativeZero, ins->snapshot()))
                return false;

            masm.ma_b(ToRegister(rhs), ToRegister(rhs), &bailNegativeZero, Assembler::Signed);
            if (!bailoutFrom(&bailNegativeZero, ins->snapshot()))
                return false;
            masm.bind(&done);
        }
    }

    return true;
}

bool
CodeGeneratorMIPS::visitDivI(LDivI *ins)
{
    // Extract the registers from this instruction
    Register lhs = ToRegister(ins->lhs());
    Register rhs = ToRegister(ins->rhs());
    Register dest = ToRegister(ins->output());
    Register temp = ToRegister(ins->getTemp(0));
    MDiv *mir = ins->mir();

    Label done;

    // Handle divide by zero.
    if (mir->canBeDivideByZero() || mir->canBeNegativeZero()) {
        if (mir->isTruncated()) {
            // Truncated division by zero is zero (Infinity|0 == 0)
            Label notzero;
            masm.ma_b(rhs, rhs, &notzero, Assembler::NonZero, ShortJump);
            masm.ma_li(dest, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&notzero);
        } else {
            MOZ_ASSERT(mir->fallible());
            Label divideByZero;
            masm.ma_b(rhs, rhs, &divideByZero, Assembler::Zero);
            if (!bailoutFrom(&divideByZero, ins->snapshot()))
                return false;
        }
    }

    // Handle an integer overflow exception from -2147483648 / -1.
    if (mir->canBeNegativeOverflow()) {
        Label notMinInt;
        masm.ma_li(temp, Imm32(INT32_MIN));
        masm.ma_b(lhs, temp, &notMinInt, Assembler::NotEqual, ShortJump);

        masm.ma_li(temp, Imm32(-1));
        if (mir->isTruncated()) {
            // (-INT32_MIN)|0 == INT32_MIN
            Label skip;
            masm.ma_b(rhs, temp, &skip, Assembler::NotEqual, ShortJump);
            masm.ma_li(dest, Imm32(INT32_MIN));
            masm.ma_b(&done, ShortJump);
            masm.bind(&skip);
        } else {
            MOZ_ASSERT(mir->fallible());
            Label negativeOverflow;
            masm.ma_b(rhs, temp, &negativeOverflow, Assembler::Equal);
            if (!bailoutFrom(&negativeOverflow, ins->snapshot()))
                return false;
        }
        masm.bind(&notMinInt);
    }

    // Handle negative 0. (0/-Y)
    if (!mir->isTruncated() && mir->canBeNegativeZero()) {
        Label nonzero;
        masm.ma_b(lhs, lhs, &nonzero, Assembler::NonZero, ShortJump);

        Label negativeZero;
        masm.ma_b(rhs, Imm32(0), &negativeZero, Assembler::LessThan);
        if (!bailoutFrom(&negativeZero, ins->snapshot()))
            return false;
        masm.bind(&nonzero);
    }
    // Note: above safety checks could not be verified as Ion seems to be
    // smarter and requires double arithmetic in such cases.

    // All regular. Lets call div.
    if (mir->isTruncated()) {
        masm.as_div(lhs, rhs);
        masm.as_mflo(dest);
    } else {
        MOZ_ASSERT(mir->fallible());

        Label remainderNonZero;
        masm.ma_div_branch_overflow(dest, lhs, rhs, &remainderNonZero);
        if (!bailoutFrom(&remainderNonZero, ins->snapshot()))
            return false;
    }

    masm.bind(&done);

    return true;
}

bool
CodeGeneratorMIPS::visitDivPowTwoI(LDivPowTwoI *ins)
{
    Register lhs = ToRegister(ins->numerator());
    Register dest = ToRegister(ins->output());
    Register tmp = ToRegister(ins->getTemp(0));
    int32_t shift = ins->shift();

    if (shift != 0) {
        MDiv *mir = ins->mir();
        if (!mir->isTruncated()) {
            // If the remainder is going to be != 0, bailout since this must
            // be a double.
            masm.ma_sll(tmp, lhs, Imm32(32 - shift));

            Label bail;
            masm.ma_b(tmp, tmp, &bail, Assembler::NonZero);
            if (!bailoutFrom(&bail, ins->snapshot()))
                return false;
        }

        if (!mir->canBeNegativeDividend()) {
            // Numerator is unsigned, so needs no adjusting. Do the shift.
            masm.ma_sra(dest, lhs, Imm32(shift));
            return true;
        }

        // Adjust the value so that shifting produces a correctly rounded result
        // when the numerator is negative. See 10-1 "Signed Division by a Known
        // Power of 2" in Henry S. Warren, Jr.'s Hacker's Delight.
        if (shift > 1) {
            masm.ma_sra(tmp, lhs, Imm32(31));
            masm.ma_srl(tmp, tmp, Imm32(32 - shift));
            masm.ma_addu(tmp, lhs);
        } else {
            masm.ma_srl(tmp, lhs, Imm32(32 - shift));
            masm.ma_addu(tmp, lhs);
        }

        // Do the shift.
        masm.ma_sra(dest, tmp, Imm32(shift));
    } else {
        masm.ma_move(dest, lhs);
    }

    return true;

}

bool
CodeGeneratorMIPS::visitModI(LModI *ins)
{
    // Extract the registers from this instruction
    Register lhs = ToRegister(ins->lhs());
    Register rhs = ToRegister(ins->rhs());
    Register dest = ToRegister(ins->output());
    Register callTemp = ToRegister(ins->callTemp());
    MMod *mir = ins->mir();
    Label done, prevent;

    masm.ma_move(callTemp, lhs);

    // Prevent INT_MIN % -1;
    // The integer division will give INT_MIN, but we want -(double)INT_MIN.
    if (mir->canBeNegativeDividend()) {
        masm.ma_b(lhs, Imm32(INT_MIN), &prevent, Assembler::NotEqual, ShortJump);
        if (mir->isTruncated()) {
            // (INT_MIN % -1)|0 == 0
            Label skip;
            masm.ma_b(rhs, Imm32(-1), &skip, Assembler::NotEqual, ShortJump);
            masm.ma_li(dest, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&skip);
        } else {
            MOZ_ASSERT(mir->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, rhs, Imm32(-1), ins->snapshot()))
                return false;
        }
        masm.bind(&prevent);
    }

    // 0/X (with X < 0) is bad because both of these values *should* be
    // doubles, and the result should be -0.0, which cannot be represented in
    // integers. X/0 is bad because it will give garbage (or abort), when it
    // should give either \infty, -\infty or NAN.

    // Prevent 0 / X (with X < 0) and X / 0
    // testing X / Y.  Compare Y with 0.
    // There are three cases: (Y < 0), (Y == 0) and (Y > 0)
    // If (Y < 0), then we compare X with 0, and bail if X == 0
    // If (Y == 0), then we simply want to bail.
    // if (Y > 0), we don't bail.

    if (mir->canBeDivideByZero()) {
        if (mir->isTruncated()) {
            Label skip;
            masm.ma_b(rhs, Imm32(0), &skip, Assembler::NotEqual, ShortJump);
            masm.ma_li(dest, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&skip);
        } else {
            MOZ_ASSERT(mir->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, rhs, Imm32(0), ins->snapshot()))
                return false;
        }
    }

    if (mir->canBeNegativeDividend()) {
        Label notNegative;
        masm.ma_b(rhs, Imm32(0), &notNegative, Assembler::GreaterThan, ShortJump);
        if (mir->isTruncated()) {
            // NaN|0 == 0 and (0 % -X)|0 == 0
            Label skip;
            masm.ma_b(lhs, Imm32(0), &skip, Assembler::NotEqual, ShortJump);
            masm.ma_li(dest, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&skip);
        } else {
            MOZ_ASSERT(mir->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, lhs, Imm32(0), ins->snapshot()))
                return false;
        }
        masm.bind(&notNegative);
    }

    masm.as_div(lhs, rhs);
    masm.as_mfhi(dest);

    // If X%Y == 0 and X < 0, then we *actually* wanted to return -0.0
    if (mir->canBeNegativeDividend()) {
        if (mir->isTruncated()) {
            // -0.0|0 == 0
        } else {
            MOZ_ASSERT(mir->fallible());
            // See if X < 0
            masm.ma_b(dest, Imm32(0), &done, Assembler::NotEqual, ShortJump);
            if (!bailoutCmpPtr(Assembler::Signed, callTemp, Imm32(0), ins->snapshot()))
                return false;
        }
    }
    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitModPowTwoI(LModPowTwoI *ins)
{
    Register in = ToRegister(ins->getOperand(0));
    Register out = ToRegister(ins->getDef(0));
    MMod *mir = ins->mir();
    Label negative, done;

    masm.ma_move(out, in);
    masm.ma_b(in, in, &done, Assembler::Zero, ShortJump);
    // Switch based on sign of the lhs.
    // Positive numbers are just a bitmask
    masm.ma_b(in, in, &negative, Assembler::Signed, ShortJump);
    {
        masm.ma_and(out, Imm32((1 << ins->shift()) - 1));
        masm.ma_b(&done, ShortJump);
    }

    // Negative numbers need a negate, bitmask, negate
    {
        masm.bind(&negative);
        masm.ma_negu(out, out);
        masm.ma_and(out, Imm32((1 << ins->shift()) - 1));
        masm.ma_negu(out, out);
    }
    if (mir->canBeNegativeDividend()) {
        if (!mir->isTruncated()) {
            MOZ_ASSERT(mir->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, out, zero, ins->snapshot()))
                return false;
        } else {
            // -0|0 == 0
        }
    }
    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitModMaskI(LModMaskI *ins)
{
    Register src = ToRegister(ins->getOperand(0));
    Register dest = ToRegister(ins->getDef(0));
    Register tmp = ToRegister(ins->getTemp(0));
    MMod *mir = ins->mir();

    if (!mir->isTruncated() && mir->canBeNegativeDividend()) {
        MOZ_ASSERT(mir->fallible());

        Label bail;
        masm.ma_mod_mask(src, dest, tmp, ins->shift(), &bail);
        if (!bailoutFrom(&bail, ins->snapshot()))
            return false;
    } else {
        masm.ma_mod_mask(src, dest, tmp, ins->shift(), nullptr);
    }
    return true;
}
bool
CodeGeneratorMIPS::visitBitNotI(LBitNotI *ins)
{
    const LAllocation *input = ins->getOperand(0);
    const LDefinition *dest = ins->getDef(0);
    MOZ_ASSERT(!input->isConstant());

    masm.ma_not(ToRegister(dest), ToRegister(input));
    return true;
}

bool
CodeGeneratorMIPS::visitBitOpI(LBitOpI *ins)
{
    const LAllocation *lhs = ins->getOperand(0);
    const LAllocation *rhs = ins->getOperand(1);
    const LDefinition *dest = ins->getDef(0);
    // all of these bitops should be either imm32's, or integer registers.
    switch (ins->bitop()) {
      case JSOP_BITOR:
        if (rhs->isConstant())
            masm.ma_or(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)));
        else
            masm.ma_or(ToRegister(dest), ToRegister(lhs), ToRegister(rhs));
        break;
      case JSOP_BITXOR:
        if (rhs->isConstant())
            masm.ma_xor(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)));
        else
            masm.ma_xor(ToRegister(dest), ToRegister(lhs), ToRegister(rhs));
        break;
      case JSOP_BITAND:
        if (rhs->isConstant())
            masm.ma_and(ToRegister(dest), ToRegister(lhs), Imm32(ToInt32(rhs)));
        else
            masm.ma_and(ToRegister(dest), ToRegister(lhs), ToRegister(rhs));
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("unexpected binary opcode");
    }

    return true;
}

bool
CodeGeneratorMIPS::visitShiftI(LShiftI *ins)
{
    Register lhs = ToRegister(ins->lhs());
    const LAllocation *rhs = ins->rhs();
    Register dest = ToRegister(ins->output());

    if (rhs->isConstant()) {
        int32_t shift = ToInt32(rhs) & 0x1F;
        switch (ins->bitop()) {
          case JSOP_LSH:
            if (shift)
                masm.ma_sll(dest, lhs, Imm32(shift));
            else
                masm.ma_move(dest, lhs);
            break;
          case JSOP_RSH:
            if (shift)
                masm.ma_sra(dest, lhs, Imm32(shift));
            else
                masm.ma_move(dest, lhs);
            break;
          case JSOP_URSH:
            if (shift) {
                masm.ma_srl(dest, lhs, Imm32(shift));
            } else {
                // x >>> 0 can overflow.
                masm.ma_move(dest, lhs);
                if (ins->mir()->toUrsh()->fallible()) {
                    Label bail;
                    masm.ma_b(dest, Imm32(0), &bail, Assembler::LessThan);
                    if (!bailoutFrom(&bail, ins->snapshot()))
                        return false;
                }
            }
            break;
          default:
            MOZ_ASSUME_UNREACHABLE("Unexpected shift op");
        }
    } else {
        // The shift amounts should be AND'ed into the 0-31 range
        masm.ma_and(dest, ToRegister(rhs), Imm32(0x1F));

        switch (ins->bitop()) {
          case JSOP_LSH:
            masm.ma_sll(dest, lhs, dest);
            break;
          case JSOP_RSH:
            masm.ma_sra(dest, lhs, dest);
            break;
          case JSOP_URSH:
            masm.ma_srl(dest, lhs, dest);
            if (ins->mir()->toUrsh()->fallible()) {
                // x >>> 0 can overflow.
                Label bail;
                masm.ma_b(dest, Imm32(0), &bail, Assembler::LessThan);
                if (!bailoutFrom(&bail, ins->snapshot()))
                    return false;
            }
            break;
          default:
            MOZ_ASSUME_UNREACHABLE("Unexpected shift op");
        }
    }

    return true;
}

bool
CodeGeneratorMIPS::visitUrshD(LUrshD *ins)
{
    Register lhs = ToRegister(ins->lhs());
    Register temp = ToRegister(ins->temp());

    const LAllocation *rhs = ins->rhs();
    FloatRegister out = ToFloatRegister(ins->output());

    if (rhs->isConstant()) {
        masm.ma_srl(temp, lhs, Imm32(ToInt32(rhs)));
    } else {
        masm.ma_srl(temp, lhs, ToRegister(rhs));
    }

    masm.convertUInt32ToDouble(temp, out);
    return true;
}

bool
CodeGeneratorMIPS::visitPowHalfD(LPowHalfD *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    FloatRegister output = ToFloatRegister(ins->output());

    Label done, skip;

    // Masm.pow(-Infinity, 0.5) == Infinity.
    masm.ma_lid(ScratchFloatReg, NegativeInfinity<double>());
    masm.ma_bc1d(input, ScratchFloatReg, &skip, Assembler::DoubleNotEqualOrUnordered, ShortJump);
    masm.as_negd(output, ScratchFloatReg);
    masm.ma_b(&done, ShortJump);

    masm.bind(&skip);
    // Math.pow(-0, 0.5) == 0 == Math.pow(0, 0.5).
    // Adding 0 converts any -0 to 0.
    masm.ma_lid(ScratchFloatReg, 0.0);
    masm.as_addd(output, input, ScratchFloatReg);
    masm.as_sqrtd(output, output);

    masm.bind(&done);
    return true;
}

MoveOperand
CodeGeneratorMIPS::toMoveOperand(const LAllocation *a) const
{
    if (a->isGeneralReg())
        return MoveOperand(ToRegister(a));
    if (a->isFloatReg()) {
        return MoveOperand(ToFloatRegister(a));
    }
    MOZ_ASSERT((ToStackOffset(a) & 3) == 0);
    int32_t offset = ToStackOffset(a);

    // The way the stack slots work, we assume that everything from
    // depth == 0 downwards is writable. However, since our frame is included
    // in this, ensure that the frame gets skipped.
    if (gen->compilingAsmJS())
        offset -= AlignmentMidPrologue;

    return MoveOperand(StackPointer, offset);
}

class js::jit::OutOfLineTableSwitch : public OutOfLineCodeBase<CodeGeneratorMIPS>
{
    MTableSwitch *mir_;
    CodeLabel jumpLabel_;

    bool accept(CodeGeneratorMIPS *codegen) {
        return codegen->visitOutOfLineTableSwitch(this);
    }

  public:
    OutOfLineTableSwitch(MTableSwitch *mir)
      : mir_(mir)
    {}

    MTableSwitch *mir() const {
        return mir_;
    }

    CodeLabel *jumpLabel() {
        return &jumpLabel_;
    }
};

bool
CodeGeneratorMIPS::visitOutOfLineTableSwitch(OutOfLineTableSwitch *ool)
{
    MTableSwitch *mir = ool->mir();

    masm.align(sizeof(void*));
    masm.bind(ool->jumpLabel()->src());
    if (!masm.addCodeLabel(*ool->jumpLabel()))
        return false;

    for (size_t i = 0; i < mir->numCases(); i++) {
        LBlock *caseblock = mir->getCase(i)->lir();
        Label *caseheader = caseblock->label();
        uint32_t caseoffset = caseheader->offset();

        // The entries of the jump table need to be absolute addresses and thus
        // must be patched after codegen is finished.
        CodeLabel cl;
        masm.ma_li(ScratchRegister, cl.dest());
        masm.as_jr(ScratchRegister);
        masm.as_nop();
        cl.src()->bind(caseoffset);
        if (!masm.addCodeLabel(cl))
            return false;
    }

    return true;
}

bool
CodeGeneratorMIPS::emitTableSwitchDispatch(MTableSwitch *mir, const Register &index,
                                           const Register &address)
{
    Label *defaultcase = mir->getDefault()->lir()->label();

    // Lower value with low value
    if (mir->low() != 0)
        masm.ma_subu(index, Imm32(mir->low()));

    // Jump to default case if input is out of range
    int32_t cases = mir->numCases();
    masm.ma_b(index, Imm32(cases), defaultcase, Assembler::AboveOrEqual);

    // To fill in the CodeLabels for the case entries, we need to first
    // generate the case entries (we don't yet know their offsets in the
    // instruction stream).
    OutOfLineTableSwitch *ool = new(alloc()) OutOfLineTableSwitch(mir);
    if (!addOutOfLineCode(ool))
        return false;

    // Compute the position where a pointer to the right case stands.
    masm.ma_li(address, ool->jumpLabel()->dest());
    masm.ma_sll(index, index, Imm32(4));
    masm.as_addu(address, address, index);

    masm.as_jr(address);
    masm.as_nop();
    return true;
}

bool
CodeGeneratorMIPS::visitMathD(LMathD *math)
{
    const LAllocation *src1 = math->getOperand(0);
    const LAllocation *src2 = math->getOperand(1);
    const LDefinition *output = math->getDef(0);

    switch (math->jsop()) {
      case JSOP_ADD:
        masm.as_addd(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_SUB:
        masm.as_subd(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_MUL:
        masm.as_muld(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_DIV:
        masm.as_divd(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("unexpected opcode");
    }
    return true;
}

bool
CodeGeneratorMIPS::visitMathF(LMathF *math)
{
    const LAllocation *src1 = math->getOperand(0);
    const LAllocation *src2 = math->getOperand(1);
    const LDefinition *output = math->getDef(0);

    switch (math->jsop()) {
      case JSOP_ADD:
        masm.as_adds(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_SUB:
        masm.as_subs(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_MUL:
        masm.as_muls(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      case JSOP_DIV:
        masm.as_divs(ToFloatRegister(output), ToFloatRegister(src1), ToFloatRegister(src2));
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("unexpected opcode");
    }
    return true;
}

bool
CodeGeneratorMIPS::visitFloor(LFloor *lir)
{
    FloatRegister input = ToFloatRegister(lir->input());
    FloatRegister scratch = ScratchFloatReg;
    Register output = ToRegister(lir->output());

    Label skipCheck, done, bail;

    // If Nan, 0 or -0 check for bailout
    masm.ma_lid(scratch, 0.0);
    masm.ma_bc1d(input, scratch, &skipCheck, Assembler::DoubleNotEqual, ShortJump);

    // If high part is not zero, it is NaN or -0, so we bail.
    masm.moveFromDoubleHi(input, SecondScratchReg);
    masm.ma_b(SecondScratchReg, Imm32(0), &bail, Assembler::NotEqual, ShortJump);
    if (!bailoutFrom(&bail, lir->snapshot()))
        return false;

    // Input was zero, so return zero.
    masm.ma_li(output, Imm32(0));
    masm.ma_b(&done, ShortJump);

    masm.bind(&skipCheck);
    masm.as_floorwd(scratch, input);
    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MAX), lir->snapshot()))
        return false;

    masm.bind(&done);

    return true;
}

bool
CodeGeneratorMIPS::visitFloorF(LFloorF *lir)
{
    FloatRegister input = ToFloatRegister(lir->input());
    FloatRegister scratch = ScratchFloatReg;
    Register output = ToRegister(lir->output());

    Label skipCheck, done, bail;

    // If Nan, 0 or -0 check for bailout
    masm.ma_lis(scratch, 0.0);
    masm.ma_bc1s(input, scratch, &skipCheck, Assembler::DoubleNotEqual, ShortJump);

    // If binary value is not zero, it is NaN or -0, so we bail.
    masm.as_mfc1(SecondScratchReg, input);
    masm.ma_b(SecondScratchReg, Imm32(0), &bail, Assembler::NotEqual, ShortJump);
    if (!bailoutFrom(&bail, lir->snapshot()))
        return false;

    // Input was zero, so return zero.
    masm.ma_li(output, Imm32(0));
    masm.ma_b(&done, ShortJump);

    masm.bind(&skipCheck);
    masm.as_floorws(scratch, input);
    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MAX), lir->snapshot()))
        return false;

    masm.bind(&done);

    return true;
}

bool
CodeGeneratorMIPS::visitRound(LRound *lir)
{
    FloatRegister input = ToFloatRegister(lir->input());
    FloatRegister temp = ToFloatRegister(lir->temp());
    FloatRegister scratch = ScratchFloatReg;
    Register output = ToRegister(lir->output());

    Label bail1, bail2, negative, end, skipCheck;

    // Load 0.5 in the temp register.
    masm.loadConstantDouble(0.5, temp);

    // Branch to a slow path for negative inputs. Doesn't catch NaN or -0.
    masm.ma_lid(scratch, 0.0);
    masm.ma_bc1d(input, scratch, &negative, Assembler::DoubleLessThan, ShortJump);

    // If Nan, 0 or -0 check for bailout
    masm.ma_bc1d(input, scratch, &skipCheck, Assembler::DoubleNotEqual, ShortJump);

    // If high part is not zero, it is NaN or -0, so we bail.
    masm.moveFromDoubleHi(input, SecondScratchReg);
    masm.ma_b(SecondScratchReg, Imm32(0), &bail1, Assembler::NotEqual, ShortJump);
    if (!bailoutFrom(&bail1, lir->snapshot()))
        return false;

    // Input was zero, so return zero.
    masm.ma_li(output, Imm32(0));
    masm.ma_b(&end, ShortJump);

    masm.bind(&skipCheck);
    masm.ma_lid(scratch, 0.5);
    masm.as_addd(scratch, input, scratch);
    masm.as_floorwd(scratch, scratch);

    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MAX), lir->snapshot()))
        return false;

    masm.jump(&end);

    // Input is negative, but isn't -0.
    masm.bind(&negative);
    masm.as_addd(temp, input, temp);

    // If input + 0.5 >= 0, input is a negative number >= -0.5 and the
    // result is -0.
    masm.ma_bc1d(temp, scratch, &bail2, Assembler::DoubleGreaterThanOrEqual);
    if (!bailoutFrom(&bail2, lir->snapshot()))
        return false;

    // Truncate and round toward zero.
    // This is off-by-one for everything but integer-valued inputs.
    masm.as_floorwd(scratch, temp);
    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    masm.bind(&end);
    return true;
}

bool
CodeGeneratorMIPS::visitRoundF(LRoundF *lir)
{
    FloatRegister input = ToFloatRegister(lir->input());
    FloatRegister temp = ToFloatRegister(lir->temp());
    FloatRegister scratch = ScratchFloatReg;
    Register output = ToRegister(lir->output());

    Label bail1, bail2, negative, end, skipCheck;

    // Load 0.5 in the temp register.
    masm.loadConstantFloat32(0.5, temp);

    // Branch to a slow path for negative inputs. Doesn't catch NaN or -0.
    masm.loadConstantFloat32(0.0, scratch);
    masm.ma_bc1s(input, scratch, &negative, Assembler::DoubleLessThan, ShortJump);

    // If Nan, 0 or -0 check for bailout
    masm.ma_bc1s(input, scratch, &skipCheck, Assembler::DoubleNotEqual, ShortJump);

    // If binary value is not zero, it is NaN or -0, so we bail.
    masm.as_mfc1(SecondScratchReg, input);
    masm.ma_b(SecondScratchReg, Imm32(0), &bail1, Assembler::NotEqual, ShortJump);
    if (!bailoutFrom(&bail1, lir->snapshot()))
        return false;

    // Input was zero, so return zero.
    masm.ma_li(output, Imm32(0));
    masm.ma_b(&end, ShortJump);

    masm.bind(&skipCheck);
    masm.ma_lis(scratch, 0.5);
    masm.as_adds(scratch, input, scratch);
    masm.as_floorws(scratch, scratch);

    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MAX), lir->snapshot()))
        return false;

    masm.jump(&end);

    // Input is negative, but isn't -0.
    masm.bind(&negative);
    masm.as_adds(temp, input, temp);

    // If input + 0.5 >= 0, input is a negative number >= -0.5 and the
    // result is -0.
    masm.ma_bc1s(temp, scratch, &bail2, Assembler::DoubleGreaterThanOrEqual);
    if (!bailoutFrom(&bail2, lir->snapshot()))
        return false;

    // Truncate and round toward zero.
    // This is off-by-one for everything but integer-valued inputs.
    masm.as_floorws(scratch, temp);
    masm.as_mfc1(output, scratch);

    if (!bailoutCmpPtr(Assembler::Equal, output, Imm32(INT_MIN), lir->snapshot()))
        return false;

    masm.bind(&end);
    return true;
}

bool
CodeGeneratorMIPS::visitTruncateDToInt32(LTruncateDToInt32 *ins)
{
    return emitTruncateDouble(ToFloatRegister(ins->input()), ToRegister(ins->output()));
}

bool
CodeGeneratorMIPS::visitTruncateFToInt32(LTruncateFToInt32 *ins)
{
    return emitTruncateFloat32(ToFloatRegister(ins->input()), ToRegister(ins->output()));
}

static const uint32_t FrameSizes[] = { 128, 256, 512, 1024 };

FrameSizeClass
FrameSizeClass::FromDepth(uint32_t frameDepth)
{
    for (uint32_t i = 0; i < JS_ARRAY_LENGTH(FrameSizes); i++) {
        if (frameDepth < FrameSizes[i])
            return FrameSizeClass(i);
    }

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
    MOZ_ASSERT(class_ != NO_FRAME_SIZE_CLASS_ID);
    MOZ_ASSERT(class_ < JS_ARRAY_LENGTH(FrameSizes));

    return FrameSizes[class_];
}

ValueOperand
CodeGeneratorMIPS::ToValue(LInstruction *ins, size_t pos)
{
    Register typeReg = ToRegister(ins->getOperand(pos + TYPE_INDEX));
    Register payloadReg = ToRegister(ins->getOperand(pos + PAYLOAD_INDEX));
    return ValueOperand(typeReg, payloadReg);
}

ValueOperand
CodeGeneratorMIPS::ToOutValue(LInstruction *ins)
{
    Register typeReg = ToRegister(ins->getDef(TYPE_INDEX));
    Register payloadReg = ToRegister(ins->getDef(PAYLOAD_INDEX));
    return ValueOperand(typeReg, payloadReg);
}

ValueOperand
CodeGeneratorMIPS::ToTempValue(LInstruction *ins, size_t pos)
{
    Register typeReg = ToRegister(ins->getTemp(pos + TYPE_INDEX));
    Register payloadReg = ToRegister(ins->getTemp(pos + PAYLOAD_INDEX));
    return ValueOperand(typeReg, payloadReg);
}

bool
CodeGeneratorMIPS::visitValue(LValue *value)
{
    const ValueOperand out = ToOutValue(value);

    masm.moveValue(value->value(), out);
    return true;
}

bool
CodeGeneratorMIPS::visitBox(LBox *box)
{
    const LDefinition *type = box->getDef(TYPE_INDEX);

    MOZ_ASSERT(!box->getOperand(0)->isConstant());

    // On x86, the input operand and the output payload have the same
    // virtual register. All that needs to be written is the type tag for
    // the type definition.
    masm.ma_li(ToRegister(type), Imm32(MIRTypeToTag(box->type())));
    return true;
}

bool
CodeGeneratorMIPS::visitBoxFloatingPoint(LBoxFloatingPoint *box)
{
    const LDefinition *payload = box->getDef(PAYLOAD_INDEX);
    const LDefinition *type = box->getDef(TYPE_INDEX);
    const LAllocation *in = box->getOperand(0);

    FloatRegister reg = ToFloatRegister(in);
    if (box->type() == MIRType_Float32) {
        masm.convertFloat32ToDouble(reg, ScratchFloatReg);
        reg = ScratchFloatReg;
    }
    masm.ma_mv(reg, ValueOperand(ToRegister(type), ToRegister(payload)));
    return true;
}

bool
CodeGeneratorMIPS::visitUnbox(LUnbox *unbox)
{
    // Note that for unbox, the type and payload indexes are switched on the
    // inputs.
    MUnbox *mir = unbox->mir();
    Register type = ToRegister(unbox->type());

    if (mir->fallible()) {
        if (!bailoutCmpPtr(Assembler::NotEqual, type, Imm32(MIRTypeToTag(mir->type())),
                           unbox->snapshot()))
            return false;
    }
    return true;
}

bool
CodeGeneratorMIPS::visitDouble(LDouble *ins)
{
    const LDefinition *out = ins->getDef(0);

    masm.ma_lid(ToFloatRegister(out), ins->getDouble());
    return true;
}

bool
CodeGeneratorMIPS::visitFloat32(LFloat32 *ins)
{
    const LDefinition *out = ins->getDef(0);
    masm.loadConstantFloat32(ins->getFloat(), ToFloatRegister(out));
    return true;
}

Register
CodeGeneratorMIPS::splitTagForTest(const ValueOperand &value)
{
    return value.typeReg();
}

bool
CodeGeneratorMIPS::visitTestDAndBranch(LTestDAndBranch *test)
{
    FloatRegister input = ToFloatRegister(test->input());

    MBasicBlock *ifTrue = test->ifTrue();
    MBasicBlock *ifFalse = test->ifFalse();

    masm.ma_lid(ScratchFloatReg, 0.0);
    // If 0, or NaN, the result is false.


    if (isNextBlock(ifFalse->lir())) {
        branchToBlock(input, ScratchFloatReg, ifTrue, Assembler::DoubleNotEqual, ShortJump);
    } else {
        branchToBlock(input, ScratchFloatReg, ifFalse, Assembler::DoubleEqualOrUnordered, ShortJump);
        jumpToBlock(ifTrue);
    }

    return true;
}

bool
CodeGeneratorMIPS::visitTestFAndBranch(LTestFAndBranch *test)
{
    FloatRegister input = ToFloatRegister(test->input());

    MBasicBlock *ifTrue = test->ifTrue();
    MBasicBlock *ifFalse = test->ifFalse();

    masm.ma_lis(ScratchFloatReg, 0.0);
    // If 0, or NaN, the result is false.


    if (isNextBlock(ifFalse->lir())) {
        branchToBlock(input, ScratchFloatReg, ifTrue, Assembler::DoubleNotEqual, false);
    } else {
        branchToBlock(input, ScratchFloatReg, ifFalse, Assembler::DoubleEqualOrUnordered, false);
        jumpToBlock(ifTrue);
    }

    return true;
}

bool
CodeGeneratorMIPS::visitCompareD(LCompareD *comp)
{
    FloatRegister lhs = ToFloatRegister(comp->left());
    FloatRegister rhs = ToFloatRegister(comp->right());
    Register dest = ToRegister(comp->output());

    Assembler::DoubleCondition cond = JSOpToDoubleCondition(comp->mir()->jsop());
    masm.ma_cmp_set_double(dest, lhs, rhs, cond);
    return true;
}

bool
CodeGeneratorMIPS::visitCompareF(LCompareF *comp)
{
    FloatRegister lhs = ToFloatRegister(comp->left());
    FloatRegister rhs = ToFloatRegister(comp->right());
    Register dest = ToRegister(comp->output());

    Assembler::DoubleCondition cond = JSOpToDoubleCondition(comp->mir()->jsop());
    masm.ma_cmp_set_float32(dest, lhs, rhs, cond);
    return true;
}


bool
CodeGeneratorMIPS::visitCompareDAndBranch(LCompareDAndBranch *comp)
{
    FloatRegister lhs = ToFloatRegister(comp->left());
    FloatRegister rhs = ToFloatRegister(comp->right());

    Assembler::DoubleCondition cond = JSOpToDoubleCondition(comp->cmpMir()->jsop());
    MBasicBlock *ifTrue = comp->ifTrue();
    MBasicBlock *ifFalse = comp->ifFalse();

    if (isNextBlock(ifFalse->lir())) {
        branchToBlock(lhs, rhs, ifTrue, cond, ShortJump);
    } else {
        branchToBlock(lhs, rhs, ifFalse, Assembler::InvertCondition(cond), ShortJump);
        jumpToBlock(ifTrue);
    }

    return true;
}

bool
CodeGeneratorMIPS::visitCompareFAndBranch(LCompareFAndBranch *comp)
{
    FloatRegister lhs = ToFloatRegister(comp->left());
    FloatRegister rhs = ToFloatRegister(comp->right());

    Assembler::DoubleCondition cond = JSOpToDoubleCondition(comp->cmpMir()->jsop());
    MBasicBlock *ifTrue = comp->ifTrue();
    MBasicBlock *ifFalse = comp->ifFalse();

    if (isNextBlock(ifFalse->lir())) {
        branchToBlock(lhs, rhs, ifTrue, cond, false);
    } else {
        branchToBlock(lhs, rhs, ifFalse, Assembler::InvertCondition(cond), false);
        jumpToBlock(ifTrue);
    }

    return true;
}

bool
CodeGeneratorMIPS::visitCompareB(LCompareB *lir)
{
    MCompare *mir = lir->mir();

    const ValueOperand lhs = ToValue(lir, LCompareB::Lhs);
    const LAllocation *rhs = lir->rhs();
    const Register output = ToRegister(lir->output());

    MOZ_ASSERT(mir->jsop() == JSOP_STRICTEQ || mir->jsop() == JSOP_STRICTNE);
    Assembler::Condition cond = JSOpToCondition(mir->compareType(), mir->jsop());

    Label notBoolean, done;
    masm.branchTestBoolean(Assembler::NotEqual, lhs, &notBoolean);
    {
        if (rhs->isConstant())
            masm.ma_cmp_set(output, lhs.payloadReg(), Imm32(rhs->toConstant()->toBoolean()), cond);
        else
            masm.ma_cmp_set(output, lhs.payloadReg(), ToRegister(rhs), cond);
        masm.jump(&done);
    }

    masm.bind(&notBoolean);
    {
        masm.move32(Imm32(mir->jsop() == JSOP_STRICTNE), output);
    }

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitCompareBAndBranch(LCompareBAndBranch *lir)
{
    MCompare *mir = lir->cmpMir();
    const ValueOperand lhs = ToValue(lir, LCompareBAndBranch::Lhs);
    const LAllocation *rhs = lir->rhs();

    MOZ_ASSERT(mir->jsop() == JSOP_STRICTEQ || mir->jsop() == JSOP_STRICTNE);

    MBasicBlock *mirNotBoolean = (mir->jsop() == JSOP_STRICTEQ) ? lir->ifFalse() : lir->ifTrue();
    branchToBlock(lhs.typeReg(), ImmType(JSVAL_TYPE_BOOLEAN), mirNotBoolean, Assembler::NotEqual);

    Assembler::Condition cond = JSOpToCondition(mir->compareType(), mir->jsop());
    if (rhs->isConstant())
        emitBranch(lhs.payloadReg(), Imm32(rhs->toConstant()->toBoolean()), cond, lir->ifTrue(),
                   lir->ifFalse());
    else
        emitBranch(lhs.payloadReg(), ToRegister(rhs), cond, lir->ifTrue(), lir->ifFalse());

    return true;
}

bool
CodeGeneratorMIPS::visitCompareV(LCompareV *lir)
{
    MCompare *mir = lir->mir();
    Assembler::Condition cond = JSOpToCondition(mir->compareType(), mir->jsop());
    const ValueOperand lhs = ToValue(lir, LCompareV::LhsInput);
    const ValueOperand rhs = ToValue(lir, LCompareV::RhsInput);
    const Register output = ToRegister(lir->output());

    MOZ_ASSERT(IsEqualityOp(mir->jsop()));

    Label notEqual, done;
    masm.ma_b(lhs.typeReg(), rhs.typeReg(), &notEqual, Assembler::NotEqual, ShortJump);
    {
        masm.ma_cmp_set(output, lhs.payloadReg(), rhs.payloadReg(), cond);
        masm.ma_b(&done, ShortJump);
    }
    masm.bind(&notEqual);
    {
        masm.move32(Imm32(cond == Assembler::NotEqual), output);
    }

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitCompareVAndBranch(LCompareVAndBranch *lir)
{
    MCompare *mir = lir->cmpMir();
    Assembler::Condition cond = JSOpToCondition(mir->compareType(), mir->jsop());
    const ValueOperand lhs = ToValue(lir, LCompareVAndBranch::LhsInput);
    const ValueOperand rhs = ToValue(lir, LCompareVAndBranch::RhsInput);

    MOZ_ASSERT(mir->jsop() == JSOP_EQ || mir->jsop() == JSOP_STRICTEQ ||
              mir->jsop() == JSOP_NE || mir->jsop() == JSOP_STRICTNE);

    MBasicBlock *notEqual = (cond == Assembler::Equal) ? lir->ifFalse() : lir->ifTrue();

    branchToBlock(lhs.typeReg(), rhs.typeReg(), notEqual, Assembler::NotEqual);
    emitBranch(lhs.payloadReg(), rhs.payloadReg(), cond, lir->ifTrue(), lir->ifFalse());

    return true;
}

bool
CodeGeneratorMIPS::visitBitAndAndBranch(LBitAndAndBranch *baab)
{
    if (baab->right()->isConstant())
        masm.ma_and(ScratchRegister, ToRegister(baab->left()), Imm32(ToInt32(baab->right())));
    else
        masm.ma_and(ScratchRegister, ToRegister(baab->left()), ToRegister(baab->right()));
    emitBranch(ScratchRegister, ScratchRegister, Assembler::NonZero, baab->ifTrue(),
               baab->ifFalse());
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSUInt32ToDouble(LAsmJSUInt32ToDouble *lir)
{
    masm.convertUInt32ToDouble(ToRegister(lir->input()), ToFloatRegister(lir->output()));
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSUInt32ToFloat32(LAsmJSUInt32ToFloat32 *lir)
{
    masm.convertUInt32ToFloat32(ToRegister(lir->input()), ToFloatRegister(lir->output()));
    return true;
}

bool
CodeGeneratorMIPS::visitNotI(LNotI *ins)
{
    masm.ma_cmp_set(ToRegister(ins->output()), ToRegister(ins->input()), Imm32(0),
                    Assembler::Equal);
    return true;
}

bool
CodeGeneratorMIPS::visitNotD(LNotD *ins)
{
    // Since this operation is not, we want to set a bit if
    // the double is falsey, which means 0.0, -0.0 or NaN.
    FloatRegister in = ToFloatRegister(ins->input());
    Register dest = ToRegister(ins->output());

    Label falsey, done;
    masm.ma_lid(ScratchFloatReg, 0.0);
    masm.ma_bc1d(in, ScratchFloatReg, &falsey, Assembler::DoubleEqualOrUnordered, ShortJump);

    masm.ma_li(dest, Imm32(0));
    masm.ma_b(&done, ShortJump);

    masm.bind(&falsey);
    masm.ma_li(dest, Imm32(1));

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitNotF(LNotF *ins)
{
    // Since this operation is not, we want to set a bit if
    // the float32 is falsey, which means 0.0, -0.0 or NaN.
    FloatRegister in = ToFloatRegister(ins->input());
    Register dest = ToRegister(ins->output());

    Label falsey, done;
    masm.ma_lis(ScratchFloatReg, 0.0);
    masm.ma_bc1s(in, ScratchFloatReg, &falsey, Assembler::DoubleEqualOrUnordered, ShortJump);

    masm.ma_li(dest, Imm32(0));
    masm.ma_b(&done, ShortJump);

    masm.bind(&falsey);
    masm.ma_li(dest, Imm32(1));

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitLoadSlotV(LLoadSlotV *load)
{
    const ValueOperand out = ToOutValue(load);
    Register base = ToRegister(load->input());
    int32_t offset = load->mir()->slot() * sizeof(js::Value);

    masm.loadValue(Address(base, offset), out);
    return true;
}

bool
CodeGeneratorMIPS::visitLoadSlotT(LLoadSlotT *load)
{
    Register base = ToRegister(load->input());
    int32_t offset = load->mir()->slot() * sizeof(js::Value);

    if (load->mir()->type() == MIRType_Double)
        masm.loadInt32OrDouble(Address(base, offset), ToFloatRegister(load->output()));
    else
        masm.ma_load(ToRegister(load->output()), Address(base, offset + NUNBOX32_PAYLOAD_OFFSET));
    return true;
}

bool
CodeGeneratorMIPS::visitStoreSlotT(LStoreSlotT *store)
{
    Register base = ToRegister(store->slots());
    int32_t offset = store->mir()->slot() * sizeof(js::Value);

    const LAllocation *value = store->value();
    MIRType valueType = store->mir()->value()->type();

    if (store->mir()->needsBarrier())
        emitPreBarrier(Address(base, offset), store->mir()->slotType());

    if (valueType == MIRType_Double) {
        masm.ma_sd(ToFloatRegister(value), Address(base, offset));
        return true;
    }

    // Store the type tag if needed.
    if (valueType != store->mir()->slotType())
        masm.storeTypeTag(ImmType(ValueTypeFromMIRType(valueType)), Address(base, offset));

    // Store the payload.
    if (value->isConstant())
        masm.storePayload(*value->toConstant(), Address(base, offset));
    else
        masm.storePayload(ToRegister(value), Address(base, offset));

    return true;
}

bool
CodeGeneratorMIPS::visitLoadElementT(LLoadElementT *load)
{
    Register base = ToRegister(load->elements());
    if (load->mir()->type() == MIRType_Double) {
        FloatRegister fpreg = ToFloatRegister(load->output());
        if (load->index()->isConstant()) {
            Address source(base, ToInt32(load->index()) * sizeof(Value));
            if (load->mir()->loadDoubles())
                masm.loadDouble(source, fpreg);
            else
                masm.loadInt32OrDouble(source, fpreg);
        } else {
            Register index = ToRegister(load->index());
            if (load->mir()->loadDoubles())
                masm.loadDouble(BaseIndex(base, index, TimesEight), fpreg);
            else
                masm.loadInt32OrDouble(base, index, fpreg);
        }
    } else {
        if (load->index()->isConstant()) {
            Address source(base, ToInt32(load->index()) * sizeof(Value));
            masm.load32(source, ToRegister(load->output()));
        } else {
            BaseIndex source(base, ToRegister(load->index()), TimesEight);
            masm.load32(source, ToRegister(load->output()));
        }
    }
    MOZ_ASSERT(!load->mir()->needsHoleCheck());
    return true;
}

void
CodeGeneratorMIPS::storeElementTyped(const LAllocation *value, MIRType valueType,
                                     MIRType elementType, const Register &elements,
                                     const LAllocation *index)
{
    if (index->isConstant()) {
        Address dest = Address(elements, ToInt32(index) * sizeof(Value));
        if (valueType == MIRType_Double) {
            masm.ma_sd(ToFloatRegister(value), Address(dest.base, dest.offset));
            return;
        }

        // Store the type tag if needed.
        if (valueType != elementType)
            masm.storeTypeTag(ImmType(ValueTypeFromMIRType(valueType)), dest);

        // Store the payload.
        if (value->isConstant())
            masm.storePayload(*value->toConstant(), dest);
        else
            masm.storePayload(ToRegister(value), dest);
    } else {
        Register indexReg = ToRegister(index);
        if (valueType == MIRType_Double) {
            masm.ma_sd(ToFloatRegister(value), BaseIndex(elements, indexReg, TimesEight));
            return;
        }

        // Store the type tag if needed.
        if (valueType != elementType)
            masm.storeTypeTag(ImmType(ValueTypeFromMIRType(valueType)), elements, indexReg);

        // Store the payload.
        if (value->isConstant())
            masm.storePayload(*value->toConstant(), elements, indexReg);
        else
            masm.storePayload(ToRegister(value), elements, indexReg);
    }
}

bool
CodeGeneratorMIPS::visitGuardShape(LGuardShape *guard)
{
    Register obj = ToRegister(guard->input());
    Register tmp = ToRegister(guard->tempInt());

    // Use second scratch. It is safe.
    masm.ma_li(SecondScratchReg, ImmGCPtr(guard->mir()->shape()));
    masm.ma_lw(tmp, Address(obj, JSObject::offsetOfShape()));
    return bailoutCmpPtr(Assembler::NotEqual, tmp, SecondScratchReg, guard->snapshot());
}

bool
CodeGeneratorMIPS::visitGuardObjectType(LGuardObjectType *guard)
{
    Register obj = ToRegister(guard->input());
    Register tmp = ToRegister(guard->tempInt());

    masm.ma_lw(tmp, Address(obj, JSObject::offsetOfType()));
    masm.ma_li(SecondScratchReg, ImmGCPtr(guard->mir()->typeObject()));
    Assembler::Condition cond =
        guard->mir()->bailOnEquality() ? Assembler::Equal : Assembler::NotEqual;
    return bailoutCmpPtr(cond, tmp, SecondScratchReg, guard->snapshot());
}

bool
CodeGeneratorMIPS::visitGuardClass(LGuardClass *guard)
{
    Register obj = ToRegister(guard->input());
    Register tmp = ToRegister(guard->tempInt());

    masm.loadObjClass(obj, tmp);
    if (!bailoutCmpPtr(Assembler::NotEqual, tmp, Imm32((uint32_t)guard->mir()->getClass()),
                       guard->snapshot()))
        return false;
    return true;
}

bool
CodeGeneratorMIPS::visitImplicitThis(LImplicitThis *lir)
{
    Register callee = ToRegister(lir->callee());
    const ValueOperand out = ToOutValue(lir);

    // The implicit |this| is always |undefined| if the function's environment
    // is the current global.
    masm.ma_lw(out.typeReg(), Address(callee, JSFunction::offsetOfEnvironment()));
    masm.ma_li(SecondScratchReg, ImmGCPtr(&gen->info().script()->global()));

    // TODO: OOL stub path.
    if (!bailoutCmpPtr(Assembler::NotEqual, out.typeReg(), SecondScratchReg, lir->snapshot()))
        return false;

    masm.moveValue(UndefinedValue(), out);
    return true;
}

bool
CodeGeneratorMIPS::visitInterruptCheck(LInterruptCheck *lir)
{
    OutOfLineCode *ool = oolCallVM(InterruptCheckInfo, lir, (ArgList()), StoreNothing());
    if (!ool)
        return false;

    void *interrupt = (void*)GetIonContext()->runtime->addressOfInterrupt();
    masm.load32(AbsoluteAddress(interrupt), t8);
    masm.ma_b(t8, Imm32(0), ool->entry(), Assembler::NotEqual);
    masm.bind(ool->rejoin());
    return true;
}

bool
CodeGeneratorMIPS::generateInvalidateEpilogue()
{
    // Ensure that there is enough space in the buffer for the OsiPoint
    // patching to occur. Otherwise, we could overwrite the invalidation
    // epilogue.
    for (size_t i = 0; i < sizeof(void *); i += Assembler::nopSize())
        masm.nop();

    masm.bind(&invalidate_);

    // Push the return address of the point that we bailed out at to the stack
    masm.Push(ra);

    // Push the Ion script onto the stack (when we determine what that
    // pointer is).
    invalidateEpilogueData_ = masm.pushWithPatch(ImmWord(uintptr_t(-1)));
    JitCode *thunk = gen->jitRuntime()->getInvalidationThunk();

    masm.branch(thunk);

    // We should never reach this point in JIT code -- the invalidation thunk
    // should pop the invalidated JS frame and return directly to its caller.
    masm.assumeUnreachable("Should have returned directly to its caller instead of here.");
    return true;
}

void
DispatchIonCache::initializeAddCacheState(LInstruction *ins, AddCacheState *addState)
{
    // Can always use the scratch register on MIPS.
    addState->dispatchScratch = ScratchRegister;
}

bool
CodeGeneratorMIPS::visitLoadTypedArrayElementStatic(LLoadTypedArrayElementStatic *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
CodeGeneratorMIPS::visitStoreTypedArrayElementStatic(LStoreTypedArrayElementStatic *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

bool
CodeGeneratorMIPS::visitAsmJSLoadHeap(LAsmJSLoadHeap *ins)
{
    const MAsmJSLoadHeap *mir = ins->mir();
    const LAllocation *ptr = ins->ptr();
    const LDefinition *out = ins->output();

    bool isSigned;
    int size;
    bool isFloat = false;
    switch (mir->viewType()) {
      case ArrayBufferView::TYPE_INT8:    isSigned = true;  size =  8; break;
      case ArrayBufferView::TYPE_UINT8:   isSigned = false; size =  8; break;
      case ArrayBufferView::TYPE_INT16:   isSigned = true;  size = 16; break;
      case ArrayBufferView::TYPE_UINT16:  isSigned = false; size = 16; break;
      case ArrayBufferView::TYPE_INT32:
      case ArrayBufferView::TYPE_UINT32:  isSigned = true;  size = 32; break;
      case ArrayBufferView::TYPE_FLOAT64: isFloat = true;   size = 64; break;
      case ArrayBufferView::TYPE_FLOAT32: isFloat = true;   size = 32; break;
      default: MOZ_ASSUME_UNREACHABLE("unexpected array type");
    }

    if (ptr->isConstant()) {
        MOZ_ASSERT(mir->skipBoundsCheck());
        int32_t ptrImm = ptr->toConstant()->toInt32();
        MOZ_ASSERT(ptrImm >= 0);
        if (isFloat) {
            if (size == 32) {
                masm.loadFloat32(Address(HeapReg, ptrImm), ToFloatRegister(out));
            } else {
                masm.loadDouble(Address(HeapReg, ptrImm), ToFloatRegister(out));
            }
        }  else {
            masm.ma_load(ToRegister(out), Address(HeapReg, ptrImm),
                         static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
        }
        return true;
    }

    Register ptrReg = ToRegister(ptr);

    if (mir->skipBoundsCheck()) {
        if (isFloat) {
            if (size == 32) {
                masm.loadFloat32(BaseIndex(HeapReg, ptrReg, TimesOne), ToFloatRegister(out));
            } else {
                masm.loadDouble(BaseIndex(HeapReg, ptrReg, TimesOne), ToFloatRegister(out));
            }
        } else {
            masm.ma_load(ToRegister(out), BaseIndex(HeapReg, ptrReg, TimesOne),
                         static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
        }
        return true;
    }

    BufferOffset bo = masm.ma_BoundsCheck(ScratchRegister);

    Label outOfRange;
    Label done;
    masm.ma_b(ptrReg, ScratchRegister, &outOfRange, Assembler::AboveOrEqual, ShortJump);
    // Offset is ok, let's load value.
    if (isFloat) {
        if (size == 32)
            masm.loadFloat32(BaseIndex(HeapReg, ptrReg, TimesOne), ToFloatRegister(out));
        else
            masm.loadDouble(BaseIndex(HeapReg, ptrReg, TimesOne), ToFloatRegister(out));
    } else {
        masm.ma_load(ToRegister(out), BaseIndex(HeapReg, ptrReg, TimesOne),
                     static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
    }
    masm.ma_b(&done, ShortJump);
    masm.bind(&outOfRange);
    // Offset is out of range. Load default values.
    if (isFloat) {
        if (size == 32)
            masm.convertDoubleToFloat32(NANReg, ToFloatRegister(out));
        else
            masm.moveDouble(NANReg, ToFloatRegister(out));
    } else {
        masm.ma_li(ToRegister(out), Imm32(0));
    }
    masm.bind(&done);

    return gen->noteHeapAccess(AsmJSHeapAccess(bo.getOffset()));
}

bool
CodeGeneratorMIPS::visitAsmJSStoreHeap(LAsmJSStoreHeap *ins)
{
    const MAsmJSStoreHeap *mir = ins->mir();
    const LAllocation *value = ins->value();
    const LAllocation *ptr = ins->ptr();

    bool isSigned;
    int size;
    bool isFloat = false;
    switch (mir->viewType()) {
      case ArrayBufferView::TYPE_INT8:
      case ArrayBufferView::TYPE_UINT8:   isSigned = false; size = 8; break;
      case ArrayBufferView::TYPE_INT16:
      case ArrayBufferView::TYPE_UINT16:  isSigned = false; size = 16; break;
      case ArrayBufferView::TYPE_INT32:
      case ArrayBufferView::TYPE_UINT32:  isSigned = true;  size = 32; break;
      case ArrayBufferView::TYPE_FLOAT64: isFloat  = true;  size = 64; break;
      case ArrayBufferView::TYPE_FLOAT32: isFloat = true;   size = 32; break;
      default: MOZ_ASSUME_UNREACHABLE("unexpected array type");
    }

    if (ptr->isConstant()) {
        MOZ_ASSERT(mir->skipBoundsCheck());
        int32_t ptrImm = ptr->toConstant()->toInt32();
        MOZ_ASSERT(ptrImm >= 0);

        if (isFloat) {
            if (size == 32) {
                masm.storeFloat32(ToFloatRegister(value), Address(HeapReg, ptrImm));
            } else {
                masm.storeDouble(ToFloatRegister(value), Address(HeapReg, ptrImm));
            }
        }  else {
            masm.ma_store(ToRegister(value), Address(HeapReg, ptrImm),
                          static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
        }
        return true;
    }

    Register ptrReg = ToRegister(ptr);
    Address dstAddr(ptrReg, 0);

    if (mir->skipBoundsCheck()) {
        if (isFloat) {
            if (size == 32) {
                masm.storeFloat32(ToFloatRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne));
            } else
                masm.storeDouble(ToFloatRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne));
        } else {
            masm.ma_store(ToRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne),
                          static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
        }
        return true;
    }

    BufferOffset bo = masm.ma_BoundsCheck(ScratchRegister);

    Label rejoin;
    masm.ma_b(ptrReg, ScratchRegister, &rejoin, Assembler::AboveOrEqual, ShortJump);

    // Offset is ok, let's store value.
    if (isFloat) {
        if (size == 32) {
            masm.storeFloat32(ToFloatRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne));
        } else
            masm.storeDouble(ToFloatRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne));
    } else {
        masm.ma_store(ToRegister(value), BaseIndex(HeapReg, ptrReg, TimesOne),
                      static_cast<LoadStoreSize>(size), isSigned ? SignExtend : ZeroExtend);
    }
    masm.bind(&rejoin);

    return gen->noteHeapAccess(AsmJSHeapAccess(bo.getOffset()));
}

bool
CodeGeneratorMIPS::visitAsmJSPassStackArg(LAsmJSPassStackArg *ins)
{
    const MAsmJSPassStackArg *mir = ins->mir();
    if (ins->arg()->isConstant()) {
        masm.ma_storeImm(Imm32(ToInt32(ins->arg())), Address(StackPointer, mir->spOffset()));
    } else {
        if (ins->arg()->isGeneralReg()) {
            masm.ma_sw(ToRegister(ins->arg()), Address(StackPointer, mir->spOffset()));
        } else {
            masm.ma_sd(ToFloatRegister(ins->arg()), Address(StackPointer, mir->spOffset()));
        }
    }

    return true;
}

bool
CodeGeneratorMIPS::visitUDiv(LUDiv *ins)
{
    Register lhs = ToRegister(ins->lhs());
    Register rhs = ToRegister(ins->rhs());
    Register output = ToRegister(ins->output());

    Label done;
    if (ins->mir()->canBeDivideByZero()) {
        if (ins->mir()->isTruncated()) {
            Label notzero;
            masm.ma_b(rhs, rhs, &notzero, Assembler::NonZero, ShortJump);
            masm.ma_li(output, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&notzero);
        } else {
            MOZ_ASSERT(ins->mir()->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, rhs, Imm32(0), ins->snapshot()))
                return false;
        }
    }

    masm.as_divu(lhs, rhs);
    masm.as_mflo(output);

    if (!ins->mir()->isTruncated()) {
        if (!bailoutCmpPtr(Assembler::LessThan, output, Imm32(0), ins->snapshot()))
            return false;
    }

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitUMod(LUMod *ins)
{
    Register lhs = ToRegister(ins->lhs());
    Register rhs = ToRegister(ins->rhs());
    Register output = ToRegister(ins->output());
    Label done;

    if (ins->mir()->canBeDivideByZero()) {
        if (ins->mir()->isTruncated()) {
            // Infinity|0 == 0
            Label notzero;
            masm.ma_b(rhs, rhs, &notzero, Assembler::NonZero, ShortJump);
            masm.ma_li(output, Imm32(0));
            masm.ma_b(&done, ShortJump);
            masm.bind(&notzero);
        } else {
            MOZ_ASSERT(ins->mir()->fallible());
            if (!bailoutCmpPtr(Assembler::Equal, rhs, Imm32(0), ins->snapshot()))
                return false;
        }
    }

    masm.as_divu(lhs, rhs);
    masm.as_mfhi(output);

    if (!ins->mir()->isTruncated()) {
        if (!bailoutCmpPtr(Assembler::LessThan, output, Imm32(0), ins->snapshot()))
            return false;
    }

    masm.bind(&done);
    return true;
}

bool
CodeGeneratorMIPS::visitEffectiveAddress(LEffectiveAddress *ins)
{
    const MEffectiveAddress *mir = ins->mir();
    Register base = ToRegister(ins->base());
    Register index = ToRegister(ins->index());
    Register output = ToRegister(ins->output());

    masm.ma_sll(output, index, Imm32(mir->scale()));
    masm.ma_addu(output, base);
    masm.ma_addu(output, Imm32(mir->displacement()));
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadGlobalVar(LAsmJSLoadGlobalVar *ins)
{
    const MAsmJSLoadGlobalVar *mir = ins->mir();
    unsigned addr = mir->globalDataOffset();
    if (mir->type() == MIRType_Int32)
        masm.ma_lw(ToRegister(ins->output()), Address(GlobalReg, addr));
    else if (mir->type() == MIRType_Float32)
        masm.ma_ls(ToFloatRegister(ins->output()), Address(GlobalReg, addr));
    else
        masm.ma_ld(ToFloatRegister(ins->output()), Address(GlobalReg, addr));
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSStoreGlobalVar(LAsmJSStoreGlobalVar *ins)
{
    const MAsmJSStoreGlobalVar *mir = ins->mir();

    MIRType type = mir->value()->type();
    MOZ_ASSERT(IsNumberType(type));
    unsigned addr = mir->globalDataOffset();
    if (mir->value()->type() == MIRType_Int32)
        masm.ma_sw(ToRegister(ins->value()), Address(GlobalReg, addr));
    else if (mir->value()->type() == MIRType_Float32)
        masm.ma_ss(ToFloatRegister(ins->value()), Address(GlobalReg, addr));
    else
        masm.ma_sd(ToFloatRegister(ins->value()), Address(GlobalReg, addr));
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadFuncPtr(LAsmJSLoadFuncPtr *ins)
{
    const MAsmJSLoadFuncPtr *mir = ins->mir();

    Register index = ToRegister(ins->index());
    Register tmp = ToRegister(ins->temp());
    Register out = ToRegister(ins->output());
    unsigned addr = mir->globalDataOffset();
    masm.ma_li(tmp, Imm32(addr));
    masm.ma_sll(index, index, Imm32(2));
    masm.ma_addu(tmp, index);
    BaseIndex source(GlobalReg, tmp, TimesOne);
    masm.load32(source, out);
    return true;
}

bool
CodeGeneratorMIPS::visitAsmJSLoadFFIFunc(LAsmJSLoadFFIFunc *ins)
{
    const MAsmJSLoadFFIFunc *mir = ins->mir();
    masm.ma_load(ToRegister(ins->output()), Address(GlobalReg, mir->globalDataOffset()));
    return true;
}

bool
CodeGeneratorMIPS::visitNegI(LNegI *ins)
{
    Register input = ToRegister(ins->input());
    masm.ma_negu(ToRegister(ins->output()), input);
    return true;
}

bool
CodeGeneratorMIPS::visitNegD(LNegD *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    FloatRegister output = ToFloatRegister(ins->output());

    masm.as_negd(output, input);
    return true;
}

bool
CodeGeneratorMIPS::visitNegF(LNegF *ins)
{
    FloatRegister input = ToFloatRegister(ins->input());
    FloatRegister output = ToFloatRegister(ins->output());

    masm.as_negs(output, input);
    return true;
}

bool
CodeGeneratorMIPS::visitForkJoinGetSlice(LForkJoinGetSlice *ins)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

JitCode *
JitRuntime::generateForkJoinGetSliceStub(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}
