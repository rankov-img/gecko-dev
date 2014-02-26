/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jscompartment.h"

#include "jit/Bailouts.h"
#include "jit/ExecutionModeInlines.h"
#include "jit/IonFrames.h"
#include "jit/IonLinker.h"
#include "jit/IonSpewer.h"
#include "jit/JitCompartment.h"
#include "jit/mips/BaselineHelpers-mips.h"
#ifdef JS_ION_PERF
# include "jit/PerfSpewer.h"
#endif
#include "jit/VMFunctions.h"

using namespace js;
using namespace js::jit;

// Size of buffer that holds non-volotile regs plus alignment slot.
static const uint32_t GENERAL_REGS_BUFF_SIZE = 10 * sizeof(intptr_t);

static const uint32_t FLOAT_REGS_BUFF_SIZE = 6 * sizeof(double);

static const uint32_t CALLEE_TOKEN_OFFSET = GENERAL_REGS_BUFF_SIZE + FLOAT_REGS_BUFF_SIZE +
                                            4 * sizeof(intptr_t);

static const uint32_t SCOPE_CHAIN_OFFSET = CALLEE_TOKEN_OFFSET + sizeof (CalleeToken);

static const uint32_t NUM_STACK_VALUES_OFFSET = SCOPE_CHAIN_OFFSET + sizeof(JSObject *);

static const uint32_t VP_OFFSET = NUM_STACK_VALUES_OFFSET + sizeof (size_t);

/*
 * Stack layout during execution of EnterJIT. This is valid for entering Ion
 * code and Baseline code without OSR. For OSR, there is additional data.
 *
 * --------- high addresses -----------
 * ...
 * -------- Arguments on stack --------
 *  Value *vp;
 *  size_t numStackValues;
 *  JSObject *scopeChain;
 *  CalleeToken calleeToken;
 * -- First 4 arguments placeholders --
 *  StackFrame *fp;
 *  Value *maxArgv;
 *  int maxArgc;
 *  void * code;     <- sp points here when function is entered.
 * ------ non-volatile registers ------
 *  void *s0;
 *  void *s1;
 *  void *s2;
 *  void *s3;
 *  void *s4;
 *  void *s5;
 *  void *s6;
 *  void *s7;
 *  void *ra;
 *  void *alignment_slot.
 * --- non-volatile float registers ---
 *  double f20;
 *  double f22;
 *  double f24;
 *  double f26;
 *  double f28;
 *  double f30;     <- sp points here while saving and restoring registers.
 * ------ JS function arguments -------
 * ... variable number of arguments
 * ------- Ion (Baseline) Frame -------
 *  uint32_t actual arguments;
 *  CalleeToken calleeToken;
 *  uint32_t frameDescriptor
 *  void * return_address
 * ------------------------------------
 * ...
 * ---------- low addresses -----------
 */

static void
GenerateReturn(MacroAssembler &masm, int returnCode)
{
    // Restore non-volatile registers
    uint32_t offset = GENERAL_REGS_BUFF_SIZE + FLOAT_REGS_BUFF_SIZE;
    masm.ma_lw(s0, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s1, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s2, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s3, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s4, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s5, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s6, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(s7, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_lw(ra, StackPointer, offset -= sizeof(intptr_t));

    // Align the stack for double values.
    offset -= sizeof(intptr_t);

    // Restore non-volatile floating point registers
    masm.ma_ld(f20, StackPointer, offset -= sizeof(double));
    masm.ma_ld(f22, StackPointer, offset -= sizeof(double));
    masm.ma_ld(f24, StackPointer, offset -= sizeof(double));
    masm.ma_ld(f26, StackPointer, offset -= sizeof(double));
    masm.ma_ld(f28, StackPointer, offset -= sizeof(double));
    masm.ma_ld(f30, StackPointer, offset -= sizeof(double));

    JS_ASSERT(offset == 0);
    masm.ma_addu(StackPointer, StackPointer,
                 Imm32(GENERAL_REGS_BUFF_SIZE + FLOAT_REGS_BUFF_SIZE));

    masm.as_jr(ra);
    masm.as_nop();
}

/*
 * This method generates a trampoline for a c++ function with the following
 * signature:
 *   void enter(void *code, int argc, Value *argv, StackFrame *fp, CalleeToken
 *              calleeToken, JSObject *scopeChain, Value *vp)
 *   ...using standard EABI calling convention
 */
JitCode *
JitRuntime::generateEnterJIT(JSContext *cx, EnterJitType type)
{
    const Register reg_code = a0;
    const Register reg_argc = a1;
    const Register reg_argv = a2;
    const Register reg_frame = a3;

    const Address slot_token(sp, CALLEE_TOKEN_OFFSET);
    const Address slot_vp(sp, VP_OFFSET);

    JS_ASSERT(OsrFrameReg == reg_frame);

    MacroAssembler masm(cx);
    AutoFlushCache afc("GenerateEnterJIT", cx->runtime()->jitRuntime());

    // Save non-volatile registers. These must be saved by the trampoline,
    // rather than the JIT'd code, because they are scanned by the conservative
    // scanner.
    uint32_t offset = GENERAL_REGS_BUFF_SIZE + FLOAT_REGS_BUFF_SIZE;
    masm.ma_subu(StackPointer, StackPointer, Imm32(offset));
    masm.ma_sw(s0, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s1, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s2, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s3, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s4, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s5, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s6, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(s7, StackPointer, offset -= sizeof(intptr_t));
    masm.ma_sw(ra, StackPointer, offset -= sizeof(intptr_t));

    // Align the stack for double values.
    offset -= sizeof(intptr_t);
    masm.ma_sd(f20, StackPointer, offset -= sizeof(double));
    masm.ma_sd(f22, StackPointer, offset -= sizeof(double));
    masm.ma_sd(f24, StackPointer, offset -= sizeof(double));
    masm.ma_sd(f26, StackPointer, offset -= sizeof(double));
    masm.ma_sd(f28, StackPointer, offset -= sizeof(double));
    masm.ma_sd(f30, StackPointer, offset -= sizeof(double));

    JS_ASSERT(offset == 0);

    // Save stack pointer into s4
    masm.movePtr(sp, s4);

    // Load calleeToken into s2.
    masm.loadPtr(slot_token, s2);

    // Save stack pointer as baseline frame.
    if (type == EnterJitBaseline)
        masm.movePtr(sp, BaselineFrameReg);

    // Load the number of actual arguments into s3.
    masm.loadPtr(slot_vp, s3);
    masm.unboxInt32(Address(s3, 0), s3);

    /***************************************************************
    Loop over argv vector, push arguments onto stack in reverse order
    ***************************************************************/

    masm.as_sll(s0, reg_argc, 3); // s0 = argc * 8
    masm.as_addu(s0, s0, reg_argv); // s0 = argv + argc * 8

    // Loop over arguments, copying them from an unknown buffer onto the Ion
    // stack so they can be accessed from JIT'ed code.
    Label header, footer;
    // If there aren't any arguments, don't do anything
    masm.ma_b(s0, reg_argv, &footer, Assembler::BelowOrEqual, ShortJump);
    {
        masm.bind(&header);

        masm.ma_subu(s0, s0, Imm32(2 * sizeof(intptr_t)));
        masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));

        //TODO: If we make stack to be always aligned, we can use l.d and s.d.
        masm.as_lw(s1, s0, sizeof(intptr_t));
        masm.as_sw(s1, StackPointer, sizeof(intptr_t));
        masm.as_lw(s1, s0, 0);
        masm.as_sw(s1, StackPointer, 0);

        masm.ma_b(s0, reg_argv, &header, Assembler::Above, ShortJump);
    }
    masm.bind(&footer);

    masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
    masm.ma_sw(s3, StackPointer, sizeof(intptr_t)); // actual arguments
    masm.ma_sw(s2, StackPointer, 0); // callee token

    masm.ma_subu(s4, s4, sp);
    masm.makeFrameDescriptor(s4, IonFrame_Entry);
    masm.push(s4); // descriptor

    CodeLabel returnLabel;
    if (type == EnterJitBaseline) {
        // Handle OSR.
        GeneralRegisterSet regs(GeneralRegisterSet::All());
        regs.take(JSReturnOperand);
        regs.take(OsrFrameReg);
        regs.take(BaselineFrameReg);
        regs.take(reg_code);

        const Address slot_numStackValues(BaselineFrameReg, NUM_STACK_VALUES_OFFSET);

        Label notOsr;
        masm.ma_b(OsrFrameReg, OsrFrameReg, &notOsr, Assembler::Zero, ShortJump);

        Register scratch = regs.takeAny();

        Register numStackValues = regs.takeAny();
        masm.load32(slot_numStackValues, numStackValues);

        // Push return address, previous frame pointer.
        masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
        masm.ma_li(scratch, returnLabel.dest());
        masm.ma_sw(scratch, StackPointer, sizeof(intptr_t));
        masm.ma_sw(BaselineFrameReg, StackPointer, 0);

        // Reserve frame.
        Register framePtr = BaselineFrameReg;
        masm.subPtr(Imm32(BaselineFrame::Size()), StackPointer);
        masm.mov(StackPointer, framePtr);

        // Reserve space for locals and stack values.
        masm.ma_sll(scratch, numStackValues, Imm32(3));
        masm.ma_subu(sp, sp, scratch);

        // Enter exit frame.
        masm.ma_addu(scratch, scratch, Imm32(BaselineFrame::Size() +
                     BaselineFrame::FramePointerOffset));
        masm.makeFrameDescriptor(scratch, IonFrame_BaselineJS);

        // Push frame descriptor and fake return address.
        masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
        masm.ma_sw(scratch, StackPointer, sizeof(intptr_t)); // Frame descriptor
        masm.ma_sw(zero, StackPointer, 0); // fake return address

        masm.enterFakeExitFrame();

        masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
        masm.ma_sw(framePtr, StackPointer, sizeof(intptr_t)); // BaselineFrame
        masm.ma_sw(reg_code, StackPointer, 0); // jitcode

        masm.setupUnalignedABICall(3, scratch);
        masm.passABIArg(BaselineFrameReg); // BaselineFrame
        masm.passABIArg(OsrFrameReg); // StackFrame
        masm.passABIArg(numStackValues);
        masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, jit::InitBaselineFrameForOsr));

        Register jitcode = regs.takeAny();
        masm.ma_lw(jitcode, StackPointer, 0);
        masm.ma_lw(framePtr, StackPointer, sizeof(intptr_t));
        masm.ma_addu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));

        JS_ASSERT(jitcode != ReturnReg);

        Label error;
        masm.addPtr(Imm32(IonExitFrameLayout::SizeWithFooter()), sp);
        masm.addPtr(Imm32(BaselineFrame::Size()), framePtr);
        masm.branchIfFalseBool(ReturnReg, &error);

        masm.jump(jitcode);

        // OOM: load error value, discard return address and previous frame
        // pointer and return.
        masm.bind(&error);
        masm.mov(framePtr, sp);
        masm.addPtr(Imm32(2 * sizeof(uintptr_t)), sp);
        masm.moveValue(MagicValue(JS_ION_ERROR), JSReturnOperand);
        masm.ma_li(scratch, returnLabel.dest());
        masm.jump(scratch);

        masm.bind(&notOsr);
        // Load the scope chain in R1.
        JS_ASSERT(R1.scratchReg() != reg_code);
        masm.loadPtr(Address(BaselineFrameReg, SCOPE_CHAIN_OFFSET), R1.scratchReg());
    }

    // Call the function with pushing return address to stack.
    masm.ma_callIonHalfPush(reg_code);

    if (type == EnterJitBaseline) {
        // Baseline OSR will return here.
        masm.bind(returnLabel.src());
        if (!masm.addCodeLabel(returnLabel))
            return nullptr;
    }

    // Pop arguments off the stack.
    // s0 <- 8*argc (size of all arguments we pushed on the stack)
    masm.pop(s0);
    masm.ma_srl(s0, s0, Imm32(4));
    masm.as_addu(sp, sp, s0);

    // Store the returned value into the slot_vp
    masm.loadPtr(slot_vp, s1);
    masm.storeValue(JSReturnOperand, Address(s1, 0));

    // Restore non-volatile registers and return.
    GenerateReturn(masm, ShortJump);

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "EnterJIT");
#endif

    return code;
}

JitCode *
JitRuntime::generateInvalidator(JSContext *cx)
{
    MacroAssembler masm(cx);

    static const uint32_t STACK_DATA_SIZE = Registers::Total * sizeof(intptr_t) +
                                            FloatRegisters::Total * sizeof(double) +
                                            2 * sizeof(intptr_t);
    static const uint32_t RETVAL_OFFSET = sizeof(intptr_t);
    static const uint32_t INV_DATA_OFFSET = RETVAL_OFFSET + sizeof(intptr_t);
    static const uint32_t GEN_REGS_OFFSET = INV_DATA_OFFSET +
                                            FloatRegisters::Total * sizeof(double);

    // Stack has to be alligned here. If not, we will have to fix it.
    masm.checkStackAlignment();

    // Make room for data on stack.
    masm.ma_subu(StackPointer, StackPointer, Imm32(STACK_DATA_SIZE));

    // Save general purpose registers
    for (uint32_t i = 0; i < Registers::Total; i++)
        masm.as_sw(Register::FromCode(i), StackPointer, GEN_REGS_OFFSET + i * sizeof(intptr_t));

    // Save floating point registers
    // We can use as_sd because stack is alligned.
    for (uint32_t i = 0; i < FloatRegisters::Total; i++)
        masm.as_sd(FloatRegister::FromCode(i), StackPointer,
                   INV_DATA_OFFSET + i * sizeof(double));

    // Pass pointer to InvalidationBailoutStack structure.
    masm.ma_addu(a0, StackPointer, Imm32(INV_DATA_OFFSET));
    // Pass pointer to return value.
    masm.ma_addu(a1, StackPointer, Imm32(RETVAL_OFFSET));
    // Pass pointer to BailoutInfo
    masm.ma_move(a2, StackPointer);

    masm.setupAlignedABICall(3);
    masm.passABIArg(a0);
    masm.passABIArg(a1);
    masm.passABIArg(a2);
    masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, InvalidationBailout));

    masm.ma_lw(a2, StackPointer, 0);
    masm.ma_lw(a1, StackPointer, RETVAL_OFFSET);
    // Remove the return address, the IonScript, the register state
    // (InvaliationBailoutStack) and the space that was allocated for the
    // return value.
    masm.ma_addu(StackPointer, Imm32(sizeof(InvalidationBailoutStack) + 2 * sizeof(intptr_t)));
    // remove the space that this frame was using before the bailout
    // (computed by InvalidationBailout)
    masm.ma_addu(StackPointer, a1);

    // Jump to shared bailout tail. The BailoutInfo pointer has to be in r2.
    JitCode *bailoutTail = cx->runtime()->jitRuntime()->getBailoutTail();
    masm.branch(bailoutTail);

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);
    IonSpew(IonSpew_Invalidate, "   invalidation thunk created at %p", (void *) code->raw());

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "Invalidator");
#endif

    return code;
}

JitCode *
JitRuntime::generateArgumentsRectifier(JSContext *cx, ExecutionMode mode, void **returnAddrOut)
{
    MacroAssembler masm(cx);

    // ArgumentsRectifierReg contains the |nargs| pushed onto the current
    // frame. Including |this|, there are (|nargs| + 1) arguments to copy.
    JS_ASSERT(ArgumentsRectifierReg == s3);

    Register numActArgsReg = t6;
    Register calleeTokenReg = t7;
    Register numArgsReg = t5;

    // Copy number of actual arguments into numActArgsReg
    masm.ma_lw(numActArgsReg, sp, IonRectifierFrameLayout::offsetOfNumActualArgs());

    // Load the number of |undefined|s to push into t1.
    masm.ma_lw(calleeTokenReg, sp, IonRectifierFrameLayout::offsetOfCalleeToken());
    masm.load16ZeroExtend(Address(calleeTokenReg, JSFunction::offsetOfNargs()), numArgsReg);

    masm.as_subu(t1, numArgsReg, s3);

    masm.moveValue(UndefinedValue(), t3, t4);

    masm.ma_move(t2, sp); // Save %sp.

    // Push undefined.
    {
        Label undefLoopTop;
        masm.bind(&undefLoopTop);

        masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
        masm.ma_sw(t3, StackPointer, sizeof(intptr_t)); // type(undefined);
        masm.ma_sw(t4, StackPointer, 0); // payload(undefined);
        masm.ma_subu(t1, t1, Imm32(1));

        masm.ma_b(t1, t1, &undefLoopTop, Assembler::NonZero, ShortJump);
    }

    // Get the topmost argument.
    masm.ma_sll(t0, s3, Imm32(3)); // t0 <- nargs * 8
    masm.as_addu(t2, t2, t0); // t2 <- t2(saved sp) + nargs * 8
    masm.ma_addu(t2, t2, Imm32(sizeof(IonRectifierFrameLayout)));

    // Push arguments, |nargs| + 1 times (to include |this|).
    {
        Label copyLoopTop, initialSkip;

        masm.ma_b(&initialSkip, ShortJump);

        masm.bind(&copyLoopTop);
        masm.ma_subu(t2, t2, Imm32(sizeof(Value)));
        masm.ma_subu(s3, s3, Imm32(1));

        masm.bind(&initialSkip);

        // Read argument and push to stack.
        masm.ma_subu(StackPointer, StackPointer, Imm32(sizeof(Value)));
        masm.ma_lw(t0, t2, sizeof(Value) / 2);
        masm.ma_sw(t0, StackPointer, sizeof(Value) / 2);
        masm.ma_lw(t0, t2, 0);
        masm.ma_sw(t0, StackPointer, 0);

        masm.ma_b(s3, s3, &copyLoopTop, Assembler::NonZero, ShortJump);
    }

    // translate the framesize from values into bytes
    masm.ma_addu(t0, numArgsReg, Imm32(1));
    masm.ma_sll(t0, t0, Imm32(3));

    // Construct sizeDescriptor.
    masm.makeFrameDescriptor(t0, IonFrame_Rectifier);

    // Construct IonJSFrameLayout.
    masm.ma_subu(StackPointer, StackPointer, Imm32(3 * sizeof(intptr_t)));
    // Push actual arguments.
    masm.ma_sw(numActArgsReg, StackPointer, 2 * sizeof(intptr_t));
    // Push callee token.
    masm.ma_sw(calleeTokenReg, StackPointer, sizeof(intptr_t));
    // Push frame descriptor.
    masm.ma_sw(t0, StackPointer, 0);

    // Call the target function.
    // Note that this code assumes the function is JITted.
    masm.ma_lw(t1, calleeTokenReg, JSFunction::offsetOfNativeOrScript());
    masm.loadBaselineOrIonRaw(t1, t1, mode, nullptr);
    masm.ma_callIonHalfPush(t1);

    uint32_t returnOffset = masm.currentOffset();

    // arg1
    //  ...
    // argN
    // num actual args
    // callee token
    // sizeDescriptor     <- sp now
    // return address

    // Remove the rectifier frame.
    // t0 <- descriptor with FrameType.
    masm.ma_lw(t0, StackPointer, 0);
    masm.ma_srl(t0, t0, Imm32(FRAMESIZE_SHIFT)); // t0 <- descriptor.

    // Discard descriptor, calleeToken and number of actual arguments.
    masm.ma_addu(StackPointer, StackPointer, Imm32(3 * sizeof(intptr_t)));

    // arg1
    //  ...
    // argN               <- sp now; t0 <- frame descriptor
    // num actual args
    // callee token
    // sizeDescriptor
    // return address

    // Discard pushed arguments.
    masm.as_addu(StackPointer, StackPointer, t0);

    masm.ret();
    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

    CodeOffsetLabel returnLabel(returnOffset);
    returnLabel.fixup(&masm);
    if (returnAddrOut)
        *returnAddrOut = (void *) (code->raw() + returnLabel.offset());

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "ArgumentsRectifier");
#endif

    return code;
}

/* There are two different stack layouts when doing bailout.
 *
 * - First case is when bailout is done trough bailout table. In this case
 * table offset is stored in $ra (look at JitRuntime::generateBailoutTable())
 * and thunk code should save it on stack. In this case frameClassId_ cannot
 * be NO_FRAME_SIZE_CLASS_ID.
 * Stack layout looks like this:
 *
 * --------- high addresses -----------
 * ...
 * ------- general purpose regs -------
 *  uintptr_t regs_[Registers::Total];
 * ------- floating-point regs --------
 *  double fpregs_[FloatRegisters::Total];
 * ------------------------------------
 * uintptr_t tableOffset_;      <- store $ra here.
 * uintptr_t frameClassId_;     <- address of BailoutStack structure.
 * -------- BailoutInfo pointer -------
 * uintptr_t align_slot;
 * BailoutInfo *bailoutInfo;    <- Bailout function updates this pointer.
 * ------------------------------------
 * ...
 * ---------- low addresses -----------
 *
 * - Other case is when bailout is done via out of line code (lazy bailout).
 * In this case frame size is stored in $ra (look at
 * CodeGeneratorMIPS::generateOutOfLineCode()) and thunk code should save it
 * on stack. Other difference is that Snapshot offset is pushed to the stack
 * by CodeGeneratorMIPS::visitOutOfLineBailout(). This method also adds one
 * empty slot so that data is properly aligned. Field frameClassId_ is forced
 * to be NO_FRAME_SIZE_CLASS_ID (See: JitRuntime::generateBailoutHandler).
 * Stack layout looks like this:
 *
 * --------- high addresses -----------
 * ...
 * --------- snapshot offset ----------
 * uintptr_t align_slot;
 * uintptr_t snapshotOffset_;
 * ------- general purpose regs -------
 *  uintptr_t regs_[Registers::Total];
 * ------- floating-point regs --------
 *  double fpregs_[FloatRegisters::Total];
 * ------------------------------------
 * uintptr_t frameSize_;
 * uintptr_t frameClassId_;     <- address of BailoutStack structure.
 * -------- BailoutInfo pointer -------
 * uintptr_t align_slot;
 * BailoutInfo *bailoutInfo;    <- pointera passed to Bailout function.
 * ------------------------------------
 * ...
 * ---------- low addresses -----------
 */
static void
GenerateBailoutThunk(JSContext *cx, MacroAssembler &masm, uint32_t frameClass)
{
    static const uint32_t STACK_DATA_SIZE = Registers::Total * sizeof(intptr_t) +
                                            FloatRegisters::Total * sizeof(double) +
                                            4 * sizeof(intptr_t);
    static const uint32_t FRAME_CLASS_OFFSET = 2 * sizeof(intptr_t);
    static const uint32_t FRAME_SIZE_OFFSET = FRAME_CLASS_OFFSET + sizeof(intptr_t);
    static const uint32_t FLOAT_REGS_OFFSET = FRAME_SIZE_OFFSET + sizeof(intptr_t);
    static const uint32_t GEN_REGS_OFFSET = FLOAT_REGS_OFFSET +
                                            FloatRegisters::Total * sizeof(double);

    // Make sure that alignment is proper.
    masm.checkStackAlignment();

    // Make room for data.
    masm.ma_subu(StackPointer, StackPointer, Imm32(STACK_DATA_SIZE));

    // Save general purpose registers.
    for (uint32_t i = 0; i < Registers::Total; i++)
        masm.as_sw(Register::FromCode(i), StackPointer, GEN_REGS_OFFSET + i * sizeof(intptr_t));

    // Save floating point registers
    // We can use as_sd because stack is alligned.
    for (uint32_t i = 0; i < FloatRegisters::Total; i++)
        masm.as_sd(FloatRegister::FromCode(i), StackPointer,
                   FLOAT_REGS_OFFSET + i * sizeof(double));

    // Store the frameSize_ or tableOffset_ stored in ra
    // See: JitRuntime::generateBailoutTable()
    // See: CodeGeneratorMIPS::generateOutOfLineCode()
    masm.ma_sw(ra, StackPointer, FRAME_SIZE_OFFSET);

    // Put frame class to stack
    masm.ma_sw(Imm32(frameClass), StackPointer, FRAME_CLASS_OFFSET);

    // Put pointer to BailoutStack as first argument to the Bailout()
    masm.ma_addu(a0, StackPointer, Imm32(FRAME_CLASS_OFFSET));
    // Put pointer to BailoutInfo
    masm.ma_move(a1, StackPointer);

    masm.setupAlignedABICall(2);
    masm.passABIArg(a0);
    masm.passABIArg(a1);
    masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, Bailout));

    // Get BailoutInfo pointer
    masm.ma_lw(a2, StackPointer, 0);

    // Remove both the bailout frame and the topmost Ion frame's stack.
    if (frameClass == NO_FRAME_SIZE_CLASS_ID) {
        // Load frameSize from stack
        masm.ma_lw(a1, StackPointer, FRAME_SIZE_OFFSET);

        // Remove bailout data from stack and add the size of
        // snapshotOffset and alignment slot
        masm.ma_addu(StackPointer, StackPointer, Imm32(STACK_DATA_SIZE + 2 * sizeof(void *)));
        // Remove frame size srom stack
        masm.as_addu(StackPointer, StackPointer, a1);
    } else {
        uint32_t frameSize = FrameSizeClass::FromClass(frameClass).frameSize();
        masm.ma_addu(StackPointer, StackPointer, Imm32(STACK_DATA_SIZE + frameSize));
    }

    // Jump to shared bailout tail. The BailoutInfo pointer has to be in a2.
    JitCode *bailoutTail = cx->runtime()->jitRuntime()->getBailoutTail();
    masm.branch(bailoutTail);
}

JitCode *
JitRuntime::generateBailoutTable(JSContext *cx, uint32_t frameClass)
{
    MacroAssembler masm(cx);

    Label bailout;
    for (size_t i = 0; i < BAILOUT_TABLE_SIZE; i++) {
        // Calculate offset to the end of table
        int32_t offset = (BAILOUT_TABLE_SIZE - i) * BAILOUT_TABLE_ENTRY_SIZE;

        // We use the 'ra' as table offset later in GenerateBailoutThunk
        masm.as_bal(BOffImm16(offset));
        masm.nop();
    }
    masm.bind(&bailout);

    GenerateBailoutThunk(cx, masm, frameClass);

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "BailoutTable");
#endif

    return code;
}

JitCode *
JitRuntime::generateBailoutHandler(JSContext *cx)
{
    MacroAssembler masm(cx);
    GenerateBailoutThunk(cx, masm, NO_FRAME_SIZE_CLASS_ID);

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "BailoutHandler");
#endif

    return code;
}

JitCode *
JitRuntime::generateVMWrapper(JSContext *cx, const VMFunction &f)
{
    JS_ASSERT(functionWrappers_);
    JS_ASSERT(functionWrappers_->initialized());
    VMWrapperMap::AddPtr p = functionWrappers_->lookupForAdd(&f);
    if (p)
        return p->value();

    MacroAssembler masm(cx);

    GeneralRegisterSet regs = GeneralRegisterSet(Register::Codes::WrapperMask);

    // Wrapper register set is a superset of Volatile register set.
    JS_STATIC_ASSERT((Register::Codes::VolatileMask & ~Register::Codes::WrapperMask) == 0);

    // The context is the first argument; a0 is the first argument register.
    Register cxreg = a0;
    regs.take(cxreg);

    // We're aligned to an exit frame, so link it up.
    masm.enterExitFrameAndLoadContext(&f, cxreg, regs.getAny(), f.executionMode);

    // Save the base of the argument set stored on the stack.
    Register argsBase = InvalidReg;
    if (f.explicitArgs) {
        argsBase = t1; // Use temporary register.
        regs.take(argsBase);
        masm.ma_addu(argsBase, sp, Imm32(IonExitFrameLayout::SizeWithFooter()));
    }

    // Reserve space for the outparameter.
    Register outReg = InvalidReg;
    switch (f.outParam) {
      case Type_Value:
        outReg = t0; // // Use temporary register.
        regs.take(outReg);
        masm.reserveStack(sizeof(Value));
        masm.ma_move(outReg, sp);
        break;

      case Type_Handle:
        outReg = t0;
        regs.take(outReg);
        masm.PushEmptyRooted(f.outParamRootType);
        masm.ma_move(outReg, sp);
        break;

      case Type_Int32:
      case Type_Pointer:
      case Type_Bool:
        outReg = t0;
        regs.take(outReg);
        masm.reserveStack(sizeof(int32_t));
        masm.ma_move(outReg, sp);
        break;

      case Type_Double:
        outReg = t0;
        regs.take(outReg);
        // Double outparam has to be 8 byte aligned because the called
        // function can use sdc1 or ldc1 instructions to access it.
        masm.ma_and(outReg, sp, Imm32(~(StackAlignment - 1)));
        masm.ma_subu(outReg, outReg, Imm32(sizeof(double)));
        masm.reserveStack(3 * sizeof(void *));
        break;

      default:
        JS_ASSERT(f.outParam == Type_Void);
        break;
    }

    masm.setupUnalignedABICall(f.argc(), regs.getAny());
    masm.passABIArg(cxreg);

    size_t argDisp = 0;

    // Copy any arguments.
    for (uint32_t explicitArg = 0; explicitArg < f.explicitArgs; explicitArg++) {
        MoveOperand from;
        switch (f.argProperties(explicitArg)) {
          case VMFunction::WordByValue:
            masm.passABIArg(MoveOperand(argsBase, argDisp), MoveOp::GENERAL);
            argDisp += sizeof(void *);
            break;
          case VMFunction::DoubleByValue:
            // Values should be passed by reference, not by value, so we
            // assert that the argument is a double-precision float.
            JS_ASSERT(f.argPassedInFloatReg(explicitArg));
            masm.passABIArg(MoveOperand(argsBase, argDisp), MoveOp::DOUBLE);
            argDisp += sizeof(double);
            break;
          case VMFunction::WordByRef:
            masm.passABIArg(MoveOperand(argsBase, argDisp, MoveOperand::EFFECTIVE_ADDRESS), MoveOp::GENERAL);
            argDisp += sizeof(void *);
            break;
          case VMFunction::DoubleByRef:
            masm.passABIArg(MoveOperand(argsBase, argDisp, MoveOperand::EFFECTIVE_ADDRESS), MoveOp::GENERAL);
            argDisp += 2 * sizeof(void *);
            break;
        }
    }

    // Copy the implicit outparam, if any.
    if (outReg != InvalidReg)
        masm.passABIArg(outReg);

    masm.callWithABI(f.wrapped);

    // Test for failure.
    switch (f.failType()) {
      case Type_Object:
        masm.branchTestPtr(Assembler::Zero, v0, v0, masm.failureLabel(f.executionMode));
        break;
      case Type_Bool:
        // Called functions return bools, which are 0/false and non-zero/true
        masm.branchIfFalseBool(v0, masm.failureLabel(f.executionMode));
        break;
      default:
        MOZ_ASSUME_UNREACHABLE("unknown failure kind");
    }

    // Load the outparam and free any allocated stack.
    switch (f.outParam) {
      case Type_Handle:
        masm.popRooted(f.outParamRootType, ReturnReg, JSReturnOperand);
        break;

      case Type_Value:
        masm.loadValue(Address(sp, 0), JSReturnOperand);
        masm.freeStack(sizeof(Value));
        break;

      case Type_Int32:
      case Type_Pointer:
        masm.load32(Address(sp, 0), ReturnReg);
        masm.freeStack(sizeof(int32_t));
        break;

      case Type_Bool:
        masm.load8ZeroExtend(Address(sp, 0), ReturnReg);
        masm.freeStack(sizeof(int32_t));
        break;

      case Type_Double:
        masm.freeStack(3 * sizeof(void *));
        if (cx->runtime()->jitSupportsFloatingPoint) {
            masm.ma_and(masm.secondScratch(), sp, Imm32(~(StackAlignment - 1)));
            masm.ma_subu(masm.secondScratch(), masm.secondScratch(), Imm32(sizeof(double)));
            masm.ma_ld(ReturnFloatReg, masm.secondScratch(), 0);
        } else {
            masm.assumeUnreachable("Unable to load into float reg, with no FP support.");
        }
        break;


      default:
        JS_ASSERT(f.outParam == Type_Void);
        break;
    }
    masm.leaveExitFrame();
    masm.retn(Imm32(sizeof(IonExitFrameLayout) +
                    f.explicitStackSlots() * sizeof(void *) +
                    f.extraValuesToPop * sizeof(Value)));

    Linker linker(masm);
    JitCode *wrapper = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);
    if (!wrapper)
        return nullptr;

    // linker.newCode may trigger a GC and sweep functionWrappers_ so we have
    // to use relookupOrAdd instead of add.
    if (!functionWrappers_->relookupOrAdd(p, &f, wrapper))
        return nullptr;

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(wrapper, "VMWrapper");
#endif

    return wrapper;
}

JitCode *
JitRuntime::generatePreBarrier(JSContext *cx, MIRType type)
{
    MacroAssembler masm(cx);

    RegisterSet save;
    if (cx->runtime()->jitSupportsFloatingPoint) {
        save = RegisterSet(GeneralRegisterSet(Registers::VolatileMask),
                           FloatRegisterSet(FloatRegisters::VolatileMask));
    } else {
        save = RegisterSet(GeneralRegisterSet(Registers::VolatileMask),
                           FloatRegisterSet());
    }
    masm.PushRegsInMask(save);

    JS_ASSERT(PreBarrierReg == a1);
    masm.movePtr(ImmPtr(cx->runtime()), a0);

    masm.setupUnalignedABICall(2, a2);
    masm.passABIArg(a0);
    masm.passABIArg(a1);

    if (type == MIRType_Value) {
        masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, MarkValueFromIon));
    } else {
        JS_ASSERT(type == MIRType_Shape);
        masm.callWithABI(JS_FUNC_TO_DATA_PTR(void *, MarkShapeFromIon));
    }

    masm.PopRegsInMask(save);
    masm.ret();

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "PreBarrier");
#endif

    return code;
}

typedef bool (*HandleDebugTrapFn)(JSContext *, BaselineFrame *, uint8_t *, bool *);
static const VMFunction HandleDebugTrapInfo = FunctionInfo<HandleDebugTrapFn>(HandleDebugTrap);

JitCode *
JitRuntime::generateDebugTrapHandler(JSContext *cx)
{
    MacroAssembler masm(cx);

    Register scratch1 = t0;
    Register scratch2 = t1;

    // Load BaselineFrame pointer in scratch1.
    masm.mov(s5, scratch1);
    masm.subPtr(Imm32(BaselineFrame::Size()), scratch1);

    // Enter a stub frame and call the HandleDebugTrap VM function. Ensure
    // the stub frame has a nullptr ICStub pointer, since this pointer is
    // marked during GC.
    masm.movePtr(ImmPtr(nullptr), BaselineStubReg);
    EmitEnterStubFrame(masm, scratch2);

    JitCode *code = cx->runtime()->jitRuntime()->getVMWrapper(HandleDebugTrapInfo);
    if (!code)
        return nullptr;

    masm.ma_subu(StackPointer, StackPointer, Imm32(2 * sizeof(intptr_t)));
    masm.ma_sw(ra, StackPointer, sizeof(intptr_t));
    masm.ma_sw(scratch1, StackPointer, 0);

    EmitCallVM(code, masm);

    EmitLeaveStubFrame(masm);

    // If the stub returns |true|, we have to perform a forced return
    // (return from the JS frame). If the stub returns |false|, just return
    // from the trap stub so that execution continues at the current pc.
    Label forcedReturn;
    masm.branchTest32(Assembler::NonZero, ReturnReg, ReturnReg, &forcedReturn);

    // ra was restored by EmitLeaveStubFrame
    masm.as_jr(ra);
    masm.nop();

    masm.bind(&forcedReturn);
    masm.loadValue(Address(s5, BaselineFrame::reverseOffsetOfReturnValue()),
                   JSReturnOperand);
    masm.ma_move(sp, s5);
    masm.pop(s5);
    masm.ret();

    Linker linker(masm);
    JitCode *codeDbg = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(codeDbg, "DebugTrapHandler");
#endif

    return codeDbg;
}


JitCode *
JitRuntime::generateExceptionTailStub(JSContext *cx)
{
    MacroAssembler masm;

    masm.handleFailureWithHandlerTail();

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "ExceptionTailStub");
#endif

    return code;
}

JitCode *
JitRuntime::generateBailoutTailStub(JSContext *cx)
{
    MacroAssembler masm;

    masm.generateBailoutTail(a1, a2);

    Linker linker(masm);
    JitCode *code = linker.newCode<NoGC>(cx, JSC::OTHER_CODE);

#ifdef JS_ION_PERF
    writePerfSpewerJitCodeProfile(code, "BailoutTailStub");
#endif

    return code;
}

