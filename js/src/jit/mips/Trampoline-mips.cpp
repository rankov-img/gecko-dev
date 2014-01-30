/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jscompartment.h"

#include "assembler/assembler/MacroAssembler.h"
#include "jit/mips/BaselineHelpers-mips.h"
#include "jit/Bailouts.h"
#include "jit/ExecutionModeInlines.h"

#include "jit/IonFrames.h"
#include "jit/IonLinker.h"
#include "jit/IonSpewer.h"
#include "jit/JitCompartment.h"
#ifdef JS_ION_PERF
# include "jit/PerfSpewer.h"
#endif
#include "jit/VMFunctions.h"

using namespace js;
using namespace js::jit;


static void
GenerateReturn(MacroAssembler &masm, int returnCode)
{
    // Restore non-volatile registers

    MOZ_ASSUME_UNREACHABLE("NYI");
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
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generateInvalidator(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generateArgumentsRectifier(JSContext *cx, ExecutionMode mode, void **returnAddrOut)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}


static void
GenerateBailoutThunk(JSContext *cx, MacroAssembler &masm, uint32_t frameClass)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

JitCode *
JitRuntime::generateBailoutTable(JSContext *cx, uint32_t frameClass)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generateBailoutHandler(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generateVMWrapper(JSContext *cx, const VMFunction &f)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generatePreBarrier(JSContext *cx, MIRType type)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

typedef bool (*HandleDebugTrapFn)(JSContext *, BaselineFrame *, uint8_t *, bool *);
static const VMFunction HandleDebugTrapInfo = FunctionInfo<HandleDebugTrapFn>(HandleDebugTrap);

JitCode *
JitRuntime::generateDebugTrapHandler(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}


JitCode *
JitRuntime::generateExceptionTailStub(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

JitCode *
JitRuntime::generateBailoutTailStub(JSContext *cx)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return nullptr;
}

