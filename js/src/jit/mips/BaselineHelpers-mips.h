/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef jit_mips_BaselineHelpers_mips_h
#define jit_mips_BaselineHelpers_mips_h

#ifdef JS_ION
#include "jit/BaselineFrame.h"
#include "jit/BaselineIC.h"
#include "jit/BaselineRegisters.h"
#include "jit/IonMacroAssembler.h"

namespace js {
namespace jit {

// Distance from sp to the top Value inside an IC stub (no return address on the stack on MIPS).
static const size_t ICStackValueOffset = 0;

inline void
EmitRestoreTailCallReg(MacroAssembler &masm)
{
    // No-op on MIPS because ra register is always holding the return address.
}

inline void
EmitRepushTailCallReg(MacroAssembler &masm)
{
    // No-op on MIPS because link register is always holding the return address.
}

inline void
EmitCallIC(CodeOffsetLabel *patchOffset, MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitEnterTypeMonitorIC(MacroAssembler &masm,
                       size_t monitorStubOffset = ICMonitoredStub::offsetOfFirstMonitorStub())
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitReturnFromIC(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitChangeICReturnAddress(MacroAssembler &masm, Register reg)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitTailCallVM(JitCode *target, MacroAssembler &masm, uint32_t argSize)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitCreateStubFrameDescriptor(MacroAssembler &masm, Register reg)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitCallVM(JitCode *target, MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

// Size of vales pushed by EmitEnterStubFrame.
static const uint32_t STUB_FRAME_SIZE = 4 * sizeof(void *);
static const uint32_t STUB_FRAME_SAVED_STUB_OFFSET = sizeof(void *);

inline void
EmitEnterStubFrame(MacroAssembler &masm, Register scratch)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitLeaveStubFrame(MacroAssembler &masm, bool calledIntoIon = false)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitStowICValues(MacroAssembler &masm, int values)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitUnstowICValues(MacroAssembler &masm, int values, bool discard = false)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitCallTypeUpdateIC(MacroAssembler &masm, JitCode *code, uint32_t objectOffset)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

template <typename AddrType>
inline void
EmitPreBarrier(MacroAssembler &masm, const AddrType &addr, MIRType type)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}

inline void
EmitStubGuardFailure(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}


} // namespace jit
} // namespace js

#endif // JS_ION

#endif /* jit_mips_BaselineHelpers_mips_h */

