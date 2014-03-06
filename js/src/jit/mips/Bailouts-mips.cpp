/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jscntxt.h"
#include "jscompartment.h"

#include "jit/Bailouts.h"
#include "jit/JitCompartment.h"

using namespace js;
using namespace js::jit;

namespace js {
namespace jit {

class BailoutStack
{
    uintptr_t frameClassId_;
    // This is pushed in the bailout handler.  Both entry points into the
    // handler inserts their own value int lr, which is then placed onto the
    // stack along with frameClassId_ above.  This should be migrated to ip.
  public:
    union {
        uintptr_t frameSize_;
        uintptr_t tableOffset_;
    };

  protected:
    mozilla::Array<double, FloatRegisters::Total> fpregs_;
    mozilla::Array<uintptr_t, Registers::Total> regs_;

    uintptr_t snapshotOffset_;
    uintptr_t padding_;

  public:
    FrameSizeClass frameClass() const {
        return FrameSizeClass::FromClass(frameClassId_);
    }
    uintptr_t tableOffset() const {
        MOZ_ASSERT(frameClass() != FrameSizeClass::None());
        return tableOffset_;
    }
    uint32_t frameSize() const {
        if (frameClass() == FrameSizeClass::None())
            return frameSize_;
        return frameClass().frameSize();
    }
    MachineState machine() {
        return MachineState::FromBailout(regs_, fpregs_);
    }
    SnapshotOffset snapshotOffset() const {
        MOZ_ASSERT(frameClass() == FrameSizeClass::None());
        return snapshotOffset_;
    }
    uint8_t *parentStackPointer() const {
        if (frameClass() == FrameSizeClass::None())
            return (uint8_t *)this + sizeof(BailoutStack);
        return (uint8_t *)this + offsetof(BailoutStack, snapshotOffset_);
    }
};

} // namespace jit
} // namespace js

IonBailoutIterator::IonBailoutIterator(const JitActivationIterator &activations,
                                       BailoutStack *bailout)
  : IonFrameIterator(activations),
    machine_(bailout->machine())
{
    uint8_t *sp = bailout->parentStackPointer();
    uint8_t *fp = sp + bailout->frameSize();

    current_ = fp;
    type_ = IonFrame_OptimizedJS;
    topFrameSize_ = current_ - sp;
    topIonScript_ = script()->ionScript();

    if (bailout->frameClass() == FrameSizeClass::None()) {
        snapshotOffset_ = bailout->snapshotOffset();
        return;
    }

    // Compute the snapshot offset from the bailout ID.
    JitActivation *activation = activations.activation()->asJit();
    JSRuntime *rt = activation->compartment()->runtimeFromMainThread();
    JitCode *code = rt->jitRuntime()->getBailoutTable(bailout->frameClass());
    uintptr_t tableOffset = bailout->tableOffset();
    uintptr_t tableStart = reinterpret_cast<uintptr_t>(code->raw());

    MOZ_ASSERT(tableOffset >= tableStart &&
              tableOffset < tableStart + code->instructionsSize());
    MOZ_ASSERT((tableOffset - tableStart) % BAILOUT_TABLE_ENTRY_SIZE == 0);

    uint32_t bailoutId = ((tableOffset - tableStart) / BAILOUT_TABLE_ENTRY_SIZE) - 1;
    MOZ_ASSERT(bailoutId < BAILOUT_TABLE_SIZE);

    snapshotOffset_ = topIonScript_->bailoutToSnapshot(bailoutId);
}

IonBailoutIterator::IonBailoutIterator(const JitActivationIterator &activations,
                                       InvalidationBailoutStack *bailout)
  : IonFrameIterator(activations),
    machine_(bailout->machine())
{
    returnAddressToFp_ = bailout->osiPointReturnAddress();
    topIonScript_ = bailout->ionScript();
    const OsiIndex *osiIndex = topIonScript_->getOsiIndex(returnAddressToFp_);

    current_ = (uint8_t*) bailout->fp();
    type_ = IonFrame_OptimizedJS;
    topFrameSize_ = current_ - bailout->sp();
    snapshotOffset_ = osiIndex->snapshotOffset();
}
