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

  private:
    mozilla::Array<double, FloatRegisters::Total> fpregs_;
    mozilla::Array<uintptr_t, Registers::Total> regs_;

    uintptr_t snapshotOffset_;

  public:
    FrameSizeClass frameClass() const {
        return FrameSizeClass::FromClass(frameClassId_);
    }
    uintptr_t tableOffset() const {
        JS_ASSERT(frameClass() != FrameSizeClass::None());
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
        JS_ASSERT(frameClass() == FrameSizeClass::None());
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
    MOZ_ASSUME_UNREACHABLE("NYI");
}

IonBailoutIterator::IonBailoutIterator(const JitActivationIterator &activations,
                                       InvalidationBailoutStack *bailout)
  : IonFrameIterator(activations),
    machine_(bailout->machine())
{
    MOZ_ASSUME_UNREACHABLE("NYI");
}
