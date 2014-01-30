/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 4 -*-
 * vim: set ts=8 sts=4 et sw=4 tw=99:
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "jit/BaselineCompiler.h"
#include "jit/BaselineHelpers.h"
#include "jit/BaselineIC.h"
#include "jit/BaselineJIT.h"
#include "jit/IonLinker.h"

#include "jsboolinlines.h"
#include "jsiter.h"

using namespace js;
using namespace js::jit;

namespace js {
namespace jit {

// ICCompare_Int32

bool
ICCompare_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
ICCompare_Double::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

// ICBinaryArith_Int32

bool
ICBinaryArith_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

bool
ICUnaryArith_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


// Mehtods that have been copied from BaselineIC.cpp because they need special
// implementation for MIPS.

bool
ICBinaryArith_BooleanWithInt32::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// Compare_Boolean
//

bool
ICCompare_Boolean::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// Compare_Int32WithBoolean
//

bool
ICCompare_Int32WithBoolean::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// ToBool_Int32
//

bool
ICToBool_Int32::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// ToBool_String
//

bool
ICToBool_String::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// ToBool_Double
//

bool
ICToBool_Double::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// ToBool_Object
//

bool
ICToBool_Object::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}

//
// IteratorMore_Native
//

bool
ICIteratorMore_Native::Compiler::generateStubCode(MacroAssembler &masm)
{
    MOZ_ASSUME_UNREACHABLE("NYI");
    return true;
}


} // namespace jit
} // namespace js
