/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* vim: set ts=8 sts=2 et sw=2 tw=80: */
/* Copyright 2013 Mozilla Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef mozilla_pkix__nssgtest_h
#define mozilla_pkix__nssgtest_h

#include "stdint.h"
#include "gtest/gtest.h"
#include "prerror.h"
#include "seccomon.h"

namespace mozilla { namespace pkix { namespace test {

class SECStatusWithPRErrorCode
{
public:
  SECStatusWithPRErrorCode(SECStatus rv, PRErrorCode errorCode)
    : mRv(rv)
    , mErrorCode(errorCode)
  {
  }

  SECStatusWithPRErrorCode(SECStatus rv)
    : mRv(rv)
    , mErrorCode(rv == SECSuccess ? 0 : PR_GetError())
  {
  }

  bool operator==(const SECStatusWithPRErrorCode& other) const
  {
    return mRv == other.mRv && mErrorCode == other.mErrorCode;
  }

private:
  const SECStatus mRv;
  const PRErrorCode mErrorCode;

  friend std::ostream& operator<<(std::ostream& os,
                                  SECStatusWithPRErrorCode const& value);

  void operator=(const SECStatusWithPRErrorCode&) /*delete*/;
};

::std::ostream& operator<<(::std::ostream&,
                           SECStatusWithPRErrorCode const&);

} } } // namespace mozilla::pkix::test

#define ASSERT_SECSuccess(rv) \
  ASSERT_EQ(::mozilla::pkix::test::SECStatusWithPRErrorCode(SECSuccess, 0), \
            ::mozilla::pkix::test::SECStatusWithPRErrorCode(rv))
#define EXPECT_SECSuccess(rv) \
  EXPECT_EQ(::mozilla::pkix::test::SECStatusWithPRErrorCode(SECSuccess, 0), \
            ::mozilla::pkix::test::SECStatusWithPRErrorCode(rv))

#define ASSERT_SECFailure(expectedError, rv) \
  ASSERT_EQ(::mozilla::pkix::test::SECStatusWithPRErrorCode(SECFailure, \
                                                            expectedError), \
            ::mozilla::pkix::test::SECStatusWithPRErrorCode(rv))
#define EXPECT_SECFailure(expectedError, rv) \
  EXPECT_EQ(::mozilla::pkix::test::SECStatusWithPRErrorCode(SECFailure, \
                                                            expectedError), \
            ::mozilla::pkix::test::SECStatusWithPRErrorCode(rv))

#endif // mozilla_pkix__nssgtest_h
