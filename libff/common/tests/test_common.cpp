/**
 *****************************************************************************
 Some tests for the functions in this directory.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cstdint>
#include <gtest/gtest.h>

#include "libff/common/utils.hpp"

using namespace libff;

TEST(Log2Test, SimpleTest) {
    // There seems to be a second log2 function that operates on floats so we added libff::.
    EXPECT_EQ(libff::log2(0), 0ull);
    EXPECT_EQ(libff::log2(1), 0ull);
    EXPECT_EQ(libff::log2(2), 1ull);
    EXPECT_EQ(libff::log2(3), 2ull);
    EXPECT_EQ(libff::log2(4), 2ull);
    EXPECT_EQ(libff::log2(5), 3ull);
    EXPECT_EQ(libff::log2(6), 3ull);
    EXPECT_EQ(libff::log2(7), 3ull);
    EXPECT_EQ(libff::log2(8), 3ull);
    EXPECT_EQ(libff::log2(9), 4ull);
}

TEST(Log2Test, PowersOfTwo) {
    for (std::size_t i = 10; i < 20; ++i)
    {
        const std::size_t k = (1ull<<i);
        EXPECT_EQ(libff::log2(k-1), i);
        EXPECT_EQ(libff::log2(k), i);
        EXPECT_EQ(libff::log2(k+1), i+1);
    }
}
