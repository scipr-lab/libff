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
#include "libff/algebra/fields/binary/gf32.hpp"
#include "libff/algebra/curves/edwards/edwards_pp.hpp"
#include "libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp"

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

template<typename FieldT>
void test_random_element()
{
    FieldT x = random_element_non_zero_one<FieldT>();
    EXPECT_NE(x, FieldT::zero());
    EXPECT_NE(x, FieldT::one());

    x = random_element_non_zero<FieldT>();
    EXPECT_NE(x, FieldT::zero());

    FieldT y = random_element_exclude(x);
    EXPECT_NE(x, y);
}

TEST(UtilsTest, RandomElementTest)
{
    init_edwards_fields();
    test_random_element<edwards_Fq3>();
    test_random_element<gf32>();
}

TEST(UtilsTest, CurveVectorSizeTest)
{
    init_edwards_params();
    init_mnt6_params();

    std::vector<edwards_G1> vec;

    vec.push_back(edwards_G1::G1_one);
    vec.push_back(edwards_G1::G1_zero);
    vec.push_back(edwards_G1::G1_one);

    EXPECT_EQ(curve_size_in_bits(vec), 552);

    std::vector<mnt6_G2> vec2;

    vec2.push_back(mnt6_G2::G2_zero);
    vec2.push_back(mnt6_G2::G2_one);
    vec2.push_back(mnt6_G2::G2_one);
    vec2.push_back(mnt6_G2::G2_zero);
    vec2.push_back(mnt6_G2::G2_zero);

    EXPECT_EQ(curve_size_in_bits(vec2), 4475);
}
