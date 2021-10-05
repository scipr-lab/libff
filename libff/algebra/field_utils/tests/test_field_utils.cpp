/**
 *****************************************************************************
 Basic tests for some of the field utils in this directory, mainly bigint
 and power.
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include "libff/algebra/field_utils/bigint.hpp"
#include "libff/algebra/field_utils/field_utils.hpp"
#include "libff/algebra/field_utils/algorithms.hpp"
#include "libff/algebra/fields/binary/gf64.hpp"
#include "libff/algebra/curves/edwards/edwards_fields.hpp"
#include <gtest/gtest.h>

using namespace libff;

using std::size_t;

template<typename FieldT>
FieldT power_naive(const FieldT &base, const std::size_t exponent)
{
    FieldT result = FieldT::one();

    for (std::size_t i = 1; i <= exponent; ++i)
    {
        result *= base;
    }

    return result;
}


TEST(ExponentiationTest, SimpleTest) {
    typedef gf64 FieldT;

    const unsigned long max_power = 3000;
    FieldT X = FieldT::random_element();

    FieldT X_i = FieldT::one();
    for (unsigned long i = 0 ; i < max_power; ++i)
    {
        const FieldT X_i_naive = power_naive<FieldT>(X, i);
        const FieldT X_i_square_and_multiply_ul = power<FieldT>(X, i);
        const FieldT X_i_square_and_multiply_ull = power<FieldT>(X, (unsigned long long) i);

        EXPECT_EQ(X_i, X_i_naive);
        EXPECT_EQ(X_i, X_i_square_and_multiply_ul);
        EXPECT_EQ(X_i, X_i_square_and_multiply_ull);

        X_i *= X;
    }
}

TEST(FieldUtilsTest, BigintTest)
{
    bigint<3> zero = bigint<3>("0");
    bigint<3> one = bigint<3>("1");
    bigint<3> x = bigint<3>("987654567895678909876876545678909876543456");
    bigint<3> y = bigint<3>("324531232345676543272920293863628304859020");

    bigint<3> zero2;
    bigint<3> z = bigint<3>("987654567895678909876876545678909876543456");
    bigint<3> w = bigint<3>("987654567895678909876876545678909876543451");

    EXPECT_EQ(zero, zero2);
    EXPECT_EQ(one, bigint<3>::one());
    EXPECT_EQ(x, x);
    EXPECT_EQ(x, z);

    EXPECT_NE(zero, one);
    EXPECT_NE(x, one);
    EXPECT_NE(x, zero);
    EXPECT_NE(x, y);

    EXPECT_FALSE(x.is_zero());
    EXPECT_FALSE(one.is_zero());
    EXPECT_TRUE(zero.is_zero());

    EXPECT_TRUE(x.is_even());
    EXPECT_TRUE(y.is_even());
    EXPECT_TRUE(zero.is_even());

    EXPECT_FALSE(one.is_even());
    EXPECT_FALSE(w.is_even());

    EXPECT_LT(zero, one);
    EXPECT_LT(zero, x);
    EXPECT_LT(zero, y);
    EXPECT_LT(one, x);
    EXPECT_LT(y, x);
    EXPECT_LT(w, x);
    EXPECT_LT(y, w);
    EXPECT_FALSE(x < y);

    x.print();
    x.print_hex();

    x.randomize();
    EXPECT_EQ(x, x);
    x.clear();
    EXPECT_EQ(x, zero);
}

TEST(FieldUtilsTest, FieldVectorConversionTest)
{
    init_edwards_fields();

    // pack_bit_vector_into_field_element_vector

    bit_vector vec;
    for (size_t i = 0; i < 12 + edwards_Fq::ceil_size_in_bits(); i++)
        vec.push_back(0);
    vec.push_back(1);
    vec.push_back(0);
    vec.push_back(1);

    std::vector<edwards_Fq> field_vec = pack_bit_vector_into_field_element_vector<edwards_Fq>(vec);

    EXPECT_EQ(field_vec.size(), 2);
    EXPECT_EQ(field_vec[0], edwards_Fq::zero());
    EXPECT_EQ(field_vec[1], edwards_Fq(40960)); // 5 * 2**13

    // convert_bit_vector_to_field_element_vector

    bit_vector vec2;
    vec2.push_back(0);
    vec2.push_back(0);
    vec2.push_back(1);
    vec2.push_back(0);
    vec2.push_back(1);

    field_vec = convert_bit_vector_to_field_element_vector<edwards_Fq>(vec2);

    EXPECT_EQ(field_vec.size(), 5);
    EXPECT_EQ(field_vec[0], edwards_Fq::zero());
    EXPECT_EQ(field_vec[1], edwards_Fq::zero());
    EXPECT_EQ(field_vec[2], edwards_Fq::one());
    EXPECT_EQ(field_vec[3], edwards_Fq::zero());
    EXPECT_EQ(field_vec[4], edwards_Fq::one());

    // convert_field_element_vector_to_bit_vector

    std::vector<edwards_Fq> field_vec2;
    field_vec2.push_back(edwards_Fq(edwards_Fq(5)));
    field_vec2.push_back(edwards_Fq::zero());

    bit_vector vec3 = convert_field_element_vector_to_bit_vector(field_vec2);

    EXPECT_EQ(vec3.size(), edwards_Fq::ceil_size_in_bits() * 2);
    EXPECT_EQ(vec3[0], 1);
    EXPECT_EQ(vec3[1], 0);
    EXPECT_EQ(vec3[2], 1);
    for (size_t i = 3; i < vec3.size(); i++)
        EXPECT_EQ(vec3[i], 0);
}
