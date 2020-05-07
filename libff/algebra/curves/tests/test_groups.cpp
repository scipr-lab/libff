/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <gtest/gtest.h>

#include <libff/algebra/curves/edwards/edwards_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>
#include <libff/common/profiling.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <sstream>

#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp>

using namespace libff;

class CurveGroupsTest: public ::testing::Test {
public:
    CurveGroupsTest()
    {
        edwards_pp::init_public_params();
        mnt4_pp::init_public_params();
        mnt6_pp::init_public_params();
        alt_bn128_pp::init_public_params();
#ifdef CURVE_BN128 // BN128 has fancy dependencies so it may be disabled
        bn128_pp::init_public_params();
#endif
    }
};

/** Returns a random element of FieldT that is not zero or one. */
template<typename GroupT>
GroupT random_element_non_zero_one()
{
    GroupT x = GroupT::random_element();
    while (x == GroupT::zero() || x == GroupT::one())
        x = GroupT::random_element();
    return x;
}

template<typename GroupT>
void test_mixed_add()
{
    GroupT base, el, result;

    base = GroupT::zero();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::zero();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = GroupT::zero();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = GroupT::random_element();
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base + el);

    base = GroupT::random_element();
    el = base;
    el.to_special();
    result = base.mixed_add(el);
    EXPECT_EQ(result, base.dbl());
}

template<typename GroupT>
void test_group()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    GroupT zero = GroupT::zero();
    EXPECT_EQ(zero, zero);
    GroupT one = GroupT::one();
    EXPECT_EQ(one, one);
    GroupT two = bigint<1>(2l) * GroupT::one();
    EXPECT_EQ(two, two);
    GroupT five = bigint<1>(5l) * GroupT::one();

    GroupT three = bigint<1>(3l) * GroupT::one();
    GroupT four = bigint<1>(4l) * GroupT::one();

    EXPECT_EQ(two+five, three+four);

    GroupT a = random_element_non_zero_one<GroupT>();
    GroupT b = random_element_non_zero_one<GroupT>();

    EXPECT_NE(one, zero);
    ASSERT_NE(a, zero);
    ASSERT_NE(a, one);
    ASSERT_NE(b, zero);
    ASSERT_NE(b, one);

    EXPECT_EQ(a.dbl(), a + a);
    EXPECT_EQ(b.dbl(), b + b);
    EXPECT_EQ(one.add(two), three);
    EXPECT_EQ(two.add(one), three);
    EXPECT_EQ(a + b, b + a);
    EXPECT_EQ(a - a, zero);
    EXPECT_EQ(a - b, a + (-b));
    EXPECT_EQ(a - b, (-b) + a);

    // handle special cases
    EXPECT_EQ(zero + (-a), -a);
    EXPECT_EQ(zero - a, -a);
    EXPECT_EQ(a - zero, a);
    EXPECT_EQ(a + zero, a);
    EXPECT_EQ(zero + a, a);

    EXPECT_EQ((a + b).dbl(), (a + b) + (b + a));
    EXPECT_EQ(bigint<1>("2") * (a + b), (a + b) + (b + a));

    EXPECT_EQ((rand1 * a) + (rand2 * a), (randsum * a));

    EXPECT_EQ(GroupT::order() * a, zero);
    EXPECT_EQ(GroupT::order() * one, zero);
    EXPECT_NE((GroupT::order() * a) - a, zero);
    EXPECT_NE((GroupT::order() * one) - one, zero);

    test_mixed_add<GroupT>();
}

template<typename GroupT>
void test_mul_by_q()
{
    GroupT a = GroupT::random_element();
    EXPECT_EQ((GroupT::field_char()*a), a.mul_by_q());
}

template<typename GroupT>
void test_output()
{
    GroupT g = GroupT::zero();

    for (size_t i = 0; i < 1000; ++i)
    {
        std::stringstream ss;
        ss << g;
        GroupT gg;
        ss >> gg;
        EXPECT_EQ(g, gg);
        /* use a random point in next iteration */
        g = GroupT::random_element();
    }
}

TEST_F(CurveGroupsTest, GroupTest)
{
    test_group<G1<edwards_pp> >();
    test_group<G2<edwards_pp> >();

    test_group<G1<mnt4_pp> >();
    test_group<G2<mnt4_pp> >();

    test_group<G1<mnt6_pp> >();
    test_group<G2<mnt6_pp> >();

    test_group<G1<alt_bn128_pp> >();
    test_group<G2<alt_bn128_pp> >();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    test_group<G1<bn128_pp> >();
    test_group<G2<bn128_pp> >();
#endif
}

TEST_F(CurveGroupsTest, OutputTest)
{
    test_output<G1<edwards_pp> >();
    test_output<G2<edwards_pp> >();

    test_output<G1<mnt4_pp> >();
    test_output<G2<mnt4_pp> >();

    test_output<G1<mnt6_pp> >();
    test_output<G2<mnt6_pp> >();

    test_output<G1<alt_bn128_pp> >();
    test_output<G2<alt_bn128_pp> >();

#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    test_output<G1<bn128_pp> >();
    test_output<G2<bn128_pp> >();
#endif
}

TEST_F(CurveGroupsTest, MulByQTest)
{
    test_mul_by_q<G2<edwards_pp> >();
    test_mul_by_q<G2<mnt4_pp> >();
    test_mul_by_q<G2<mnt6_pp> >();
    test_mul_by_q<G2<alt_bn128_pp> >();
}
