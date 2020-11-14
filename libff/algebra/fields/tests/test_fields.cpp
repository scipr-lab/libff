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
#include <libff/common/utils.hpp>
#ifdef CURVE_BN128
#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#endif
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/curves/bls12_381/bls12_381_pp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

using namespace libff;

class FpnFieldsTest: public ::testing::Test {
public:
    FpnFieldsTest()
    {
        init_edwards_fields();
        init_mnt4_fields();
        init_mnt6_fields();
        init_alt_bn128_fields();
        init_bls12_381_fields();
#ifdef CURVE_BN128 // BN128 has fancy dependencies so it may be disabled
        init_bn128_fields();
#endif
    }
};

template<typename FieldT>
void test_Frobenius()
{
    FieldT a = FieldT::random_element();
    EXPECT_EQ(a.Frobenius_map(0), a);
    FieldT a_q = a ^ FieldT::field_char();
    for (size_t power = 1; power < 10; ++power)
    {
        const FieldT a_qi = a.Frobenius_map(power);
        EXPECT_EQ(a_qi, a_q);

        a_q = a_q ^ FieldT::field_char();
    }
}

template<typename FieldT>
void test_sqrt()
{
    // sqrt() is very slow for bls12_381 fields, so this is only run 2 times.
    for (size_t i = 0; i < 2; i++)
    {
        FieldT a = FieldT::random_element();
        FieldT a_sq = a.squared();
        FieldT a_sq_sqrt = a_sq.sqrt();
        EXPECT_TRUE(a_sq_sqrt == a || a_sq_sqrt == -a);
    }
}

template<typename FieldT>
void test_field()
{
    bigint<1> rand1 = bigint<1>("76749407");
    bigint<1> rand2 = bigint<1>("44410867");
    bigint<1> randsum = bigint<1>("121160274");

    FieldT zero = FieldT::zero();
    FieldT one = FieldT::one();
    FieldT a = FieldT::random_element();
    FieldT a_ser;
    a_ser = reserialize(a);
    EXPECT_EQ(a_ser, a);

    FieldT b = FieldT::random_element();
    FieldT c = FieldT::random_element();
    FieldT d = FieldT::random_element();

    EXPECT_EQ(a * a, a.squared());
    EXPECT_EQ((a + b).squared(), a.squared() + a * b + b * a + b.squared());
    EXPECT_EQ((a + b) * (c + d), a * c + a * d + b * c + b * d);
    EXPECT_EQ(a - b, a + (-b));
    EXPECT_EQ(a - b, (-b) + a);

    EXPECT_EQ((a^rand1) * (a^rand2), a^randsum);

    a = random_element_non_zero<FieldT>();
    b = random_element_non_zero<FieldT>();
    c = random_element_non_zero<FieldT>();
    ASSERT_NE(a, zero);
    ASSERT_NE(b, zero);
    ASSERT_NE(c, zero);

    EXPECT_EQ(a * a.inverse(), one);
    EXPECT_EQ((a + b) * c.inverse(), a * c.inverse() + (b.inverse() * c).inverse());

    test_sqrt<FieldT>();
    test_Frobenius<FieldT>();
}

template<typename FieldT>
void test_two_squarings()
{
    FieldT a = FieldT::random_element();
    FieldT a_squared = a * a;
    EXPECT_EQ(a_squared, a.squared());
    EXPECT_EQ(a_squared, a.squared_complex());
    EXPECT_EQ(a_squared, a.squared_karatsuba());
}

template<typename FieldT>
void test_unitary_inverse()
{
    ASSERT_EQ(FieldT::extension_degree() % 2, 0);
    FieldT a = random_element_non_zero<FieldT>();
    ASSERT_NE(a, FieldT::zero());
    FieldT aqcubed_minus1 = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    ASSERT_NE(aqcubed_minus1, FieldT::zero());
    EXPECT_EQ(aqcubed_minus1.inverse(), aqcubed_minus1.unitary_inverse());
}

template<typename FieldT>
void test_cyclotomic_squaring();

template<>
void test_cyclotomic_squaring<Fqk<edwards_pp> >()
{
    typedef Fqk<edwards_pp> FieldT;
    EXPECT_EQ(FieldT::extension_degree() % 2, 0);
    FieldT a = random_element_non_zero<FieldT>();
    ASSERT_NE(a, FieldT::zero());
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^((q^(k/2)-1)*(q+1))
    FieldT beta = a_unitary.Frobenius_map(1) * a_unitary;
    EXPECT_EQ(beta.cyclotomic_squared(), beta.squared());
}

template<>
void test_cyclotomic_squaring<Fqk<mnt4_pp> >()
{
    typedef Fqk<mnt4_pp> FieldT;
    ASSERT_EQ(FieldT::extension_degree() % 2, 0);
    FieldT a = random_element_non_zero<FieldT>();
    ASSERT_NE(a, FieldT::zero());
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^(q^(k/2)-1)
    FieldT beta = a_unitary;
    EXPECT_EQ(beta.cyclotomic_squared(), beta.squared());
}

template<>
void test_cyclotomic_squaring<Fqk<mnt6_pp> >()
{
    typedef Fqk<mnt6_pp> FieldT;
    ASSERT_EQ(FieldT::extension_degree() % 2, 0);
    FieldT a = random_element_non_zero<FieldT>();
    ASSERT_NE(a, FieldT::zero());
    FieldT a_unitary = a.Frobenius_map(FieldT::extension_degree()/2) * a.inverse();
    // beta = a^((q^(k/2)-1)*(q+1))
    FieldT beta = a_unitary.Frobenius_map(1) * a_unitary;
    EXPECT_EQ(beta.cyclotomic_squared(), beta.squared());
}

template<typename ppT>
void test_all_fields()
{
    test_field<Fr<ppT> >();
    test_field<Fq<ppT> >();
    test_field<Fqe<ppT> >();
    test_field<Fqk<ppT> >();

    test_unitary_inverse<Fqk<ppT> >();
}

template<typename Fp4T>
void test_Fp4_toom_cook()
{
    typedef typename Fp4T::my_Fp FieldT;
    for (size_t i = 0; i < 100; ++i)
    {
        const Fp4T a = Fp4T::random_element();
        const Fp4T b = Fp4T::random_element();
        const Fp4T correct_res = a * b;

        Fp4T res;

        const FieldT
            &a0 = a.c0.c0,
            &a1 = a.c1.c0,
            &a2 = a.c0.c1,
            &a3 = a.c1.c1;

        const FieldT
            &b0 = b.c0.c0,
            &b1 = b.c1.c0,
            &b2 = b.c0.c1,
            &b3 = b.c1.c1;

        FieldT
            &c0 = res.c0.c0,
            &c1 = res.c1.c0,
            &c2 = res.c0.c1,
            &c3 = res.c1.c1;

        const FieldT v0 = a0 * b0;
        const FieldT v1 = (a0 + a1 + a2 + a3) * (b0 + b1 + b2 + b3);
        const FieldT v2 = (a0 - a1 + a2 - a3) * (b0 - b1 + b2 - b3);
        const FieldT v3 = (a0 + FieldT(2)*a1 + FieldT(4)*a2 + FieldT(8)*a3) * (b0 + FieldT(2)*b1 + FieldT(4)*b2 + FieldT(8)*b3);
        const FieldT v4 = (a0 - FieldT(2)*a1 + FieldT(4)*a2 - FieldT(8)*a3) * (b0 - FieldT(2)*b1 + FieldT(4)*b2 - FieldT(8)*b3);
        const FieldT v5 = (a0 + FieldT(3)*a1 + FieldT(9)*a2 + FieldT(27)*a3) * (b0 + FieldT(3)*b1 + FieldT(9)*b2 + FieldT(27)*b3);
        const FieldT v6 = a3 * b3;

        const FieldT beta = Fp4T::non_residue;

        c0 = v0 + beta*(FieldT(4).inverse()*v0 - FieldT(6).inverse()*(v1 + v2) + FieldT(24).inverse() * (v3 + v4) - FieldT(5) * v6);
        c1 = - FieldT(3).inverse()*v0 + v1 - FieldT(2).inverse()*v2 - FieldT(4).inverse()*v3 + FieldT(20).inverse() * v4 + FieldT(30).inverse() * v5 - FieldT(12) * v6 + beta * ( - FieldT(12).inverse() * (v0 - v1) + FieldT(24).inverse()*(v2 - v3) - FieldT(120).inverse() * (v4 - v5) - FieldT(3) * v6);
        c2 = - (FieldT(5)*(FieldT(4).inverse()))* v0 + (FieldT(2)*(FieldT(3).inverse()))*(v1 + v2) - FieldT(24).inverse()*(v3 + v4) + FieldT(4)*v6 + beta*v6;
        c3 = FieldT(12).inverse() * (FieldT(5)*v0 - FieldT(7)*v1) - FieldT(24).inverse()*(v2 - FieldT(7)*v3 + v4 + v5) + FieldT(15)*v6;

        EXPECT_EQ(res, correct_res);

        // {v0, v3, v4, v5}
        const FieldT u = (FieldT::one() - beta).inverse();
        EXPECT_EQ(v0, u * c0 + beta * u * c2 - beta * u * FieldT(2).inverse() * v1 - beta * u * FieldT(2).inverse() * v2 + beta * v6);
        EXPECT_EQ(v3, - FieldT(15) * u * c0 - FieldT(30) * u * c1 - FieldT(3) * (FieldT(4) + beta) * u * c2 - FieldT(6) * (FieldT(4) + beta) * u * c3 + (FieldT(24) - FieldT(3) * beta * FieldT(2).inverse()) * u * v1 + (-FieldT(8) + beta * FieldT(2).inverse()) * u * v2
               - FieldT(3) * (-FieldT(16) + beta) * v6);
        EXPECT_EQ(v4, - FieldT(15) * u * c0 + FieldT(30) * u * c1 - FieldT(3) * (FieldT(4) + beta) * u * c2 + FieldT(6) * (FieldT(4) + beta) * u * c3 + (FieldT(24) - FieldT(3) * beta * FieldT(2).inverse()) * u * v2 + (-FieldT(8) + beta * FieldT(2).inverse()) * u * v1
               - FieldT(3) * (-FieldT(16) + beta) * v6);
        EXPECT_EQ(v5, - FieldT(80) * u * c0 - FieldT(240) * u * c1 - FieldT(8) * (FieldT(9) + beta) * u * c2 - FieldT(24) * (FieldT(9) + beta) * u * c3 - FieldT(2) * (-FieldT(81) + beta) * u * v1 + (-FieldT(81) + beta) * u * v2
               - FieldT(8) * (-FieldT(81) + beta) * v6);

        // c0 + beta c2 - (beta v1)/2 - (beta v2)/ 2 - (-1 + beta) beta v6,
        // -15 c0 - 30 c1 - 3 (4 + beta) c2 - 6 (4 + beta) c3 + (24 - (3 beta)/2) v1 + (-8 + beta/2) v2 + 3 (-16 + beta) (-1 + beta) v6,
        // -15 c0 + 30 c1 - 3 (4 + beta) c2 + 6 (4 + beta) c3 + (-8 + beta/2) v1 + (24 - (3 beta)/2) v2 + 3 (-16 + beta) (-1 + beta) v6,
        // -80 c0 - 240 c1 - 8 (9 + beta) c2 - 24 (9 + beta) c3 - 2 (-81 + beta) v1 + (-81 + beta) v2 + 8 (-81 + beta) (-1 + beta) v6
    }
}

TEST_F(FpnFieldsTest, GeneralTest)
{
    test_all_fields<edwards_pp>();

    test_all_fields<mnt4_pp>();

    test_all_fields<mnt6_pp>();

    test_all_fields<alt_bn128_pp>();
    test_field<alt_bn128_Fq6>(); // alt_bn128_Fq6 is not included in test_all_fields.

    test_all_fields<bls12_381_pp>();
    test_field<bls12_381_Fq6>(); // bls12_381_Fq6 is not included in test_all_fields.
#ifdef CURVE_BN128       // BN128 has fancy dependencies so it may be disabled
    test_field<Fr<bn128_pp> >();
    test_field<Fq<bn128_pp> >();
#endif
}

TEST_F(FpnFieldsTest, SquareTest)
{
    test_cyclotomic_squaring<Fqk<edwards_pp> >();

    test_two_squarings<Fqe<mnt4_pp> >();
    test_cyclotomic_squaring<Fqk<mnt4_pp> >();
}

TEST_F(FpnFieldsTest, TomCookTest)
{
    test_Fp4_toom_cook<mnt4_Fq4>();
}
