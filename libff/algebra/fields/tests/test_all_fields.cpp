/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
#include <gtest/gtest.h>
#include <set>

#include "libff/algebra/curves/mnt/mnt4/mnt4_fields.hpp"
#include "libff/algebra/curves/mnt/mnt6/mnt6_fields.hpp"
#include "libff/algebra/curves/alt_bn128/alt_bn128_fields.hpp"
#include "libff/algebra/fields/binary/gf32.hpp"
#include "libff/algebra/fields/binary/gf64.hpp"
#include "libff/algebra/fields/binary/gf128.hpp"
#include "libff/algebra/fields/binary/gf192.hpp"
#include "libff/algebra/fields/binary/gf256.hpp"

using namespace libff;

class AllFieldsTest: public ::testing::Test {
public:
    // p, q, r are three different primes
    typedef alt_bn128_Fq Fp;
    typedef alt_bn128_Fq2 Fp2;
    typedef alt_bn128_Fq6 Fp6_3_2;
    typedef alt_bn128_Fq12 Fp12_2_3_2;

    typedef mnt4_Fq Fq;
    typedef mnt4_Fq2 Fq2;
    typedef mnt4_Fq4 Fq4;

    typedef mnt6_Fq Fr;
    typedef mnt6_Fq3 Fr3;
    typedef mnt6_Fq6 Fr6_2_3;
    
    AllFieldsTest()
    {
        init_alt_bn128_fields();
        init_mnt4_fields();
        init_mnt6_fields();
    }
};

/** Returns a random element of FieldT that is not equal to y. */
template<typename FieldT>
FieldT random_element_exclude(FieldT y)
{
    FieldT x = FieldT::random_element();
    while (x == y)
        x = FieldT::random_element();
    return x;
}

template<typename FieldT>
void expect_equal_or_negative(FieldT x, FieldT y)
{
    EXPECT_TRUE(x == y || x == -y);
}

template<typename FieldT>
void test_field()
{
    // constants
    const FieldT one = FieldT::one();
    const FieldT zero = FieldT::zero();
    const FieldT two = one + one;

    /******************* Test standard field axioms and properties. *******************/

    FieldT x = FieldT::random_element();
    FieldT y = FieldT::random_element();
    FieldT z = FieldT::random_element();
    FieldT w = random_element_exclude(zero);
    
    EXPECT_EQ(x + y, y + x); // commutative law of addition
    EXPECT_EQ((x + y) + z, x + (y + z)); // associative law of addition
    EXPECT_EQ(x + zero, x); // additive identity
    EXPECT_EQ(x + (-x), zero); // additive inverse
    EXPECT_EQ(x * y, y * x); // commutative law of multiplication
    EXPECT_EQ((x * y) * z, x * (y * z)); // associative law of multiplication
    EXPECT_EQ(x * one, x); // multiplicative identity
    EXPECT_EQ(w * w.inverse(), one); // multiplicative inverse
    EXPECT_EQ((x + y) * z, x * z + y * z); // distributive law

    EXPECT_EQ(-zero, zero);
    EXPECT_EQ(one.inverse(), one);
    EXPECT_EQ(-(-x), x);
    EXPECT_EQ(w.inverse().inverse(), w);
    EXPECT_EQ(x * zero, zero);
    EXPECT_EQ(x * (-y), -(x * y));

    /*********************** Test +, -, *, zero(), and one(). ***********************/

    x.randomize();
    y = x + x;
    EXPECT_EQ(x * two, y);
    EXPECT_EQ(x + y, x + x + x);
    EXPECT_EQ(two + two, two * two);

    EXPECT_EQ(x - x, zero);
    EXPECT_EQ(x - y, -x);
    y.randomize();
    EXPECT_EQ(x + (-y), x - y);
    EXPECT_EQ(x - y, -(y - x));

    EXPECT_EQ(x * (-one), -x);
    EXPECT_EQ((-x) * (-y), x * y);
    z.randomize();
    w.randomize();
    EXPECT_EQ((x + y) * (z + w), x * z + x * w + y * z + y * w);

    EXPECT_NE(zero, one);
    EXPECT_NE(one, two);
    EXPECT_NE(x + one, x);
    x = random_element_exclude(zero);
    y = random_element_exclude(zero);
    z = random_element_exclude(one);
    if (two == zero)
        EXPECT_EQ(x, -x);
    else
        EXPECT_NE(x, -x);
    EXPECT_NE(x + y, x);
    EXPECT_NE(x * z, x);
    y = random_element_exclude(x);
    z = random_element_exclude(two);
    EXPECT_NE(x - y, zero);
    EXPECT_NE(x * z, x + x);

    // test assignment
    x.randomize();
    y = x;
    x += one;
    EXPECT_NE(x, y);
    y += one;
    EXPECT_EQ(x, y);
    x.square();
    EXPECT_EQ(x, y * y);

    /******************** Test squared(), inverse(), sqrt(), and ^. ********************/

    x.randomize();
    FieldT x2 = x * x;
    y = x;
    EXPECT_EQ(x.squared(), x2);
    EXPECT_EQ(x, y);
    x += x;
    EXPECT_EQ(x, y + y);
    y.randomize();
    EXPECT_EQ(x.squared() + two * x * y + y.squared(), (x + y).squared());
    EXPECT_EQ((x * y).squared(), x.squared() * y.squared());

    expect_equal_or_negative((x * x).sqrt(), x);
    expect_equal_or_negative((x * (x + x + x + x)).sqrt(), two * x);
    expect_equal_or_negative(one.sqrt(), one);
    expect_equal_or_negative((two + two).sqrt(), two);
    expect_equal_or_negative(zero.sqrt(), zero);

    x = random_element_exclude(zero);
    y = random_element_exclude(zero);
    EXPECT_EQ(x.squared().inverse(), x.inverse().squared());
    EXPECT_EQ((x * y).inverse(), x.inverse() * y.inverse());
    EXPECT_EQ((x * y.inverse()).inverse(), x.inverse() * y);

    x.randomize();
    y.randomize();
    EXPECT_EQ(x.squared(), x^2);
    EXPECT_EQ(x * x * x, x^3);
    EXPECT_EQ(x.squared().squared().squared() * x, x^9);

    // I tried using random bigints, but it was too much work
    const bigint<1> pow1 = bigint<1>("64703871");
    const bigint<1> pow2 = bigint<1>("42796681");
    const bigint<1> sum = bigint<1>("107500552");
    const bigint<1> diff = bigint<1>("21907190");
    EXPECT_EQ((x^pow1) * (x^pow2), x^sum);
    EXPECT_EQ((x * y)^pow1, (x^pow1) * (y^pow1));

    x = random_element_exclude(zero);
    EXPECT_EQ(x.inverse()^pow1, (x^pow1).inverse());
    EXPECT_EQ((x^pow1) * (x.inverse()^pow2), x^diff);

    /******************** Test +=, -=, *=, square(), inverse(), ^=. ********************/

    x.randomize();
    y.randomize();
    z = x + y;
    x += y;
    EXPECT_EQ(x, z);
    x.randomize();
    z = x - y;
    x -= y;
    EXPECT_EQ(x, z);
    x.randomize();
    z = x * y;
    x *= y;
    EXPECT_EQ(x, z);

    x.randomize();
    z = x.squared();
    x.square();
    EXPECT_EQ(x, z);
    x = random_element_exclude(zero);
    z = x.inverse();
    x.invert();
    EXPECT_EQ(x, z);
    x.randomize();
    z = x^82;
    x ^= 82;
    EXPECT_EQ(x, z);
    x.randomize();
    z = x^pow1;
    x ^= pow1;
    EXPECT_EQ(x, z);

    /****************** Test is_zero(), print(), clear(), <<, and >>. ******************/

    EXPECT_TRUE(zero.is_zero());
    EXPECT_FALSE(one.is_zero());
    EXPECT_FALSE(random_element_exclude(zero).is_zero());

    EXPECT_NO_THROW(x.print());

    x.clear();
    EXPECT_EQ(x, zero);
    EXPECT_TRUE(x.is_zero());

    x.randomize();
    EXPECT_EQ(reserialize<FieldT>(x), x);

    /****************** Test extension_degree() and size_in_bits(). ******************/

    EXPECT_GE(FieldT::extension_degree(), 1);
    EXPECT_GE(FieldT::size_in_bits(), 1);
}

template<typename FieldT>
void test_op_profiling()
{
    FieldT::add_cnt = 0;
    FieldT::sub_cnt = 0;
    FieldT::mul_cnt = 0;
    FieldT::sqr_cnt = 0;
    FieldT::inv_cnt = 0;

    FieldT x = FieldT::random_element();
    FieldT y = FieldT::random_element();
    FieldT one = FieldT::one();
    FieldT zero = FieldT::zero();
    x += y;
    y = x + one;
    x *= x + y;
    y = x * y + one.inverse() - (x * one);
    y.square();
    x -= y.squared();
    x = random_element_exclude(zero);
    x.invert();
    x = x - one;
    x *= x * x.squared();

    EXPECT_EQ(FieldT::add_cnt, 4);
    EXPECT_EQ(FieldT::sub_cnt, 3);
    EXPECT_EQ(FieldT::mul_cnt, 5);
    EXPECT_EQ(FieldT::sqr_cnt, 3);
    EXPECT_EQ(FieldT::inv_cnt, 2);
}

template<typename FieldT>
void test_prime_field()
{
    EXPECT_GE(FieldT::field_char().num_bits(), 2); // characteristic is at least 2

    FieldT x = FieldT::random_element();
    FieldT x_q = x;
    for (size_t power = 0; power < 10; ++power)
    {
        const FieldT x_qi = x.Frobenius_map(power);
        EXPECT_EQ(x_qi, x_q);

        x_q = x_q^FieldT::field_char();
    }
}

template<typename FieldT>
void test_binary_field()
{
    FieldT zero = FieldT::zero();
    FieldT one = FieldT::one();
    FieldT x = FieldT::random_element();
    FieldT y = random_element_exclude<FieldT>(x);
    FieldT z = x;

    std::vector<uint64_t> x_words = x.as_words();
    std::vector<uint64_t> y_words = y.as_words();
    std::vector<uint64_t> z_words = z.as_words();
    std::vector<uint64_t> zero_words = zero.as_words();
    EXPECT_NE(x_words, y_words);
    EXPECT_EQ(x_words, z_words);
    for (uint64_t word : zero_words)
        EXPECT_EQ(word, 0);

    EXPECT_GE(FieldT::modulus_, 1);
    const uint64_t bits = FieldT::num_bits;
    EXPECT_GE(bits, 1);
    const bigint<1> characteristic = FieldT::template field_char<1>();
    EXPECT_EQ(characteristic, bigint<1>(2));

    FieldT generator = FieldT::multiplicative_generator;
    x = generator;
    EXPECT_NE(generator, zero);
    EXPECT_NE(generator, one);
    std::set<std::vector<uint64_t> > values;
    for (uint16_t i = 0; i < 10000; i++)
    {
        if (x == one)
            break;
        EXPECT_EQ(values.find(x.as_words()), values.end()); // generator^n never repeats
        values.insert(x.as_words());
        x *= generator;
    }
}

TEST_F(AllFieldsTest, AllFieldsApiTest)
{
    test_field<AllFieldsTest::Fp>();
    test_field<AllFieldsTest::Fp2>();
    test_field<AllFieldsTest::Fp6_3_2>();
    test_field<AllFieldsTest::Fp12_2_3_2>();

    test_field<AllFieldsTest::Fq4>();

    test_field<AllFieldsTest::Fr3>();
    test_field<AllFieldsTest::Fr6_2_3>();

    // test_field<gf32>();
    // test_field<gf64>();
    // test_field<gf128>();
    // test_field<gf192>();
    // test_field<gf256>();
}

#ifdef PROFILE_OP_COUNTS
TEST_F(AllFieldsTest, AllFieldsOpCountTest)
{
    test_op_profiling<AllFieldsTest::Fp>();
    test_op_profiling<AllFieldsTest::Fp2>();
    test_op_profiling<AllFieldsTest::Fp6_3_2>();
    test_op_profiling<AllFieldsTest::Fp12_2_3_2>();

    test_op_profiling<AllFieldsTest::Fq4>();

    test_op_profiling<AllFieldsTest::Fr3>();
    test_op_profiling<AllFieldsTest::Fr6_2_3>();

    test_op_profiling<gf32>();
    test_op_profiling<gf64>();
    test_op_profiling<gf128>();
    test_op_profiling<gf192>();
    test_op_profiling<gf256>();
}
#endif

TEST_F(AllFieldsTest, PrimeFieldsApiTest)
{
    test_prime_field<AllFieldsTest::Fp>();
    test_prime_field<AllFieldsTest::Fp2>();
    test_prime_field<AllFieldsTest::Fp6_3_2>();
    test_prime_field<AllFieldsTest::Fp12_2_3_2>();

    test_prime_field<AllFieldsTest::Fq4>();

    test_prime_field<AllFieldsTest::Fr3>();
    test_prime_field<AllFieldsTest::Fr6_2_3>();
}

TEST_F(AllFieldsTest, BinaryFieldsApiTest)
{
    test_binary_field<gf32>();
    test_binary_field<gf64>();
    test_binary_field<gf128>();
    test_binary_field<gf192>();
    test_binary_field<gf256>();
}
