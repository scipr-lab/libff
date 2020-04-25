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
#include <libff/algebra/curves/alt_bn128/alt_bn128_pp.hpp>
#include <libff/algebra/fields/prime_extension/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime_extension/fp6_3over2.hpp>

namespace libff {

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

    // EXPECT_EQ((x * x).sqrt(), x);
    // EXPECT_EQ((x * (x + x + x + x)).sqrt(), two * x);
    // EXPECT_EQ(one.sqrt(), one);
    // EXPECT_EQ(zero.sqrt(), zero);
    // EXPECT_EQ((two + two).sqrt(), two);

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
    bigint<1> pow1 = bigint<1>("64703871");
    bigint<1> pow2 = bigint<1>("42796681");
    bigint<1> sum = bigint<1>("107500552");
    bigint<1> diff = bigint<1>("21907190");
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
}

TEST_F(AllFieldsTest, AllFieldsApiTest)
{
    test_field<AllFieldsTest::Fp>();
    test_field<AllFieldsTest::Fp2>();
}

}
