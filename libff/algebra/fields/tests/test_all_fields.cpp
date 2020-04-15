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
#include <libff/algebra/fields/prime/fp12_2over3over2.hpp>
#include <libff/algebra/fields/prime/fp6_3over2.hpp>

namespace libff {

class AllFieldsTest: public ::testing::Test { 
public:
    static const mp_size_t bitcount = 254;
    static const mp_size_t limbs = (bitcount+GMP_NUMB_BITS-1)/GMP_NUMB_BITS;

    static bigint<limbs> p;

    typedef Fp_model<limbs, p> Fp;
    typedef Fp2_model<limbs, p> Fp2;
    typedef Fp6_3over2_model<limbs, p> Fp6_3_2;
    typedef Fp12_2over3over2_model<limbs, p> Fp12_2_3_2;
    
    AllFieldsTest()
    { 
        init_Fp();
        init_Fp2();
    }

private:
    /** Same as alt_bn128_Fq */
    void init_Fp()
    {
        p = bigint<limbs>("21888242871839275222246405745257275088696311157297823662689037894645226208583");
        assert(Fp::modulus_is_valid());
        if (sizeof(mp_limb_t) == 8)
        {
            Fp::Rsquared = bigint<limbs>("3096616502983703923843567936837374451735540968419076528771170197431451843209");
            Fp::Rcubed = bigint<limbs>("14921786541159648185948152738563080959093619838510245177710943249661917737183");
            Fp::inv = 0x87d20782e4866389;
        }
        else if (sizeof(mp_limb_t) == 4)
        {
            Fp::Rsquared = bigint<limbs>("3096616502983703923843567936837374451735540968419076528771170197431451843209");
            Fp::Rcubed = bigint<limbs>("14921786541159648185948152738563080959093619838510245177710943249661917737183");
            Fp::inv = 0xe4866389;
        }
        Fp::num_bits = 254;
        Fp::euler = bigint<limbs>("10944121435919637611123202872628637544348155578648911831344518947322613104291");
        Fp::s = 1;
        Fp::t = bigint<limbs>("10944121435919637611123202872628637544348155578648911831344518947322613104291");
        Fp::t_minus_1_over_2 = bigint<limbs>("5472060717959818805561601436314318772174077789324455915672259473661306552145");
        Fp::multiplicative_generator = Fp("3");
        Fp::root_of_unity = Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582");
        Fp::nqr = Fp("3");
        Fp::nqr_to_t = Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582");
    }

    /** Same as alt_bn128_Fq2 */
    void init_Fp2()
    {
        Fp2::euler = bigint<2*limbs>("239547588008311421220994022608339370399626158265550411218223901127035046843189118723920525909718935985594116157406550130918127817069793474323196511433944");
        Fp2::s = 4;
        Fp2::t = bigint<2*limbs>("29943448501038927652624252826042421299953269783193801402277987640879380855398639840490065738714866998199264519675818766364765977133724184290399563929243");
        Fp2::t_minus_1_over_2 = bigint<2*limbs>("14971724250519463826312126413021210649976634891596900701138993820439690427699319920245032869357433499099632259837909383182382988566862092145199781964621");
        Fp2::non_residue = Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582");
        Fp2::nqr = Fp2(Fp("2"),Fp("1"));
        Fp2::nqr_to_t = Fp2(Fp("5033503716262624267312492558379982687175200734934877598599011485707452665730"), Fp("314498342015008975724433667930697407966947188435857772134235984660852259084"));
        Fp2::Frobenius_coeffs_c1[0] = Fp("1");
        Fp2::Frobenius_coeffs_c1[1] = Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582");
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
    const FieldT one = FieldT::one();
    const FieldT zero = FieldT::zero();
    const FieldT two = one + one;

    /******************* Test standard field axioms and properties. *******************/

    FieldT x = random_element_exclude(zero);
    FieldT y = FieldT::random_element();
    FieldT z = FieldT::random_element();
    EXPECT_EQ(x + y, y + x); // commutative law of addition
    EXPECT_EQ((x + y) + z, x + (y + z)); // associative law of addition
    EXPECT_EQ(x + zero, x); // additive identity
    EXPECT_EQ(x + (-x), zero); // additive inverse
    EXPECT_EQ(x * y, y * x); // commutative law of multiplication
    EXPECT_EQ((x * y) * z, x * (y * z)); // associative law of multiplication
    EXPECT_EQ(x * one, x); // multiplicative identity
    EXPECT_EQ(x * x.inverse(), one); // multiplicative inverse
    EXPECT_EQ((x + y) * z, x * z + y * z); // distributive law

    EXPECT_EQ(-zero, zero);
    EXPECT_EQ(one.inverse(), one);
    EXPECT_EQ(-(-x), x);
    EXPECT_EQ(x.inverse().inverse(), x);
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

    // test square()
    x = FieldT::random_element();
    FieldT x2 = x * x;
    x.square();
    EXPECT_EQ(x, x2);

    // test squared()
    x.randomize();
    x2 = x * x;
    y = x;
    EXPECT_EQ(x.squared(), x2);
    EXPECT_EQ(x, y);
    x += x;
    EXPECT_EQ(x, y + y);

    // test +=
    z.randomize();
    y = x;
    x += z;
    EXPECT_EQ(x, y + z);
    x += zero;
    EXPECT_EQ(x, y + z);
    x += one;
    EXPECT_NE(x, y + z);
}

TEST(AllFieldsTest, AllFieldsApiTest)
{
    // test_field<AllFieldsTest::Fp>();
    // test_field<AllFieldsTest::Fp2>();

    alt_bn128_pp::init_public_params();
    test_field<alt_bn128_Fq>();
    test_field<alt_bn128_Fq2>();
}

}
