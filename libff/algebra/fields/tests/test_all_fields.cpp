/**
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/
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

using namespace libff;

template<typename FieldT>
void test_field()
{
    // test square()
    FieldT x = FieldT::random_element();
    FieldT x2 = x * x;
    x.square();
    assert(x == x2);

    // test squared()
    x.randomize();
    x2 = x * x;
    FieldT y = x;
    assert(x.squared() == x2);
    assert(x == y);
    x += x;
    assert(x == y + y);

    // test +=
    FieldT z = FieldT::random_element();
    y = x;
    x += z;
    assert(x == y + z);
    x += FieldT::zero();
    assert(x == y + z);
    x += FieldT::one();
    assert(x != y + z);
}

int main(void)
{
    alt_bn128_pp::init_public_params();
    test_field<alt_bn128_Fq>();
    test_field<alt_bn128_Fq2>();

    std::cout << "All tests passed!\n";
}
