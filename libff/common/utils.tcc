/** @file
 *****************************************************************************
 Implementation of templatized utility functions
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef UTILS_TCC_
#define UTILS_TCC_

namespace libff {

template<typename T>
size_t ceil_size_in_bits(const std::vector<T> &v)
{
    return v.size() * T::ceil_size_in_bits();
}

template<typename T>
T random_element_non_zero_one()
{
    T x = T::random_element();
    while (x == T::zero() || x == T::one())
        x = T::random_element();
    return x;
}

template<typename FieldT>
FieldT random_element_non_zero()
{
    FieldT x = FieldT::random_element();
    while (x.is_zero())
        x = FieldT::random_element();
    return x;
}

} // libff

#endif // UTILS_TCC_
