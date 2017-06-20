/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef BN_UTILS_HPP_
#define BN_UTILS_HPP_
#include <vector>
#include "third_party/ate-pairing/include/bn.h"

namespace libff {

template<typename FieldT>
void bn_batch_invert(std::vector<FieldT> &vec);

} // libff

#include "bn_utils.tcc"

#endif // BN_UTILS_HPP_
