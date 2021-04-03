#include "libff/algebra/field_utils/algorithms.hpp"

namespace libff {

template<mp_size_t m>
gf128& gf128::operator^=(const bigint<m> &pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t m>
gf128 gf128::operator^(const bigint<m> &pow) const
{
    return power<gf128>(*this, pow);
}

} // namespace libff
