#include "libff/algebra/field_utils/algorithms.hpp"

namespace libff {

template<mp_size_t m>
gf64& gf64::operator^=(const bigint<m> &pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t m>
gf64 gf64::operator^(const bigint<m> &pow) const
{
    return power<gf64>(*this, pow);
}

} // namespace libff
