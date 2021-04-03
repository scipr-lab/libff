#include "libff/algebra/field_utils/algorithms.hpp"

namespace libff {

template<mp_size_t m>
gf192& gf192::operator^=(const bigint<m> &pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t m>
gf192 gf192::operator^(const bigint<m> &pow) const
{
    return power<gf192>(*this, pow);
}

} // namespace libff
