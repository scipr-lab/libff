/** @file
 *****************************************************************************

 Implementation of interfaces for the (extension) field Fp4.

 See fp4.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP4_TCC_
#define FP4_TCC_

#include <libff/algebra/field_utils/field_utils.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>

namespace libff {

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> Fp4_model<n, modulus>::mul_by_non_residue(const Fp2_model<n, modulus> &elt)
{
    return Fp2_model<n, modulus>(non_residue * elt.c1, elt.c0);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp4_model<n, modulus>::zero()
{
    return Fp4_model<n,modulus>(my_Fp2::zero(),
                                my_Fp2::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> Fp4_model<n, modulus>::one()
{
    return Fp4_model<n,modulus>(my_Fp2::one(),
                                my_Fp2::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::random_element()
{
    Fp4_model<n, modulus> r;
    r.c0 = my_Fp2::random_element();
    r.c1 = my_Fp2::random_element();

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp4_model<n,modulus>::randomize()
{
    (*this) = Fp4_model<n, modulus>::random_element();
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp4_model<n,modulus>::operator==(const Fp4_model<n,modulus> &other) const
{
    return (this->c0 == other.c0 && this->c1 == other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp4_model<n,modulus>::operator!=(const Fp4_model<n,modulus> &other) const
{
    return !(operator==(other));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::operator+(const Fp4_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    return Fp4_model<n,modulus>(this->c0 + other.c0,
                                this->c1 + other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::operator-(const Fp4_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->sub_cnt++;
#endif
    return Fp4_model<n,modulus>(this->c0 - other.c0,
                                this->c1 - other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp4_model<n, modulus> &rhs)
{
#ifdef PROFILE_OP_COUNTS
    rhs.mul_cnt++;
#endif
    return Fp4_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n, modulus> operator*(const Fp2_model<n, modulus> &lhs, const Fp4_model<n, modulus> &rhs)
{
#ifdef PROFILE_OP_COUNTS
    rhs.mul_cnt++;
#endif
    return Fp4_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::operator*(const Fp4_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->mul_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */

    const my_Fp2 &B = other.c1, &A = other.c0,
        &b = this->c1, &a = this->c0;
    const my_Fp2 aA = a*A;
    const my_Fp2 bB = b*B;

    const my_Fp2 beta_bB = Fp4_model<n,modulus>::mul_by_non_residue(bB);
    return Fp4_model<n,modulus>(aA + beta_bB,
                                (a+b)*(A+B) - aA  - bB);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::mul_by_023(const Fp4_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->mul_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */
    assert(other.c0.c1.is_zero());

    const my_Fp2 &B = other.c1, &A = other.c0,
        &b = this->c1, &a = this->c0;
    const my_Fp2 aA = my_Fp2(a.c0 * A.c0, a.c1 * A.c0);
    const my_Fp2 bB = b*B;

    const my_Fp2 beta_bB = Fp4_model<n,modulus>::mul_by_non_residue(bB);
    return Fp4_model<n,modulus>(aA + beta_bB,
                                (a+b)*(A+B) - aA  - bB);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::operator-() const
{
    return Fp4_model<n,modulus>(-this->c0,
                                -this->c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::operator^(const unsigned long pow) const
{
    return power<Fp4_model<n, modulus> >(*this, pow);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp4_model<n, modulus> Fp4_model<n,modulus>::operator^(const bigint<m> &exponent) const
{
    return power<Fp4_model<n, modulus> >(*this, exponent);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m, const bigint<m>& modulus_p>
Fp4_model<n, modulus> Fp4_model<n,modulus>::operator^(const Fp_model<m, modulus_p> &exponent) const
{
    return (*this)^(exponent.as_bigint());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::operator+=(const Fp4_model<n,modulus>& other)
{
    (*this) = *this + other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::operator-=(const Fp4_model<n,modulus>& other)
{
    (*this) = *this - other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::operator*=(const Fp4_model<n,modulus>& other)
{
    (*this) = *this * other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::operator^=(const unsigned long pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::operator^=(const bigint<m> &pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::squared() const
{
#ifdef PROFILE_OP_COUNTS
    this->sqr_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Complex) */

    const my_Fp2 &b = this->c1, &a = this->c0;
    const my_Fp2 ab = a * b;

    return Fp4_model<n,modulus>((a+b)*(a+Fp4_model<n,modulus>::mul_by_non_residue(b))-ab-Fp4_model<n,modulus>::mul_by_non_residue(ab),
                                ab + ab);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::square()
{
    (*this) = squared();
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::inverse() const
{
#ifdef PROFILE_OP_COUNTS
    this->inv_cnt++;
#endif
    /* From "High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves"; Algorithm 8 */
    const my_Fp2 &b = this->c1, &a = this->c0;
    const my_Fp2 t1 = b.squared();
    const my_Fp2 t0 = a.squared() - Fp4_model<n,modulus>::mul_by_non_residue(t1);
    const my_Fp2 new_t1 = t0.inverse();

    return Fp4_model<n,modulus>(a * new_t1, - (b * new_t1));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus>& Fp4_model<n,modulus>::invert()
{
    (*this) = inverse();
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    return Fp4_model<n,modulus>(c0.Frobenius_map(power),
                                Frobenius_coeffs_c1[power % 4] * c1.Frobenius_map(power));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::unitary_inverse() const
{
    return Fp4_model<n,modulus>(this->c0,
                                -this->c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::cyclotomic_squared() const
{
    const my_Fp2 A = this->c1.squared();
    const my_Fp2 B = this->c1 + this->c0;
    const my_Fp2 C = B.squared() - A;
    const my_Fp2 D = Fp4_model<n,modulus>::mul_by_non_residue(A); // Fp2(A.c1 * non_residue, A.c0)
    const my_Fp2 E = C - D;
    const my_Fp2 F = D + D + my_Fp2::one();
    const my_Fp2 G = E - my_Fp2::one();

    return Fp4_model<n,modulus>(F, G);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp4_model<n, modulus> Fp4_model<n,modulus>::cyclotomic_exp(const bigint<m> &exponent) const
{
    Fp4_model<n,modulus> res = Fp4_model<n,modulus>::one();
    Fp4_model<n,modulus> this_inverse = this->unitary_inverse();

    bool found_nonzero = false;
    std::vector<long> NAF = find_wnaf(1, exponent);

    for (long i = static_cast<long>(NAF.size() - 1); i >= 0; --i)
    {
        if (found_nonzero)
        {
            res = res.cyclotomic_squared();
        }

        if (NAF[i] != 0)
        {
            found_nonzero = true;

            if (NAF[i] > 0)
            {
                res = res * (*this);
            }
            else
            {
                res = res * this_inverse;
            }
        }
    }

    return res;
}

template<mp_size_t n, const bigint<n>& modulus>
Fp4_model<n,modulus> Fp4_model<n,modulus>::sqrt() const
{
    return tonelli_shanks_sqrt(*this);
}

template<mp_size_t n, const bigint<n>& modulus>
std::vector<uint64_t> Fp4_model<n,modulus>::to_words() const
{
    std::vector<uint64_t> words = c0.to_words();
    std::vector<uint64_t> words1 = c1.to_words();
    words.insert(words.end(), words1.begin(), words1.end());
    return words;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp4_model<n,modulus>::from_words(std::vector<uint64_t> words)
{
    std::vector<uint64_t>::const_iterator vec_start = words.begin();
    std::vector<uint64_t>::const_iterator vec_center = words.begin() + words.size() / 2;
    std::vector<uint64_t>::const_iterator vec_end = words.end();
    std::vector<uint64_t> words0(vec_start, vec_center);
    std::vector<uint64_t> words1(vec_center, vec_end);
    // Fp_model's from_words() takes care of asserts about vector length.
    return c0.from_words(words0) && c1.from_words(words1);
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream &out, const Fp4_model<n, modulus> &el)
{
    out << el.c0 << OUTPUT_SEPARATOR << el.c1;
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp4_model<n, modulus> &el)
{
    in >> el.c0 >> el.c1;
    return in;
}

} // libff

#endif // FP4_TCC_
