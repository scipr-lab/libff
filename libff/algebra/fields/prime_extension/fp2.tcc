/** @file
 *****************************************************************************
 Implementation of arithmetic in the finite field F[p^2].
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef FP2_TCC_
#define FP2_TCC_

#include <libff/algebra/field_utils/field_utils.hpp>

namespace libff {

using std::size_t;

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::zero()
{
    return Fp2_model<n, modulus>(my_Fp::zero(), my_Fp::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::one()
{
    return Fp2_model<n, modulus>(my_Fp::one(), my_Fp::zero());
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::random_element()
{
    Fp2_model<n, modulus> r;
    r.c0 = my_Fp::random_element();
    r.c1 = my_Fp::random_element();

    return r;
}

template<mp_size_t n, const bigint<n>& modulus>
void Fp2_model<n,modulus>::randomize()
{
    (*this) = Fp2_model<n, modulus>::random_element();
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp2_model<n,modulus>::operator==(const Fp2_model<n,modulus> &other) const
{
    return (this->c0 == other.c0 && this->c1 == other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp2_model<n,modulus>::operator!=(const Fp2_model<n,modulus> &other) const
{
    return !(operator==(other));
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator+(const Fp2_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    return Fp2_model<n,modulus>(this->c0 + other.c0,
                                this->c1 + other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator-(const Fp2_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->sub_cnt++;
#endif
    return Fp2_model<n,modulus>(this->c0 - other.c0,
                                this->c1 - other.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n, modulus> operator*(const Fp_model<n, modulus> &lhs, const Fp2_model<n, modulus> &rhs)
{
#ifdef PROFILE_OP_COUNTS
    rhs.mul_cnt++;
#endif
    return Fp2_model<n,modulus>(lhs*rhs.c0,
                                lhs*rhs.c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator*(const Fp2_model<n,modulus> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->mul_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba) */
    const my_Fp
        &A = other.c0, &B = other.c1,
        &a = this->c0, &b = this->c1;
    const my_Fp aA = a * A;
    const my_Fp bB = b * B;

    return Fp2_model<n,modulus>(aA + non_residue * bB,
                                (a + b)*(A+B) - aA - bB);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator-() const
{
    return Fp2_model<n,modulus>(-this->c0,
                                -this->c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator^(const unsigned long pow) const
{
    return power<Fp2_model<n, modulus>>(*this, pow);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp2_model<n,modulus> Fp2_model<n,modulus>::operator^(const bigint<m> &pow) const
{
    return power<Fp2_model<n, modulus>, m>(*this, pow);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::operator+=(const Fp2_model<n,modulus>& other)
{
    (*this) = *this + other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::operator-=(const Fp2_model<n,modulus>& other)
{
    (*this) = *this - other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::operator*=(const Fp2_model<n,modulus>& other)
{
    (*this) = *this * other;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::operator^=(const unsigned long pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
template<mp_size_t m>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::operator^=(const bigint<m> &pow)
{
    (*this) = *this ^ pow;
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::squared() const
{
    return squared_complex();
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::square()
{
    (*this) = squared();
    return (*this);
}


template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::squared_karatsuba() const
{
#ifdef PROFILE_OP_COUNTS
    this->sqr_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Karatsuba squaring) */
    const my_Fp &a = this->c0, &b = this->c1;
    const my_Fp asq = a.squared();
    const my_Fp bsq = b.squared();

    return Fp2_model<n,modulus>(asq + non_residue * bsq,
                                (a + b).squared() - asq - bsq);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::squared_complex() const
{
#ifdef PROFILE_OP_COUNTS
    this->sqr_cnt++;
#endif
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on Pairing-Friendly Fields.pdf; Section 3 (Complex squaring) */
    const my_Fp &a = this->c0, &b = this->c1;
    const my_Fp ab = a * b;

    return Fp2_model<n,modulus>((a + b) * (a + non_residue * b) - ab - non_residue * ab,
                                ab + ab);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::inverse() const
{
#ifdef PROFILE_OP_COUNTS
    this->inv_cnt++;
#endif
    const my_Fp &a = this->c0, &b = this->c1;

    /* From "High-Speed Software Implementation of the Optimal Ate Pairing over Barreto-Naehrig Curves"; Algorithm 8 */
    const my_Fp t0 = a.squared();
    const my_Fp t1 = b.squared();
    const my_Fp t2 = t0 - non_residue * t1;
    const my_Fp t3 = t2.inverse();
    const my_Fp c0 = a * t3;
    const my_Fp c1 = - (b * t3);

    return Fp2_model<n,modulus>(c0, c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus>& Fp2_model<n,modulus>::invert()
{
    (*this) = inverse();
    return (*this);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::Frobenius_map(unsigned long power) const
{
    return Fp2_model<n,modulus>(c0,
                                Frobenius_coeffs_c1[power % 2] * c1);
}

template<mp_size_t n, const bigint<n>& modulus>
Fp2_model<n,modulus> Fp2_model<n,modulus>::sqrt() const
{
    return tonelli_shanks_sqrt(*this);
}

template<mp_size_t n, const bigint<n>& modulus>
std::vector<uint64_t> Fp2_model<n,modulus>::to_words() const
{
    std::vector<uint64_t> words = c0.to_words();
    std::vector<uint64_t> words1 = c1.to_words();
    words.insert(words.end(), words1.begin(), words1.end());
    return words;
}

template<mp_size_t n, const bigint<n>& modulus>
bool Fp2_model<n,modulus>::from_words(std::vector<uint64_t> words)
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
std::ostream& operator<<(std::ostream &out, const Fp2_model<n, modulus> &el)
{
    out << el.c0 << OUTPUT_SEPARATOR << el.c1;
    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream &in, Fp2_model<n, modulus> &el)
{
    in >> el.c0 >> el.c1;
    return in;
}

template<mp_size_t n, const bigint<n>& modulus>
std::ostream& operator<<(std::ostream& out, const std::vector<Fp2_model<n, modulus> > &v)
{
    out << v.size() << "\n";
    for (const Fp2_model<n, modulus>& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

template<mp_size_t n, const bigint<n>& modulus>
std::istream& operator>>(std::istream& in, std::vector<Fp2_model<n, modulus> > &v)
{
    v.clear();

    size_t s;
    in >> s;

    char b;
    in.read(&b, 1);

    v.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        Fp2_model<n, modulus> el;
        in >> el;
        v.emplace_back(el);
    }

    return in;
}

} // libff
#endif // FP2_TCC_
