/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/bn128/bn128_g1.hpp>
#include <libff/algebra/curves/bn128/bn_utils.hpp>

namespace libff {

#ifdef PROFILE_OP_COUNTS
long long bn128_G1::add_cnt = 0;
long long bn128_G1::dbl_cnt = 0;
#endif

std::vector<size_t> bn128_G1::wnaf_window_table;
std::vector<size_t> bn128_G1::fixed_base_exp_window_table;
bn128_G1 bn128_G1::G1_zero = {};
bn128_G1 bn128_G1::G1_one = {};
bool bn128_G1::initialized = false;
bigint<bn128_G1::h_limbs> bn128_G1::h;

bn::Fp bn128_G1::sqrt(const bn::Fp &el)
{
    size_t v = bn128_Fq_s;
    bn::Fp z = bn128_Fq_nqr_to_t;
    bn::Fp w = mie::power(el, bn128_Fq_t_minus_1_over_2);
    bn::Fp x = el * w;
    bn::Fp b = x * w;

#if DEBUG
    // check if square with Euler's criterion
    bn::Fp check = b;
    for (size_t i = 0; i < v-1; ++i)
    {
        bn::Fp::square(check, check);
    }

    assert(check == bn::Fp(1));
#endif

    // compute square root with Tonelli--Shanks
    // (does not terminate if not a square!)

    while (b != bn::Fp(1))
    {
        size_t m = 0;
        bn::Fp b2m = b;
        while (b2m != bn::Fp(1))
        {
            // invariant: b2m = b^(2^m) after entering this loop
            bn::Fp::square(b2m, b2m);
            m += 1;
        }

        int j = v-m-1;
        w = z;
        while (j > 0)
        {
            bn::Fp::square(w, w);
            --j;
        } // w = z^2^(v-m-1)

        z = w * w;
        b = b * z;
        x = x * w;
        v = m;
    }

    return x;
}

bn128_G1::bn128_G1()
{
    if (bn128_G1::initialized)
    {
        this->X = G1_zero.X;
        this->Y = G1_zero.Y;
        this->Z = G1_zero.Z;
    }
}

void bn128_G1::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        bn128_G1 copy(*this);
        copy.to_affine_coordinates();
        std::cout << "(" << copy.X.toString(10) << " : " << copy.Y.toString(10) << " : " << copy.Z.toString(10) << ")\n";
    }
}

void bn128_G1::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        std::cout << "(" << X.toString(10) << " : " << Y.toString(10) << " : " << Z.toString(10) << ")\n";
    }
}

void bn128_G1::to_affine_coordinates()
{
    if (this->is_zero())
    {
        X = 0;
        Y = 1;
        Z = 0;
    }
    else
    {
        bn::Fp r;
        r = Z;
        r.inverse();
        bn::Fp::square(Z, r);
        X *= Z;
        r *= Z;
        Y *= r;
        Z = 1;
    }
}

void bn128_G1::to_special()
{
    this->to_affine_coordinates();
}

bool bn128_G1::is_special() const
{
    return (this->is_zero() || this->Z == 1);
}

bool bn128_G1::is_zero() const
{
    return Z.isZero();
}

bool bn128_G1::operator==(const bn128_G1 &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    bn::Fp Z1sq, Z2sq, lhs, rhs;
    bn::Fp::square(Z1sq, this->Z);
    bn::Fp::square(Z2sq, other.Z);
    bn::Fp::mul(lhs, Z2sq, this->X);
    bn::Fp::mul(rhs, Z1sq, other.X);

    if (lhs != rhs)
    {
        return false;
    }

    bn::Fp Z1cubed, Z2cubed;
    bn::Fp::mul(Z1cubed, Z1sq, this->Z);
    bn::Fp::mul(Z2cubed, Z2sq, other.Z);
    bn::Fp::mul(lhs, Z2cubed, this->Y);
    bn::Fp::mul(rhs, Z1cubed, other.Y);

    return (lhs == rhs);
}

bool bn128_G1::operator!=(const bn128_G1& other) const
{
    return !(operator==(other));
}

bn128_G1 bn128_G1::operator+(const bn128_G1 &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case, and then all the rest
    if (this->operator==(other))
    {
        return this->dbl();
    }
    else
    {
        return this->add(other);
    }
}

bn128_G1 bn128_G1::operator-() const
{
    bn128_G1 result(*this);
    bn::Fp::neg(result.Y, result.Y);
    return result;
}

bn128_G1 bn128_G1::operator-(const bn128_G1 &other) const
{
    return (*this) + (-other);
}

bn128_G1 bn128_G1::add(const bn128_G1 &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    bn::Fp this_coord[3], other_coord[3], result_coord[3];
    this->fill_coord(this_coord);
    other.fill_coord(other_coord);
    bn::ecop::ECAdd(result_coord, this_coord, other_coord);

    bn128_G1 result(result_coord);
    return result;
}

bn128_G1 bn128_G1::mixed_add(const bn128_G1 &other) const
{
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

#ifdef DEBUG
    assert(other.is_special());
#endif

    // check for doubling case

    // using Jacobian coordinates so:
    // (X1:Y1:Z1) = (X2:Y2:Z2)
    // iff
    // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
    // iff
    // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

    // we know that Z2 = 1

    bn::Fp Z1Z1;
    bn::Fp::square(Z1Z1, this->Z);
    const bn::Fp &U1 = this->X;
    bn::Fp U2;
    bn::Fp::mul(U2, other.X, Z1Z1);
    bn::Fp Z1_cubed;
    bn::Fp::mul(Z1_cubed, this->Z, Z1Z1);

    const bn::Fp &S1 = this->Y;
    bn::Fp S2;
    bn::Fp::mul(S2, other.Y, Z1_cubed); // S2 = Y2*Z1*Z1Z1

    if (U1 == U2 && S1 == S2)
    {
        // dbl case; nothing of above can be reused
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif

    bn128_G1 result;
    bn::Fp H, HH, I, J, r, V, tmp;
    // H = U2-X1
    bn::Fp::sub(H, U2, this->X);
    // HH = H^2
    bn::Fp::square(HH, H);
    // I = 4*HH
    bn::Fp::add(tmp, HH, HH);
    bn::Fp::add(I, tmp, tmp);
    // J = H*I
    bn::Fp::mul(J, H, I);
    // r = 2*(S2-Y1)
    bn::Fp::sub(tmp, S2, this->Y);
    bn::Fp::add(r, tmp, tmp);
    // V = X1*I
    bn::Fp::mul(V, this->X, I);
    // X3 = r^2-J-2*V
    bn::Fp::square(result.X, r);
    bn::Fp::sub(result.X, result.X, J);
    bn::Fp::sub(result.X, result.X, V);
    bn::Fp::sub(result.X, result.X, V);
    // Y3 = r*(V-X3)-2*Y1*J
    bn::Fp::sub(tmp, V, result.X);
    bn::Fp::mul(result.Y, r, tmp);
    bn::Fp::mul(tmp, this->Y, J);
    bn::Fp::sub(result.Y, result.Y, tmp);
    bn::Fp::sub(result.Y, result.Y, tmp);
    // Z3 = (Z1+H)^2-Z1Z1-HH
    bn::Fp::add(tmp, this->Z, H);
    bn::Fp::square(result.Z, tmp);
    bn::Fp::sub(result.Z, result.Z, Z1Z1);
    bn::Fp::sub(result.Z, result.Z, HH);
    return result;
}

bn128_G1 bn128_G1::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif

    bn::Fp this_coord[3], result_coord[3];
    this->fill_coord(this_coord);
    bn::ecop::ECDouble(result_coord, this_coord);

    bn128_G1 result(result_coord);
    return result;
}

bn128_G1 bn128_G1::mul_by_cofactor() const
{
    // Cofactor = 1
    return (*this);
}

bn128_G1 bn128_G1::zero()
{
    return G1_zero;
}

bn128_G1 bn128_G1::one()
{
    return G1_one;
}

bn128_G1 bn128_G1::random_element()
{
    return bn128_Fr::random_element().as_bigint() * G1_one;
}

std::ostream& operator<<(std::ostream &out, const bn128_G1 &g)
{
    bn128_G1 gcopy(g);
    gcopy.to_affine_coordinates();

    out << (gcopy.is_zero() ? '1' : '0') << OUTPUT_SEPARATOR;

#ifdef NO_PT_COMPRESSION
    /* no point compression case */
#ifndef BINARY_OUTPUT
    out << gcopy.X << OUTPUT_SEPARATOR << gcopy.Y;
#else
    out.write((char*) &gcopy.X, sizeof(gcopy.X));
    out.write((char*) &gcopy.Y, sizeof(gcopy.Y));
#endif

#else
    /* point compression case */
#ifndef BINARY_OUTPUT
    out << gcopy.X;
#else
    out.write((char*) &gcopy.X, sizeof(gcopy.X));
#endif
    out << OUTPUT_SEPARATOR << (((unsigned char*)&gcopy.Y)[0] & 1 ? '1' : '0');
#endif

    return out;
}

bool bn128_G1::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          y^2 = x^3 + b

          We are using Jacobian coordinates, so equation we need to check is actually

          (y/z^3)^2 = (x/z^2)^3 + b
          y^2 / z^6 = x^3 / z^6 + b
          y^2 = x^3 + b z^6
        */
        bn::Fp X2, Y2, Z2;
        bn::Fp::square(X2, this->X);
        bn::Fp::square(Y2, this->Y);
        bn::Fp::square(Z2, this->Z);

        bn::Fp X3, Z3, Z6;
        bn::Fp::mul(X3, X2, this->X);
        bn::Fp::mul(Z3, Z2, this->Z);
        bn::Fp::square(Z6, Z3);

        return (Y2 == X3 + bn128_coeff_b * Z6);
    }
}

std::istream& operator>>(std::istream &in, bn128_G1 &g)
{
    char is_zero;
    in.read((char*)&is_zero, 1); // this reads is_zero;
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

#ifdef NO_PT_COMPRESSION
    /* no point compression case */
#ifndef BINARY_OUTPUT
    in >> g.X;
    consume_OUTPUT_SEPARATOR(in);
    in >> g.Y;
#else
    in.read((char*) &g.X, sizeof(g.X));
    in.read((char*) &g.Y, sizeof(g.Y));
#endif

#else
    /* point compression case */
    bn::Fp tX;
#ifndef BINARY_OUTPUT
    in >> tX;
#else
    in.read((char*)&tX, sizeof(tX));
#endif
    consume_OUTPUT_SEPARATOR(in);
    unsigned char Y_lsb;
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + b)
    if (!is_zero)
    {
        g.X = tX;
        bn::Fp tX2, tY2;
        bn::Fp::square(tX2, tX);
        bn::Fp::mul(tY2, tX2, tX);
        bn::Fp::add(tY2, tY2, bn128_coeff_b);

        g.Y = bn128_G1::sqrt(tY2);
        if ((((unsigned char*)&g.Y)[0] & 1) != Y_lsb)
        {
            bn::Fp::neg(g.Y, g.Y);
        }
    }
#endif

    /* finalize */
    if (!is_zero)
    {
        g.Z = bn::Fp(1);
    }
    else
    {
        g = bn128_G1::zero();
    }

    return in;
}

std::ostream& operator<<(std::ostream& out, const std::vector<bn128_G1> &v)
{
    out << v.size() << "\n";
    for (const bn128_G1& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }
    return out;
}

std::istream& operator>>(std::istream& in, std::vector<bn128_G1> &v)
{
    v.clear();

    size_t s;
    in >> s;
    consume_newline(in);
    v.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        bn128_G1 g;
        in >> g;
        consume_OUTPUT_NEWLINE(in);
        v.emplace_back(g);
    }
    return in;
}

void bn128_G1::batch_to_special_all_non_zeros(std::vector<bn128_G1> &vec)
{
    std::vector<bn::Fp> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z);
    }
    bn_batch_invert<bn::Fp>(Z_vec);

    const bn::Fp one = 1;

    for (size_t i = 0; i < vec.size(); ++i)
    {
        bn::Fp Z2, Z3;
        bn::Fp::square(Z2, Z_vec[i]);
        bn::Fp::mul(Z3, Z2, Z_vec[i]);

        bn::Fp::mul(vec[i].X, vec[i].X, Z2);
        bn::Fp::mul(vec[i].Y, vec[i].Y, Z3);
        vec[i].Z = one;
    }
}

} // libff
