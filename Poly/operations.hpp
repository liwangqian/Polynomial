#ifndef OPERATIONS_HPP_INCLUDED
#define OPERATIONS_HPP_INCLUDED

namespace poly
{

template<typename eT>
arma_inline
const bool
operator==(const Polynomial<eT>& x, const Polynomial<eT>& y)
{
    return x.equal(y);
}

template<typename T1>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value,
    const T1&
>::type
operator+(const T1& X)
{
    return X;
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T1, T2, glue_plus>
>::type
operator+(const T1& X, const T2& Y)
{
    return glue_mt<T1, T2, glue_plus>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_supported_elem_type<T2>::value,
    const glue_mt<T1, T2, glue_plus>
>::type
operator+(const T1& X, const T2 Y)
{
    return glue_mt<T1, T2, glue_plus>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T2>::value && is_supported_elem_type<T1>::value,
    const glue_mt<T2, T1, glue_plus>
>::type
operator+(const T1 X, const T2& Y)
{
    return glue_mt<T2, T1, glue_plus>(Y, X);
}

template<typename T1>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value,
    glue_mt<T1, typename T1::elem_type, glue_negate>
>::type
operator-(const T1& X)
{
    typedef typename T1::elem_type eT;
    return glue_mt<T1, eT, glue_negate>(X, eT(0));
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T1, T2, glue_subs>
>::type
operator-(const T1& X, const T2& Y)
{
    return glue_mt<T1, T2, glue_subs>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_supported_elem_type<T2>::value,
    const glue_mt<T1, T2, glue_subs>
>::type
operator-(const T1& X, const T2 Y)
{
    return glue_mt<T1, T2, glue_subs>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_supported_elem_type<T1>::value && is_PolyType<T2>::value,
    const glue_mt<glue_mt<T2, typename T2::elem_type, glue_negate>, T1, glue_plus>
>::type
operator-(const T1 X, const T2& Y)
{
    typedef typename T2::elem_type eT2;
    typedef glue_mt<T2, eT2, glue_negate> n_T;
    const  n_T n_Y(Y, eT2(0));
    return glue_mt<n_T, T1, glue_plus>(n_Y, X);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T1, T2, glue_mul>
>::type
operator*(const T1& X, const T2& Y)
{
    return glue_mt<T1, T2, glue_mul>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_supported_elem_type<T2>::value,
    const glue_mt<T1, T2, glue_mul>
>::type
operator*(const T1& X, const T2 Y)
{
    return glue_mt<T1, T2, glue_mul>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_supported_elem_type<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T2, T1, glue_mul>
>::type
operator*(const T1 X, const T2& Y)
{
    return glue_mt<T2, T1, glue_mul>(Y, X);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T1, T2, glue_div>
>::type
operator/(const T1& X, const T2& Y)
{
    return glue_mt<T1, T2, glue_div>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_supported_elem_type<T2>::value,
    const glue_mt<T1, T2, glue_div>
>::type
operator/(const T1& X, const T2 Y)
{
    return glue_mt<T1, T2, glue_div>(X, Y);
}

template<typename T1, typename T2>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value && is_PolyType<T2>::value,
    const glue_mt<T1, T2, glue_rem>
>::type
operator%(const T1& X, const T2& Y)
{
    return glue_mt<T1, T2, glue_rem>(X, Y);
}

template<typename T1>
arma_inline typename
std::enable_if
<
    is_PolyType<T1>::value,
    const glue_mt<T1, size_t, glue_pow>
>::type
operator^(const T1& X, const size_t Y)
{
    return glue_mt<T1, size_t, glue_pow>(X, Y);
}

///////////////////////////////////

/*
template<typename eT1, typename eT2>
arma_inline
const glue_st<Polynomial<eT1>, Polynomial<eT2>, glue_plus>
operator+(const Polynomial<eT1>& X, const Polynomial<eT2>& Y)
{
    return glue_st<Polynomial<eT1>, Polynomial<eT2>, glue_plus>(X, Y);
}


template<typename T1, typename U1, typename T2, typename U2>
arma_inline
const glue_st<glue_st<T1, U1, glue_plus>, glue_st<T2, U2, glue_plus>, glue_plus>
operator+(const glue_st<T1, U1, glue_plus>& X, const glue_st<T2, U2, glue_plus>& Y)
{
    return glue_st<glue_st<T1, U1, glue_plus>, glue_st<T2, U2, glue_plus>, glue_plus>(X, Y);
}
*/


}

#endif // OPERATIONS_HPP_INCLUDED
