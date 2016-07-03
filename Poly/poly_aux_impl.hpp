#ifndef POLY_AUX_IMPL_HPP_INCLUDED
#define POLY_AUX_IMPL_HPP_INCLUDED

namespace poly
{

template<typename T1, typename T2>
arma_hot arma_inline void
aux::copy(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem)
{
    if(n_elem == 1)
    {
        dest[0] = src[0];
        return;
    }
    uword i, j;
    for(i = 0, j = 1; j < n_elem; i+=2, j+=2)
    {
        dest[i] = src[i];
        dest[j] = src[j];
    }
    if(i < n_elem)
        dest[i] = src[i];
}

template<typename T1, typename T2>
arma_hot arma_inline void
aux::copy_backward(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem)
{
    const uword endn = n_elem - 1;
    if(n_elem == 1)
    {
        dest[0] = src[0];
        return;
    }
    uword i, j;
    for(i = 0, j = 1; j < n_elem; i+=2, j+=2)
    {
        dest[i] = src[endn - i];
        dest[j] = src[endn - j];
    }
    if(i < n_elem)
        dest[i] = src[endn - i];

}

template<typename eT>
arma_hot arma_inline void
aux::inplace_plus(eT* dest, const eT* src, const uword n_elem)
{
    arrayops::inplace_plus(dest, src, n_elem);
}

template<typename eT>
arma_hot arma_inline void
aux::inplace_minus(eT* dest, const eT* src, const uword n_elem)
{
    arrayops::inplace_minus(dest, src, n_elem);
}

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline
R
aux::conv_impl_p1(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size)
{
    R acc = R(0);
    for(uword x_i = 0, h_i = ii; x_i <= ii; ++x_i, --h_i)
    {
        acc += X[x_i] * H[h_i];
    }
    return acc;
}

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline
R
aux::conv_impl_p2(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size)
{
    R acc = R(0);
    uword h_i = h_size - 1;
    for(uword x_i = ii - h_i; x_i <= ii; ++x_i, --h_i)
    {
        acc += X[x_i] * H[h_i];
    }
    return acc;
}

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline
R
aux::conv_impl_p3(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size)
{
    R acc = R(0);
    uword h_i = h_size - 1;
    for(uword x_i = ii - h_size + 1; x_i < x_size; ++x_i, --h_i)
    {
        acc += X[x_i] * H[h_i];
    }
    return acc;
}

template<typename eT>
arma_hot inline void
aux::conv(eT* r_mem, const eT* A, const uword A_size, const eT* B, const uword B_size)
{
    const uword h_size = std::min(A_size, B_size);
    const uword x_size = std::max(A_size, B_size);
    const uword r_size = h_size + x_size - 1;

    const eT*   h_mem = (A_size <= B_size) ? A : B;
    const eT*   x_mem = (A_size <= B_size) ? B : A;

    for(uword r_i = 0; r_i < h_size - 1; ++r_i)
    {
        r_mem[r_i] = conv_impl_p1<const eT*, const eT*, eT>(r_i, h_mem, h_size, x_mem, x_size);
    }

    for(uword r_i = h_size-1; r_i < r_size - (h_size-1); ++r_i)
    {
        r_mem[r_i] = conv_impl_p2<const eT*, const eT*, eT>(r_i, h_mem, h_size, x_mem, x_size);
    }

    for(uword r_i = r_size - (h_size-1); r_i < r_size; ++r_i)
    {
        r_mem[r_i] = conv_impl_p3<const eT*, const eT*, eT>(r_i, h_mem, h_size, x_mem, x_size);
    }

}

template<typename eT>
arma_hot arma_inline void
aux::inplace_mul(eT* dest, const eT val, const uword n_elem)
{
    arrayops::inplace_mul(dest, val, n_elem);
}


template<typename T1, typename T2>
arma_hot inline void
aux::inplace_plus(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem)
{
    uword i,j;
    for(i= 0, j = 1; j < n_elem; i += 2, j += 2)
    {
        dest[i] += src[i];
        dest[j] += src[j];
    }
    if(i < n_elem)
        dest[i] += src[i];
}

template<typename T1, typename T2>
arma_hot inline void
aux::inplace_minus(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem)
{
    uword i,j;
    for(i= 0, j = 1; j < n_elem; i += 2, j += 2)
    {
        dest[i] -= src[i];
        dest[j] -= src[j];
    }
    if(i < n_elem)
        dest[i] -= src[i];
}


}

#endif // POLY_AUX_IMPL_HPP_INCLUDED
