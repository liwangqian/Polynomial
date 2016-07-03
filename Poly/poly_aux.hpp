#ifndef POLY_AUX_HPP_INCLUDED
#define POLY_AUX_HPP_INCLUDED

namespace poly
{

struct aux
{

template<typename T1, typename T2>
arma_hot arma_inline static void
copy(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem);

template<typename T1, typename T2>
arma_hot arma_inline static void
copy_backward(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem);

template<typename eT>
arma_hot arma_inline static
void
inplace_plus(eT* dest, const eT* src, const uword n_elem);

template<typename eT>
arma_hot arma_inline static
void
inplace_minus(eT* dest, const eT* src, const uword n_elem);

template<typename eT>
arma_hot arma_inline static void
inplace_mul(eT* dest, const eT val, const uword n_elem);

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline static
R
conv_impl_p1(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size);

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline static
R
conv_impl_p2(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size);

template
<
    typename T1,
    typename T2,
    typename R
>
arma_hot arma_inline static
R
conv_impl_p3(const uword ii, typename boost::call_traits<T1>::param_type H, const uword h_size, typename boost::call_traits<T2>::param_type X, const uword x_size);

template<typename eT>
arma_hot arma_inline static
void
conv(eT* out_mem, const eT* X, const uword x_size, const eT* Y, const uword y_size);

template<typename T1, typename T2>
arma_hot inline static
void
inplace_plus(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem);

template<typename T1, typename T2>
arma_hot inline static
void
inplace_minus(typename boost::call_traits<T1>::param_type dest, typename boost::call_traits<T2>::param_type src, const uword n_elem);


};

}

#endif // POLY_AUX_HPP_INCLUDED
