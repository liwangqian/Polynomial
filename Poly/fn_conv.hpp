#ifndef FN_CONV_HPP_INCLUDED
#define FN_CONV_HPP_INCLUDED

namespace poly
{

template
<
    typename eT1, typename T1,
    typename eT2, typename T2,
    typename eT = typename common_type<eT1, eT2>::type
>
arma_inline Polynomial<eT>
conv(const PolyBase<eT1, T1>& p1, const PolyBase<eT2, T2>& p2)
{
    typedef typename common_type<eT1, eT2>::type elem_type;
    const Proxy<T1> px(p1.get_ref());
    const Proxy<T2> py(p2.get_ref());
    Polynomial<elem_type> r = px.Q * py.Q;
    return std::move(r);
}


}

#endif // FN_CONV_HPP_INCLUDED
