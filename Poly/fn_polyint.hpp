#ifndef FN_POLYINT_HPP_INCLUDED
#define FN_POLYINT_HPP_INCLUDED

namespace poly
{

//! Integration of polynomial
template<typename eT, typename T>
arma_inline Polynomial<eT>
polyint(const PolyBase<eT, T>& p, const uword n = 1)
{
    const Expression<T> ep(p.get_ref());
    int deg = ep.size() - 1;
    while(deg >= 0 && ep[deg] == eT(0)) --deg;
    if(deg < 0) return Polynomial<eT>();

    Polynomial<eT> r(deg+1+n);
    for(uword i = 0; i < n; ++i)  r[i] = eT(0);
    eT in = 1;
    for(uword i = 2; i <= n; ++i)  in *= i;
    for(uword i = n; i < r.m_size; ++i)
    {
        r[i] = ep[i-n] / in;
        in /= (i+1-n);
        in *= (i+1);
    }

    return std::move(r);
}


}

#endif // FN_POLYINT_HPP_INCLUDED
