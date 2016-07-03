#ifndef FN_POLYDER_HPP_INCLUDED
#define FN_POLYDER_HPP_INCLUDED

namespace poly
{

//! Difference of the polynomial
template<typename eT, typename T>
arma_inline Polynomial<eT>
polyder(const PolyBase<eT, T>& p, const uword n = 1)
{
    const Expression<T> ep(p.get_ref());
    int deg = ep.size() - 1;
    while(deg >= 0 && ep[deg] == eT(0)) --deg;
    if(deg < 0) return Polynomial<eT>();

    Polynomial<eT> r(deg-n+1);
    eT in = 1;
    for(uword i = 2; i <= n; ++i) in *= i;
    for(uword i = 0; i < r.m_size; ++i)
    {
        r[i] = ep[i+n] * in;
        in /= (i+1);
        in *= (i+n+1);
    }

    return std::move(r);
}

}

#endif // FN_POLYDER_HPP_INCLUDED
