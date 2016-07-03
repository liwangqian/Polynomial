#ifndef FN_LCM_HPP_INCLUDED
#define FN_LCM_HPP_INCLUDED

namespace poly
{

template<typename eT>
arma_inline Polynomial<eT>
lcm(const Polynomial<eT>& x, const Polynomial<eT>& y)
{
    Polynomial<eT> r;
    if(!(x.is_empty() || y.is_empty()))
        r = x*y/gcd(x,y);
    return std::move(r);
}

}

#endif // FN_LCM_HPP_INCLUDED
