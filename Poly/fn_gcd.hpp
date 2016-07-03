#ifndef FN_GCD_HPP_INCLUDED
#define FN_GCD_HPP_INCLUDED

namespace poly
{

//! Compute the greatest-common-divisor using Euclidean algorithm.
template<typename eT>
inline Polynomial<eT>
gcd(const Polynomial<eT>& x, const Polynomial<eT>& y)
{
    Polynomial<eT> r;
    Polynomial<eT> u;
    Polynomial<eT> v;

    const int xdeg = x.degree();
    const int ydeg = y.degree();

    if(xdeg < 0 || ydeg < 0)
        return std::move(u);

    if((xdeg == 0 && x[0] == 0) || (ydeg == 0 && y[0] == 0))
    {
        u.resize(1);
        u[0] = eT(0);
        return std::move(u);
    }

    const bool maxo = (xdeg >= ydeg);
    u = maxo ? x : y;
    v = maxo ? y : x;

    while(1)
    {
        r = u % v;
        u = v;
        v = r;
        if(v.degree() == 0 && v[0] == 0)
            break;
    }

    if(u.degree() == 0) //x and y are co-primes
        u[0] = eT(1);
    return std::move(u);
}


}

#endif // FN_GCD_HPP_INCLUDED
