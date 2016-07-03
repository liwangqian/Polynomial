#ifndef FN_POW_HPP_INCLUDED
#define FN_POW_HPP_INCLUDED

namespace poly
{

template<typename eT, typename T>
arma_inline Polynomial<eT>
pow(const PolyBase<eT, T>& expr, const uword n = 2)
{
    Polynomial<eT>  R(expr.get_ref());
    R ^= n;
    return std::move(R);
}


}

#endif // FN_POW_HPP_INCLUDED
