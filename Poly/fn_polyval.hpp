#ifndef FN_POLYVAL_HPP_INCLUDED
#define FN_POLYVAL_HPP_INCLUDED

namespace poly
{

template<typename eT>
arma_inline eT
polyval(const Polynomial<eT>& p, const eT& x)
{
    return p.eval(x);
}

template<typename eT, typename T>
arma_inline complex<typename common_type<eT, T>::type >
polyval(const Polynomial<eT>& p, const complex<T>& x)
{
    return p.eval(x);
}

//! Matrix/vector polynomial evaluation
template<typename eT, typename T, typename R = typename std::common_type<eT, T>::type >
arma_inline Mat<R>
polyvalm(const Polynomial<eT>& p, const Mat<T>& m)
{
    const uword p_size = p.m_size;
    Mat<R> res(m.n_rows, m.n_cols);
    if(m.is_vec())
    {
        for(uword i = 0; i < res.n_elem; ++i)
            res[i] = p.eval(m[i]);
    }
    else if(m.is_square())
    {
        Mat<R> I = eye(m.n_rows, m.n_cols);
        res = p[p_size-1] * I;
        for(int i = p.m_size - 1; i > 0; --i)
            res = res * m + p[i-1] * I;
    }
    else
    {
        arma_debug_warn(true, "error in polyvalm(p, m): the m must be vector or square matrix!");
        res.reset();
    }

    return std::move(res);
}



}

#endif // FN_POLYVAL_HPP_INCLUDED
