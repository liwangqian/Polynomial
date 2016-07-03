#ifndef FN_DECONV_HPP_INCLUDED
#define FN_DECONV_HPP_INCLUDED

namespace poly
{


template
<
    typename eT1, typename T1,
    typename eT2, typename T2,
    typename eT = typename common_type<eT1, eT2>::type
>
arma_inline Polynomial<eT>
deconv(const PolyBase<eT1, T1>& U, const PolyBase<eT2, T2>& V, Polynomial<eT>& R)
{
                    R = U.get_ref(); //reminder
    const Proxy<T2> B_(V.get_ref());
    const int       A_size = R.size();

    Polynomial<eT> Q(A_size); //quoter
    Q.fill(eT(0));

    int k, j, n = A_size - 1, nv = B_.Q.degree();

    arma_debug_check(nv < 0, "error in deconv(): poly divide/mod by empty polynomial!");
    arma_warn(nv == 0 && B_[nv] == eT(0), "warning in deconv(): poly divide/mod by empty polynomial!");

    for(k = n-nv; k >= 0; --k)
    {
        Q[k] = R[nv+k] / B_[nv];
        for(j = nv+k-1; j>= k; --j)
            R[j] -= Q[k] * B_[j-k];
    }
    for(j = nv; j <= n; ++j) R[j] = 0.0;

    return std::move(Q);
}

}

#endif // FN_DECONV_HPP_INCLUDED
