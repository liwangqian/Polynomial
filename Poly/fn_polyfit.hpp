#ifndef FN_POLYFIT_HPP_INCLUDED
#define FN_POLYFIT_HPP_INCLUDED

namespace poly
{


namespace detail
{

//! phi = [1 x x^2 ... x^n]
template<typename eT>
arma_inline void
gen_phi(Col<eT>& phi, const eT x)
{
    phi[0] = eT(1);
    for(uword i = 1; i < phi.n_rows; ++i)
        phi[i] = phi[i-1] * x;
}


}

//! Polynomial data fit using least-square method
template<typename eT>
inline Polynomial<eT>
polyfit(const eT* x, const eT* y, const uword length, const uword order)
{
    const uword size = order + 1;
    Col<eT> Q(size);   //parameter to be estimated.
    Q.fill(eT(0));        //initialized with zeros.
    Col<eT> phi(size); //
    Mat<eT> P = eye(size, size) * 1e6;
    Col<eT> K;
    Mat<eT> I = eye(size, size);

    for(uword i = 0; i < length; ++i)
    {
        detail::gen_phi(phi, x[i]);
        K = P * phi / (1 + as_scalar(phi.st() * P * phi)); //K(t) = P(t-1)*phi(t)*inv(I+phi'*P(t-1)*phi).
        P = (I - K * phi.st()) * P; //P(t) = (I-K(t)*phi')*P(t-1).
        Q = Q + K * (y[i] - as_scalar(phi.st() * Q)); //Q(t) = Q(t-1) + K(t)*(y(t)-phi'Q(t-1)).
    }
    return std::move(Polynomial<eT>(Q));
}

}

#endif // FN_POLYFIT_HPP_INCLUDED
