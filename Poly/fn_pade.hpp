#ifndef FN_PADE_HPP_INCLUDED
#define FN_PADE_HPP_INCLUDED

namespace poly
{

//! Generate the hankel matrix.
template<typename eT>
void hankel(Mat<eT>& H, const Row<eT>& c, const Row<eT>& r)
{
    if(c.n_elem < 1)
        return;
    if(c.n_elem == 1)
    {
        H.resize(1);
        H(0) = c(0);
        return;
    }
    const Row<eT> P = join_rows(c, r.subvec(1, r.n_elem - 1)); // P = [c r(2:end)]
    H.resize(c.n_cols, r.n_cols);
    for(uword i = 0; i < H.n_rows; ++i)
        for(uword j = 0; j < H.n_cols; ++j)
            H(i,j) = P(i+j);//
}

//! Rotate the matrix 90 degree anti-clockwise
template<typename eT>
void rot90(Mat<eT>& out_A, const Mat<eT>& in_A)
{
    out_A.resize(in_A.n_cols, in_A.n_rows);
    const uword endc = in_A.n_cols - 1;
    for(uword i = 0; i < out_A.n_rows; ++i)
        out_A.row(i) = in_A.col(endc - i).st();
}

/**
 * \brief Compute the pade approximation of the given coefficients of taylor series.
 * @param c	        The coefficients of function's taylor series,
 *                  NOTE that the coefficients stored from low order to high order,
 *                  and the container c must support operator[] for element access.
 * @param c_size    The size of the coefficients, c_size > r+m.
 * @param N	        The coefficients of the numerator Polynomial.
 * @param D	        The coefficients of the denominator Polynomial.
 * @param r     	The order of the numerator.
 * @param m     	The order of the denominator, m >= r.
 * @return		    true if no error occurred, else false.
 */
template<typename C, typename eT>
bool pade(typename boost::call_traits<C>::param_type c, const uword c_size,
          Polynomial<eT>& N, Polynomial<eT>& D,
          const uword r, const uword m)
{
    if(r>m || r+m >= c_size)
    {
        arma_warn(true, "Warning: not proper parameters!");
        return false;
    }

    const uword w_size = m;
    const uword v_size = m;
	Col<eT> w(w_size); // w = -c(r+2:m+r+1)'
	Row<eT> v(v_size); // v = [c(r+1:-1:1) zeros(1,m-r-1)]
	const uword k = r + 1;
	for(uword i = k; i <= r + m; ++i)
        w(i-k) = -c(i);
    for(uword i = 0; i <= r; ++i)
        v(i) = c(r-i);
    for(uword i = r+1; i < m; ++i)
        v(i) = eT(0);
    Mat<eT> Wh, W;
    Row<eT> cw(m);
    for(uword i = 0; i < m; ++i)
        cw(i) = c(m+r-1-i);
    hankel(Wh, cw, v);
    rot90(W, Wh);
    Mat<eT> Vh, V;
    Row<eT> cv(r);
    if(r > 0)
    {
        for(uword i = 0; i < r; ++i)
            cv(i) = c(r-i-1);
        hankel(Vh, cv, cv);
        rot90(V, Vh);
    }
    Row<eT> x(1+m); //x=[1 (W/w)']
    Row<eT> y(1+r); //y=[1 x(2:r+1)*V'+c(2:r+1)]
    x(0) = 1; y(0) = 1;
    x.subvec(1, m) = solve(W, w).st();
    if(r > 0)
        y.subvec(1, r) = x.subvec(1,r)*strans(V);
    for(uword i = 1; i <= r; ++i)
        y(i) += c(i);
    N.resize(r+1);
    D.resize(m+1);
    for(uword i = 0; i <= m; ++i)
        D(i) = x(m-i) / x(m);
    for(uword i = 0; i <= r; ++i)
        N(i) = y(r-i) / x(m);
    return true;
}




}

#endif // FN_PADE_HPP_INCLUDED
