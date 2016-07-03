#ifndef FN_POLY_HPP_INCLUDED
#define FN_POLY_HPP_INCLUDED

namespace poly
{

//! computes the characteristic polynomial of the square matrix when m is square matrix,
//! when m is a vector, computes the polynomial which roots is given by the vector.
//! the program uses the Leverrier Feddeev algorithm.
template<typename eT>
inline
Polynomial<eT> poly(const Mat<eT>& m)
{
    const uword nr = m.n_rows;
    const uword nc = m.n_cols;
    Polynomial<eT> P;

    if(nr == nc) // m is square matrix
    {
        Mat<eT> I = eye(nr, nc);
        Mat<eT> R = I;
        P.resize(1+nc);
        P[nc] = eT(1);

        for(uword i = 1; i <= nc; i++)
        {
            R = m * R;
            P[nc-i] = -1.0/i * trace(R);
            R += P[nc-i] * I;
        }
    }
    else if(nr == 1 || nc == 1) // m is vector
    {
        const uword n = m.n_elem;
        P.resize(n+1);
        P[n] = eT(1);

        for(uword i = 0; i < n; i++)
        {
            for(int j = i; j >= 0; j--)
            {
                P[n-(j+1)] -= m[i] * P[n-j];
            }
        }
    }
    else
        arma_debug_warn(true, "Error in calling poly(A): Argument must be a vector or a square matrix!");

    return std::move(P);
}

}

#endif // FN_POLY_HPP_INCLUDED
