#ifndef FN_MAT2POLY_HPP_INCLUDED
#define FN_MAT2POLY_HPP_INCLUDED

namespace poly
{

//! The parameter n specify the slices of the matrix m,
//! that's to say matrix m is constructed by n same size matrix: size = [rows(m), cols(m)/n].
template<typename eT>
inline void
mat2poly(const Mat<eT>& m, field<Polynomial<eT> >& pa, const uword n)
{
    if(n == 0) return;
    if((m.n_cols % n) != 0)
    {
        arma_warn(true, "warning: mat2poly(m,p,n) -> the matrix m is not proper constructed!");
        return;
    }

    const uword c = m.n_cols / n; //number of columns of each matrix slice.
    const uword r = m.n_rows;

    pa.set_size(r, c);  //size of poly array.

    Polynomial<eT> p(n);
    uword k;
    for(uword h = 0; h < r; ++h)
    {
        for(uword i = 0; i < c; ++i)
        {
            k = 0;
            for(uword j = 0; j < n; ++j)
            {
                p[k++] = m(h, i+j*c);
            }
            pa(h,i) = p;
        }
    }
}

template<typename eT>
arma_inline void
mat2poly(const Mat<eT>& m, Polynomial<eT>& pa)
{
    pa = poly(m);
}


}

#endif // FN_MAT2POLY_HPP_INCLUDED
