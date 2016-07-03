#ifndef FN_POLY2MAT_HPP_INCLUDED
#define FN_POLY2MAT_HPP_INCLUDED

namespace poly
{

template<typename eT>
inline void
poly2mat(const field<Polynomial<eT> >& pa, Mat<eT>& m)
{
    const uword r = pa.n_rows;
    const uword c = pa.n_cols;
    uword n = 0;
    uword k;

    for(uword i = 0; i < r; ++i)
    {
        for(uword j = 0; j < c; ++j)
        {
            k = pa(i,j).degree();
            if(k > n)   n = k;
        }
    }

    n++;
    m.set_size(r, n*c);

    for(uword i = 0; i < r; ++i)
    {
        for(uword j = 0; j < c; ++j)
        {
            for(uword k = 0; k < n; ++k)
                m(i, k*c+j) = pa(i,j)[k];
        }
    }
}

template<typename eT>
inline void
poly2mat(const Polynomial<eT>& p, Mat<eT>& m)
{
    const uword n = p.degree();
    m.set_size(n,n);

    eT head = p[n];

    for(uword i = 0; i < n-1; ++i)
    {
        m(i, i+1) = 1.0;
    }
    for(uword i = 0; i < n; ++i)
    {
        m(n-1, i) = -p[i]/head;
    }
}

}

#endif // FN_POLY2MAT_HPP_INCLUDED
