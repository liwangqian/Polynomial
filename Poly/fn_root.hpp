#ifndef FN_ROOT_HPP_INCLUDED
#define FN_ROOT_HPP_INCLUDED

namespace poly
{

template<typename eT, typename T>
inline
std::vector<std::complex<double> >
root(const PolyBase<eT, T>& expr)
{
    typedef std::complex<double>   root_type;
    typedef std::vector<root_type> roots_holder;
    const Expression<T> EX(expr.get_ref());
    const int e_size = EX.size();
    roots_holder result;
    if(e_size <= 1)
        return std::move(result);

    double* coefs = new double[e_size];
    aux::copy_backward<double*, Expression<T> >(coefs, EX, e_size);

    const uword endn = e_size - 1;
    int pb = 0, pe = endn;
    while(coefs[pb] == 0.0 && pb < e_size) ++pb; //strip the begin zeros
    while(coefs[pe] == 0.0 && pe > pb    ) --pe; //strip the end zeros

    result.resize(endn - pe); //set the zero roots

    if(pe != pb) //compute the nonzero roots
    {
        const uword nz_roots = pe - pb; //number of nonzero roots
        double* zeror = new double[nz_roots];
        double* zeroi = new double[nz_roots];
        int     info[nz_roots + 1];
        int nz = rpoly_impl(coefs+pb, nz_roots, zeror, zeroi, info);
        for(int i = 0; i < nz; ++i)
            result.push_back(root_type(zeror[i], zeroi[i]));

        delete zeroi;
        delete zeror;
    }
    delete coefs;

    return std::move(result);
}



}

#endif // FN_ROOT_HPP_INCLUDED
