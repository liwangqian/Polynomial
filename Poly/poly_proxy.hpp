#ifndef POLY_PROXY_HPP_INCLUDED
#define POLY_PROXY_HPP_INCLUDED

namespace poly
{

template<typename T>
struct Proxy;

template<typename eT>
struct Proxy< Polynomial<eT> >
{
    typedef eT                                       elem_type;
    typedef typename get_pod_type<eT>::result        pod_type;
    typedef Polynomial<eT>                           stored_type;
    typedef const eT*                                ea_type;

    arma_aligned const Polynomial<eT>& Q;

    arma_inline explicit
    Proxy(const Polynomial<eT>& A)
        : Q(A)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return Q.m_size; }
    arma_inline elem_type operator[] (const uword i) const { return Q[i];           }
    arma_inline ea_type                     get_ea() const { return Q.m_data; }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
};

template<typename T1, typename T2, typename glue_type>
struct Proxy< glue_mt<T1,T2,glue_type> >
{
    typedef typename glue_mt<T1,T2,glue_type>::elem_type    elem_type;
    typedef typename get_pod_type<elem_type>::result        pod_type;
    typedef Polynomial<elem_type>                           stored_type;
    typedef const elem_type*                                ea_type;

    arma_aligned const Polynomial<elem_type> Q;//direct evaluation

    arma_inline explicit
    Proxy(const glue_mt<T1,T2,glue_type>& A)
        : Q(A)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return Q.m_size; }
    arma_inline elem_type operator[] (const uword i) const { return Q[i];     }
    arma_inline ea_type                     get_ea() const { return Q.m_data; }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return false; }
};




}

#endif // POLY_PROXY_HPP_INCLUDED
