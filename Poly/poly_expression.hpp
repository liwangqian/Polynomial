#ifndef POLY_EXPRESSION_HPP_INCLUDED
#define POLY_EXPRESSION_HPP_INCLUDED

namespace poly
{

//Default type to Scalar
template<typename T>
struct Expression
{
    typedef T                                        elem_type;
    typedef T                                        pod_type;
    typedef T                                        stored_type;

    arma_aligned const T Q;

    arma_inline explicit
    Expression(const T A)
        : Q(A)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return 1;                         }
    arma_inline elem_type operator[] (const uword i) const { return i == 0 ? Q : elem_type(0); }

    template<typename T1>
    arma_inline bool is_alias(const Polynomial<T1>& X)                   const { return false; }
};

template<typename eT>
struct Expression< Polynomial<eT> >
{
    typedef eT                                       elem_type;
    typedef typename get_pod_type<eT>::result        pod_type;
    typedef Polynomial<eT>                           stored_type;

    arma_aligned const Polynomial<eT>& Q;

    arma_inline explicit
    Expression(const Polynomial<eT>& A)
        : Q(A)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return Q.m_size; }
    arma_inline elem_type operator[] (const uword i) const { return Q[i];     }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_plus> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;


    arma_aligned const expr1_type Q;
    arma_aligned const expr2_type R;

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_plus>& A)
        : Q(A.A), R(A.B)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return std::max(Q.size(), R.size()); }
    arma_inline elem_type operator[] (const uword i) const { return Q[i] + R[i];                  }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return Q.is_alias(X) || R.is_alias(X); }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_subs> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned const expr1_type Q;
    arma_aligned const expr2_type R;

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_subs>& A)
        : Q(A.A), R(A.B)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return std::max(Q.size(), R.size()); }
    arma_inline elem_type operator[] (const uword i) const { return Q[i] - R[i];                  }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return Q.is_alias(X) || R.is_alias(X); }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_mul> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned const Proxy<T1>        Q;
    arma_aligned Polynomial<elem_type>  R;

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_mul>& A)
        : Q(A.A), R(A.B)
    {
        arma_extra_debug_sigprint();
        R *= Q.Q;
    }

    arma_inline uword                    size() const { return Q.size() + R.size() - 1; }
    inline elem_type operator[] (const uword i) const
    {
        return R[i];
    }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return void_ptr(&Q.Q) == void_ptr(&X); }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_div> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned Polynomial<elem_type>  Q;//direct evaluation of the left expression.
    arma_aligned Polynomial<elem_type>  R;//direct evaluation of the right expression.
    arma_aligned Polynomial<elem_type>  T;//direct evaluation of the divide expression.

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_div>& A)
        : Q(A.A), R(A.B)
    {
        arma_extra_debug_sigprint();
        Polynomial<elem_type> M;
        T = deconv(Q, R, M);
    }

    arma_inline uword                    size() const { return Q.size() + R.size() - 1; }
    inline elem_type operator[] (const uword i) const
    {
        return T[i];
    }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return false; }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_rem> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned Polynomial<elem_type>  Q;//direct evaluation of the left expression.
    arma_aligned Polynomial<elem_type>  R;//direct evaluation of the right expression.
    arma_aligned Polynomial<elem_type>  T;//direct evaluation of the divide expression.

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_rem>& A)
        : Q(A.A), R(A.B)
    {
        arma_extra_debug_sigprint();
        deconv(Q, R, T);
    }

    arma_inline uword                    size() const { return Q.size() + R.size() - 1; }
    inline elem_type operator[] (const uword i) const
    {
        return T[i];
    }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return false; }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_pow> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned Polynomial<elem_type>  Q;//direct evaluation of the left expression.
    arma_aligned const std::size_t      R;//exponent
    arma_aligned const std::size_t      lhs_size; //size of the left operand.

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_pow>& A)
        : Q(A.A), R(A.B), lhs_size(Q.m_size)
    {
        arma_extra_debug_sigprint();
        Q ^= R; //redirect to Polynomial::operator^(size_t).
    }

    arma_inline uword                    size() const { return (lhs_size-1) * R + 1; }
    inline elem_type operator[] (const uword i) const
    {
        return Q[i];
    }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return false; }
};

template<typename T1, typename T2>
struct Expression< glue_mt<T1, T2, glue_negate> >
{
    typedef Expression<T1>                           expr1_type;
    typedef Expression<T2>                           expr2_type;
    typedef typename expr1_type::elem_type           eT1;
    typedef typename expr2_type::elem_type           eT2;
    typedef typename common_type<eT1, eT2>::type     elem_type;
    typedef typename get_pod_type<elem_type>::result pod_type;

    arma_aligned const expr1_type Q;

    arma_inline explicit
    Expression(const glue_mt<T1, T2, glue_negate>& A)
        : Q(A.A)
    {
        arma_extra_debug_sigprint();
    }

    arma_inline uword                         size() const { return Q.size(); }
    arma_inline elem_type operator[] (const uword i) const { return -Q[i];    }

    template<typename eT2>
    arma_inline bool is_alias(const Polynomial<eT2>& X) const { return Q.is_alias(X); }
};



}


#endif // POLY_EXPRESSION_HPP_INCLUDED
