#ifndef FUNCTORS_HPP_INCLUDED
#define FUNCTORS_HPP_INCLUDED

namespace poly
{

struct glue_plus
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_plus>& Expr);

};

struct glue_subs
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_subs>& Expr);
};

struct glue_mul
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_mul>& Expr);
};

struct glue_div
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_div>& Expr);
};

struct glue_rem
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_rem>& Expr);
};

struct glue_negate
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_negate>& Expr);
};

struct glue_pow
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_pow>& Expr);
};

struct glue_int
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_int>& Expr);
};

struct glue_diff
{
    template<typename eT, typename T1, typename T2>
    static void
    apply(Polynomial<eT>& R, const glue_st<T1, T2, glue_diff>& Expr);
};

}

#endif // FUNCTORS_HPP_INCLUDED
