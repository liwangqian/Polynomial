#ifndef GLUE_ST_IMPL_HPP_INCLUDED
#define GLUE_ST_IMPL_HPP_INCLUDED

namespace poly
{

template<typename T1, typename T2, typename glue_type>
glue_st<T1, T2, glue_type>::glue_st(lhs_stored_type in_A, rhs_stored_type in_B)
    : A(in_A), B(in_B)
{
    arma_extra_debug_sigprint();
}

template<typename T1, typename T2, typename glue_type>
glue_st<T1, T2, glue_type>::~glue_st()
{
    arma_extra_debug_sigprint();
}



}

#endif // GLUE_ST_IMPL_HPP_INCLUDED
