#ifndef GLUE_MT_IMPL_HPP_INCLUDED
#define GLUE_MT_IMPL_HPP_INCLUDED

namespace poly
{

template<typename T1, typename T2, typename glue_type>
glue_mt<T1, T2, glue_type>::glue_mt(lhs_stored_type in_A, rhs_stored_type in_B)
    : A(in_A), B(in_B)
{
    arma_extra_debug_sigprint();
}

template<typename T1, typename T2, typename glue_type>
glue_mt<T1, T2, glue_type>::~glue_mt()
{
    arma_extra_debug_sigprint();
}


}

#endif // GLUE_MT_IMPL_HPP_INCLUDED
