#ifndef GLUE_MT_HPP_INCLUDED
#define GLUE_MT_HPP_INCLUDED

namespace poly
{


//! class glue_mt is the holder for expressions with mixed operator types.
template
<
    typename T1,
    typename T2,
    typename glue_type
>
struct glue_mt
: public PolyBase<typename get_common_elem_type<T1, T2>::type, glue_mt<T1, T2, glue_type> >
{
    typedef typename get_common_elem_type<T1, T2>::type elem_type;
    typedef typename get_pod_type<elem_type>::result    pod_type;
    typedef T1                                          lhs_type;
    typedef T2                                          rhs_type;
    typedef const T1&                                   lhs_stored_type; //! Must be PolyBase and it's derived class.
    typedef typename boost::call_traits<T2>::param_type rhs_stored_type; //! Can be Scalar and Polynomial, PolyGlue,...

    glue_mt(lhs_stored_type in_A, rhs_stored_type in_B);
    ~glue_mt();

    lhs_stored_type A;  //!< first operand
    rhs_stored_type B;  //!< second operand
};


}

#endif // GLUE_MT_HPP_INCLUDED
