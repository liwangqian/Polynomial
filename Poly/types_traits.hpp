#ifndef TYPES_TRAITS_HPP_INCLUDED
#define TYPES_TRAITS_HPP_INCLUDED

namespace poly
{

template<typename T>
struct is_supported_elem_type
{
    static const bool value = std::is_arithmetic<T>::value;
};

template<typename T, typename R = void>
struct supported_elem_type_only
    : public std::enable_if<is_supported_elem_type<T>::value, R>
{

};

template<typename S, typename D, typename R = void>
struct convertible_type_only
    : public std::enable_if<std::is_convertible<S, D>::value, R>
{

};

template<typename T>
struct is_polynomial : public std::integral_constant<bool, false> {};

template<typename eT>
struct is_polynomial<Polynomial<eT> > : public std::integral_constant<bool, true> {};

template<typename T>
struct is_glue_mt :  public std::integral_constant<bool, false> {};

template<typename T1, typename T2, typename glue_type>
struct is_glue_mt<glue_mt<T1, T2, glue_type> > :  public std::integral_constant<bool, true> {};

template<typename T>
struct is_glue_st :  public std::integral_constant<bool, false> {};

template<typename T1, typename T2, typename glue_type>
struct is_glue_st<glue_st<T1, T2, glue_type> > :  public std::integral_constant<bool, true> {};

template<typename T>
struct is_PolyType
{
    static const bool value = is_polynomial<T>::value ||
                              is_glue_mt<T>::value   ||
                              is_glue_st<T>::value ;
};

template<typename T>
struct get_elem_type
{
    typedef T   type;
};

template<typename eT, typename T>
struct get_elem_type<PolyBase<eT, T> >
{
    typedef eT  type;
};

template<typename eT>
struct get_elem_type<Polynomial<eT> >
{
    typedef eT  type;
};

template<typename T1, typename T2, typename g>
struct get_elem_type<glue_mt<T1, T2, g> >
{
    typedef typename glue_mt<T1, T2, g>::elem_type   type;
};

template<typename T1, typename T2, typename g>
struct get_elem_type<glue_st<T1, T2, g> >
{
    typedef typename glue_st<T1, T2, g>::elem_type   type;
};

template<typename T1, typename T2>
struct get_common_elem_type
{
    typedef typename get_elem_type<T1>::type            eT1;
    typedef typename get_elem_type<T2>::type            eT2;
    typedef typename std::common_type<eT1, eT2>::type   type;
};

}

#endif // TYPES_TRAITS_HPP_INCLUDED
