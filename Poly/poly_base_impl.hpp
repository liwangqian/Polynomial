#ifndef POLY_BASE_IMPL_HPP_INCLUDED
#define POLY_BASE_IMPL_HPP_INCLUDED

namespace poly
{

template<typename eT, typename derived>
arma_inline const derived&
PolyBase<eT, derived>::get_ref() const
{
    return static_cast<const derived&>(*this);
}

template<typename eT, typename derived>
arma_inline void
PolyBase<eT, derived>::print(const char* extra_text) const
{
    const Proxy<derived> P(get_ref());
    P.Q.impl_print(extra_text);
}

}

#endif // POLY_BASE_IMPL_HPP_INCLUDED
