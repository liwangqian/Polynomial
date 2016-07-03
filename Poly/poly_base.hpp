#ifndef POLY_BASE_HPP_INCLUDED
#define POLY_BASE_HPP_INCLUDED

namespace poly
{

template<typename eT, typename derived>
struct PolyBase
{
    typedef eT                  elem_type;
    typedef derived             derived_type;
    typedef std::size_t         size_type;

    const derived& get_ref() const;

    //! Print the Polynomial to std::cout;
    void print(const char* = 0) const;
};


}

#endif // POLY_BASE_HPP_INCLUDED
