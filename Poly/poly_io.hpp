#ifndef POLY_IO_HPP_INCLUDED
#define POLY_IO_HPP_INCLUDED

namespace poly
{

struct poly_ostream
{

template<typename eT>
arma_inline static void
print(std::ostream& o, const Polynomial<eT>& p, const bool modify);

template<typename eT>
arma_inline static bool
save(const Polynomial<eT>& p, const std::string& file_name, const file_type = raw_ascii);

template<typename eT>
arma_inline static bool
save(const Polynomial<eT>& p, std::ostream& o);


};

struct poly_istream
{



};


}

#endif // POLY_IO_HPP_INCLUDED
