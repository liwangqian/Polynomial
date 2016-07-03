#ifndef POLY_IO_IMPL_HPP_INCLUDED
#define POLY_IO_IMPL_HPP_INCLUDED

namespace poly
{

template<typename eT>
arma_inline void
poly_ostream::print(std::ostream& o, const Polynomial<eT>& p, const bool modify)
{
    arma_extra_debug_sigprint();
    const arma_ostream_state stream_state(o);
    const std::streamsize cell_width = modify ? arma_ostream::modify_stream(o, p.m_data, p.m_size) : o.width();

    const uword p_size = p.size();
    const uword n_elem = p.degree() + 1;

    if(n_elem == 0) //all zero elements or empty
    {
        if(p_size > 0)
        {
            if(cell_width > 0)
            {
                for(uword i = 0; i < p_size; ++i)
                {
                    o.width(cell_width);
                    arma_ostream::print_elem(o, p[i], modify);
                }
            }
            else
            {
                for(uword i = 0; i < p_size; ++i)
                {
                    arma_ostream::print_elem(o, p[i], modify);
                    o << ' ';
                }
            }
        }
        else
        {
            o << "[]";
        }
    }
    else // with nonzero elements
    {
        if(cell_width > 0)
        {
            for(uword i = 0; i < n_elem; ++i)
            {
                o.width(cell_width);
                arma_ostream::print_elem(o, p[i], modify);
            }
        }
        else
        {
            for(uword i = 0; i < n_elem; ++i)
            {
                arma_ostream::print_elem(o, p[i], modify);
                o << ' ';
            }
        }
    }
    o << '\n';
    o.flush();
    stream_state.restore(o);
}



}

#endif // POLY_IO_IMPL_HPP_INCLUDED
