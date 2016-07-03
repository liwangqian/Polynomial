#ifndef FN_POLY2SYM_HPP_INCLUDED
#define FN_POLY2SYM_HPP_INCLUDED

namespace poly
{


template<typename eT>
arma_inline std::string
poly2sym(const Polynomial<eT>& p, const char* sym = "x")
{
    const int deg = p.degree();
    if(deg < 0)
        return std::string("[]");

    const char* symbol = (sym == 0) ? "x" : sym;

    std::ostringstream in;

    int ix = 0;
    for(; ix <= deg; ++ix)
        if(p[ix] != eT(0))
            break;

    eT tmp = p[ix];
    if(ix == 0)
    {
        in << tmp;
        ++ix;
        tmp = p[ix];
        if(tmp != eT(0))
        {
            in << ((tmp > 0) ? " + " : " - ");
            if(tmp != eT(1) && tmp != eT(-1))
                in << std::fabs(tmp) << " ";
            in << symbol;
        }
    }
    else if(ix == 1)
    {
        if(tmp != eT(1) && tmp != eT(-1))
            in << tmp << " ";
        in << symbol;
    }
    else
    {
        if(tmp != eT(1) && tmp != eT(-1))
            in << tmp << " ";
        in << symbol << "^" << ix;
    }

    for(int i = ix+1; i <= deg; ++i)
    {
        if(p[i] == eT(0))
            continue;
        tmp = p[i];
        in << ((tmp > 0) ? " + " : " - ");
        if(tmp != eT(1) && tmp != eT(-1))
            in << std::fabs(tmp) << " ";
        in << symbol << "^" << i;
    }

    return in.str();
}


}

#endif // FN_POLY2SYM_HPP_INCLUDED
