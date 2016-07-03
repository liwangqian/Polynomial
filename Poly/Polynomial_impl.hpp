#ifndef POLYNOMIAL_IMPL_HPP_INCLUDED
#define POLYNOMIAL_IMPL_HPP_INCLUDED

namespace poly
{


template<typename eT, typename T>
arma_inline std::ostream&
operator<<(std::ostream& o, const PolyBase<eT, T>& expr)
{
    const Proxy<T> P(expr.get_ref());
    poly_ostream::print(o, P.Q, true);
    return o;
}

template<typename eT>
arma_inline void
Polynomial<eT>::init()
{
    arma_extra_debug_sigprint();
    access::rwp(m_data) = new elem_type[m_size];
}

template<typename eT>
arma_inline void
Polynomial<eT>::init(const eT* c, bool high_end)
{
    arma_extra_debug_sigprint();
    eT* d_ptr = new elem_type[m_size];
    access::rwp(m_data) = d_ptr;
    if(high_end)
        aux::copy<eT*, const eT*>(d_ptr, c, m_size);
    else
        aux::copy_backward<eT*, const eT*>(d_ptr, c, m_size);
}

template<typename eT>
arma_inline void
Polynomial<eT>::init(const char* expr, bool high_end)
{
    arma_extra_debug_sigprint();
    std::istringstream is(expr);
    eT elem;
    vector<eT> vec;
    while(is >> elem)
    {
        vec.push_back(elem);
    }
    access::rw(m_size) = vec.size();
    init(vec.data(), high_end);
}

template<typename eT>
template<typename T>
arma_inline void
Polynomial<eT>::init(typename boost::call_traits<T>::param_type c, bool high_end)
{
    arma_extra_debug_sigprint();
    typedef typename boost::call_traits<T>::param_type c_type;
    eT* d_ptr = new elem_type[m_size];
    access::rwp(m_data) = d_ptr;
    if(high_end)
        aux::copy<eT*, c_type>(d_ptr, c, m_size);
    else
        aux::copy_backward<eT*, c_type>(d_ptr, c, m_size);
}

template<typename eT>
arma_inline
Polynomial<eT>::~Polynomial()
{
    arma_extra_debug_sigprint_this(this);
    static_assert(is_supported_elem_type<eT>::value, "The element type is not supported!");
    if(m_data)
        delete m_data;
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial()
    : m_size(0), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const size_type size)
    : m_size(size), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    init();
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const vector<eT>& v, const bool high_end)
    : m_size(v.size()), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(v.data(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const valarray<eT>& v, const bool high_end)
    : m_size(v.size()), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(&v[0], high_end);
}

template<typename eT>
template<std::size_t N>
arma_inline
Polynomial<eT>::Polynomial(const array<eT, N>& v, const bool high_end)
    : m_size(N), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(v.data(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const Polynomial& rhs)
    : m_size(rhs.m_size), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(rhs.m_data, true);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(Polynomial&& temp)
    : m_size(temp.m_size), m_data(0)
{
//    cout << "move constructor" << endl;
    arma_extra_debug_sigprint_this(this);
    std::swap(access::rwp(m_data), access::rwp(temp.m_data));
    access::rw(temp.m_size) = 0;
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const initializer_list<eT> l, const bool high_end)
    : m_size(l.size()), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(l.begin(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const string& expr, const bool high_end)
    : m_size(0), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    init(expr.c_str(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const char* expr, const bool high_end)
    : m_size(0), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(expr == 0) return;
    init(expr, high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const Row<eT>& r, const bool high_end)
    : m_size(r.n_cols), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(r.memptr(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const Col<eT>& c, const bool high_end)
    : m_size(c.n_rows), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init(c.memptr(), high_end);
}

template<typename eT>
arma_inline
Polynomial<eT>::Polynomial(const eT* x, const size_type size, const bool high_end)
    : m_size(size), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(x == 0) return;
    init(x, high_end);
}

template<typename eT>
template<typename T>
arma_inline
Polynomial<eT>::Polynomial(const vector<T>& v, const bool high_end, typename convertible_type_only<T, eT>::type*)
    : m_size(v.size()), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init<const T*>(v.data(), high_end);
}

template<typename eT>
template<typename T>
arma_inline
Polynomial<eT>::Polynomial(const valarray<T>& v, const bool high_end, typename convertible_type_only<T, eT>::type*)
    : m_size(v.size()), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init<const T*>(&v[0], high_end);
}

template<typename eT>
template<typename T, std::size_t N>
arma_inline
Polynomial<eT>::Polynomial(const array<T, N>& v, const bool high_end, typename convertible_type_only<T, eT>::type*)
    : m_size(N), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init<const T*>(&v[0], high_end);
}

template<typename eT>
template<typename T>
arma_inline
Polynomial<eT>::Polynomial(const Row<T>& r, const bool high_end, typename convertible_type_only<T, eT>::type*)
    : m_size(r.n_cols), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init<const T*>(r.memptr(), high_end);
}

template<typename eT>
template<typename T>
arma_inline
Polynomial<eT>::Polynomial(const Col<T>& c, const bool high_end, typename convertible_type_only<T, eT>::type*)
    : m_size(c.n_rows), m_data(0)
{
    arma_extra_debug_sigprint_this(this);
    if(m_size == 0) return;
    init<const T*>(c.memptr(), high_end);
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(const Polynomial& rhs)
{
    arma_extra_debug_sigprint();
    if(this != &rhs && !rhs.is_empty())
    {
        access::rw(m_size) = rhs.m_size;
        if(m_data) delete m_data;
        init(rhs.m_data, true);
    }
    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(const initializer_list<eT> rhs)
{
    arma_extra_debug_sigprint();
    access::rw(m_size) = rhs.size();

    if(m_size == 0) return *this;

    if(m_data) delete m_data;
    init(rhs.begin(), true);
    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(const string& rhs)
{
    arma_extra_debug_sigprint();
    if(rhs.empty()) return;
    if(m_data) delete m_data;
    init(rhs.c_str(), true);
    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(const char* rhs)
{
    arma_extra_debug_sigprint();
    if(rhs == 0)    return;
    if(m_data) delete m_data;
    init(rhs, true);
    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(Polynomial&& rhs)
{
    arma_extra_debug_sigprint();
    if(&rhs != this)
    {
        swap(rhs);
    }

    return *this;
}

template<typename eT>
template<typename T>
arma_inline
Polynomial<eT>::Polynomial(const Polynomial<T>& rhs, const typename convertible_type_only<T, eT>::type*)
    : m_size(rhs.m_size), m_data(0)
{
    arma_extra_debug_sigprint();
    if(m_size == 0) return;
    init<const T*>(rhs.m_data, true);
}

template<typename eT>
template<typename T>
arma_inline const typename convertible_type_only<T, eT, Polynomial<eT> >::type&
Polynomial<eT>::operator=(const Polynomial<T>& other)
{
    arma_extra_debug_sigprint();
    access::rw(m_size) = other.size();
    if(m_size == 0) return;
    if(m_data) delete m_data;
    init<const T*>(other.m_data, true);
    return *this;
}

template<typename eT>
template<typename eT1, typename T>
arma_inline
Polynomial<eT>::Polynomial(const PolyBase<eT1, T>& rhs, const typename convertible_type_only<eT1, eT>::type*)
    : m_data(0), m_size(0)
{
    arma_extra_debug_sigprint();

    const Expression<T> expr(rhs.get_ref());
    access::rw(m_size) = expr.size();
    if(m_size == 0) return;
    init<Expression<T> >(expr);
}

template<typename eT>
template<typename eT1, typename T>
arma_inline const typename convertible_type_only<eT1, eT, Polynomial<eT> >::type&
Polynomial<eT>::operator=(const PolyBase<eT1, T>& rhs)
{
    arma_extra_debug_sigprint();
    const Expression<T> expr(rhs.get_ref());
    if(expr.is_alias(*this) == false)
    {
        access::rw(m_size) = expr.size();
        if(m_size > 0)
        {
            if(m_data) delete m_data;
            init<Expression<T> >(expr);
        }
    }

    return *this;
}

template<typename eT>
template<typename T1, typename T2, typename glue_type>
arma_inline
Polynomial<eT>::Polynomial(const glue_mt<T1, T2, glue_type>& rhs)
    : m_size(0), m_data(0)
{
    arma_extra_debug_sigprint();
    *this = rhs;
}

template<typename eT>
template<typename T1, typename T2, typename glue_type>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator=(const glue_mt<T1, T2, glue_type>& rhs)
{
    arma_extra_debug_sigprint();

    typedef Expression<glue_mt<T1, T2, glue_type> > expression_type;
    const expression_type expr(rhs);
    access::rw(m_size)  = expr.size();
    if(m_size == 0) return *this;

    const bool is_alias = expr.is_alias(*this);//check for alias
    eT* t_mem = 0;

    if(is_alias == true)
    {
        t_mem = new eT[m_size];
        aux::copy<eT*, const expression_type>(t_mem, expr, m_size);
    }

    if(m_data)  delete m_data;

    if(is_alias == true)
        access::rwp(m_data) = t_mem;
    else
        init<expression_type>(expr, true);

    return *this;
}

template<typename eT>
arma_inline const bool
Polynomial<eT>::equal(const Polynomial<eT>& y) const
{
    if(this == &y)
        return true;

    const int deg = degree();
    if(deg != y.degree())
        return false;

    for(int i = 0; i<= deg; ++i)
    {
        if(m_data[i] != y[i])
            return false;
    }

    return true;
}

template<typename eT>
arma_inline void
Polynomial<eT>::swap(Polynomial& x)
{
    std::swap(access::rw(m_size), access::rw(x.m_size));
    std::swap(access::rwp(m_data), access::rwp(x.m_data));
}

template<typename eT>
arma_inline eT
Polynomial<eT>::monic()
{
    const size_type end = degree();
    eT head = m_data[end];

    if(head != eT(1))
        for(size_type i = 0; i <= end; ++i)
            access::rw(m_data[i]) /= head;
    return head;
}

template<typename eT>
arma_inline const int
Polynomial<eT>::degree() const
{
    arma_extra_debug_sigprint();
    int d = m_size - 1;
    while(d >= 0 && m_data[d] == eT(0)) --d;
    d = (m_size > 0 && d < 0) ? 0 : d;  //all zero polynomial is not empty polynomial
    return d;
}

template<typename eT>
arma_inline const typename Polynomial<eT>::size_type
Polynomial<eT>::size() const
{
    arma_extra_debug_sigprint();
    return m_size;
}

template<typename eT>
arma_inline bool
Polynomial<eT>::is_empty() const
{
    arma_extra_debug_sigprint();
    return (m_size == 0);
}

template<typename eT>
arma_inline bool
Polynomial<eT>::is_zero() const
{
    arma_extra_debug_sigprint();
    const bool f = (m_size == 0) || (degree() == 0 && m_data[0] == elem_type(0));
    return f;
}

template<typename eT>
arma_inline void
Polynomial<eT>::fill(const eT fill_value)
{
    arma_extra_debug_sigprint();
    eT* d_ptr = access::rwp(m_data);
    for(size_type i = 0; i < m_size; i++)
        d_ptr[i] = fill_value;
}

template<typename eT>
arma_inline void
Polynomial<eT>::reset()
{
    arma_extra_debug_sigprint();
    if(m_data)
    {
        delete m_data;
        access::rwp(m_data) = 0;
        access::rw(m_size) = 0;
    }
}

template<typename eT>
arma_inline void
Polynomial<eT>::resize(const size_type size, bool copy)
{
    arma_extra_debug_sigprint();
    if(m_size < size)
    {
        eT* t_mem = new eT[size];
        if(copy)
            aux::copy<eT*, const eT*>(t_mem, m_data, m_size);
        for(size_type i = (copy ? m_size : 0); i < size; ++i)
            t_mem[i] = eT(0);

        delete m_data;
        access::rwp(m_data) = t_mem;

    }
    access::rw(m_size)  = size;
}

template<typename eT>
arma_inline const eT
Polynomial<eT>::coef(const size_type order) const
{
    return order < m_size ? m_data[order] : eT(0);
}

template<typename eT>
arma_inline eT&
Polynomial<eT>::coef(const size_type order)
{
    arma_debug_check(order >= m_size, "Polynomial<eT>::coef() : Index out of bounds!");
    return access::rw(m_data[order]);
}

template<typename eT>
arma_inline eT
Polynomial<eT>::eval(const eT x) const
{
    const int deg = degree();
    if(deg < 0 )
    {
        arma_warn(true, "Warning: evaluate an empty polynomial!");
        return numeric_limits<eT>::quiet_NaN();
    }
    if(deg == 0 || x == eT(0)) return m_data[0];

    eT acc = m_data[deg];
    for(int i = deg - 1; i >= 0; --i)
        acc = acc * x + m_data[i];
    return acc;
}

template<typename eT>
template<typename T>
arma_inline complex<typename common_type<T,eT>::type>
Polynomial<eT>::eval(const complex<T>& x) const
{
    typedef complex<typename common_type<T,eT>::type> result_type;

    const int deg = degree();
    if(deg < 0 )
    {
        arma_warn(true, "Warning: evaluate of an empty polynomial!");
        return numeric_limits<result_type>::quiet_NaN();
    }
    if(deg == 0 || x == eT(0)) return m_data[0];

    result_type acc(m_data[deg]);
    for(int i = deg - 1; i >= 0; --i)
        acc = acc * x + m_data[i];
    return acc;
}

template<typename eT>
arma_inline const eT
Polynomial<eT>::at(const size_type index) const
{
    arma_extra_debug_sigprint();
    arma_debug_check(index >= m_size, "Polynomial::at() : index out of bound!");
    return m_data[index];
}

template<typename eT>
arma_inline eT&
Polynomial<eT>::at(const size_type index)
{
    arma_extra_debug_sigprint();
    arma_debug_check(index >= m_size, "Polynomial::at() : index out of bound!");
    return access::rw(m_data[index]);
}

template<typename eT>
arma_inline const eT
Polynomial<eT>::operator[](const size_type index) const
{
    arma_extra_debug_sigprint();
    return index < m_size ? m_data[index] : eT(0);
}

template<typename eT>
arma_inline eT&
Polynomial<eT>::operator[](const size_type index)
{
    arma_extra_debug_sigprint();
    return access::rw(m_data[index]);
}

template<typename eT>
arma_inline const eT
Polynomial<eT>::operator()(const size_type index) const
{
    arma_extra_debug_sigprint();
    return index < m_size ? m_data[index] : eT(0);
}

template<typename eT>
arma_inline eT&
Polynomial<eT>::operator()(const size_type index)
{
    arma_extra_debug_sigprint();
    arma_debug_check(index >= m_size, "Polynomial::operator() : index out of bound!");
    return access::rw(m_data[index]);
}

template<typename eT>
arma_inline const eT*
Polynomial<eT>::begin() const
{
    arma_extra_debug_sigprint();
    return m_data;
}

template<typename eT>
arma_inline eT*
Polynomial<eT>::begin()
{
    arma_extra_debug_sigprint();
    return access::rwp(m_data);
}

template<typename eT>
arma_inline const eT*
Polynomial<eT>::end() const
{
    arma_extra_debug_sigprint();
    return m_data + m_size;
}

template<typename eT>
arma_inline eT*
Polynomial<eT>::end()
{
    arma_extra_debug_sigprint();
    return access::rwp(m_data) + m_size;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator+=(const Polynomial& rhs)
{
    arma_extra_debug_sigprint();

    if(rhs.m_size == 0) return *this;

    const size_type r_size = std::max(m_size, rhs.m_size);

    resize(r_size, true);
    eT* d_ptr = access::rwp(m_data);
    if(r_size == 1)
        d_ptr[0] += rhs.m_data[0];
    else
        aux::inplace_plus(d_ptr, rhs.m_data, r_size);

    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator-=(const Polynomial& rhs)
{
    arma_extra_debug_sigprint();

    if(rhs.m_size == 0) return *this;

    const size_type r_size = std::max(m_size, rhs.m_size);

    resize(r_size, true);
    eT* d_ptr = access::rwp(m_data);
    if(r_size == 1)
        d_ptr[0] -= rhs.m_data[0];
    else
        aux::inplace_minus(d_ptr, rhs.m_data, r_size);

    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator*=(const Polynomial& rhs)
{
    arma_extra_debug_sigprint();
    if(m_size == 0 || rhs.m_size == 0)
    {
        reset();
        return *this;
    }

    if(rhs.m_size == 1)
    {
        return (*this*=rhs[0]);
    }

    const size_type r_size = m_size + rhs.m_size - 1;
    eT*             r_mem  = new eT[r_size];

    aux::conv(r_mem, m_data, m_size, rhs.m_data, rhs.m_size);

    delete m_data;
    access::rwp(m_data) = r_mem;
    access::rw(m_size)  = r_size;

    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator+=(const elem_type rhs)
{
    arma_extra_debug_sigprint();
    if(m_size == 0)
    {
        access::rw(m_size) = 1;
        init();
        access::rw(m_data[0]) = eT(0);
    }
    access::rw(m_data[0]) += rhs;

    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator-=(const elem_type rhs)
{
    arma_extra_debug_sigprint();
    return (*this += -rhs);
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator*=(const elem_type rhs)
{
    arma_extra_debug_sigprint();
    if(m_size > 0 && rhs != eT(1))
    {
        eT* d_ptr = access::rwp(m_data);
        aux::inplace_mul(d_ptr, rhs, m_size);
    }
    return *this;
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator/=(const elem_type rhs)
{
    arma_extra_debug_sigprint();
    return (*this *= (eT(1) / rhs));
}

template<typename eT>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator^=(const size_type n)
{
    arma_extra_debug_sigprint();
    if(m_size == 0) return *this;
    if(n == 0 && m_size > 0)
    {
        resize(1);
        access::rw(m_data[0]) = eT(1);
        return *this;
    }
    if(n == 1)  return *this;

    const size_type size   = m_size - 1;
    const size_type r_size = n * size + 1;
    const size_type t_size = (n-1) * size + 1;
    eT*             r_mem  = new eT[r_size];
    eT*             t_mem  = new eT[t_size];

    aux::copy<eT*, const eT*>(t_mem, m_data, m_size);

    size_type rs;
    size_type ts;
    for(size_type i = 2; i < n; ++i)
    {
        rs = i * size + 1;
        ts = (i - 1) * size + 1;
        aux::conv(r_mem, t_mem, ts, m_data, m_size);
        aux::copy<eT*, const eT*>(t_mem, r_mem, rs);
    }
    aux::conv(r_mem, t_mem, t_size, m_data, m_size);

    delete m_data;
    delete t_mem;
    access::rwp(m_data) = r_mem;
    access::rw(m_size)  = r_size;

    return *this;
}

template<typename eT>
template<typename T>
arma_inline const typename convertible_type_only<T, eT, Polynomial<eT> >::type&
Polynomial<eT>::operator+=(const Polynomial<T>& rhs)
{
    arma_extra_debug_sigprint();

    if(rhs.m_size == 0) return *this;

    const size_type r_size = std::max(m_size, rhs.m_size);

    resize(r_size, true);
    eT* d_ptr = access::rwp(m_data);
    if(r_size == 1)
        d_ptr[0] += rhs.m_data[0];
    else
        aux::inplace_plus<eT*, const T*>(d_ptr, rhs.m_data, r_size);

    return *this;
}

template<typename eT>
template<typename T>
arma_inline const typename convertible_type_only<T, eT, Polynomial<eT> >::type&
Polynomial<eT>::operator-=(const Polynomial<T>& rhs)
{
    arma_extra_debug_sigprint();

    if(rhs.m_size == 0) return *this;

    const size_type r_size = std::max(m_size, rhs.m_size);

    resize(r_size, true);
    eT* d_ptr = access::rwp(m_data);
    if(r_size == 1)
        d_ptr[0] -= rhs.m_data[0];
    else
        aux::inplace_minus<eT*, const T*>(d_ptr, rhs.m_data, r_size);

    return *this;
}

template<typename eT>
template<typename T>
arma_inline const typename convertible_type_only<T, eT, Polynomial<eT> >::type&
Polynomial<eT>::operator*=(const Polynomial<T>& rhs)
{
    arma_extra_debug_sigprint();
    if(m_size == 0 || rhs.m_size == 0)
    {
        reset();
        return *this;
    }

    if(rhs.m_size == 1)
    {
        return (*this*=rhs[0]);
    }

    const size_type r_size = m_size + rhs.m_size - 1;
    eT*             r_mem  = new eT[r_size];
    eT*             t_mem  = new eT[rhs.m_size];

    aux::copy<eT*, const T*>(t_mem, rhs.m_data, rhs.m_size);
    aux::conv(r_mem, m_data, m_size, t_mem, rhs.m_size);

    delete m_data;
    delete t_mem;
    access::rwp(m_data) = r_mem;
    access::rw(m_size)  = r_size;

    return *this;
}

template<typename eT>
template<typename T1, typename T2, typename glue_type>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator+=(const glue_mt<T1, T2, glue_type>& rhs)
{
    typedef glue_mt<T1, T2, glue_type> ExprType;

    const Expression<ExprType> expr(rhs);
    const uword e_size = expr.size();

    if(e_size == 0) return *this;

    const size_type r_size = std::max<uword>(e_size, m_size);
    resize(r_size, true);

    eT* d_ptr = access::rwp(m_data);

    if(r_size == 1)
        d_ptr[0] += expr[0];
    else
        aux::inplace_plus<eT*, const Expression<ExprType> >(d_ptr, expr, r_size);

    return *this;
}

template<typename eT>
template<typename T1, typename T2, typename glue_type>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator-=(const glue_mt<T1, T2, glue_type>& rhs)
{
    typedef glue_mt<T1, T2, glue_type> ExprType;
    const Expression<ExprType> expr(rhs);
    const uword e_size = expr.size();

    if(e_size == 0) return *this;

    const size_type r_size = std::max<uword>(e_size, m_size);
    resize(r_size, true);

    eT* d_ptr = access::rwp(m_data);

    if(r_size == 1)
        d_ptr[0] -= expr[0];
    else
        aux::inplace_minus<eT*, const Expression<ExprType> >(d_ptr, expr, r_size);

    return *this;
}

template<typename eT>
template<typename T1, typename T2, typename glue_type>
arma_inline const Polynomial<eT>&
Polynomial<eT>::operator*=(const glue_mt<T1, T2, glue_type>& rhs)
{
    arma_extra_debug_sigprint();

    typedef glue_mt<T1, T2, glue_type> ExprType;
    const Expression<ExprType> expr(rhs);
    const size_type e_size = expr.size();

    if(m_size == 0 || e_size == 0)
    {
        reset();
        return *this;
    }

    if(e_size == 1)
    {
        return (*this*=expr[0]);
    }

    const size_type r_size = m_size + e_size - 1;
    eT* r_mem = new eT[r_size];
    eT* t_mem = new eT[e_size];

    aux::copy<eT*, const Expression<ExprType> >(t_mem, expr, e_size);
    aux::conv(r_mem, m_data, m_size, t_mem, e_size);

    delete m_data;
    delete t_mem;
    access::rwp(m_data) = r_mem;
    access::rw(m_size)  = r_size;

    return *this;
}

template<typename eT>
inline void
Polynomial<eT>::impl_print(const char* extra_text) const
{
    if(extra_text)  cout << extra_text << " = ";
    poly_ostream::print(cout, *this, true);
}



}

#endif // POLYNOMIAL_IMPL_HPP_INCLUDED
