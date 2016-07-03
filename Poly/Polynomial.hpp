#ifndef POLYNOMIAL_HPP_INCLUDED
#define POLYNOMIAL_HPP_INCLUDED

namespace poly
{

using namespace std;

template<typename eT>
class Polynomial : public PolyBase<eT, Polynomial<eT> >
{

public:
    typedef PolyBase<eT, Polynomial<eT>>        base_type;

    typedef typename base_type::elem_type       elem_type;
    typedef typename get_pod_type<eT>::result   pod_type;

    typedef typename base_type::size_type       size_type;

    //properties, read-only
    arma_aligned const size_type        m_size;
    arma_aligned const elem_type* const m_data;

    //destructor
    arma_inline
    ~Polynomial();

    //constructors
    arma_inline
    Polynomial();

    arma_inline
    explicit Polynomial(const size_type size);

    arma_inline
    explicit Polynomial(const vector<eT>& v, const bool high_end = true);

    arma_inline
    explicit Polynomial(const valarray<eT>& x, const bool high_end = true);

    template<std::size_t N>
    arma_inline
    explicit Polynomial(const array<eT, N>& x, const bool high_end = true);

    arma_inline
    Polynomial(const Polynomial&          other);
    arma_inline
    Polynomial(Polynomial&&               temp); //std::move
    arma_inline
    Polynomial(const initializer_list<eT> list, const bool high_end = true);
    arma_inline
    Polynomial(const string&              expr, const bool high_end = true);
    arma_inline
    Polynomial(const char*                expr, const bool high_end = true);
    arma_inline
    explicit Polynomial(const Row<eT>& r,       const bool high_end = true);
    arma_inline
    explicit Polynomial(const Col<eT>& c,       const bool high_end = true);
    arma_inline
    Polynomial(const eT* array, const size_type size, const bool high_end = true);

    template<typename T>
    arma_inline
    explicit Polynomial(const vector<T>&   v, const bool high_end = true, typename convertible_type_only<T, eT>::type* junk = 0);
    template<typename T>
    arma_inline
    explicit Polynomial(const valarray<T>& x, const bool high_end = true, typename convertible_type_only<T, eT>::type* junk = 0);
    template<typename T, size_t N>
    arma_inline
    explicit Polynomial(const array<T, N>& x, const bool high_end = true, typename convertible_type_only<T, eT>::type* junk = 0);

    template<typename T>
    arma_inline
    explicit Polynomial(const Row<T>& v, const bool high_end = true,      typename convertible_type_only<T, eT>::type* junk = 0);
    template<typename T>
    arma_inline
    explicit Polynomial(const Col<T>& v, const bool high_end = true,      typename convertible_type_only<T, eT>::type* junk = 0);

    arma_inline
    const Polynomial& operator=(const Polynomial&          other);
    arma_inline
    const Polynomial& operator=(const initializer_list<eT> list);
    arma_inline
    const Polynomial& operator=(const string&              expr);
    arma_inline
    const Polynomial& operator=(const char*                expr);
    arma_inline
    const Polynomial& operator=(      Polynomial&&         rhs); //std::move

    template<typename T>
    arma_inline
    Polynomial(const Polynomial<T>& rhs, const typename convertible_type_only<T, eT>::type* = 0);

    template<typename T>
    arma_inline const typename convertible_type_only<T, eT, Polynomial>::type&
    operator=(const Polynomial<T>& rhs);

    template<typename eT1, typename T>
    arma_inline
    Polynomial(const PolyBase<eT1, T>& rhs, const typename convertible_type_only<eT1, eT>::type* = 0);

    template<typename eT1, typename T>
    arma_inline const typename convertible_type_only<eT1, eT, Polynomial>::type&
    operator= (const PolyBase<eT1, T>& rhs);

    template<typename T1, typename T2, typename glue_type>
    arma_inline
    Polynomial(const glue_mt<T1, T2, glue_type>& rhs);

    template<typename T1, typename T2, typename glue_type>
    arma_inline const Polynomial&
    operator= (const glue_mt<T1, T2, glue_type>& rhs);

    arma_inline
    const bool
    equal(const Polynomial&) const;

    arma_inline
    void swap(Polynomial&);

    arma_inline
    eT monic();

    //methods
    arma_inline
    const int       degree() const;
    arma_inline
    const size_type size() const;
    arma_inline
    bool            is_empty() const;

    bool            is_zero() const;

    arma_inline
    void            fill(const eT fill_value);
    arma_inline
    void            reset();
    arma_inline
    void            resize(const size_type size, bool copy = false);
    arma_inline
    const eT        coef(const size_type order) const;
    arma_inline
    eT&             coef(const size_type order);

    arma_inline
    eT              eval(const eT x) const;
    template<typename T>
    arma_inline
    complex<typename common_type<T,eT>::type>
    eval(const complex<T>& x) const;

    arma_inline
    const eT        at(const size_type index) const;
    arma_inline
    eT&             at(const size_type index);
    arma_inline
    const eT        operator[](const size_type index) const;
    arma_inline
    eT&             operator[](const size_type index);
    arma_inline
    const eT        operator()(const size_type index) const;
    arma_inline
    eT&             operator()(const size_type index);

    arma_inline
    const eT*       begin() const;
    arma_inline
    eT*             begin();
    arma_inline
    const eT*       end() const;
    arma_inline
    eT*             end();

    arma_inline const Polynomial& operator+=(const Polynomial& rhs);
    arma_inline const Polynomial& operator-=(const Polynomial& rhs);
    arma_inline const Polynomial& operator*=(const Polynomial& rhs);

    arma_inline const Polynomial& operator+=(const elem_type   rhs);
    arma_inline const Polynomial& operator-=(const elem_type   rhs);
    arma_inline const Polynomial& operator*=(const elem_type   rhs);
    arma_inline const Polynomial& operator/=(const elem_type   rhs);
    arma_inline const Polynomial& operator^=(const size_type   ord); //! power of the Polynomial


    template<typename T>
    arma_inline const typename convertible_type_only<T, eT, Polynomial>::type&
    operator+=(const Polynomial<T>& rhs);
    template<typename T>
    arma_inline const typename convertible_type_only<T, eT, Polynomial>::type&
    operator-=(const Polynomial<T>& rhs);
    template<typename T>
    arma_inline const typename convertible_type_only<T, eT, Polynomial>::type&
    operator*=(const Polynomial<T>& rhs);

    template<typename T1, typename T2, typename glue_type>
    arma_inline const Polynomial&
    operator+=(const glue_mt<T1, T2, glue_type>& expr);
    template<typename T1, typename T2, typename glue_type>
    arma_inline const Polynomial&
    operator-=(const glue_mt<T1, T2, glue_type>& expr);
    template<typename T1, typename T2, typename glue_type>
    arma_inline const Polynomial&
    operator*=(const glue_mt<T1, T2, glue_type>& expr);

    inline
    void impl_print(const char* extra_text) const;

    inline bool
    save(const std::string   name, const file_type type = arma_binary, const bool print_status = true) const;
    inline bool
    save(      std::ostream& os,   const file_type type = arma_binary, const bool print_status = true) const;
    inline bool
    load(const std::string   name, const file_type type = auto_detect, const bool print_status = true);
    inline bool
    load(      std::istream& is,   const file_type type = auto_detect, const bool print_status = true);


private:
    void init();
    void init(const eT* c, bool high_end);
    void init(const char* expr, bool high_end);
    template<typename T>
    void init(typename boost::call_traits<T>::param_type c, bool high_end);


};



}

#endif // POLYNOMIAL_HPP_INCLUDED
