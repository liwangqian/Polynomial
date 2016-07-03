#include <iostream>
#include "polynomial"
#include <boost/call_traits.hpp>

using namespace std;
using namespace arma::poly;

typedef Polynomial<double> dpoly_t;
typedef Polynomial<int> ipoly_t;
typedef vector<double> dvector;
typedef vector<complex<double> > root_holder;

template<typename Container>
void print(const Container& c, const string& text)
{
    typedef typename Container::elem_type eT;
    cout << "[" << text << "] = " ;
    cout << setprecision(4) << setiosflags(std::ios_base::fixed);
    for(auto x : c)
    {
        cout << setw(9) << x << " ";
    }

    cout << endl;
}

template<typename eT>
/*const*/Polynomial<eT> test_move()
{
    Polynomial<eT> p({1,2,3});
    return std::move(p);
}

void test0()
{
    cout << "test0..." << endl;

    valarray<double> va = {1,2,3,4};
    dvector v1 = {1,2,3};

    dpoly_t p1(v1);
    cout << p1.degree() << endl;
    cout << p1.is_empty() << endl;
    cout << p1[1] << endl;
    cout << p1(1) << endl;
    cout << p1.at(1) << endl;
    dpoly_t px = test_move<double>();
    px.print("px");

    dpoly_t p2 = {1,2,3,4,5,6,7};
    p2.print("p2");

    dpoly_t p3(9);
    p3.fill(12.9);
    p3.print("p3");

    p2 += p3;
    p2.print("p2 += p3");

    p2 -= p3;
    p2.print("p2 -= p3");

    p2 *= p3;
    p2.print("p2 *= p3");

    dpoly_t y1 = {1.0, 2, 3};
    y1 += 3.0;
    y1.print("y1 += 3.0");

    y1 -= 3.0;
    y1.print("y1 -= 3.0");

    y1 *= 3.0;
    y1.print("y1 *= 3.0");

    y1 /= 3.0;
    y1.print("y1 /= 3.0");

    y1 /= 0.0;
    y1.print("y1 /= 0.0");

    dpoly_t ps = "1 2 3 4 5.0";
    ps.print("ps");

    dpoly_t p4 = {1,2,3};
    p4 ^= 6;
    p4.print("p4");

    ipoly_t ip1 = {1,2,3,4};
    dpoly_t dp1 = {3,4,5};

    dp1 += ip1;
    dp1.print("dp1 += ip1");

    dp1 -= ip1;
    dp1.print("dp1 -= ip1");

    ip1 += dp1;
    ip1.print("ip1 += dp1");

    ip1 -= dp1;
    ip1.print("ip1 -= dp1");

    ip1 *= dp1;
    ip1.print("ip1 *= dp1");

    ipoly_t ip2 = {1,2,3};
    ipoly_t ip3 = {1,2};
    dpoly_t dp4 = 1.2 + 2 + ip1 + ip2 * dp1 - ip3 + 1;
    dp4.print("dp4");

    ipoly_t ip4 = {1,2,3,4};
    ipoly_t ip5 = -ip4 + ip2 * dp1 - ip3 + 1;
    ip5.print("ip5");

    ip4 += -ip2 * ip3;
    ip4.print("ip4");

    ip4 -= -ip2 * ip3;
    ip4.print("ip4");

    ip4 *= -ip2 * ip3;
    ip4.print("ip4");

    ipoly_t ip6 = {1,2,1,2,1};
    root_holder roots = root(ip6);
    cout << "root(ip6) = ";
    for(complex<double> x : roots)
        cout << x << " ";
    cout << endl;

    arma::Mat<double> A;
    A << 1.0+1e-5 << 1.0      << arma::endr
      << 0.0      << 1000.0-1e-5 << arma::endr;
    A.print("A");
    arma::Mat<double> f_A = arma::exp(A);
    f_A.print("exp(A)");
    arma::exp(A).print("Error");


    dpoly_t dp5({1.0,1.0,1.0/2,1.0/6,1.0/24,1.0/120});
    dp5.print("dp5");
    dpoly_t N, D;
    pade<dpoly_t>(dp5, 6, N, D, 0, 5);
    N.print("N");
    D.print("D");
    cout << "dp5^2 = " << (dp5^2) << endl;


    arma::Row<double> c = {1,2,3};
    arma::Row<double> r = {1,2,3,4};
    arma::Mat<double> H;
    arma::Mat<double> W;
    hankel(H, c, c);
    rot90(W, H);
    H.print("H");
    W.print("W");

    dpoly_t dp6 = {1,2,1,1}, dp7 = {1,1,1};
    dpoly_t dp9 = dp6 % dp7 + dp6;
    dp9.print("dp9");

    dpoly_t dp10 = {1,1};
    dpoly_t dp11 = dp10 ^ 3;
    dp11.print("dp11");

    cout << dp7.degree() << endl;
    cout << dp7.eval(1) << endl;

    dpoly_t dp13 = {1,2,3,4,5,6};
    print(polyint(dp13+dp10,1), "int(dp13)");
    print(polyder(dp13+dp10,1), "diff(dp13)");

    dp13.print("dp13");
    (dp13+dp10).print("dp13+dp10");

    dpoly_t dp0 = {1,2,3,4};
    dp0.print("dp0");
    dp9.print("dp9");
    dp0 = dp0 * dp7 + dp9;
    dp0.print("dp0 * dp7 + dp9");

    cout << pow(dp7+2) << endl;

    double x[201] = {0}, y[201] = {0};
    for(int i = 1; i < 201; ++i)
    {
        x[i] = x[i-1] + 0.02;
        y[i] = sin(x[i]);
    }

    cout << polyfit(x, y, 201, 8) << endl;

    dpoly_t R;
    cout << dp0 << endl << dp1 << endl;
    dpoly_t Q = deconv(dp0, dp1, R);
    cout << Q << endl;
    cout << R << endl;
    cout << (conv(Q, dp1) + R) << endl;


    arma::Mat<double> A1("1 2 3;4 5 6;7 8 9");
    arma::Row<double> V1("1 2 3");
    A1.print("A1");
    V1.print("V1");

    print(poly(V1), "poly(V1)");

    dpoly_t p_a = poly(A1);
    print(p_a, "poly(A1)");

    arma::Row<double> rv = {1, 2, 3};
    print(polyvalm(p_a, rv), "polyvalm(v)");

    arma::Mat<double> mm = polyvalm(p_a, A1);
    cout << mm << endl;
    cout << arma::norm(mm, 2) << endl;

    arma::Mat<double> A3 = std::move(A1);
    A3.print("A3");
    A1.print("A1");

    dpoly_t p_s = {0, 1, 0, 3, 5};
    cout << poly2sym(p_s) << endl;

    arma::Mat<double> me("1 2 3 4 5 6;2 3 4 5 6 7;3 4 5 6 7 8");
    arma::field<Polynomial<double> > pa;
    mat2poly(me, pa, 3);
    cout << pa << endl;

    arma::Mat<double> ms;
    poly2mat(pa, ms);
    cout << ms << endl;
}

void test1()
{
    cout << "test1..." << endl;

    dpoly_t p1 = {1, 2, 1};
    dpoly_t p2 = {1, 1};
    dpoly_t p3;

    cout << (p1/p2) << endl;
    cout << gcd(p2, p1) << endl;
    cout << lcm(p1, p2) << endl;


    dpoly_t p4 = p1*p2;

    cout << p4 << endl;

    p4 = lcm(p1, p3);
    cout << p4 << endl;

    p4 = gcd(p1, p3);
    cout << p4 << endl;
}

int main()
{
    std::cout << "Hello world!" << std::endl;

    test0();

    ////////
    test1();
    ////////
    return 0;
}
