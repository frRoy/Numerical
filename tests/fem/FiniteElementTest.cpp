/**
 *  @file    FiniteElementTest.cpp
 *  @brief   Test the FiniteElement class.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fem/FiniteElement.hpp"


namespace
{
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
/**
*/
double test_lagrange_polynomial(int d, int i) {
    Vec points = Vec::LinSpaced(d+1, 0.0, 1.0);
    double out = numerical::fem::lagrange_polynomial<double>(
        0.25, i, points);
    return out;
}
/**
*/
double test_basis(int d, int i, int j, double x, double h=0.25){
    return numerical::fem::phi(i, j, x, h, d);
}
/**
*  Check the integration of f(x) = 3*x^2+2*x+3 for x in [-1, 1]. The 
*  integration gives: x^3+x^2+3x which evaluates to 8.0.
*/
double test_quadrature(int n){
    std::vector<std::vector<double>> out = 
        numerical::fem::GaussLegendre<double>(n);
    std::vector<double> points = out[0];
    std::vector<double> weigths = out[1];
    double res = 0.0;
    for(int i=0; i<points.size(); i++){
        res += weigths[i]*(3*points[i]*points[i]+2*points[i]+3.0);
    }
    return res;
}

}  // namespace

TEST_CASE( "FEM functions tests are computed", "[mesh]" )
{
    CHECK(test_lagrange_polynomial(1, 0) == Approx(0.75));
    CHECK(test_lagrange_polynomial(1, 1) == Approx(0.25));
    CHECK(test_lagrange_polynomial(2, 0) == Approx(0.375));
    CHECK(test_lagrange_polynomial(2, 1) == Approx(0.75));
    CHECK(test_lagrange_polynomial(2, 2) == Approx(-0.125));
    // test basis
    // d = 0
    CHECK(test_basis(0, 0, 0, 0.0) == Approx(1.0));
    // d = 1
    CHECK(test_basis(1, 0, 0, 0.0) == Approx(0.5));
    CHECK(test_basis(1, 0, 1, 0.0) == Approx(0.5));
    CHECK(test_basis(1, 0, 0, 0.25) == Approx(0.375));
    CHECK(test_basis(1, 0, 1, 0.25) == Approx(0.625));
    // d = 2
    CHECK(test_basis(2, 0, 0, 0.0) == Approx(-0.0));
    CHECK(test_basis(2, 0, 1, 0.0) == Approx(1.0));
    CHECK(test_basis(2, 0, 2, 0.0) == Approx(0.0));
    CHECK(test_basis(2, 0, 0, 0.25) == Approx(-0.09375));
    CHECK(test_basis(2, 0, 1, 0.25) == Approx(0.9375));
    CHECK(test_basis(2, 0, 2, 0.25) == Approx(0.15625));
    // test basis derivatives
    // d = 0
    CHECK(test_basis(0, 1, 0, 0.0, 0.5) == Approx(0.0));
    // d = 1
    CHECK(test_basis(1, 1, 0, 0.0, 0.5) == Approx(-2.0));
    CHECK(test_basis(1, 1, 1, 0.0, 0.5) == Approx(2.0));
    CHECK(test_basis(1, 1, 0, 0.25, 0.5) == Approx(-2.0));
    CHECK(test_basis(1, 1, 1, 0.25, 0.5) == Approx(2.0));
    // d = 2
    CHECK(test_basis(2, 1, 0, 0.0, 0.5) == Approx(-2.0));
    CHECK(test_basis(2, 1, 1, 0.0, 0.5) == Approx(0.0));
    CHECK(test_basis(2, 1, 2, 0.0, 0.5) == Approx(2.0));
    CHECK(test_basis(2, 1, 0, 0.25, 0.5) == Approx(-1.0));
    CHECK(test_basis(2, 1, 1, 0.25, 0.5) == Approx(-2.0));
    CHECK(test_basis(2, 1, 2, 0.25, 0.5) == Approx(3.0));
    // quadrature
    CHECK(test_quadrature(2) == Approx(8.0));
}
