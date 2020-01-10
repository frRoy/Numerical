/**
 *  @file    FiniteElement.hpp
 *  @brief   Defines basis functions.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#ifndef FINITEELEMENT_H
#define FINITEELEMENT_H

#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Core>
#include "spdlog/spdlog.h"

namespace numerical {
namespace fem {
/**
*  A linear mapping from \f$X\f$ to \f$x\f$ can be written as:
*
*  \f[
*    x = \frac{1}{2}\left(x_L+x_R\right)+\frac{1}{2}\left(x_R-x_L\right)X
*  \f]
*
*  @param X The local coordinate in the reference cell.
*  @param omega_e The left and right global coordinates of the element.
*  @return The global coordinate.
*/
template <typename T>
T affine_mapping(T X, const T (&omega_e)[2]){
    T x_L = omega_e[0];
    T x_R = omega_e[1];
    return 0.5*(x_L + x_R) + 0.5*(x_R - x_L)*X;
}
/**
*  Returns the \f$i\f$-th Lagrange polynomial. 
*
*  \f[
*      \phi_i(x) = \prod_{j=0, j\neq i}^N \frac{x-x_j}{x_i-x_j}
*  \f]
*
*  @param x The local coordinate in the reference cell.
*  @param i The Lagrange polynomial index.
*  @param nodes The interpolation nodes in the reference cell.
*  @return The \f$i\f$-th Lagrange polynomial.
*/
template <typename T>
T lagrange_polynomial(T x, int i, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& nodes){
    T p = 1.0;
    for(int j=0; j<nodes.size(); j++){
        if(j != i){
            p *= (x - nodes[j])/(nodes[i] - nodes[j]);
          }
    }
    return p;
};
/**
*  Returns the first derivative of the \f$i\f$-th Lagrange polynomial. 
*
*  \f[
*    \begin{align}
*      \frac{\partial}{\partial x}\phi_i(x) &= 
*      \frac{\partial}{\partial X}\phi_i(X)\frac{\partial X}{\partial x}\\
*      &= \frac{\partial}{\partial X}\phi_i(X)\frac{2}{h}
*    \end{align}
*  \f]
*
*  For details see Langtangen2019, p.188 Cellwise computation.
*
*  The first derivative of the Lagrange polynomial with respect to local 
*  coordinates \f$X\f$ can be obtained using the logarithm of the Lagrange 
*  polynamial, i.e.
*
*  \f[
*    \textrm{ln}\left(\prod_{j=0, j\neq i}^N \frac{x-x_j}{x_i-x_j}\right) =
*      \sum_{j=0, j\neq i}^N \textrm{ln}\left(\frac{x-x_j}{x_i-x_j}\right)
*  \f]
*
*  The derivative of the logarithm is
*
*  \f[
*    \frac{\partial}{\partial X}\left(\textrm{ln}(f(x))\right) = \frac{1}{f(X)}
*      \frac{\partial}{\partial X}f(X)
*  \f]
*
*  with \f$f(x)=\frac{x-x_j}{x_i-x_j}\f$, we get
*
*  \f[
*    \frac{\partial}{\partial X}\phi_i(X) = \phi_i(X)\sum_{j=0, j\neq i}^N 
*      \frac{1}{x-x_j}
*  \f]
*
*  Note that this result is only valid for \f$(x-x_j)\neq 0\f$.
*  
*  @param x The local coordinate.
*  @param i The polynomial index.
*  @param nodes The interpolation nodes in the reference cell.
*  @param h The element size.
*  @param d The polynomial order, default=1.
*  @return The \f$i\f$-th Lagrange polynomial.
*/
template <typename T>
T lagrange_polynomial_first_derivative(T x, int i, 
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& nodes, T h, int d=1){
    T p = 0.0;
    T denom = 1.0;
    T lp = lagrange_polynomial(x, i, nodes);
    for(int j=0; j<nodes.size(); j++){
        if(j != i){ // d > 1
            denom *= (nodes[i] - nodes[j]);
            if(d == 1){
                p += 1.0;
            } else if(d ==2){
                p += (x - nodes[j]);
            }
          }
    }
    p *= 2.0 / h / denom;
    return p;
};
/**
*  Returns all local basis function \f$\phi\f$ and their derivatives,
*  in physical coordinates, as functions of the local point \f$X\f$ in a 
*  reference cell with \f$d+1\f$ nodes.
*
*  Example:
*  ~~~~~~~~~~~~~~~{.cpp}
*  #include "numerical/fem/FiniteElement.hpp"
*
*  using namespace numerical;
*
*  int main()
*  {
*    double f0 = fem::phi(0, 0, 0.0);  // basis func 0 at X=0
*    double fp0 = fem::phi(0, 0, 0.0);  // 1st x-derivative at X=0
*    spdlog::info("f0: {}, fp0: {}", f0, fp0);
*    return 0;
*  }
*  ~~~~~~~~~~~~~~~
*  Output:
*  ~~~~~~~~~~~~~~~{.sh}
*  $ f0: 0.5, fp0: -2.0
*  ~~~~~~~~~~~~~~~
*
*  @param i 0 if function, 1 if first derivative.
*  @param j The index of the basis function.
*  @param x The local coordinate.
*  @param h The element size.
*  @param d The polynomial order.
*  @return The value of the basis function/derivative at the local 
*    coordinate.
*/
template <typename T>
T phi(int i, int j, T x, T h=2.0, int d=1){
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
    Vec nodes = Vec::LinSpaced(d+1, -1.0, 1.0);  // local nodes
    if(j>d){
        spdlog::error("Invalid value, index j must be in the range [0,d].");
        return 0.0;
    }
    if(std::abs(x)>1){
        spdlog::error("Invalid value, x must be in the range [-1,1].");
        return 0.0;
    }
    if(i==0){  // function
        return lagrange_polynomial(x, j, nodes);
    } else if(i==1){
        return lagrange_polynomial_first_derivative(x, j, nodes, h, d);
    } else{
        spdlog::error("Invalid value, 'i' is either 0 or 1.");
        return 0.0;
    }
}
/**
*  Return points and weights for Gauss-Legendre rules on [-1,1]. The number of 
*  points implemented are 1-4.
*
*  Numerical integration rules can be expressed in a common form:
*
*  \f[
*    \int_{-1}^{1} g(X)dX\approx \sum_{j=0}^M w_jg(\bar{X}_j)
*  \f]
*
*  where \f$\bar{X}_j\f$ are integration points, and \f$w_j\f$ are integration
*  weights, for \f$j=0,\ldots ,M\f$.  
*
*  More accurate rules, for a given \f$M\f$, arise if the location of the 
*  integration points are optimized for polynomial integrands. The
*  Gauss-Legendre quadrature constitute one such class of integration methods. 
*  Two widely applied Gauss-Legendre rules in this family have the choice:
*
*  for \f$M=1\f$:
*
*  \f[
*    \bar{X}_0 = -\frac{1}{\sqrt{3}},\quad \bar{X}_1 = \frac{1}{\sqrt{3}},
*      \quad w_0=w_1=1
*  \f]
*
*  For \f$M=2\f$:
*
*  \f[
*    \bar{X}_0 = -\sqrt{\frac{3}{5}},\quad \bar{X}_1 = 0,
*      \quad \bar{X}_2 = \sqrt{\frac{3}{5}}, \quad w_0=w_2=\frac{5}{9},
*      \qaud w_1 = \frac{8}{9}
*  \f]
*
*  These rules integrate 3rd and 5th degree polynomials exactly. In general,
*  an \f$M\f$-point Gauss-Legendre rule integrates a polynomial of degree 
*  \f$2M + 1\f$ exactly.
*
*  @param n The number of integration points.
*  @return The points and weights.
*/
template <typename T>
std::vector<std::vector<T>> GaussLegendre(int n){
    std::vector<T> points;
    std::vector<T> weights;
    if(n == 1){
        points = {0.0};
        weights = {2.0};
    } else if(n == 2){
        points = {-0.5773502691896257645091488, 
                   0.5773502691896257645091488};
        weights = {1.0, 1.0};
    } else if(n == 3){
        points = {-0.7745966692414833770358531, 
                  0.0,
                  0.7745966692414833770358531};
        weights = {0.5555555555555555555555556,
                   0.8888888888888888888888889,
                   0.5555555555555555555555556};
    } else if(n == 4){
        points = {-0.8611363115940525752239465,
                  -0.3399810435848562648026658,
                  0.3399810435848562648026658,
                  0.8611363115940525752239465};
        weights = {0.3478548451374538573730639,
                   0.6521451548625461426269361,
                   0.6521451548625461426269361,
                   0.3478548451374538573730639};
    } else {
        spdlog::error("Gauss-Legendre: The number of points is not supported");
        return {{}};
    }
    return {points, weights};
}

}  // namespace fem
} // namespace numerical

#endif  // FINITEELEMENT_H
