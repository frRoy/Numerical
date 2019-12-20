/**
 *  @file    diffusion.hpp
 *  @brief   Collections of diffusion problem, definitions.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <string>
#include <Eigen/Core>
#include "spdlog/spdlog.h"
#include "numerical/fdm/fdm.hpp"

namespace bench{
namespace diffusion{

typedef Eigen::VectorXd Vec;
typedef Eigen::ArrayXd Arr;

/**
* Computes the 1D diffusion problem defined in @ref reference_diffusion_a() 
* using the finite difference method. Here we solve the problem using the
* \f$\Theta\f$-scheme in time and a centered difference scheme in space, i.e.
*
* \f[
*    u_i^{n+1} - \Theta\left(\alpha\left(
*        \frac{u_{i-1}^{n+1}-2u_{i}^{n+1}+u_{i+1}^{n+1}}{\Delta x^2}
*        \right)\right) = \left(1-\Theta\right)\left(\alpha\left(
*        \frac{u_{i-1}^{n}-2u_{i}^{n}+u_{i+1}^{n}}{\Delta x^2}
*        \right)\right) + \Theta\Delta t f_i^{n+1}+\left(1-\Theta\right)f_i^n
*        + u_i^n
* \f]
*
* For \f$\Theta=0\f$ we have the forward Euler scheme in time, for 
* \f$\Theta=1\f$ we have the backward Euler scheme in time and for 
* \f$\Theta=1/2\f$ we have the Crank-Nicolson scheme in time.
*
* In matrix form we have:
*
* \f[
*    \mathbf{A}\mathbf{u} = \mathbf{b}
* \f]
*
* Using the Fourrier coefficient \f$F=\alpha\Delta t/\Delta x^2\f$, the matrix
* entries are:
*
* \f[
*    A_{i,i-1}=-F\Theta,~A{i, i}=1+2F\Theta,~A_{i, i+1}=-F\Theta
* \f].
*
* The vector \f$\mathbf{u}\f$ represents the unknown at the discrete time 
* \f$t=t^{n+1}\f$.
*
* The RHS can be written as:
*
* \f[
*    b_i = u_i^n + F\left(1-\Theta\right)u_{i+1}^n-2u_i^n+u_{i-1}^n +
*        \Delta t \Theta f_i^{n+1} + \Delta t \left(1-\Theta\right)f_i^n
* \f]
*
* For more information see pp. 181-183 of:
*
* R. J. Leveque, Finite Difference Methods for Ordinary and Partial 
* Differential Equations: Steady-State and Time-Dependent Problems, SIAM, 2007.
*
* @return The solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
bool fdm_diffusion_a();

class FDMDiffusionA : public numerical::fdm::Problem<double>{
public:
    FDMDiffusionA(numerical::fdm::Parameters<double>* p): 
        numerical::fdm::Problem<double>(p){
        }
    double left(double y, double z, double t);
    double right(double y, double z, double t);
    double source(const Eigen::Matrix<double, 3, 1>& x, double t);
    /**
    * Reference solution for the 1D diffusion problem a, defined as:
    *
    * \f[
    *    \frac{\partial u}{\partial t} = -\alpha \frac{\partial^2 u}{\partial x^2}
    *                                    + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with Dirichlet boundary condition \f$u(0, t) = 0\f$ and  \f$ u(L, t) = 0\f$
    * for \f$t>0\f$.
    *
    * The diffusion coefficient, \f$ \alpha(L, t) = 1\f$ is constant and uniform.
    *
    * We define the source term and initial condition with a reference 
    * (manufactured) solution that respects the boundary conditions at the 
    * extremities of the interval:
    *
    * \f[
    *    u(x, t) = 5tx\left(L-x\right)
    * \f]
    *
    * We get the source term substituting the manufactured solution in the 
    * diffusion equation, i.e.:
    *
    * \f[ 
    *    f(x, t) = 5x\left(L-x\right)+10\alpha t
    * \f]
    *
    * The initial condition is simply set to the manufactured solution at 
    * \f$t=0\f$:
    *
    * \f[ 
    *    u(x, 0) = u_0 = 0
    * \f] 
    *
    * @param x The spatial mesh node coordinates.
    * @param t The temporal mesh node.
    * @return The reference solution.
    */
    double reference(const Eigen::Matrix<double, 3, 1>& x, double t);
};

}  // namespace diffusion
}  // namespace bench

#endif  // DIFFUSION_H
