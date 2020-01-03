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
* Computes the 1D diffusion problem defined in @ref FDMDiffusionA::reference() 
* using the finite difference method and computes the L2-norm of the error
* between the exact and computed solution.
*
* @return The L2-norm of the error between the computed and exact solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
double fdm_diffusion_a();
/**
* Computes the 1D diffusion problem defined in @ref FDMDiffusionB::reference() 
* using the finite difference method and computes the L2-norm of the error
* between the exact and computed solution.
*
* @return The L2-norm of the error between the computed and exact solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
double fdm_diffusion_b();
/**
* Computes the 1D diffusion problem defined in @ref FDMDiffusionC::reference() 
* using the finite difference method and computes the L2-norm of the error
* between the exact and computed solution.
*
* @return The L2-norm of the error between the computed and exact solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
double fdm_diffusion_c();
/**
* Computes the 1D diffusion problem defined in @ref FDMDiffusionD::reference() 
* using the finite difference method and computes the L2-norm of the error
* between the exact and computed solution.
*
* @return The L2-norm of the error between the computed and exact solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
double fdm_diffusion_d();
/**
* Computes the 1D diffusion problem defined in @ref FDMDiffusionE::reference() 
* using the finite difference method and computes the L2-norm of the error
* between the exact and computed solution.
*
* @return The L2-norm of the error between the computed and exact solution.
* @see numerical::fdm::Parameters
* @see numerical::fdm::Mesh
* @see numerical::fdm::SparseSolver
*/
double fdm_diffusion_e();

/**
* This class provides tools to compute the finite difference problem 
* defined in @ref FDMDiffusionA::reference().
*
* The class is derived form the @ref numerical::fdm::Problem class and 
* overrides the member function @ref numerical::fdm::Problem::left(),
* @ref numerical::fdm::Problem::right(), 
* @ref numerical::fdm::Problem::initial_value(),
* @ref numerical::fdm::Problem::source(), and 
* @ref numerical::fdm::Problem::reference().
*/
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
    *    \frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}
    *                                 {\partial x^2}
    *                               + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with homogeneous Dirichlet boundary condition \f$u(0, t) = 0\f$ and  
    * \f$ u(L, t) = 0\f$ for \f$t>0\f$.
    *
    * The diffusion coefficient, \f$ \alpha(L, t)\f$ is constant and 
    * uniform over the interval (line).
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

/**
* This class provides tools to compute the finite difference problem 
* defined in @ref FDMDiffusionB::reference().
*
* The class is derived form the @ref numerical::fdm::Problem class and 
* overrides the member function @ref numerical::fdm::Problem::left(),
* @ref numerical::fdm::Problem::right(), 
* @ref numerical::fdm::Problem::initial_value(),
* @ref numerical::fdm::Problem::source(), and 
* @ref numerical::fdm::Problem::reference().
*/
class FDMDiffusionB : public numerical::fdm::Problem<double>{
public:
    FDMDiffusionB(numerical::fdm::Parameters<double>* p): 
        numerical::fdm::Problem<double>(p){
        }
    double left(double y, double z, double t);
    double right(double y, double z, double t);
    double source(const Eigen::Matrix<double, 3, 1>& x, double t);
    /**
    * Reference solution for the 1D diffusion problem b, defined as:
    *
    * \f[
    *    \frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}
    *                                 {\partial x^2}
    *                                 + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with homogeneous Neumann boundary condition 
    * \f$\left. \frac{\partial}{\partial x}u(x, t)\right|_{x=0} = 0\f$ and  
    * \f$\left. \frac{\partial}{\partial x}u(x, t)\right|_{x=L} = 0\f$.
    *
    * The diffusion coefficient, \f$ \alpha(L, t)\f$ is constant and 
    * uniform over the interval (line).
    *
    * We define the source term and initial condition with a reference 
    * (manufactured) solution that respects the boundary conditions at the 
    * extremities of the interval:
    *
    * \f[
    *    \frac{\partial}{\partial x}u(x, t) = 5tx\left(L-x\right)
    * \f]
    *
    * which leads to
    *
    * \f[
    *    u(x, t) = 5tx\left(\frac{Lx}{2}-\frac{x^2}{3}\right)
    * \f]
    *
    * We get the source term substituting the manufactured solution in the 
    * diffusion equation, i.e.:
    *
    * \f[ 
    *    f(x, t) = 5x\left(\frac{Lx}{2}-\frac{x^2}{3}\right)-
    *      5\alpha t\left(L-2x\right)
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

/**
* This class provides tools to compute the finite difference problem 
* defined in @ref FDMDiffusionC::reference().
*
* The class is derived form the @ref numerical::fdm::Problem class and 
* overrides the member function @ref numerical::fdm::Problem::left(),
* @ref numerical::fdm::Problem::right(), 
* @ref numerical::fdm::Problem::initial_value(),
* @ref numerical::fdm::Problem::source(), and 
* @ref numerical::fdm::Problem::reference().
*/
class FDMDiffusionC : public numerical::fdm::Problem<double>{
public:
    FDMDiffusionC(numerical::fdm::Parameters<double>* p): 
        numerical::fdm::Problem<double>(p){
        }
    double left(double y, double z, double t);
    double right(double y, double z, double t);
    double source(const Eigen::Matrix<double, 3, 1>& x, double t);
    /**
    * Reference solution for the 1D diffusion problem c, defined as:
    *
    * \f[
    *    \frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}
    *                                 {\partial x^2}
    *                                 + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with non-homogeneous Neumann boundary condition 
    * \f$\mathbf{\hat{n}}\cdot\left. \frac{\partial}{\partial x}
    *  u(x, t)\right|_{x=0} = -5tL\f$ and  
    * \f$\mathbf{\hat{n}}\cdot\left. \frac{\partial}{\partial x}
    *  u(x, t)\right|_{x=L} = -5tL\f$.
    *
    * The diffusion coefficient, \f$ \alpha(L, t)\f$ is constant and 
    * uniform over the interval (line).
    *
    * We define the source term and initial condition with a reference 
    * (manufactured) solution that respects the boundary conditions at the 
    * extremities of the interval:
    *
    * \f[
    *    \frac{\partial}{\partial x}u(x, t) = 5t\left(L-2x\right)
    * \f]
    *
    * which leads to
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

/**
* This class provides tools to compute the finite difference problem 
* defined in @ref FDMDiffusionD::reference().
*
* The class is derived form the @ref numerical::fdm::Problem class and 
* overrides the member function @ref numerical::fdm::Problem::left(),
* @ref numerical::fdm::Problem::right(), 
* @ref numerical::fdm::Problem::initial_value(),
* @ref numerical::fdm::Problem::source(), and 
* @ref numerical::fdm::Problem::reference().
*/
class FDMDiffusionD : public numerical::fdm::Problem<double>{
public:
    FDMDiffusionD(numerical::fdm::Parameters<double>* p): 
        numerical::fdm::Problem<double>(p){
        }
    double left(double y, double z, double t);
    double right(double y, double z, double t);
    double source(const Eigen::Matrix<double, 3, 1>& x, double t);
    /**
    * Reference solution for the 1D diffusion problem d, defined as:
    *
    * \f[
    *    \frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}
    *                                 {\partial x^2}
    *                                 + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with non-homogeneous Dirichlet boundary condition 
    * \f$u(0, t) = 5tL^2/4\f$ and 
    * \f$\mathbf{\hat{n}}\cdot\left. \alpha(u)\frac{\partial}{\partial x}
    *  u(x, t)\right|_{x=L} = g\f$.
    *
    * The diffusion coefficient, \f$ \alpha(u)=u\f$.
    *
    * We define the source term and initial condition with a reference 
    * (manufactured) solution that respects the boundary conditions at the 
    * extremities of the interval:
    *
    * \f[
    *    u(x, t) = 5t\left(x-\frac{L}{2}\right)^2
    * \f]
    *
    * We get the source term substituting the manufactured solution in the 
    * diffusion equation, i.e.:
    *
    * \f[ 
    *    f(x, t) = 5\left(x-\frac{L}{2}\right)^2-
    *      10\alpha t
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

/**
* This class provides tools to compute the finite difference problem 
* defined in @ref FDMDiffusionE::reference().
*
* The class is derived form the @ref numerical::fdm::Problem class and 
* overrides the member function @ref numerical::fdm::Problem::left(),
* @ref numerical::fdm::Problem::right(), 
* @ref numerical::fdm::Problem::initial_value(),
* @ref numerical::fdm::Problem::source(), and 
* @ref numerical::fdm::Problem::reference().
*/
class FDMDiffusionE : public numerical::fdm::Problem<double>{
public:
    FDMDiffusionE(numerical::fdm::Parameters<double>* p): 
        numerical::fdm::Problem<double>(p){
        }
    double left(double y, double z, double t);
    double right(double y, double z, double t);
    double source(const Eigen::Matrix<double, 3, 1>& x, double t);
    /**
    * Reference solution for the nonlinear 1D diffusion problem e, defined as:
    *
    * \f[
    *    \frac{\partial u}{\partial t} = \frac{\partial}{\partial x}
    *        \left(\alpha(u)\frac{\partial u}{\partial x}\right)
    *        + f(x,t)\quad x \in (0, L),~t \in (0, T]
    * \f]
    *
    * with non-homogeneous mixed Dirichlet/Neumann boundary conditions 
    * \f$u(0, t) = u_D\f$ and \f$u(L, t)= 5tL^2/4\f$.
    *
    * The diffusion coefficient, \f$ \alpha(L, t)\f$ is constant and 
    * uniform over the interval (line).
    *
    * We define the source term and initial condition with a reference 
    * (manufactured) solution that respects the boundary conditions at the 
    * extremities of the interval:
    *
    * \f[
    *    u(x, t) = 5t\left(x-\frac{L}{2}\right)^2
    * \f]
    *
    * We get the source term substituting the manufactured solution in the 
    * diffusion equation, i.e.:
    *
    * \f[ 
    *    f(x, t) = 5\left(x-\frac{L}{2}\right)^2-
    *      10\alpha t
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
