/**
 *  @file    ProblemTest.cpp
 *  @brief   Test the Problem class @see numerical::fdm::Problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/FDProblem.hpp"


namespace
{

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
typedef numerical::fdm::FDProblem<double> Problem;
typedef numerical::fdm::Parameters<double> Parameters;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

/**
* Get the RHS vector for the homogeneous Nuumann boundary conditions in 1D.
*
* @return The RHS vector.
*/
Vec test_rhs_homogeneous_neumann(){
	Parameters* params = new Parameters();
    params->lengths = {{0.0, 1.0}};
    params->n = {10};
    params->nt = {100};
    params->tend = 1.0;

    Problem p = Problem(params);
    int types[6] = {1, 1, 0, 0, 0, 0};
    p.bc_types(types);
    Vec rhs = Vec::Constant(p.n(), 0.0);
    Vec u_n = Vec::LinSpaced(p.n(), 0, 90); // [0, 9, 18, ...]
    Vec alpha = Vec::Constant(p.n(), 1.0);
    Vec f_n = Vec::Constant(p.n(), 0.0);
    Vec f = Vec::Constant(p.n(), 0.0);
    double dx = p.dx()[0], dy = p.dx()[1], dz = p.dx()[2], dt = p.dt();
    double theta = p.theta(), t=dt;
    delete params;
    p.rhs_bc(rhs, u_n, alpha, f_n, f, dx, dy, dz, dt, theta, t);
    return rhs;
}

}  // end namespace

TEST_CASE( "FDM:  problem tests", "[FDProblem, RHS, Neumann]" )
{
  Vec out = test_rhs_homogeneous_neumann();
  CHECK(1 == 1); 
}
