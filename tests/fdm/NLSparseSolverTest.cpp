/**
 *  @file    NLSparseSolverTest.cpp
 *  @brief   Test the NL SparseSolve class @see numerical::fdm::NLSparseSolver.
 *  @author  Francois Roy
 *  @date    01/03/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/NLSparseSolver.hpp"


namespace
{

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef Eigen::Triplet<double> T;

T trip = T(1,1,3.0);

void func(numerical::fdm::SparseSolver<double> *xyz) { 
    xyz->solve(); 
}

bool test_func(){
    numerical::fdm::Parameters<double>* params = 
        new numerical::fdm::Parameters<double>();
    params->lengths = {{0.0, 1.0}};
    params->n = {10};
    numerical::fdm::Problem<double>* p = 
        new numerical::fdm::Problem<double>(params);
    numerical::fdm::NLSparseSolver<double>* s = new
        numerical::fdm::NLSparseSolver<double>(p);
    func(s);  //  should run NL version
    delete params;
    delete p;
    delete s;
    return true;
}


}  // end namespace

TEST_CASE( "NLSparseSolver are computed", "[NL sparse solver]" )
{
  REQUIRE( trip.value() == Approx(3.0) );
  CHECK(test_func());
}
