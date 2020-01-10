/**
 *  @file    FVProblemTest.cpp
 *  @brief   Test the FVProblem class.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fvm/FVProblem.hpp"


namespace 
{
typedef numerical::fvm::FVProblem<double> Problem;

/**
* Pass the FVProblem to a function that gets the base class as argument. The
* function should use the member function of the FVProblem.
*
* @param test_case The thest case.
* @return The actual value of the tested member function.
*/
std::vector<double> test_derived(int test_case){
    std::vector<double> lengths = {0.0, 1.0};
    double t0 = 0.0, tend = 1.0;
    int nx = 10, nt = 10;
    Problem* s = new Problem();
    std::vector<double> out = {999.0};
    if(test_case == 0){
        // out = test_vertices(s);
    }
    delete s;
    return out;
}

}  // namespace

TEST_CASE( "FVM: Problem tests", "[FVProblem]" )
{
    test_derived(0);
    CHECK(1 == 1);
}
