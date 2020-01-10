/**
 *  @file    ParametersTest.cpp
 *  @brief   Test the Parameters class.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/Parameters.hpp"


namespace
{

typedef numerical::fdm::Parameters<double> Params;  

/**
* Set length in 1D.
*
* @return The spatial dimension through size.
*/
int test_set_lengths_1d() {
    Params params = Params();
    params.lengths = {{0.0, 1.0}}; 
    return params.lengths.size();
}

}  // end namespace

TEST_CASE( "PARAMETERS: parameters tests.", "[Parameters]" )
{
  CHECK( test_set_lengths_1d() == 1 );
}
