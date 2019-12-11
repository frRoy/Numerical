/**
 *  @file    MeshTest.cpp
 *  @brief   Test the @c Mesh class.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/Mesh.hpp"


namespace
{

typedef std::vector<std::vector<double>> Vect;
typedef numerical::fdm::Mesh<double> Mesh;
/**
*/
double test_lengths_1d( std::vector<std::vector<double>> lengths) {
    Mesh m = Mesh(lengths, 0.1, 1.1, {20}, 22);
    return m.dx()[0];
}
/**
*/
double test_lengths_2d( std::vector<std::vector<double>> lengths) {
    Mesh m = Mesh(lengths, 0.1, 1.1, {20, 30}, 22);
    return m.dx()[1];
}

}  // namespace

TEST_CASE( "Mesh tests are computed", "[mesh]" )
{
  REQUIRE( test_lengths_1d({{0.1, 0.2}}) == Approx(0.1/20.0) );
  REQUIRE( test_lengths_2d({{0.1, 0.2}, {0.3, 0.4}}) == Approx(0.1/30.0) );
  REQUIRE_THROWS_WITH( test_lengths_1d({{0.1, 0.2}, {0.3, 0.4}}), 
    Catch::Contains( "n and lengths" ) && Catch::Contains( "different" ) );
  REQUIRE_THROWS_WITH( test_lengths_1d({{0.1}}), 
    Catch::Contains( "item" ) && Catch::Contains( "not a 2-list" ) );
}
