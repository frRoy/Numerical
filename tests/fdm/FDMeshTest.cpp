/**
 *  @file    FDMeshTest.cpp
 *  @brief   Test the finite difference mesh.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include <catch2/catch.hpp>
#include "numerical/fdm/FDMesh.hpp"


namespace
{
typedef std::vector<std::vector<double>> Vect;
typedef numerical::fdm::FDMesh<double> Mesh;
/**
* Get the space increment in 1D.
*
* @param lengths The minimum and maximum global x-coordinates.
* @return The space increment alongthe x-direction.
*/
double test_lengths_1d( std::vector<std::vector<double>> lengths) {
    Mesh m = Mesh(lengths, 0.1, 1.1, {20}, 22);
    // m.coordinates();
    return m.dx()[0];
}
/**
* Get the space increment in 2D.
*
* @param lengths The minimum and maximum global x- and y-coordinates.
* @return The space increment alongthe y-direction.
*/
double test_lengths_2d( std::vector<std::vector<double>> lengths) {
    Mesh m = Mesh(lengths, 0.1, 1.1, {2, 3}, 22);
    // m.coordinates();
    return m.dx()[1];
}
/**
* Get the space increment in 3D.
*
* @param lengths The minimum and maximum global x- and y- and z-coordinates.
* @return The space increment alongthe z-direction.
*/
double test_lengths_3d( std::vector<std::vector<double>> lengths) {
    Mesh m = Mesh(lengths, 0.1, 1.1, {2, 3, 4}, 22);
    // m.coordinates();
    return m.dx()[2];
}
/**
* Get the chosen boundary in 1D.
*
* @param dir If dir = 0, left otherwise right.
* @return the left or right boundary node index.
*/
std::vector<int> test_bnd_1d(int dir=0) {
    Mesh m = Mesh({{0.0, 1.0}}, 0.1, 1.1, {10}, 22);
    std::vector<int> out = m.left();
    if (dir == 1){  // right
      out = m.right();
    }
    return out;
}
/**
* Get the chosen boundary in in 2D.
*
* @param dir If dir = 0, left, dir = 1, right, dir = 2, bottom otherwise top.
* @return the left, right, bottom or top boundary node index.
*/
std::vector<int> test_bnd_2d(int dir=0) {
    Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, {2, 3}, 22);
    std::vector<int> out;
    if (dir == 0) {  // left
        out = m.left();   
    }
    else if (dir == 1){  // right
      out = m.right();
    }
    else if (dir == 2){  // bottom
        out = m.bottom();
    }
    else{  // top
        out = m.top();
    } 
    return out;
}
/**
* Get the chosen boundary in in 3D.
*
* @param dir If dir = 0, left, dir = 1, right, dir = 2, bottom, dir = 3, top, 
*  dir = 4, front otherwise back.
* @return the left, right, bottom, top, front, or back boundary node index.
*/
std::vector<int> test_bnd_3d(int dir=0) {
    Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, 
        {2, 3, 4}, 22);
    std::vector<int> out;
    if (dir == 0) {  // left
        out = m.left();   
    }
    else if (dir == 1){  // right
      out = m.right();
    }
    else if (dir == 2){  // bottom
        out = m.bottom();
    }
    else if (dir == 3){  // top
        out = m.top();
    }
    else if (dir ==4){  // front
       out = m.front();
    }
    else{  // back
        out = m.back();
    }
    return out;
}
/**
* Get the boundary indices.
*
* @param dim The spatial dimension.
* @return The boundary indices.
*/
std::vector<int> test_boundaries(int dim=1){
    std::vector<int> out;
    if(dim == 1){
        Mesh m = Mesh({{0.0, 1.0}}, 0.1, 1.1, {10}, 22);
        out = m.boundaries();
    }
    else if(dim == 2){
        Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, {2, 3}, 22);
        out = m.boundaries();
    }
    else{  // dim =3
        Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, 
            {2, 3, 4}, 22);
        out = m.boundaries();
    }
    return out; 
}
/**
* Get the spatial coordinates.
*
* @param dim The spatial dimension.
* @return The spatial coordinates.
*/
std::vector<Eigen::Matrix<double, 3, 1>> test_coordinates(int dim){
    std::vector<Eigen::Matrix<double, 3, 1>> out;
    if(dim == 1){
        Mesh m = Mesh({{0.0, 1.0}}, 0.1, 1.1, {10}, 22);
        out = m.coordinates();
    }
    else if(dim == 2){
        Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, {2, 3}, 22);
        out = m.coordinates();
    }
    else{  // dim =3
        Mesh m = Mesh({{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}}, 0.1, 1.1, 
            {2, 3, 4}, 22);
        out = m.coordinates();
    }
    return out;
}

}  // namespace

TEST_CASE( "FDM: Mesh lengths tests.", "[FDMesh, lengths]" )
{
  REQUIRE( test_lengths_1d({{0.1, 0.2}}) == Approx(0.1/20.0) );
  REQUIRE( test_lengths_2d({{0.1, 0.2}, {0.3, 0.4}}) == Approx(0.1/3.0) );
  REQUIRE( test_lengths_3d({{0.0, 0.1}, {0.0, 0.1}, {0.0, 0.1}}) == 
    Approx(0.1/4.0) );
  REQUIRE_THROWS_WITH( test_lengths_1d({{0.1, 0.2}, {0.3, 0.4}}), 
    Catch::Contains( "n and lengths" ) && Catch::Contains( "different" ) );
  REQUIRE_THROWS_WITH( test_lengths_1d({{0.1}}), 
    Catch::Contains( "item" ) && Catch::Contains( "not a 2-list" ) );
}

TEST_CASE("FDM: Mesh boundary tests.", "[FDMesh, boundaries]")
{
  // TODO the order doesn't matter, sorting the vector or transforming to sets
  // would give a correct comparison.
  std::vector<int> left_1d = {0};
  std::vector<int> right_1d = {10};
  std::vector<int> left_2d = {0, 3, 6, 9};
  std::vector<int> right_2d = {2, 5, 8, 11};
  std::vector<int> bottom_2d = {0, 1, 2};
  std::vector<int> top_2d = {9, 10, 11};
  std::vector<int> left_3d = {0, 12, 24, 36, 48,
                              3, 15, 27, 39, 51,
                              6, 18, 30, 42, 54,
                              9, 21, 33, 45, 57};
  std::vector<int> right_3d = {2,  14, 26, 38, 50,
                               5,  17, 29, 41, 53,
                               8,  20, 32, 44, 56,
                               11, 23, 35, 47, 59};
  std::vector<int> bottom_3d = {0, 12, 24, 36, 48,
                                1, 13, 25, 37, 49,
                                2, 14, 26, 38, 50};
  std::vector<int> top_3d = { 9, 21, 33, 45, 57,
                             10, 22, 34, 46, 58,
                             11, 23, 35, 47, 59};
  std::vector<int> front_3d = { 48, 51, 54, 57,
                                49, 52, 55, 58, 
                                50, 53, 56, 59};
  std::vector<int> back_3d = { 0, 3, 6, 9,
                               1, 4, 7, 10, 
                               2, 5, 8, 11};
  std::vector<int> bnd_1 = {0, 10};
  std::vector<int> bnd_2 = {0, 1, 2, 3, 5, 6, 8, 9, 10, 11};
  std::vector<int> bnd_3 = {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                             11, 12, 13, 14, 15, 17, 18, 20, 21, 22, 23, 
                             24, 25, 26, 27, 29, 30, 32, 33, 34, 35, 36, 
                             37, 38, 39, 41, 42, 44, 45, 46, 47, 48, 49, 
                             50, 51, 52, 53, 54, 55, 56, 57, 58, 59 };
  CHECK(test_bnd_1d() == left_1d);
  CHECK(test_bnd_1d(1) == right_1d);
  CHECK(test_bnd_2d() == left_2d);
  CHECK(test_bnd_2d(1) == right_2d);
  CHECK(test_bnd_2d(2) == bottom_2d);
  CHECK(test_bnd_2d(3) == top_2d);
  CHECK(test_bnd_3d() == left_3d);
  CHECK(test_bnd_3d(1) == right_3d);
  CHECK(test_bnd_3d(2) == bottom_3d);
  CHECK(test_bnd_3d(3) == top_3d);
  CHECK(test_bnd_3d(4) == front_3d);
  CHECK(test_bnd_3d(5) == back_3d);
  CHECK(test_boundaries(1) == bnd_1);
  CHECK(test_boundaries(2) == bnd_2);
  CHECK(test_boundaries(3) == bnd_3);
}
TEST_CASE("FDM: Mesh coordinates tests.", "[FDMesh, coordinates]")
{  
    std::vector<Eigen::Matrix<double, 3, 1>> coords;
    std::vector<double> c1 = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 
                              0.6, 0.7, 0.8, 0.9, 1.0};
    coords = test_coordinates(1);
    // spdlog::info("{}", coords[0][1]);
    for (int i=0; i<coords.size(); i++){
        CHECK(coords[i][0] == Approx(c1[i]));
    }
}
