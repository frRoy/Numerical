/**
 *  @file    FEMeshTest.cpp
 *  @brief   Test the finite element Mesh class.
 *  @author  Francois Roy
 *  @date    01/06/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fem/FEMesh.hpp"


namespace
{
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef numerical::fem::FEMesh<double> Mesh;
typedef  double (numerical::fem::FEMesh<double>::*fctptr)(double x);

/**
* Get the vertices.
*
* @return The vertices list.
*/
Vec test_vertices() {
    Mesh m = Mesh({0.0, 1.0}, 0.1, 1.1, 10, 22);
    return m.vertices();
}
/**
* Get the cells.
*
* @return The cells list.
*/
std::vector<std::vector<int>> test_cells() {
    Mesh m = Mesh({0.0, 1.0}, 0.1, 1.1, 10, 22);
    return m.cells();
}
/**
* Get the dof maps.
*
* @param d Polynomial order.
* @return The dof maps list.
*/
std::vector<std::vector<int>> test_dof_map(int d) {
    Mesh m = Mesh({0.0, 1.0}, 0.1, 1.1, 10, 22, d);
    return m.dof_map();
}

}  // namespace

TEST_CASE( "FEM: Mesh vertices tests.", "[FEMesh, vertices]" )
{
    Vec vertices = test_vertices();
    std::vector<double> ref_v = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 
                                 0.9, 1.0};
    for(int i=0; i<vertices.size(); i++){
        CHECK(vertices[i] == Approx(ref_v[i]));
    }
}
TEST_CASE( "FEM: Mesh cells tests.", "[FEMesh, cells]" )
{   
    std::vector<std::vector<int>> cells = test_cells();
    std::vector<std::vector<int>> ref = {{0, 1}, {1, 2}, {2, 3}, {3, 4},
                                         {4, 5}, {5, 6}, {6, 7}, {7, 8},
                                         {8, 9}, {9, 10}};
    for(int i=0; i<cells.size(); i++){
        CHECK(cells[i] == ref[i]);
    }
}
/*
TEST_CASE( "FEM: dof map tests", "[FEMesh, dof_map]" )
{
    // degree = 0
    std::vector<std::vector<int>> dof_map = test_dof_map(0);
    std::vector<std::vector<int>> ref = {{0}, {1}, {2}, {3}, {4}, {5}, {6}, 
                                         {7}, {8}, {9}, {10}, {11}};
    for(int i=0; i<dof_map.size(); i++){
        CHECK(dof_map[i] == ref[i]);
    }
    // degree = 1
    dof_map = test_dof_map(1);
    ref = {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8},
           {8, 9}, {9, 10}};
    for(int i=0; i<dof_map.size(); i++){
        CHECK(dof_map[i] == ref[i]);
    }
    // degree = 2
    dof_map = test_dof_map(2);
    ref = {{0, 1, 2}, {2, 3, 4}, {4, 5, 6}, {6, 7, 8}, {9, 10, 11},
           {11, 12, 13}, {13, 14, 15}, {15, 16, 17}, {17, 18, 19}, 
           {19, 20, 21}};
    for(int i=0; i<dof_map.size(); i++){
        CHECK(dof_map[i] == ref[i]);
    }
}*/
