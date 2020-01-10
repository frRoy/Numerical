/**
 *  @file    FVMeshTest.cpp
 *  @brief   Test the FVMesh class.
 *  @author  Francois Roy
 *  @date    01/10/2020
 */
#include <catch2/catch.hpp>
#include "numerical/fvm/FVMesh.hpp"


namespace 
{
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;
typedef numerical::fvm::FVMesh<double> Mesh;

/**
* Return the space increment.
*
* @param mesh The mesh instance.
* @return The space increment.
*/
std::vector<double> dx(numerical::Mesh<double> *mesh) { 
    return {mesh->dx()}; 
}
/**
* Get the vertices.
*
* @return The vertices list.
*/
std::vector<double> test_vertices(numerical::Mesh<double> *mesh) {
    std::vector<double> out;
    Vec v = mesh->vertices();
    for(int i=0; i<v.size(); i++){
        out.push_back(v[i]);
    }
    return out;
}
/**
* Pass the FVMesh to a function that gets the base class as argument. The
* function should use the member function of the FVMesh.
*
* The default discretization scheme is node centered.
*
* @param test_case The thest case.
* @return The actual value of the tested member function.
*/
std::vector<double> test_derived(int test_case){
    std::vector<double> lengths = {0.0, 1.0};
    double t0 = 0.0, tend = 1.0;
    int nx = 10, nt = 10;
    // node centered by default
    Mesh* s = new Mesh(lengths, t0, tend, nx, nt);
    std::vector<double> out = {999.0};
    if(test_case == 0){
        out = test_vertices(s);
    }
    delete s;
    return out;
}
/**
* Same as above but for the cell centered scheme.
*
* @param test_case The thest case.
* @return The actual value of the tested member function.
*/
std::vector<double> test_derived_cell_centered(int test_case){
    std::vector<double> lengths = {0.0, 1.0};
    double t0 = 0.0, tend = 1.0;
    int nx = 10, nt = 10;
    Mesh* s = new Mesh(lengths, t0, tend, nx, nt, 1);
    std::vector<double> out = {999.0};
    if(test_case == 0){
        out = test_vertices(s);
    }
    delete s;
    return out;
}
/**
* Throw esception for non-existent scheme.
*
* @return The actual value of the tested member function.
*/
void test_non_existing_scheme(){
    std::vector<double> lengths = {0.0, 1.0};
    double t0 = 0.0, tend = 1.0;
    int nx = 10, nt = 10;
    Mesh* s = new Mesh(lengths, t0, tend, nx, nt, 3);
    delete s;
}

}  // namespace

TEST_CASE( "FVM: Mesh scheme tests", "[mesh, cells]" )
{
    CHECK_THROWS_WITH( test_non_existing_scheme(), 
      Catch::Contains( "Not a valid" ) && 
      Catch::Contains( "discretization scheme" ) );
}
TEST_CASE( "FVM: Mesh vertices tests.", "[FEMesh, vertices]" )
{
    std::vector<double> vertices = test_derived(0);
    std::vector<double> ref_v = {0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 
        0.75, 0.85, 0.95};
    for(int i=0; i<vertices.size(); i++){
        CHECK(vertices[i] == Approx(ref_v[i]));
    }
    // cell centered
    vertices = test_derived_cell_centered(0);
    ref_v = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 
        0.7, 0.8, 0.9, 1.0};
    for(int i=0; i<vertices.size(); i++){
        CHECK(vertices[i] == Approx(ref_v[i]));
    }
}
