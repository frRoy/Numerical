/**
 *  @file    Problem.hpp
 *  @brief   Define a finite difference problem.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef PROBLEM_H
#define PROBLEM_H

// #include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/SparseCore>
#include "spdlog/spdlog.h"
#include "Parameters.hpp"
#include "Mesh.hpp"

namespace numerical {

namespace fdm {

/**
 * Defines the finite difference problem over the hypercube.
 * This class only define the 1D diffusion problem with Dirichlet/Neumann 
 * boundary conditions and heterogenous diffusion coefficient (for now).
 */
template <typename T>
class Problem {
typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vec;
typedef std::vector<Eigen::Matrix<T, 3, 1>> Coord;
// function pointer to member functions
typedef  T (Problem<T>::*fctptr)(T x1, T x2, T time);
private:
    int m_dim; // space dimenions of the problem 
    Vec m_t;
    std::vector<T> m_dx;
    T m_dt;
protected:
    Parameters<T>* m_params;  // default parameters
    Mesh<T>* m_mesh;  // mesh
    int m_bc_types[6];  // boundary types (default = Dirichlet)
    Coord m_coordinates;
    // array of function pointers to member functions
    fctptr m_functions[6] = {&Problem<T>::left, &Problem<T>::right,
                             &Problem<T>::bottom, &Problem<T>::top,
                             &Problem<T>::back, &Problem<T>::front};
public:
    Problem(Parameters<T>* params): 
      m_params(params),
      m_bc_types{0, 0, 0, 0, 0, 0} {
        // define other variables variables
        m_dim = m_params->lengths.size();
  	    // define mesh
  	    m_mesh = new Mesh<T>(m_params->lengths, m_params->t0, m_params->tend, 
          m_params->n, m_params->nt);
        m_coordinates = m_mesh->coordinates();  // constant
        m_t = m_mesh->t();
        m_dx = m_mesh->dx();
        m_dt = m_mesh->dt();
    }
    virtual ~Problem(){
        delete m_mesh;
    }

    /**
    * Diffusion coefficient value. The default is a constant obtained from 
    * m_params.
    *
    * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
    *    node.
    * @param t The discrete time.
    * @return The diffusion coefficient value at a specified mesh location.
    */
    virtual T alpha(const Eigen::Matrix<T, 3, 1>& x, T t){
  	    return m_params->alpha;
    }
    /**
    *  @return The mesh coordinates.
    */
    std::vector<Eigen::Matrix<T, 3, 1>> coordinates(){
        return m_coordinates;
    }
    /**
    * Set the boundary types, with the following code:
    *
    *  0 = Dirichlet
    *  1 = Neumann
    *
    * For the six boundaries in the following order:
    *
    * 0- letf
    * 1- right
    * 2- bottom
    * 3- top
    * 4- back
    * 5- front
    *
    *  @param types The type of boundary conditions for all boundaries.
    */
    void bc_types(const int (&types)[6]){
        for (int i=0; i<6; i++){
            m_bc_types[i] = types[i];
            // spdlog::info("{}", m_bc_types[i]);
        }
    }
    /**
    *
    * @param boundary The boundary selection.
    * @return The boundary type of the selected boundary,
    */
    int bc_type(int boundary){
        return m_bc_types[boundary];
    }
    /**
    * User defined function: back boundary.
    * 
    * Back boundary value, i.e. at z = length[2][0]. Only for 3D models.
    *
    * @param x The \f$x\f$-coordinate.
    * @param y The \f$y\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T back(T x, T y, T t){
        return 0.0;
    }
    /**
    * User defined function: bottom boundary.
    *
    * Bottom boundary value, i.e. at y = lengths[1][0]. Only for 2D and 
    * 3D models.
    *
    * @param x The \f$x\f$-coordinate.
    * @param z The \f$z\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T bottom(T x, T z, T t){
        return 0.0;
    }
    /**
    * Set the value of the diagonals of the coefficient matrix at the indices 
    * of the boundary nodes.
    *
    *  @param dia The main diagonal.
    *  @param lower The lower diagonal.
    *  @param upper The upper diagonal.
    *  @param lower_a The second lower diagonal.
    *  @param upper_a The second upper diagonal.
    *  @param lower_b The third lower diagonal.
    *  @param upper_b The third upper diagonal.
    *  @param t The discrete time.
    */
    void coeffs_bc(Vec& dia, Vec& lower, Vec& upper, Vec& lower_a, 
        Vec& upper_a, Vec& lower_b, Vec& upper_b, T t=0.0){
    }
    /**
    * @return The spatial dimension of the problem.
    */
    int dim(){
        return m_dim;
    }
    /**
    * @return the constant time increment.
    */
    T dt() {
        return m_dt;
    }
    /**
    * @return the spatial increment.
    */
    std::vector<T> dx() {
        return m_dx;
    }
    /**
    * User defined function: front boundary.
    *
    * Front boundary value, i.e. at z = length[2][0]. Only for 3D models.
    *
    * @param x The \f$x\f$-coordinate.
    * @param y The \f$y\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T front(T x, T y, T t){
        return 0.0;
    }
    /**
    * User defined function: initial  value.
    *
    * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
    *    node.
    * @return The initial value at a specified mesh location.
    */
    virtual T initial_value(const Eigen::Matrix<T, 3, 1>& x){
        return 0.0;
    }
    /**
    * User defined function: left boundary.
    *
    * Left boundary value, i.e. for x = lengths[0][0].
    *
    * @param y The \f$y\f$-coordinate.
    * @param z The \f$z\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T left(T y, T z, T t){
        return 0.0;
    }
    /**
    * @return The total number of spatial mesh nodes.
    */
    int n() {
        int out = m_params->n[0] + 1;
        if (m_dim == 2){
            out = (m_params->n[0] + 1)*(m_params->n[1] + 1);
        }
        if (m_dim == 3){
            out = (m_params->n[0] + 1)*(m_params->n[1] + 1)*
                (m_params->n[2] + 1);
        }
        return out;
    }
    /**
    *  @return The number of division along the \f$x\f$-axis.
    */
    int n_x(){
        return m_params->n[0]; 
    }
    /**
    * @return The number of division along the \f$y\f$-axis.
    */
    int n_y(){
        return m_params->n[1]; 
    }
    /**
    * @return The number of division along the \f$z\f$-axis.
    */
    int n_z(){
        return m_params->n[2]; 
    }
    /**
    * @return The total number of discrete time steps.
    */
    int n_t(){
        return m_params->nt;
    }
    /**
    *  Get the reference solution at a specified mesh location.
    *
    * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
    *    node.
    * @param t The discrete time.
    * @return The reference solution at a specified mesh location.
    */
    virtual T reference(const Eigen::Matrix<T, 3, 1>& x, T t){
        return 0.0;
    }
    /**
    *  Set the value of the RHS vector at the indices of the boundary
    *  nodes.
    *
    *  @param rhs The rhs vector.
    *  @param u_n The solution vector at previous time step.
    *  @param alpha The rhs vector.
    *  @param f_n The source vector at privious time step.
    *  @param f The source vector.
    *  @param dx The space increment in the \f$x\f$-direction.
    *  @param dy The space increment in the \f$y\f$-direction.
    *  @param dz The space increment in the \f$z\f$-direction.
    *  @param dt The time increment.
    *  @param theta The scheme coefficient.
    *  @param t The discrete time.
    */
    void rhs_bc(Vec& rhs, Vec& u_n, const Vec& alpha, const Vec& f_n, 
                const Vec& f, const T dx, const T dy, const T dz,
                const T dt, const T theta, const T t){
        // boundary node addresses
        const std::vector<int>& left_nodes = m_mesh->left();
        const std::vector<int>& right_nodes = m_mesh->right();
        const std::vector<int>& bottom_nodes = m_mesh->bottom();
        const std::vector<int>& top_nodes = m_mesh->top();
        const std::vector<int>& back_nodes = m_mesh->back();
        const std::vector<int>& front_nodes = m_mesh->front();
        // array of addresses to specific boundry nodes 
        std::vector<int> nodes[6] = {
            left_nodes, right_nodes,
            bottom_nodes, top_nodes,
            back_nodes, front_nodes
            };
        // indices of the coordinates for each function
        // 0 for x, 1 for y, 2 for z
        int ind[6][2] = {{1,2},{1,2},{0,2},{0,2},{0,1},{0,1}};
        T d[6] = {dt/dx/dx, dt/dx/dx, dt/dy/dy, dt/dy/dy, dt/dz/dz, dt/dz/dz};
        T d_i[6] = {dx, dx, dy, dy, dz, dz};
        // the number of boundary to set up depend on the dimension
        int n_bnd = 2 * m_dim;
        const Coord& c = coordinates();
        for(int bnd=0; bnd<n_bnd; bnd++){
            if (m_bc_types[bnd] == 0){  // Dirichlet
                for(int i=0; i<nodes[bnd].size(); i++){
                    // the rhs value is equal to the value returned by the
                    // boundary function
                    rhs[nodes[bnd][i]] = (this->*m_functions[bnd]) (
                        c[nodes[bnd][i]][ind[bnd][0]], 
                        c[nodes[bnd][i]][ind[bnd][1]], 
                        t);
                }
            } else if(m_bc_types[bnd] == 1){  // Neumann
                // TODO: fix corners for 2D and 3D and fix edges for 3D
                T sign = 1.0, g, g_n, alpha_m, alpha_p;
                int k, side = 1;
                for(int i=0; i<nodes[bnd].size(); i++){  // scaled by 1/2
                    k = nodes[bnd][i];
                    g = (this->*m_functions[bnd]) (
                            c[k][ind[bnd][0]], 
                            c[k][ind[bnd][1]], 
                            t + dt);
                    g_n = (this->*m_functions[bnd]) (
                            c[k][ind[bnd][0]], 
                            c[k][ind[bnd][1]], 
                            t);
                    if (bnd % 2){  // right, top, front
                        sign = -1.0;
                        side = -1;
                        // alpha_m = 0.5 * (alpha[k] + alpha[k-1]);
                        // get value on the right of the boundary
                        // alpha_p = 0.5 * (alpha[k] + alpha[k]);
                    } else{  // left, bottom, back
                        sign = 1.0;
                        side = 1;
                        // get value on the left of the boundary
                        // alpha_m = 0.5 * (alpha[k] + alpha[k]);
                        // alpha_p = 0.5 * (alpha[k] + alpha[k+1]);
                    }
                    // TODO fix value of alpha
                    rhs[k] = 0.5 * (
                        u_n[k] + 
                        theta*(sign*d[bnd]*alpha[k]*g*2.0*d_i[bnd] + dt*f[k]) +
                        (1.0-theta)*(2*alpha[k]*d_i[bnd]*(u_n[k+side]-u_n[k])
                            +sign*d[bnd]*alpha[k]*g_n*2.0*d_i[bnd] + dt*f_n[k])
                        );
                }
            } else{
                throw std::invalid_argument(
                "Not a valid boundary condition type.");
            }
        }
    }
    /**
    * User defined function: right boundary.
    *
    * Right boundary value, i.e. for x = lengths[0][1].
    *
    * @param y The \f$y\f$-coordinate.
    * @param z The \f$z\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T right(T y, T z, T t){
        return 0.0;
    }
    
    /**
    *  Get the computed solution at a specified mesh location.
    *
    * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
    *    node.
    * @param t The discrete time.
    * @return The computed solution at a specified mesh location.
    */
    virtual T solution(const Eigen::Matrix<T, 3, 1>& x, T t){
        return 0.0;
    }
    /**
    * User defined function: source term.
    *
    * @param x The \f$x\f$-, \f$y\f$-, and \f$z\f$-coordinates of the mesh 
    *    node.
    * @param t The discrete time.
    * @return The source term at a specified mesh location.
    */
    virtual T source(const Eigen::Matrix<T, 3, 1>& x, T t){
        return 0.0;
    }
    /**
    * @return The theta parameter for the finite difference scheme. 
    */
    T theta(){
        return m_params->theta;
    }
    /**
    *
    */
    Vec t(){
        return m_t;
    }
    /**
    * User defined function: top boundary.
    *
    * Top boundary value, i.e. at y = lengths[1][1]. Only for 2D and 3D models.
    *
    * @param x The \f$x\f$-coordinate.
    * @param z The \f$z\f$-coordinate.
    * @param t The discrete time.
    * @return The boundary value at a specified mesh location.
    */
    virtual T top(T x, T z, T t){
        return 0.0;
    }
    /**
    *
    *  @return The initial value vector.
    */
    Vec u_0(){
        int num = n();
        Vec u_0(num), x(num), y(num), z(num);
        Coord coords = coordinates();
        for(int i=0; i<num; i++) {
            x[i] = coords[i][0];
            y[i] = coords[i][1];
            z[i] = coords[i][2];
            u_0[i] = initial_value(coords[i]);
        } 
        return u_0;
    }
};

}  // namespace fdm

}  // namespace numerical

#endif  // PROBLEM_H
