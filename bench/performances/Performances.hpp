/**
 *  @file    Performances.hpp
 *  @brief   Check performances for different methods.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#ifndef PERFORMANCES_H
#define PERFORMANCES_H

#include <vector>
#include <set>
#include <array>
#include <algorithm>
#include <iterator>
#include <Eigen/SparseCore>
#include "spdlog/spdlog.h"

namespace bench {

namespace performances {

const unsigned int n = 10000;
const Eigen::VectorXd u = Eigen::VectorXd::Random(n + 1);

/*!
* Computes a simple looped difference.
* @return The computed vector.
*/
Eigen::VectorXd loop();
/*!
* Computes a vectorized looped difference.
* @return The computed vector.
*/
Eigen::VectorXd vectorized_loop();
/*!
*  Initialize vector: left boundary of a rectangular grid
*
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The vector containing the grid nodes belonging to the left boundary.
*/
std::vector<int> left_v(int nx, int ny);
/*!
*  Initialize vector: bottom boundary of a rectangular grid
*
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The vector containing the grid nodes belonging to the bottom 
*  boundary.
*/
std::vector<int> bottom_v(int nx, int ny);
/*!
*  Intersection of vectors: bottom left corner
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return the bottom-left corner (the node at the intersection bewteen the 
*  left and bottom lines)
*/
int bottom_left_corner_v(int nx, int ny);
/*!
*  Same as above but with a different algorithm.
*/
int bottom_left_corner_v_1(int nx, int ny);
/*!
*  Union of vectors: bottom and left lines.
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return the nodes at the union bewteen the left and bottom lines.
*/
std::vector<int> union_v(int nx, int ny);
/*!
*  Initialize set: left boundary of a rectangular grid
*
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The vector containing the grid nodes belonging to the left boundary.
*/
std::set<int> left_s(int nx, int ny);
/*!
*  Initialize set: bottom boundary of a rectangular grid
*
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The vector containing the grid nodes belonging to the bottom 
*  boundary.
*/
std::set<int> bottom_s(int nx, int ny);
/*!
*  Intersection of sets: bottom left corner
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return the bottom-left corner (the node at the intersection bewteen the 
*  left and bottom lines)
*/
int bottom_left_corner_s(int nx, int ny);
/*!
*  Union of sets: bottom and left lines.
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return the nodes at the union bewteen the left and bottom lines.
*/
std::vector<int> union_s(int nx, int ny);
/*!
*  Change value of a vector u at given indices with a set.
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The new vector.
*/
Eigen::VectorXd uleft_s(int nx, int ny);

/*!
*  Change value of a vector u at given indices with a vector.
*  @param nx The number of division along the x-axis.
*  @param ny The number of dicision along the y-axis.
*  @return The new vector.
*/
Eigen::VectorXd uleft_v(int nx, int ny);

}  // namespace performances
}  // namespace bench

#endif  // PERFORMANCES_H