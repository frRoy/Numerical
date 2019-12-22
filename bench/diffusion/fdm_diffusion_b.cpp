/**
 *  @file fdm_diffusion_b.cpp
 *  @brief Run a 1D diffusion problem with Neumann boundary conditions.
 *  @author Francois Roy
 *  @date 12/20/2019
 */
#include "spdlog/spdlog.h"
#include "numerical/fdm/Parameters.hpp"
#include "diffusion.hpp"

namespace bench {
namespace diffusion {
typedef numerical::fdm::Parameters<double> Params;
typedef numerical::fdm::Problem<double> Problem;

double FDMDiffusionB::left(double y, double z, double t){
     return 0.0;
}
double FDMDiffusionB::right(double y, double z, double t){
     return 0.0;
}
double FDMDiffusionB::source(const Eigen::Matrix<double, 3, 1>& x, double t){
    std::vector<double> lx = FDMDiffusionB::m_params->lengths[0];
    double alpha = FDMDiffusionB::alpha(x, t);
    double length = lx[1] - lx[0];
    return 5*x[0]*(length - x[0]) + 10.0*alpha*t;
}
double FDMDiffusionB::reference(const Eigen::Matrix<double, 3, 1>& x, 
        double t){
    // TODO: Define this
    std::vector<double> lx = FDMDiffusionB::m_params->lengths[0];
    double length = lx[1] - lx[0];
    // spdlog::info("length: {}, x: {}, t: {}", length, x[0], t);
    return 5.0 * t * x[0] * (length - x[0]);
}

double fdm_diffusion_b(){
    Params* p = new Params();
    p->lengths = {{0.0, 1.0}};
    p->n = {10};
    p->nt = {20};
    p->tend = 0.01;
    p->alpha = 0.1;
    FDMDiffusionB *prob = new FDMDiffusionB(p);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(prob);
    int types[6] = {1, 1, 0, 0, 0, 0};
    prob->bc_types(types);
    s.solve();
    delete p;
    delete prob;
    // TODO return the L2-norm
	return 1.0;
}

}  // namespace diffusion

}  // namespace bench
