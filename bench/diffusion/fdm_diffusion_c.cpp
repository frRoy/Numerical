/**
 *  @file fdm_diffusion_c.cpp
 *  @brief Run a 1D diffusion problem with Neumann boundary conditions.
 *  @author Francois Roy
 *  @date 12/01/2019
 */
#include "spdlog/spdlog.h"
#include "numerical/fdm/Parameters.hpp"
#include "diffusion.hpp"

namespace bench {
namespace diffusion {
typedef numerical::fdm::Parameters<double> Params;
typedef numerical::fdm::Problem<double> Problem;

double FDMDiffusionC::left(double y, double z, double t){
    std::vector<double> lx = FDMDiffusionC::m_params->lengths[0];
    double length = lx[1] - lx[0];
     return 5.0*t*length;  // outward direction
}
double FDMDiffusionC::right(double y, double z, double t){
    std::vector<double> lx = FDMDiffusionC::m_params->lengths[0];
    double length = lx[1] - lx[0];
     return -5.0*t*length;  // outward direction
}
double FDMDiffusionC::source(const Eigen::Matrix<double, 3, 1>& x, double t){
    std::vector<double> lx = FDMDiffusionC::m_params->lengths[0];
    double alpha = FDMDiffusionC::alpha(x, t);
    double length = lx[1] - lx[0];
    return 5*x[0]*(length - x[0]) + 10.0*alpha*t;
}
double FDMDiffusionC::reference(const Eigen::Matrix<double, 3, 1>& x, 
        double t){
    std::vector<double> lx = FDMDiffusionC::m_params->lengths[0];
    double length = lx[1] - lx[0];
    // spdlog::info("length: {}, x: {}, t: {}", length, x[0], t);
    return 5.0 * t * x[0] * (length - x[0]);
}

double fdm_diffusion_c(){
    Params* p = new Params();
    p->lengths = {{0.0, 1.0}};
    p->n = {10};
    p->nt = {10};
    p->tend = 0.4;
    p->alpha = 0.1;
    p->theta = 1.0;
    FDMDiffusionC *prob = new FDMDiffusionC(p);
    numerical::fdm::SparseSolver<double> s = 
        numerical::fdm::SparseSolver<double>(prob);
    int types[6] = {1, 1, 0, 0, 0, 0};
    prob->bc_types(types);
    double l2norm = s.solve();
    delete p;
    delete prob;
    return l2norm;
}

}  // namespace diffusion

}  // namespace bench
