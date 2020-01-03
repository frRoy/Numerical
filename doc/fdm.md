# Finite Difference Method {#fdm}

### Table of Contents
* [Diffusion Equation](#diffusion_fdm)
* [Convection-Diffusion Equation](#conv_diff_fdm)  --> TODO
* [Wave Equation](#wave_fdm)  --> TODO
* [Dynamical Systems](#dynamics_fdm)  --> TODO
___

# Diffusion Equation  {#diffusion_fdm}

The Laplace/Poisson eqautions are just a steady-state variation of the
time-dependent diffusion equation.

For the diffusion problem, we are using two different finite difference schemes
for the temporal and spatial mesh.

## Centered Difference Scheme (space)

Describe the concept here... see [[1]](#references)

## \f$\Theta\f$-Scheme (time)

The important feature of this time discretization scheme is that we can
implement one formula and then generate a family of well-known and
widely used schemes, namely:

* \f$\Theta\f$=0 gives the Forward Euler scheme in time.
* \f$\Theta\f$=1 gives the Backward Euler scheme in time.
* \f$\Theta\f$=1/2 gives the Crank-Nicolson scheme in time.

For a general 3D problem, the unifying \f$\Theta\f$-scheme is defined as

\f[
    \frac{u^{n+1}_{i,j,k} - u^n_{i,j,k}} {\Delta t} = 
      \Theta G\left(u^{n+1}_{i,j,k}\right) +
      \left(1-\Theta\right) G\left(u^{n}_{i,j,k}\right)
\f]

where \f$G\f$ corresponds to the centered difference scheme, i.e.

\f[
    G(u^{n}_{i,j,k}) = 
      \mathcal{F}_x\left(+u^n_{i-1,j,k}-2u^n_{i,j,k}+u^n_{i+1,j,k}\right) +
      \mathcal{F}_y\left(+u^n_{i,j-1,k}-2u^n_{i,j,k}+u^n_{i,j+1,k}\right) + 
      \mathcal{F}_z\left(+u^n_{i,j,k-1}-2u^n_{i,j,k}+u^n_{i,j,k+1}\right) +
      \Delta t f^n_{i,j,k}.
\f]

Where the Fourier coefficient are defined as:

\f[
    \mathcal{F}_x=\frac{\alpha \Delta t}{\Delta x^2},\quad
    \mathcal{F}_y=\frac{\alpha \Delta t}{\Delta y^2},\quad
    \mathcal{F}_z=\frac{\alpha \Delta t}{\Delta z^2}
\f]

A similar expression can be obtained for \f$G(u^{n+1}_{0, j})\f$, by replacing
\f$n\f$ by \f$n+1\f$.

Mapping the mesh nodes to a specific index \f$i\f$, we can define the scheme 
on a matrix system of the form

\f[
    \mathbf{A}\mathbf{u^{n+1}} = \mathbf{b}
\f]

Where the coefficient of the coefficient matrix \f$\mathbf{A}\f$ can be
defined collecting all the unknown on the left side of the scheme equation.

For interior mesh nodes, the scheme leads to a coefficient matrix with entries:

\f[
    A_{i, i-1} = -\mathcal{F}_{[]}\Theta,\quad 
    A_{i,i}=1+2\mathcal{F}_{[]}\Theta, \quad
    A_{i, i+1} = -\mathcal{F}_{[]}\Theta
\f]

where \f$\mathcal{F}_{[]}\f$ represents the Fourier coefficient along one of
the grid axis, i.e. \f$\{x,y,z\}\f$.

For interior mesh nodes, the right-hand side (RHS) vector is then

\f[
    b_i = u^n_i + \mathcal{F}_x\left(1-\Theta\right)
      \left(u^n_{i-1,j,k}-2u^n_{i,j,k}+u^n_{i+1,j,k}\right)+
      \mathcal{F}_y\left(1-\Theta\right)
      \left(u^n_{i,j-1,k}-2u^n_{i,j,k}+u^n_{i,j+1,k}\right) + 
      \mathcal{F}_z\left(1-\Theta\right)
      \left(u^n_{i,j,k-1}-2u^n_{i,j,k}+u^n_{i,j,k+1}\right) + 
      \Delta t\Theta f^{n+1}_{i,j,k}+
      \Delta t\left(1-\Theta\right) f^{n}_{i,j,k}
\f]

### Mapping

The mapping \f$p=m(i,j,k)\f$ from a mesh node with indices \f$(i,j,k)\f$ to
the corresponding unknown indice \f$p\f$ in the equation system can be defined
using

\f[
    p = m(i, j, k) = i+ j(n_x + 1) + k((n_x + 1)(n_y+1)))
\f]

This notation is obatined by numbering the mesh nodes from 
left (\f$x=x_\textrm{min}\f$) to right (\f$x=x_\textrm{max}\f$), to
bottom (\f$y=y_\textrm{min}\f$) to top (\f$y=y_\textrm{max}\f$), and to 
back (\f$z=z_\textrm{min}\f$) to front (\f$z=z_\textrm{max}\f$).

\anchor fig_mesh_fdm
\image html bloc.png <"Figure 1: A typical mesh for a bloc. The scene shows 
the bottom (orange), left (magenta) and front (blue) boundaries.">
\image latex bloc.eps "bloc" width=10cm

### Heterogeneous Media

For a spatially-dependent diffusion coefficient \f$\alpha=\alpha(x,y,z)\f$, 
the diffusion equation becomes

\f[
    \frac{\partial u}{\partial t} = \nabla\left(\alpha(x,y,z)\nabla u\right) 
      + f(x,y,z,t)
\f]

Which give the following discretized version for a line (1D)

\f[
    \frac{u_i^{n+1}-u^n_i}{\Delta t} = \Theta\frac{1}{\Delta x^2}\left(
      \alpha_{i+1/2}\left(u^{n+1}_{i+1}-u^{n+1}_i\right)-
      \alpha_{i-1/2}\left(u^{n+1}_{i}-u^{n+1}_{i+1}\right)\right)+
      \left(1-\Theta\right)\frac{1}{\Delta x^2}\left(
      \alpha_{i+1/2}\left(u^{n}_{i+1}-u^{n}_i\right)-
      \alpha_{i-1/2}\left(u^{n}_{i}-u^{n}_{i+1}\right)\right)+
      \Theta f^{n+1}_i + \left(1-\Theta\right)f^n_i
\f]

where 

\f[
    \alpha_{i+1/2} = \frac{1}{2}\left(\alpha_i+\alpha_{i+1}\right)
\f]

and

\f[
    \alpha_{i-1/2} = \frac{1}{2}\left(\alpha_i+\alpha_{i-1}\right)
\f]

are the arithmetic means of \f$\alpha(x,y,z)\f$.

## Boundary Conditions

### Dirichlet

At the Dirichelt boundary condition we fix \f$u\f$ at the boundary. For 
example, on the left boundary of a line (1D), we have the following discretized
Dirichlet boundary condition:

\f[
    u^n_{0} = u_d
\f]

For a constant diffusion coefficient, the coefficient matrix and RHS vector 
for the first two-rows would then be updated such that:

\f[
  \left(
   \begin{array}{cccc}
     1 & 0 & 0 & \ldots\\
     0 & 1+2\Theta\mathcal{F}_x & -\Theta\mathcal{F}_x & \ldots\\
     \vdots & \vdots & \vdots & \ddots
   \end{array} 
  \right) 
  \left(
    \begin{array}{c}
      u^{n+1}_0\\u^{n+1}_1\\ \vdots
    \end{array}
  \right) = 
  \left(
    \begin{array}{c}
      u_d\\b_1\\ \vdots
    \end{array}
  \right)
\f]

where \f$b_1\f$ is

\f[
    b_1 = u^n_1 + \mathcal{F}_x
    \left(1-\Theta\right)\left(u^n_{2}-2u^n_1+u^n_0\right) + 
    \Delta t\Theta f^{n+1}_1 + \Delta t\left(1-\Theta\right)f^n_i
\f]


### Neumann

A natural approximation to the normal derivative is a one sided difference, 
for example, on the left boundary of a rectangle (2D), and using a constant 
diffusion coefficient \f$\alpha\f$ we have the following discretized 
Neumann boundary condition:

\f[
    \frac{\partial u}{\partial x} \approx 
        \frac{u^n_{-1, j} - u^n_{1, j}}{2\Delta x} = g^n_{0,j}
\f]

where \f$g^n_{0,j,k}\f$ is the set value at the spatial node \f$\{i,j,k\}\f$.

Note that the value \f$u^n_{-1, j, k}\f$ is not well defined and need to be 
eliminated from the system of equations. This can be achieved substituting 
the Neumann boundary condition in the centered difference scheme, i.e.

\f[
    \frac{u^{n+1}_{0,j}-u^{n}_{0,j}}{\Delta t} = 
      \Theta G(u^{n+1}_{0, j}) + \left(1-\Theta\right)G(u^{n}_{0, j})
\f]

where

\f[
    G(u^{n}_{0, j}) = 
      \mathcal{F}_x\left(\left(2g^n_{0,j}\Delta x+u^n_{1,j}\right)
      -2u^n_{0,j}+u^n_{1,j}\right) +
      \mathcal{F}_y\left(u^n_{0,j-1}-2u^n_{0,j}+u^n_{0,j+1}\right)+
      \Delta t f^n_{0, j}.
\f]

A similar expression can be obtained for \f$G(u^{n+1}_{0, j})\f$, by replacing
\f$n\f$ by \f$n+1\f$.

The coefficient matrix and RHS vector for the fifth row of a 4 by 3 
rectangular grid would then be updated such that:

\f[
 [2 \Theta\mathcal{F}_x, -2\Theta Fx, 0, 0, \dots] [u^{n+1}_5] = b_5
\f]

where \f$b_5\f$ is

\f[
    b_5 = 
    u^n_5+\left(1-\Theta\right)\left(\mathcal{F}_x\left(2u^n_6 -2u^n_5\right)+
    \mathcal{F}_y\left(u^n_0 - 2u^n_5 + u^n_{10}\right)\right) +
    \Theta\left(\Delta tf^{n+1}_5+2\mathcal{F}_x\Delta xg^{n+1}_5\right)+
    \left(1-\Theta\right)\left(\Delta tf^{n}_5+
    2\mathcal{F}_x\Delta xg^{n}_5\right)
\f]

In order to keep the coefficient matrix symmetric, one may choose to scale
the fifth row of the augmented system by 0.5, which gives:

\f[
 [\Theta\mathcal{F}_x, -\Theta Fx, 0, 0, \dots] [u^{n+1}_5] = b_5
\f]

where \f$b_5\f$ is

\f[
    b_5 = \frac{1}{2}\left(
    u^n_5+\left(1-\Theta\right)\left(\mathcal{F}_x\left(2u^n_6 -2u^n_5\right)+
    \mathcal{F}_y\left(u^n_0 - 2u^n_5 + u^n_{10}\right)\right) +
    \Theta\left(\Delta tf^{n+1}_5+2\mathcal{F}_x\Delta xg^{n+1}_5\right)+
    \left(1-\Theta\right)\left(\Delta tf^{n}_5+
    2\mathcal{F}_x\Delta xg^{n}_5\right)\right)
\f]

**Note that the coefficient matrix may not be symmetric for heterogeous 
media.**

At corner points, the norm vector is not well defined. We use average of two
directional derivatives to get an approximation. 

Taking (0, 0) as an example, we have ...

\f[
    \begin{align}
    u_{-1, 0} - u_{1, 0} &= 2\Delta x g_{1,1}\\
    u_{0, -1} - u_{0, 1} &= 2\Delta x g_{1,1}
    \end{align}
\f]

## Nonlinear Problems

There three possible approaches:

1. Discretize the problem in time and linearize the time-discrete (stationary) 
  PDE problem at each time level into a sequence of linear PDE problems.
2. Discretize the problem in time and space and solve the resulting nonlinear
  algebraic equations at each time level.
3. Operator splitting.

### First approach (PDE level)
Linearizing time-discrete PDEs directly prior to discretization 
in space. In this case we can use either the Picard or Newton method for
implicit shemes (backward Euler of Crank-Nicolson). For the forward Euler
method, we can directly introduce the nonlinear term in the discretized
equation. The disadvantage of the explicit method is the strict stability
criterion.

### Second approach (Algebraic level)
Carry out the discretization in space of the time-discrete 
nonlinear PDE problem and get a system of nonlinear algebraic equations.

### Third approach

### Picard Method

The Picard method iterates on the spatialy discretized equation for a given
time step until the difference between the last and new solution match at a
given tolerence value.

that is iterates until

\f[
  u^{n, k+1} = u^{, k} + \textrm{TOL}
\f]

### Newton Method

**PDE level**:

\f[
  u^{n, k+1} = u^{n, k} + \delta u
\f]

**Algebraic level**:

### Generalized method (PDE level)

\f[
  \frac{u^{n}-u^{n-1}}{\Delta t} = \nabla\cdot\left(\alpha(u^{n, k})
    \nabla u\right) + f(u^{n ,k}) + \gamma\left(
    \nabla\cdot\left(\alpha'(u^{n, k})\left(u^{n}-u^{n, k}\right)
    \nabla u^{n,k}\right)+f'(u^{n, k})\left(u^n-u^{n,k}\right)\right)
\f]

where the primed quantities represent the derivative with respect to \f$u\f$.

Note that if \f$\gamma = 0\f$, we have the Picard method, and if 
\f$\gamma=1\f$, we have the Newton method.

## References        {#references}

1. H. P. Langtangen and S. Linge. Finite Difference Computing with PDEs: 
  A Modern Software Approach. Springer, 2017.