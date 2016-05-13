/* TO DO:
 * - AQUI: time stepping works, do results make sense?
 *   + De momento no cambia mucho y de repente se des estabiliza. Por que?
 *   + probar con elastico puro
 * - Use another time step? (maybe try relaxation time?)
 * - Compare solutions with analytical
 * - Remove comparason with turcotte and savage
 * - Weakening mechanism
 *    + The functions ElasticModulus and Viscosity will be useless after the first time step
 */

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2007 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Martin Kronbichler, Uppsala University,
 *          Wolfgang Bangerth, Texas A&M University 2007, 2008
 */
/*! \mainpage
 * \section intro Introduction
 * This code simulates the deformation of a viscoelastic medium in presence
 * of a discrete fault under the anti-plane shear approximation.
 * Under this assumption there is no displacement in the modeled plane
 * \f$(v_{x}=v_{y}=0)\f$ and there is no spatial variation in the direction
 * perpendicular to the modeled plane \f$(\partial z=0)\f$. If the case of a
 * uniform shear modulus and shear viscosity, the problem is reduced to the
 * Poisson equation where the right-hand side depends on the velocity values
 * from the previous time step.
 *
 * \note This code has been modified from deall.II tutorial
 * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_31.html">step-31</a>.
 *
 * \subsection equations Viscoelastic approximation
 * We are interested in a system in which large deformations occur. In addition,
 * we would like to take advantage of the methods that have already been studied.
 * Therefore, we will try to adapt the equations that drive our problem so they
 * resemble fluid flow equations, for which extensive research has been done
 * regarding finite element methods and solving algorithms. For such purpose,
 * we will start with the equation of conservation of momentum for an
 * incompressible fluid:
 * \f[\partial_i{(\sigma_{ij})}=f_j\f]
 * where \f$i\f$ and \f$j\f$ run from 1 to 3 and designate spatial coordinate,
 * \f$\sigma\f$ is the stress tensor and \f$f\f$ is the total external force, in
 * this case, only the gravitational force
 * \f$f_j=f^g_j=\rho g \left(1-\alpha T \right) \delta_{iy}\f$. The total stress can
 * be separated in the deviatoric component \f$\tau_{ij}\f$ and pressure \f$P\f$:
 * \f[ \sigma_{ij} = \tau_{ij}-P \delta_{ij} \f]
 *
 * and the equation of conservation of motion will be:
 * \f[\partial_i \tau_{ij} - \partial_j P =
 * \rho g \left( 1 - \alpha T \right) \delta_{jy} \f]
 *
 *
 * Since we are going to use a semi-viscous approach we need a constitutive
 * equation to relate the deviatoric stress and strain rate. For that, we will
 * follow the approach from \cite moresi_et_al_02 and \cite moresi_et_al_03
 * and consider that the total deformation is acomodated by a linear combination
 * of elastic and viscous deformation:
 * \f[ \dot \varepsilon = \dot \varepsilon^e + \dot \varepsilon^v =
 * \frac{\check{\tau}}{2\mu} + \frac{\tau}{2\eta} \f]
 *
 * where \f$\dot\varepsilon\f$ is the strain rate, \f$\check{\tau}\f$
 * is the Jaumann corotatonal stress rate for an element of the continuum,
 * \f$\mu\f$ is the shear modulus and \f$\eta\f$ is shear viscosity. The
 * strain rate can be expressed in terms of the velocity of the continuum
 * \f$u\f$:
 * \f[ \dot \varepsilon_{ij} =
 * \frac{1}{2} \left( \partial_j v_i + \partial_i v_j \right) \f]
 *
 * The Jaumann derivative is an observer independent objective rate (see
 * /cite harder_91 for a discussion about the effect of using different objective
 * derivatives on convection solutions). Jaumann stress can be expressed as:
 * \f[ \check{\tau} = \frac{D\tau}{Dt}+\tau\omega-\omega\tau\f]
 *
 * where \f$\omega\f$ is the material spin tensor, which can be expressed
 * in terms of the velocity:
 * \f[\omega_{ij} =
 * \frac{1}{2} \left( \partial_j v_i - \partial_i v_j \right) \f]
 *
 * \subsection num_approach Numerical approach
 * To use a semi-viscous approach we need to manipulate the constitutive equation
 * to write the deviatoric stress in terms of the strain rate (therefore in terms
 * of derivatives of the velocity) and the substitute it in the momentum equation.
 * This will give us an PDE of the velocities, which we will have to solve.To achieve
 * this, the first step is to express the Jaumann stress rate as a first order
 * difference form:
 * \f[\check{\tau} = \frac{\tau^t - \tau ^{t-1}}{\Delta t} + \tau^{t-1}\omega^{t-1}
 * - \omega^{t-1}\tau^{t-1} \f]
 *
 * where the superscript \f$t\f$ indicates the current time step and \f$t-1\f$ the
 * previous one, and \f$\Delta t\f$ is the time increment during the last step.
 * The deviatoric stress tensor for the current time step can therefore be expressed
 * in terms of the current strain rate and past deviatoric stress and material spin
 * tensors:
 * \f[ \tau^t = 2 \eta_{ef} \dot\varepsilon^t +
 * \eta_{ef} \left[\frac{1}{\mu\Delta t} \tau ^{t-1} +
 * \frac {1}{\mu} \left( \omega^{t-1}\tau^{t-1} -
 * \tau^{t-1}\omega^{t-1} \right) \right] \f]
 *
 * where \f$\eta_{ef}\f$ is an effective viscosity that depends on both the elastic
 * and viscous coefficients:
 * \f[ \eta_{ef} = \eta \frac{\mu \Delta t}{\mu \Delta t + \eta} \f]
 *
 * Using this in the momentum equation we obtain:
 * \f[ \partial_i \left(2 \eta_{ef} \dot\varepsilon_{ij} \right) - \partial_j P =
 * f^g_j
 * - \partial_i \left\{ \eta_{ef} \left[\frac{1}{\mu\Delta t} \tau ^{t-1}_{ij} +
 * \frac {1}{\mu} \left( \omega^{t-1}_{ik}\tau^{t-1}_{kj} -
 * \tau^{t-1}_{ik}\omega^{t-1}_{kj} \right) \right] \right\} \f]
 *
 * This equation is very similar to the equation of conservation of momentum
 * for an incompressible fluid, but the discretization of the stress rate has
 * included an additional internal elastic force force that depends on the
 * stress and material spin from the previous time steps. If we define this
 * internal elastic force the divergence of the second-order tensor \f$ g^e_{ij}\f$:
 * \f[ f^e_j = -\partial_i g^e_{ij}  = -\partial_i \left\{ \eta_{ef}
 * \left[\frac{1}{\mu\Delta t} \tau ^{t-1}_{ij} +
 * \frac {1}{\mu} \left( \omega^{t-1}_{ik}\tau^{t-1}_{kj} -
 * \tau^{t-1}_{ik}\omega^{t-1}_{kj} \right) \right] \right\} \f]
 *
 * the PDE we need to solve is simply:
 * \f[ \partial_i \left(2 \eta_{ef} \dot\varepsilon_{ij} \right) - \partial_j P =
 * f^g_j + f^e_j \f]
 *
 * Therefore, we will follow the next scheme to solve the problem:
 *   -# Input initial conditions for the stress \f$\tau^0\f$ and velocity \f$v^0\f$.
 *   -# Compute the initial material spin \f$\omega^0(v^0)\f$.
 *   -# Compute the internal elastic force that will be used in the firts time step
 *   \f$f^{e^1}_j = \partial_i g^{e^1}_{ij}(\tau^0,\omega^0)\f$.
 *   -# Solve the momentum equation to obtain the velocity of the first time step
 *   \f$v^1\f$.
 *   -# Compute the strain rate of the first time step \f$\dot\varepsilon^1(v^1)\f$,
 *   stress \f$\tau^1(\dot{\varepsilon}^1,\tau^0,\omega^0)\f$ and the material spin
 *   \f$\omega^1(v^1)\f$.
 *   -# Repeat steps 3 to 5.
 *
 *
 *
 * \subsubsection anti_plane Anti-plane shear approximation
 * The last equation is valid for elastic problems in any dimensions. However,
 * solving it can be computationally demanding. For certain problems it is enough
 * to solve the equiation under certain approximations. In this case we are going
 * to model the displacement in the \f$x-y\f$ plane produced by a fault that runs
 * parallel to the \f$y-z\f$ plane and is sliding parallel to the \f$z\f$ direction.
 * \image html aps.png "Figure 1: Schematic representation of a fault and the modeled plane"
 * \image latex aps.eps "Figure 1: Schematic representation of a fault and the modeled plane"
 *
 * In this case, we can considered that displacement only occurs in the \f$z\f$
 * direction \f$(v_{x}=v_{y}=0)\f$. Moreover, under this approximation no
 * magnitude can vary along the \f$z\f$ direction \f$(\frac{\partial}{\partial z}=0)\f$.
 * Therefore, under this approximation
 * \f$\varepsilon_{xx}=\varepsilon_{yy}=\varepsilon_{zz}=\varepsilon_{xy}=\varepsilon_{yx}=0\f$
 * and \f$\omega_{xx}=\omega_{yy}=\omega_{zz}=\omega_{xy}=\omega_{yx}=0\f$
 * and the governing equations are reduced to:
 * \f[-\partial_x P = f^e_x \f]
 * \f[-\partial_y P = f^g_y + f^e_y \f]
 * \f[\partial_x(2\eta_{ef}\varepsilon_{xz})+\partial_y(2\eta_{ef}\varepsilon_{yz})=f^e_z\f]
 *
 * where
 * \f[f^e_x = -\partial_x \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{xx} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{xz}\tau^{t-1}_{zx} - \tau^{t-1}_{xz}\omega^{t-1}_{zx} \right)
 * \right] \right\} -\partial_y \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{yx} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{yz}\tau^{t-1}_{zx} - \tau^{t-1}_{yz}\omega^{t-1}_{zx} \right)
 * \right] \right\} \f]
 *
 * \f[f^e_y = -\partial_x \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{xy} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{xz}\tau^{t-1}_{zy} - \tau^{t-1}_{xz}\omega^{t-1}_{zy} \right)
 * \right] \right\} -\partial_y \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{yy} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{yz}\tau^{t-1}_{zy} - \tau^{t-1}_{yz}\omega^{t-1}_{zy} \right)
 * \right] \right\} \f]
 *
 * \f[f^e_z = -\partial_x \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{xz} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{xz}\tau^{t-1}_{zz} - \tau^{t-1}_{xx}\omega^{t-1}_{xz} -
 * \tau^{t-1}_{xy}\omega^{t-1}_{yz} \right) \right] \right\}
 * -\partial_y \left\{\eta_{ef}\left[ \frac{1}{\mu\Delta t} \tau^{t-1}_{yz} +
 * \frac{1}{\mu} \left(\omega^{t-1}_{yz}\tau^{t-1}_{zz} - \tau^{t-1}_{yx}\omega^{t-1}_{xz} -
 * \tau^{t-1}_{yy}\omega^{t-1}_{yz} \right) \right] \right\}\f]
 *
 * The velocity at time \f$t\f$ only appears in the \f$z\f$ component equation. Therefore,
 * we will only solve said equation. In terms of \f$v_z\f$ the equation is:
 * \f[\partial_x \left(\eta_{ef}\partial_x v_z \right) +
 * \partial_y \left(\eta_{ef}\partial_y v_z \right) = f^e_z \f]
 *
 * and for uniform effective viscosity we obtain the Poisson equation:
 *
 * \f[\eta_{ef}\Delta v_z = f^e_z \f]
 *
 *
 * \section model Modeling
 * \subsection setup Setup
 * Here we will solve the visco-elastic equation in a two-dimensional domain
 * of dimensions \f$a\f$ by \f$b\f$. We will solve for the velocity in the
 * direction perpendicular to the plane \f$(v_z)\f$ caused by a fault also
 * perpendicular to the modeled plane.
 *
 * There are four possible boundary conditions:
 *
 * - \f$\Gamma_0\f$ - The locked part above the fault, with imposed zero velocity, \f$v_z=0\f$
 * (homogeneous Dirichlet boundary conditions).
 *
 * - \f$\Gamma_1\f$ - Regions of imposed velocity, \f$v_z=V\f$ (non-homogeneous Dirichlet
 * boundary conditions).
 *
 * - \f$\Gamma_2\f$ - Zero tangential stress boundaries, \f$\boldsymbol{n} \cdot \tau - P =0\f$
 * (homogeneous Neumann boundary conditions).
 *
 * - \f$\Gamma_3\f$ - Imposed tangential stress \f$\boldsymbol{n} \cdot \tau - P =S\f$
 * (non-homogeneous Neumann boundary conditions).
 *
 *
 * \subsection weak_form The Weak Form
 *
 * In order to solve the equation using the Finite Element Method we need to derive
 * the weak form. We will begin with the partial differential equation that governs
 * the velocity in a visco-elastic medium under the anti-plane shear approximation:
 * \f[ \nabla \left(\eta_{ef}\nabla v_z\right)=f^e_z\f]
 *
 * where i=1,2. First we multiply the equation from the left by a test function \f$v\f$,
 * and then we integrate over the whole domain:
 * \f[\int_{\Omega} v \cdot \nabla\left(\eta_{ef}\nabla v_z \right) =
 * \int_{\Omega} v \cdot f^e_z \phi \f]
 *
 * Integrating by parts and using the Gauss theorem
 * \f[ \int_{\Omega} \nabla v \cdot \eta_{ef}\nabla v_z =
 * -\int_{\Omega} v \cdot f^e_z + \int_{\partial \Omega} v \cdot \eta_{ef} \boldsymbol{n} \nabla v_z \f]
 *
 * We can integrate by parts the term corresponding to the internal elastic force
 * if we take into account that it is in fact the divergence of another function
 * \f$ f^e_z = -\partial_i g^e_{iz} \f$. Therefore, we obtain:
 * \f[ \int_{\Omega} \nabla v \cdot \eta_{ef}\nabla v_z =
 * -\int_{\Omega} \nabla v \cdot \boldsymbol{g^e_z} +
 * \int_{\partial \Omega} v \cdot \eta_{ef} \boldsymbol{n} \nabla v_z +
 * \int_{\partial \Omega} v \cdot \boldsymbol{n} \boldsymbol{g^e_z} \f]
 *
 * Considering that the boundary terms are related to the stress tensor
 * \f$ \tau_{ij}\f$ we obtain:
 * \f[ \int_{\Omega} \nabla v \cdot \eta_{ef}\nabla v_z =
 * -\int_{\Omega} \nabla v \cdot \boldsymbol{g^e_z} +
 * \int_{\partial \Omega} v \cdot \boldsymbol{n} \tau_{iz} \f]
 *
 * or:
 * \f[ \left( \nabla v, \eta_{ef}\nabla v_z \right)_{\Omega} =
 * - \left( \nabla v, \boldsymbol{g^e_z} \right)_{\Omega} +
 * \left(v, \boldsymbol{n} \boldsymbol{\tau_{z}} \right)_{\partial \Omega} \f]
 *
 * where \f$\boldsymbol{\tau_{z}} = \tau_{iz}\f$.
 *
 * In those boundaries where the velocity is imposed (Dirichlet boundary conditions), this is
 * \f$\Gamma_0\f$ and \f$\Gamma_1\f$,the test function \f$v = 0\f$. The contribution of those
 * boundaries in which the imposed tangential stress is zero (i.e., \f$\Gamma_2\f$) will not
 * contribute to the right hand side either. On the other hand,those boundaries in which we
 * impose a non-zero tangential stress (non-homogeneous Neumann boundary conditions), this is
 * \f$\Gamma_3\f$, the test functions \f$v \neq 0\f$, and that part of the boundary integral
 * has to be included in the weak form.
  * \f[ \left( \nabla v, \eta_{ef}\nabla v_z \right)_{\Omega} =
 * - \left( \nabla v, \boldsymbol{g^e_z} \right)_{\Omega} +
 * \left(v, S \right)_{\Gamma^3} \f]
 *
 * where the operator \f$(a, b)\f$ is \f$ \int a \cdot b \f$ and \f$S\f$ is the imposed tangential
 * stress at the boundary.
 *
 * We can now  replace the exact
 * solution \f$ v_z \f$ for an approximate solution in terms of the shape functions,
 * \f$ \sum_{j} U_j v_j(\boldsymbol{x}) \f$. Here \f$U_j\f$ are the expansion coefficients we need to
 * determine to find the approximate solution to the laplace equation. Then, the problem
 * is reduced to solving
 * \f[AU=B\f]
 * where the matrix A and the vector B are defined as:
 * \f[A_{ij} = (v_i,v_j)_{\Omega}\f]
 * \f[B_i = -(\nabla v_i, \boldsymbol{g^e_z})_\Omega
 * + (v_i, S)_{\Gamma_3} \f]
 *
 *
 * \subsection adaptive_refinement Adaptive Mesh Refinement in Time Dependent Problems
 * Every time step we need to use the solution from the two last time steps and the visco-elastic
 * stress from the previous one. The old solutions are stored in a
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classVector.html">Vector</a> with as
 * many elements as the degrees of freedom the problem has. The old stress is stored at each
 * quadrature and can be accessed through a user pointer that each cell holds. Every time the
 * mesh is refined, the number of cells and degrees of freedom changes, therefore, as part of
 * the mesh refinement we need to transfer the solution vectors and the old stress from the
 * old mesh to the new one. For the solution vectors we just need to use
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classSolutionTransfer.html">SolutionTransfer</a>
 * (see for example
 * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_26.html#codeHeatEquationrefine_meshcode">Step-26</a>).
 *
 * Transferring the old stress is slightly more complicated, we need to transfer the data sotred
 * in the quadrature points to a finite element field that is defined everywhere so that we can
 * later transfer it to the new mesh (using
 * <a href="https://www.dealii.org/8.4.0/doxygen/deal.II/classSolutionTransfer.html">SolutionTransfer</a>
 * and then interpolate it to the new quadrature points. We need a discontinuous field that
 * matches the values in the quadrature points (we will use a Discontinuous Galerking finite element
 * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFE__DGQ.html">FE_DGQ</a>).
 */


/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth and Ralf Hartmann, University of Heidelberg, 2000
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include <typeinfo>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <math.h>

#define PI 3.14159265

namespace vsf
{
  using namespace dealii;

  /**
   * \struct PointHistory
   * We first define a structure where we will store the old data that we will need
   * to compute the right hand side of the equation. In this case it is enough if we save
   * the elastic stress from the previous time step.
   */
  template <int dim>
  struct PointHistory
  {
      Tensor<2,3>       old_stress;
      Tensor<1,2>       elastic_stress;
  };
  /**
   * We declare and define some function classes that represent the body force, the effective
   * viscosity and the boundary values. For simplicity we have chosen to define all relevant parameters
   * in a base class, which later will serve as the base class for the rest of the classes.
   * These parameters are the width and height of the model (FunctionBase::width and
   * FunctionBase::height), lock depth \f$d\f$ (FunctionBase::locked_depth),
   * the imposed velocity at the boundary \f$V\f$ (FunctionBase::boundary_velocity) and the imposed
   * stress at the boundary \f$S\f$ (FunctionBase::boundary_stress).
   *
   * \note These parameters have to be hardcoded here because the solution and boundary
   * conditions classes are derived from
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFunction.html">Function</a>
   * and are used directly by other deal.II class
   * functions. Adding any extra input parameters in a class that inherits from
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/classFunction.html">Function</a>
   * would require modifying the source code. Another option would be to input these parameter
   * via an input file (see deal.II tutorial
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_33.html">step-33</a>).
   *
   * All the parameters used in the code are set here to keep things in order.
   */

  template <int dim>
  class FunctionBase
  {
  protected:
    static const double width; /**< Model's width. Will not be used here, but we want to keep input parameters together.*/
    static const double height; /**< Model's height.*/
    static const double locked_depth; /**< Depth of the locked region.*/
    static const double boundary_velocity; /**< Velocity imposed in some boundary.*/
    static const double boundary_stress; /**< Stress imposed in some boundaries.*/
    static const double elastic_modulus; /**< Reference elastic modulus.*/
    static const double viscosity; /**< Reference viscosity.*/
  };

  template <int dim>
  const double FunctionBase<dim>::width = 2.0e5;

  template <int dim>
  const double FunctionBase<dim>::height = 1.0e5;

  template <int dim>
  const double FunctionBase<dim>::locked_depth = 0.25e5;

  template <int dim>
  const double FunctionBase<dim>::boundary_velocity = 10/(3.6e3*24.0*365.0);

  template <int dim>
  const double FunctionBase<dim>::boundary_stress = 1e5;

  template <int dim>
  const double FunctionBase<dim>::elastic_modulus = 1.0e11;

  template <int dim>
  const double FunctionBase<dim>::viscosity = 1e20;



  /** This class is used to compute the analytical solution of the surface displacement caused
   * by a fault in a completely elastic medium using \cite turcotte_spence_74.
   */
  template <int dim>
  class Solution_Turcotte : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Solution_Turcotte () : Function<dim>() {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Solution::value extracts the value of the analytical solution in one point ("p").
   */
  template <int dim>
  double Solution_Turcotte<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
  {
    /** This solution is only valid at the surface, therefore, if the point at which it
     * will be evaluated is not in the surface, an exception will be thrown.
     */
    Assert (std::fabs(p[1] - 0)<1e-12, ExcNotImplemented ());

    const double c_1 = PI/(2*this->height);
    const double c_2 = sin (c_1 * this->locked_depth);
    return (this->boundary_stress/c_1)*log((sinh(c_1*p[0])+sqrt(pow(sinh(c_1*p[0]),2)+pow(c_2,2)))/c_2);
  }

  /** Solution::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Solution_Turcotte<dim>::value_list (const std::vector<Point<dim> >   &points,
                                           std::vector<double>              &values,
                                           const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Solution_Turcotte::value(points[i]);
  }



  /** This class is used to compute the analytical solution of the surface displacement caused
   * by a fault in a completely elastic medium using \cite savage_burford_73.
   */
  template <int dim>
  class Solution_Savage : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Solution_Savage () : Function<dim>() {}
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Solution::value extracts the value of the analytical solution in one point ("p").
   */
  template <int dim>
  double Solution_Savage<dim>::value (const Point<dim>   &p,
                             const unsigned int) const
  {
    /** This solution is only valid at the surface, therefore, if the point at which it
     * will be evaluated is not in the surface, an exception will be thrown.
     */
    Assert (std::fabs(p[1] - 0)<1e-12, ExcNotImplemented ());

    return (this->boundary_velocity/PI)*atan(p[0]/(this->locked_depth));
  }

  /** Solution::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Solution_Savage<dim>::value_list (const std::vector<Point<dim> >   &points,
                                  std::vector<double>              &values,
                                  const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Solution_Savage::value(points[i]);
  }



  /**
  * The template class Elastic_Modulus gives the value of the elastic modulus in the requested point/s
  */
  template <int dim>
  class Elastic_Modulus : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Elastic_Modulus () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Elastic_Modulus::value extracts the value in one point ("p"). In the current implementation
   * we consider a uniform elastic modulus. To use a spatially variable elastic modulus modify this
   * method accordingly.
   */
  template <int dim>
  double Elastic_Modulus<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
     return this->elastic_modulus;
  }

  /** Elastic_Modulus::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Elastic_Modulus<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Elastic_Modulus::value(points[i]);
  }



  /**
  * The template class Viscosity gives the value of the viscosity in the requested point/s
  */
  template <int dim>
  class Viscosity : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Viscosity () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Viscosity::value extracts the value in one point ("p"). In the current implementation
   * we consider a uniform elastic modulus. To use a spatially variable elastic modulus modify this
   * method accordingly.
   */
  template <int dim>
  double Viscosity<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
     return this->viscosity;
  }

  /** Viscosity::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Viscosity<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Viscosity::value(points[i]);
  }



   /**
   * The template class BodyForce gives the value of body force (on the right hand side of the
   * elastic equation) in the requested point/s
   */
   template <int dim>
   class BodyForce : public Function<dim>,
   protected FunctionBase<dim>
   {
   public:
     BodyForce () : Function<dim>() {}

     virtual double value (const Point<dim>   &p,
                           const unsigned int  component = 0) const;
     virtual void value_list (const std::vector<Point<dim> >      &points,
                              std::vector<double>                 &values,
                              const unsigned int                  component = 0) const;
   };

   /** RightHandSide::value extracts the value in one point ("p"). In the current implementation
   * we consider a zero body force. To use a spatially variable elastic modulus modify this
   * method accordingly.
   */
   template <int dim>
   double BodyForce<dim>::value (const Point<dim>   &p,
                                     const unsigned int) const
   {
      return 0.0;
   }

   /** BodyForce::value_list allows to extract values for several points ("points") at once and
    * outputs them in a vector ("values").
    */
   template <int dim>
   void BodyForce<dim>::value_list (const std::vector<Point<dim> >   &points,
                                std::vector<double>              &values,
                                const unsigned int               component) const
   {
     const unsigned int n_points = points.size();
     Assert (values.size() == n_points,
             ExcDimensionMismatch (values.size(), n_points));
     Assert (component == 0,
             ExcIndexRange (component, 0, 1));

     for (unsigned int i=0; i<n_points; ++i)
       values[i] = BodyForce::value(points[i]);
   }



  /**
  * The template class Dirichlet_BC computes the values of the non-homogeneous Dirichlet BC
  * in each point.
  */
  template <int dim>
  class Dirichlet_BC : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
      Dirichlet_BC () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Dirichlet_BC::value extracts the value in one point ("p"). To use a spatially variable BC
   * modify this method accordingly.
   */
  template <int dim>
  double Dirichlet_BC<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
    /** The displacement imposed at the boundary is only half the total displacement because we are only
     * modeling half the domain.
     */
    return this->boundary_velocity/2;
  }

  /** Dirichlet_BC::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Dirichlet_BC<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Dirichlet_BC::value(points[i]);
  }



  /**
   * The template class Neumann_BC computes the values of the non-homogeneous BC in each point. These boundary
   * conditions are no used in the current implementation of this model, but this class included here for
   * completeness.
   */
  template <int dim>
  class Neumann_BC : public Function<dim>,
  protected FunctionBase<dim>
  {
  public:
    Neumann_BC () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
    virtual void value_list (const std::vector<Point<dim> >      &points,
                             std::vector<double>                 &values,
                             const unsigned int                  component = 0) const;
  };

  /** Neumann_BC::value extracts the value in one point ("p"). To use a spatially variable BC
   * modify this method accordingly.
   */
  template <int dim>
  double Neumann_BC<dim>::value (const Point<dim>   &p,
                                    const unsigned int) const
  {
    return this->boundary_stress;
  }

  /** Neumann_BC::value_list allows to extract values for several points ("points") at once and
   * outputs them in a vector ("values").
   */
  template <int dim>
  void Neumann_BC<dim>::value_list (const std::vector<Point<dim> >   &points,
                               std::vector<double>              &values,
                               const unsigned int               component) const
  {
    const unsigned int n_points = points.size();
    Assert (values.size() == n_points,
            ExcDimensionMismatch (values.size(), n_points));
    Assert (component == 0,
            ExcIndexRange (component, 0, 1));

    for (unsigned int i=0; i<n_points; ++i)
      values[i] = Neumann_BC::value(points[i]);
  }



  /**
   * The main class of the code, contains all the functions that do the job of setting up and
   * solving the model.
   *
   * \note Since the model parameters have been set in FunctionBase, this new class
   * will also inherit from it.
   */
  template <int dim>
  class ApShear : protected FunctionBase<dim>
  {
  public:
    enum RefinementMode
    {
      global_refinement, adaptive_refinement
    };

    enum Model
    {
      turcotte, savage
    };

    ApShear (const FiniteElement<dim>  &fe, /**< The code is written so different elements can be used and are defined as input parameters.*/
             const FiniteElement<dim>  &history_fe, /**< FE to store history data*/
             const unsigned int        degree, /**< To be consistent, we define the degree of the finite elements as an input parameter. */
             const RefinementMode      refinement_mode,/**< The mesh can be refined either globally (global_refinement) or adaptively (adaptive_refinement), which one to use is defined as an input parameters */
             const Model               model /**< The code can be used to simmulate the deformation caused by a fault using two models.*/
                          );

    ~ApShear ();

    void run ();

  private:
    void do_first_time_step (); /**< Create the initial mesh and perform preliminary refinements and solve for the initial time step*/
    void do_time_step (); /**< Solve for successive time steps (includes mesh refinement every several time steps).*/

    void make_grid_turcotte (); /**<  Create initial mesh by refining globally several times. Prepare boundaries for Turcotte and Spence Model */
    void make_grid_savage (); /**<  Create initial mesh by refining globally several times. Prepare boundaries for Savage and Burford Model */
    void refine_grid (const unsigned int min_grid_level,
                      const unsigned int max_grid_level); /**<  Refines the mesh adaptively. */
    void qpoints_to_DG (); /**< The information stored in the quadrature points is moved to a Discontinuous Galerkin space*/
    void DG_to_qpoints (); /**< The information stored in a Discontinious Galerking space is moved to the quadrature points*/

    void setup_system (); /**<  Setup DoFs, renumbering and constraints. */
    void setup_quadrature_point_history (); /**< Initialize the data stored at quadrature points. */
    void update_quadrature_point_history (); /**< Update the data stored at quadrature points and sore it in a DG space. */
    void assemble_system (); /**<  Assemble stiffness matrix and RHS. */

    void solve (); /**< The solver an refinement steps are defined here.*/
    void get_time_step_and_displacement (); /**< Obtain the maximal velocity to determine the time step increase and compute the total diplacement*/

    void compare_solutions_turcotte (); /**< The numerical solution is compared against Turcotte and Spence Model*/
    void compare_solutions_savage (); /**< The numerical solution is compared against Savage and Burford Model*/
    void output_results (); /* Results are output in vtk format.*/

    const unsigned int                            degree;

    Triangulation<dim>                            triangulation;
    DoFHandler<dim>                               dof_handler,
                                                  history_dof_handler;
    const QGauss<dim>                             quadrature_formula;

    SmartPointer<const FiniteElement<dim> >       fe,
                                                  history_fe;

    ConstraintMatrix                              constraints;

    SparsityPattern                               sparsity_pattern;
    SparseMatrix<double>                          system_matrix;
    
    Vector<double>                                system_rhs,
                                                  solution,
                                                  old_solution,
                                                  displacement;


    std::vector<PointHistory<dim> >               quadrature_point_history;
    std::vector< std::vector< Vector<double> > >  history_field ;
    std::vector< Vector<double> >                 elastic_stress;

    const RefinementMode                          refinement_mode;

    const Model                                   model;

    const unsigned int                            n_initial_global_refinement,
                                                  n_pre_refinment;
    unsigned int                                  step;

    double                                        time_step,
                                                  total_time;

    std::vector<double>                           surface_q_points,
                                                  surface_analytical_solution,
                                                  surface_numerical_solution;
    double                                        error;
  };


  template <int dim>
  ApShear<dim>::ApShear (const FiniteElement<dim> &fe,
                         const FiniteElement<dim> &history_fe,
                         const unsigned int       degree,
                         const RefinementMode     refinement_mode,
                         const                    Model model) :
    degree (degree),
    dof_handler (triangulation),
    history_dof_handler (triangulation),
    quadrature_formula(degree+1),
    fe (&fe),
    history_fe (&history_fe),
    history_field (3, std::vector< Vector<double> >(3)),
    elastic_stress (dim),
    refinement_mode (refinement_mode),
    model (model),
    n_initial_global_refinement (2),
    n_pre_refinment ((refinement_mode==global_refinement)?2:5),
    step (0),
    time_step (1e-6),
    total_time (0),
    error (0)
  {}


  template <int dim>
  ApShear<dim>::~ApShear ()
  {
    dof_handler.clear ();
  }


  /**
   * The problem is solved in several steps:
   *   - Initial time step (do_first_time_step()):
   *     -# Create an initial uniform grid and set boundary conditions
   *        (make_grid_turcotte() and make_grid_savage()).
   *     -# Initialize the history values stored at quadrature points with zeros
   *        (setup_quadrature_point_history()).
   *     -# Setup DoFs and constraints (setup_system()).
   *     -# Assemble the stiffness matrix and RHS (assemble_system()).
   *     -# Solve (solve()).
   *     -# Refine the mesh adaptively (refine_grid()) and repeat 2-5.
   *
   *   - Successive time steps (do_time_step()):
   *     -# Calculate the time step and compute the total displacement (get_time_step_andisplacement()).
   *     -# Update the values stored in the quadrature points using the new solution
   *        (update_quadrature_point_history()).
   *     -# Assemble the stiffness matrix and RHS (setup_system()).
   *     -# Solve (solve()).
   *     -# Every several time steps refine the grid.
   *       + a) Transfer the quadrature point history to a Discontinuous Galerkin field
   *            (qpoints_to_DG()).
   *       + b) Refine the mesh and interpolate the old solution and stress in the new mesh
   *            (refine_grid()).
   *       + c) Transfer the data stored in the DG field to the new quadrature points
   *          (DG_to_qpoints()).
   */
  template <int dim>
  void ApShear<dim>::run ()
  {
    do_first_time_step ();

    for (step=1; step<21; ++step)
      {
        do_time_step ();
     }
    output_results ();
  }


  /**
   * In the first time step we will create a uniform grid and set the boundary conditions depending
   * on the model that we want to solve. Then we will solve the problem and refine the mesh adaptively
   * several times
   */
  template <int dim>
  void ApShear<dim>::do_first_time_step ()
  {
    std::cout << "Step 0, t = 0" << std::endl;

    std::cout << "Refining initial mesh...";
    for (unsigned int refinement_step=0; refinement_step<=n_pre_refinment; ++refinement_step)
      {
        if (refinement_step==0)
          switch (model)
          {
          case turcotte:
            make_grid_turcotte ();
            break;
          case savage:
            make_grid_savage ();
            break;
          default:
            Assert (false, ExcNotImplemented());
          }
        else
          refine_grid (n_initial_global_refinement,
                       n_initial_global_refinement + n_pre_refinment);

        assemble_system ();
        solve ();
      }
    std::cout << triangulation.n_active_cells() << " cells." << std::endl;
    output_results ();
  }


  /**
   * In successive time steps we will update the visco-elastic stress (using the
   * solution from the last two time steps) and we will solve again. The mesh will
   * be refined adaptively every several time steps.
   *
   */
  template <int dim>
  void ApShear<dim>::do_time_step ()
  {
    get_time_step_and_displacement ();
    total_time += time_step;
    std::cout << "Step " << step << ", t = " << total_time << std::endl;
    update_quadrature_point_history ();
    assemble_system ();
    solve ();
    qpoints_to_DG ();
    output_results ();

    if (step % 5 == 0)
      {

        std::cout << "Refining grid...";
        qpoints_to_DG ();
        refine_grid (n_initial_global_refinement,
                             n_initial_global_refinement + n_pre_refinment);


        DG_to_qpoints ();
        std::cout << triangulation.n_active_cells() << " cells." << std::endl;
      }
 }


  /**
   * The initial mesh will consist of regular rectangular elements. This mesh will be refined
   * in successive steps using refine_grid(). First a mesh is created with
   * elements of the highest possible quality (as close to square as possible). Then this mesh
   * is refined globally a given number of times.
   *
   * \note If the fault is very small it can be ignored (if it's size is smaller than the
   * element size). If this happens, include more global refinement steps.
   */
  template <int dim>
  void ApShear<dim>::make_grid_turcotte ()
  {
    /* The domain can be rectangle, but we want the cells go be as square as possible. */
    Point<dim> corner_low_left(0.0,-this->height);
    Point<dim> corner_up_right(this->width,0.0);

    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = std::max(1.0,round(this->width/this->height));
    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               repetitions,
                                               corner_low_left,
                                               corner_up_right);

    triangulation.refine_global (n_initial_global_refinement);

    /**
     * Once the mesh has been created and refined, we assign an indicator to each part of the
     * boundary, depending on what kind of boundary it is. There are four possible boundary
     * conditions:
     *
     * - 0) The locked part above the fault, with imposed zero displacement (homogeneous
     * Dirichlet boundary conditions).
     *
     * - 1) Regions of imposed displacement (non-homogeneous Dirichlet
     * boundary conditions).
     *
     * - 2) Free boundaries, with zero imposed stress (homogeneous Neumann boundary
     * conditions).
     *
     * - 3) Imposed stress (non-homogeneous boundary conditions.
     */
    typename Triangulation<dim>::cell_iterator cell = triangulation.begin (),
                                 endc = triangulation.end();
      for (; cell!=endc; ++cell)
        for (unsigned int face_number=0;
             face_number<GeometryInfo<dim>::faces_per_cell;
             ++face_number)
          /** For Turcotte and Spence's model \cite turcotte_spence_74, the
           * locked part (0) will be on the left boundary, from the top to a certain
           * depth \f$d\f$. From that depth until the bottom of the model the
           * boundary will be free (2). The top and bottom boundaries will also
           * be free (2). Finally, the leftmost boundary will have an imposed
           * stress \f$S\f$ (3).
           */
          if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Locked part
              &&
              (cell->face(face_number)->center()(1) > -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (0);
          else if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Fault
              &&
              (cell->face(face_number)->center()(1) <= -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (2);
          else if (std::fabs(cell->face(face_number)->center()(1) - (0)) < 1e-12) // Top
            cell->face(face_number)->set_boundary_indicator (2);
          else if  (std::fabs(cell->face(face_number)->center()(1) - (-this->height)) < 1e-12) // Bottom
            cell->face(face_number)->set_boundary_indicator (2);
          else if (std::fabs(cell->face(face_number)->center()(0) - (this->width)) < 1e-12) // Left
            cell->face(face_number)->set_boundary_indicator (3);

      setup_system ();
      setup_quadrature_point_history ();
  }


  /**
   * The initial mesh will consist of regular rectangular elements. This mesh will be refined
   * in successive steps using refine_grid(). First a mesh is created with
   * elements of the highest possible quality (as close to square as possible). Then this mesh
   * is refined globally a given number of times.
   *
   * \note If the fault is very small it can be ignored (if it's size is smaller than the
   * element size). If this happens, include more global refinement steps.
   */
  template <int dim>
  void ApShear<dim>::make_grid_savage ()
  {
    /* The domain can be rectangle, but we want the cells go be as square as possible. */
    Point<dim> corner_low_left(0.0,-this->height);
    Point<dim> corner_up_right(this->width,0.0);

    std::vector<unsigned int> repetitions(dim,1);
    repetitions[0] = std::max(1.0,round(this->width/this->height));
    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               repetitions,
                                               corner_low_left,
                                               corner_up_right);

    triangulation.refine_global (n_initial_global_refinement);

    /**
     * Once the mesh has been created and refined, we assign an indicator to each part of the
     * boundary, depending on what kind of boundary it is. There are four possible boundary
     * conditions:
     *
     * - 0) The locked part above the fault, with imposed zero displacement (homogeneous
     * Dirichlet boundary conditions).
     *
     * - 1) Regions of imposed displacement (non-homogeneous Dirichlet
     * boundary conditions).
     *
     * - 2) Free boundaries, with zero imposed stress (homogeneous Neumann boundary
     * conditions).
     *
     * - 3) Imposed stress (non-homogeneous boundary conditions.
     */
    typename Triangulation<dim>::cell_iterator cell = triangulation.begin (),
                                 endc = triangulation.end();
      for (; cell!=endc; ++cell)
        for (unsigned int face_number=0;
             face_number<GeometryInfo<dim>::faces_per_cell;
             ++face_number)

          /** For Savage and Burford model \cite savage_burford_73, the
           * locked part (0) will be on the left boundary, from the top to a certain
           * depth \f$d\f$. From that depth until the bottom of the model we impose
           * a displacement \f$D\f$ (1). The bottom boundary will have the same imposed
           * displacement (1). Finally, the top and left boundaries will be free (2).
           */
          if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Locked part
              &&
              (cell->face(face_number)->center()(1) > -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (0);
          else if ((std::fabs(cell->face(face_number)->center()(0) - (0)) < 1e-12) // Fault
              &&
              (cell->face(face_number)->center()(1) <= -this->locked_depth))
            cell->face(face_number)->set_boundary_indicator (1);
          else if (std::fabs(cell->face(face_number)->center()(1) - (0)) < 1e-12) // Top
            cell->face(face_number)->set_boundary_indicator (2);
          else if  (std::fabs(cell->face(face_number)->center()(1) - (-this->height)) < 1e-12) // Bottom
            cell->face(face_number)->set_boundary_indicator (1);
          else if (std::fabs(cell->face(face_number)->center()(0) - (this->width)) < 1e-12) // Left
            cell->face(face_number)->set_boundary_indicator (2);

      setup_system ();
      setup_quadrature_point_history ();
  }


  template <int dim>
  void ApShear<dim>::refine_grid (const unsigned int min_grid_level,
                                  const unsigned int max_grid_level)
  {
    switch (refinement_mode)
      {
      case global_refinement:
      {
        triangulation.refine_global (1);
        break;
      }

      /**
       * The first step is to identify which cells need to be refined, for which we use we use a
       * Kelly Error Estimator. We will refine 60% of the cells that have the highest error and coarsen
       * the rest).
       */
      case adaptive_refinement:
      {
        Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

        KellyErrorEstimator<dim>::estimate (dof_handler,
                                            QGauss<dim-1>(3),
                                            typename FunctionMap<dim>::type(),
                                            solution,
                                            estimated_error_per_cell);

        GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                         estimated_error_per_cell,
                                                         0.6, 0.4);
        /**
         * We do not want to make any cell too small, because that would mean (according to
         * the Courant–Friedrichs–Lewy condition) that we have to make the time steps very small.
         * Therefore, we will mark those cell that have already been refined a certain number of
         * times so they are not refined any further. a similar thing happens with cells that
         * are too big, therefore, we will also mark the cells that have reached a certain maximum
         * size so they are not coarsened any more.
         */

        if (triangulation.n_levels() > max_grid_level)
          for (typename Triangulation<dim>::active_cell_iterator
               cell = triangulation.begin_active(max_grid_level);
               cell != triangulation.end(); ++cell)
            cell->clear_refine_flag ();
        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active(min_grid_level);
             cell != triangulation.end_active(min_grid_level); ++cell)
          cell->clear_coarsen_flag ();

        /**
         * Here we also need to transfer the #old_solution, the #solution and the #history_field
         * (where the old stress is stored) to the new mesh. First we have to prepare the
         * vectors that will be transferred to the new grid (the old grid will disappear
         * once we have done the refinement so we have to transfer the data at the same
         * time as the refinement).
         */

        SolutionTransfer<dim> solution_trans(dof_handler);
        std::vector<Vector<double> > unrefined_solution (3, Vector<double>(dof_handler.n_dofs()));
        unrefined_solution[0] = solution;
        unrefined_solution[1] = old_solution;
        unrefined_solution[2] = displacement;

        SolutionTransfer<dim> history_trans(history_dof_handler);
        std::vector<Vector<double> > unrefined_history_field (9+dim, Vector<double>(history_dof_handler.n_dofs()));
        for (unsigned int i=0; i<3; ++i)
          for (unsigned int j=0; j<3; ++j)
            unrefined_history_field [3 * i + j] = history_field[i][j];
        for (unsigned int i=0; i<dim; ++i)
          unrefined_history_field [9+i] =  elastic_stress[i];

        triangulation.prepare_coarsening_and_refinement();
        solution_trans.prepare_for_coarsening_and_refinement(unrefined_solution);
        history_trans.prepare_for_coarsening_and_refinement(unrefined_history_field);

        /**
         * Once everything is ready we can actually refine and coarsen the mesh.
         */
        triangulation.execute_coarsening_and_refinement ();

        /**
         * Once the new mesh is ready, we need to resent the system and the data at the quadrature
         * points.
         */
        setup_system ();
        setup_quadrature_point_history ();

        /**
         * Finally, we can interpolate the #old_solution, the #solution and the
         * #history_field in the new mesh.
         */
        std::vector<Vector<double> > refined_solution (3, Vector<double>(dof_handler.n_dofs()));
        solution_trans.interpolate(unrefined_solution, refined_solution);
        solution = refined_solution[0];
        old_solution = refined_solution[1];
        displacement = refined_solution[2];

        std::vector<Vector<double> > refined_history_field (9+dim, Vector<double>(history_dof_handler.n_dofs()));
        history_trans.interpolate(unrefined_history_field, refined_history_field);
        for (unsigned int i=0; i<3; ++i)
          for (unsigned int j=0; j<3; ++j)
            history_field[i][j] = refined_history_field[3 * i + j];
        for (unsigned int i=0; i<dim; ++i)
          elastic_stress[i] = refined_history_field[9+i];

        break;
      }

      default:
      {
        Assert (false, ExcNotImplemented());
      }
      }
  }


  /**
   * As part of the refinement process, we need to store the old stress in a FE field defined
   * everywhere. Once that has been done, we will be able to transfer the old stress (or any
   * data stored in the quadrature point using PointHistory) from the old mesh to the refined
   * mesh. For more details see
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_18.html#Refinementduringtimesteps">Step-18</a>
   */
  template <int dim>
  void ApShear<dim>::qpoints_to_DG ()
  {
    std::vector< std::vector< Vector<double> > >
                 local_history_values_at_qpoints (3, std::vector< Vector<double> >(3)),
                 local_history_fe_values (3, std::vector< Vector<double> >(3));
    for (unsigned int i=0; i<3; i++)
      for (unsigned int j=0; j<3; j++)
      {
        local_history_values_at_qpoints[i][j].reinit(quadrature_formula.size());
        local_history_fe_values[i][j].reinit(history_fe->dofs_per_cell);
      }
    std::vector< Vector<double> >
                 local_elastic_stress_at_qpoints (dim),
                 local_elastic_stress_fe_values (dim);
    for (unsigned int i=0; i<dim; i++)
      {
        local_elastic_stress_at_qpoints[i].reinit(quadrature_formula.size());
        local_elastic_stress_fe_values[i].reinit(history_fe->dofs_per_cell);
      }
    FullMatrix<double> qpoint_to_dof_matrix (history_fe->dofs_per_cell,
                                             quadrature_formula.size());
    FETools::compute_projection_from_quadrature_points_matrix
              (*history_fe,
               quadrature_formula, quadrature_formula,
               qpoint_to_dof_matrix);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();
    for (; cell!=endc; ++cell, ++dg_cell)
      {
        PointHistory<dim> *local_quadrature_points_history
               = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
        Assert (local_quadrature_points_history >=
                    &quadrature_point_history.front(),
                    ExcInternalError());
        Assert (local_quadrature_points_history <
                    &quadrature_point_history.back(),
                    ExcInternalError());
        for (unsigned int i=0; i<3; i++)
          for (unsigned int j=0; j<3; j++)
          {
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              local_history_values_at_qpoints[i][j](q)
                       = local_quadrature_points_history[q].old_stress[i][j];
            qpoint_to_dof_matrix.vmult (local_history_fe_values[i][j],
                                        local_history_values_at_qpoints[i][j]);
            dg_cell->set_dof_values (local_history_fe_values[i][j],
                                     history_field[i][j]);
          }
        for (unsigned int i=0; i<dim; i++)
          {
            for (unsigned int q=0; q<quadrature_formula.size(); ++q)
              local_elastic_stress_at_qpoints[i](q)
                   = local_quadrature_points_history[q].elastic_stress[i];
            qpoint_to_dof_matrix.vmult (local_elastic_stress_fe_values[i],
                                        local_elastic_stress_at_qpoints[i]);
            dg_cell->set_dof_values (local_elastic_stress_fe_values[i],
                                     elastic_stress[i]);
          }
      }
  }


  /**
   * After the mesh has been refined, we need to transfer the old stress back to the
   * quadrature points. For more details see
   * <a href="https://www.dealii.org/8.2.0/doxygen/deal.II/step_18.html#Refinementduringtimesteps">Step-18</a>
   */
  template <int dim>
  void ApShear<dim>::DG_to_qpoints ()
  {
    std::vector< std::vector< Vector<double> > >
                   local_history_values_at_qpoints (3, std::vector< Vector<double> >(3)),
                   local_history_fe_values (3, std::vector< Vector<double> >(3));
      for (unsigned int i=0; i<3; i++)
        for (unsigned int j=0; j<3; j++)
        {
          local_history_values_at_qpoints[i][j].reinit(quadrature_formula.size());
          local_history_fe_values[i][j].reinit(history_fe->dofs_per_cell);
        }

    std::vector<Vector<double> >
                       local_elastic_stress_at_qpoints (dim),
                       local_elastic_stress_fe_values (dim);
          for (unsigned int i=0; i<dim; i++)
            {
              local_elastic_stress_at_qpoints[i].reinit(quadrature_formula.size());
              local_elastic_stress_fe_values[i].reinit(history_fe->dofs_per_cell);
            }

    FullMatrix<double> dof_to_qpoint_matrix (quadrature_formula.size(),
                                             history_fe->dofs_per_cell);
    FETools::compute_interpolation_to_quadrature_points_matrix
              (*history_fe,
               quadrature_formula,
               dof_to_qpoint_matrix);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();
    for (; cell != endc; ++cell, ++dg_cell)
    {
     PointHistory<dim> *local_quadrature_points_history
             = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
      Assert (local_quadrature_points_history >=
                  &quadrature_point_history.front(),
                  ExcInternalError());
      Assert (local_quadrature_points_history <
                  &quadrature_point_history.back(),
                  ExcInternalError());
      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
        {
          dg_cell->get_dof_values (history_field[i][j],
                                   local_history_fe_values[i][j]);
          dof_to_qpoint_matrix.vmult (local_history_values_at_qpoints[i][j],
                                      local_history_fe_values[i][j]);
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            local_quadrature_points_history[q].old_stress[i][j]
                       = local_history_values_at_qpoints[i][j](q);
        }
      for (unsigned int i=0; i<dim; ++i)
        {
          dg_cell->get_dof_values (elastic_stress[i],
                                   local_elastic_stress_fe_values[i]);
          dof_to_qpoint_matrix.vmult (local_elastic_stress_at_qpoints[i],
                                      local_elastic_stress_fe_values[i]);
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            local_quadrature_points_history[q].elastic_stress[i]
                       = local_elastic_stress_at_qpoints[i](q);
        }
    }
  }


  /**
   * Here we set up DoFs, renumber them for efficiency and set required constraints.
   */
  template <int dim>
  void ApShear<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (*fe);
    history_dof_handler.distribute_dofs (*history_fe);

    /**
     * Renumber DoFs to increase efficiency of the solver. In this case, it is not really necessary
     * for the solver and preconditioners that are being used
     */
    DoFRenumbering::Cuthill_McKee (dof_handler);

    /**
     * Set constraints for hanging nodes (due to adaptive refinement) and for homogeneous and
     * non-homogeneous BC (indicators 0 and 1 respectively).
     */
    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              1,
                                              Dirichlet_BC<dim>(),
                                              constraints);
    constraints.close ();


    CompressedSparsityPattern c_sparsity(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                    c_sparsity,
                                    constraints,
                                    /*keep_constrained_dofs = */ false);

    sparsity_pattern.copy_from(c_sparsity);

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    old_solution.reinit (dof_handler.n_dofs());
    displacement.reinit (dof_handler.n_dofs());
    for (unsigned int i=0; i<3; i++)
          for (unsigned int j=0; j<3; j++)
            history_field[i][j].reinit(history_dof_handler.n_dofs());
    for (unsigned int i=0; i<dim;  i++)
      elastic_stress[i].reinit(history_dof_handler.n_dofs());

    system_rhs.reinit (dof_handler.n_dofs());
  }


  template <int dim>
  void ApShear<dim>::setup_quadrature_point_history ()
  {

    triangulation.clear_user_data();

    /* Only necesary if implementing parallel code (the number of quadrature points
     * might change between processors
    {
      std::vector<PointHistory<dim> > tmp;
      tmp.swap (quadrature_point_history);
    }
    */
    quadrature_point_history.resize (triangulation.n_active_cells() *
                                     quadrature_formula.size());
    unsigned int history_index = 0;
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        cell->set_user_pointer (&quadrature_point_history[history_index]);
        history_index += quadrature_formula.size();
      }
    Assert (history_index == quadrature_point_history.size(),
            ExcInternalError());
  }

  template <int dim>
  void ApShear<dim>::update_quadrature_point_history ()
  {
    FEValues<dim>                 fe_values (*fe, quadrature_formula,
                                             update_quadrature_points | update_values |
                                             update_gradients);

    const Elastic_Modulus<dim>    elastic_modulus;
    std::vector<double>           em_values (quadrature_formula.size());
    const Viscosity<dim>          viscosity;
    std::vector<double>           visc_values (quadrature_formula.size());

    std::vector<Tensor<1,dim> >   new_solution_gradients (quadrature_formula.size());
    std::vector<Tensor<1,dim> >   old_solution_gradients (quadrature_formula.size());

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
    {
      PointHistory<dim> *local_quadrature_points_history
        = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

      fe_values.reinit (cell);

      elastic_modulus.value_list(fe_values.get_quadrature_points(),em_values);
      viscosity.value_list(fe_values.get_quadrature_points(),visc_values);

      fe_values.get_function_gradients (solution,
                                        new_solution_gradients);
      fe_values.get_function_gradients (old_solution,
                                        old_solution_gradients);

      for (unsigned int q=0; q<quadrature_formula.size(); ++q)
        {
          const double effective_viscosity = visc_values[q] *
                                             em_values[q] *
                                             time_step /
                                             (em_values[q] * time_step + visc_values[q]);


          Tensor<2,3> new_stress;
          Tensor<1,dim> g_e_z;

          Tensor<2,3> old_stress = local_quadrature_points_history[q].old_stress;

          //The values of the strain rate and the spin are only valid under the anti-plane shear approximation
          Tensor<2,3> new_strain_rate;
          new_strain_rate[0][2] = new_solution_gradients[q][0];
          new_strain_rate[1][2] = new_solution_gradients[q][1];
          new_strain_rate[2][0] = new_solution_gradients[q][0];
          new_strain_rate[2][1] = new_solution_gradients[q][1];

          Tensor<2,3> old_spin;
          old_spin[0][2] = -old_solution_gradients[q][0];
          old_spin[1][2] = -old_solution_gradients[q][1];
          old_spin[2][0] = old_solution_gradients[q][0];
          old_spin[2][1] = old_solution_gradients[q][1];

          Tensor<2,3> new_spin;
          new_spin[0][2] = -new_solution_gradients[q][0];
          new_spin[1][2] = -new_solution_gradients[q][1];
          new_spin[2][0] = new_solution_gradients[q][0];
          new_spin[2][1] = new_solution_gradients[q][1];


          for (unsigned int i=0; i < 3; ++i)
            for (unsigned int j=0; j < 3; ++j)
              {
                new_stress[i][j] = 2 * effective_viscosity * new_strain_rate[i][j] +
                                       effective_viscosity * 1/(em_values[q] * time_step) * old_stress[i][j];
                for (unsigned int k=0; k < 3; ++k)
                  {
                    new_stress[i][j] += effective_viscosity * 1/em_values[q] *
                                        (old_spin[i][k] * old_stress[k][j] -
                                         old_stress[i][k] * old_spin[k][j]);
                  }
              }
          for (unsigned int l=0; l<dim; ++l)
              {
                g_e_z[l] = (effective_viscosity *
                            new_stress[l][2] /
                            em_values[q] /
                            time_step);
                for (unsigned int k=0; k<3; ++k)
                  g_e_z[l] += (effective_viscosity *
                               ((new_spin[l][k] *
                                 new_stress[k][2]) -
                                (new_stress[l][k] *
                                 new_spin[k][2])) /
                               em_values[q]);
              }
          local_quadrature_points_history[q].old_stress = new_stress;
          local_quadrature_points_history[q].elastic_stress = g_e_z;
        }
    }
  }


  /**
   * The stiffness matrix and RHS are assembled here.
   */
  template <int dim>
  void ApShear<dim>::assemble_system ()
  {
    QGauss<dim-1> face_quadrature_formula(degree+1);

    const unsigned int n_q_points    = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const unsigned int dofs_per_cell = fe->dofs_per_cell;

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    FEValues<dim>  fe_values (*fe, quadrature_formula,
                              update_values   | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                      update_values         | update_quadrature_points  |
                                      update_normal_vectors | update_JxW_values);

    const Elastic_Modulus<dim>   elastic_modulus;
    std::vector<double>          em_values (n_q_points);
    const Viscosity<dim>         viscosity;
    std::vector<double>          visc_values (n_q_points);
    const BodyForce<dim>         body_force;
    std::vector<double>          bf_values (n_q_points);
    const Neumann_BC<dim>        neumann_bc;
    std::vector<double>          nbc_values (n_face_q_points);
    std::vector<Tensor<1,dim> >  old_solution_gradients (n_q_points);

    old_solution = solution;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    /**
     * To assemble the stiffness matrix and RHS for the given problem we assemble
     * the local matrices for each element and the transfer it to the global
     * matrix.
     */
    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit (cell);

        PointHistory<dim> *local_quadrature_points_history
                          = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());

        elastic_modulus.value_list (fe_values.get_quadrature_points(),em_values);
        viscosity.value_list (fe_values.get_quadrature_points(),visc_values);
        body_force.value_list (fe_values.get_quadrature_points(),bf_values);
        fe_values.get_function_gradients (old_solution, old_solution_gradients);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            /*
             * To compute the contribution of the internal elastic force in the right hand
             * side of the weak form, we need the effective viscosity, the deviatoric stress
             * from the previous time step (stored at each quadrature point) and the spin
             * from the previous time step (depends on the partial derivatives of the solution).
             */
            const double effective_viscosity = visc_values[q] *
                                               em_values[q] *
                                               time_step /
                                               (em_values[q] * time_step + visc_values[q]);

            Tensor<2,3> old_stress = local_quadrature_points_history[q].old_stress;

            Tensor<2,3> old_spin;
            old_spin[0][2] = -old_solution_gradients[q][0];
            old_spin[1][2] = -old_solution_gradients[q][1];
            old_spin[2][0] = old_solution_gradients[q][0];
            old_spin[2][1] = old_solution_gradients[q][1];

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (em_values [q] *
                                       fe_values.shape_grad(i,q) *
                                       fe_values.shape_grad(j,q)*
                                       fe_values.JxW(q));

                Tensor <1,dim> g_e_z;
                for (unsigned int l=0; l<dim; ++l)
                  {
                    g_e_z[l] = (effective_viscosity *
                                old_stress[l][2] /
                                em_values[q] /
                                time_step);
                    for (unsigned int k=0; k<3; ++k)
                      g_e_z[l] += (effective_viscosity *
                                   ((old_spin[l][k] *
                                     old_stress[k][2]) -
                                    (old_stress[l][k] *
                                     old_spin[k][2])) /
                                   em_values[q]);
                  }
                cell_rhs(i) += -((fe_values.shape_grad(i,q)) *
                                  g_e_z *
                                 fe_values.JxW(q)) -
                               (fe_values.shape_value(i,q) *
                                bf_values [q] *
                                fe_values.JxW(q));
              }
          }

        for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
          /**
           * Apply homogeneous Neumann boundary conditions. We loop over all the faces in
           * each cell. We check if that face is in a boundary and if it is a boundary of type 3
           * (non-homogeneous Neumann boundary conditions), to those faces we add the corresponing
           * contribution.
           */
          if (cell->face(face_number)->at_boundary()
              &&
              (cell->face(face_number)->boundary_indicator() == 3))
            {
              fe_face_values.reinit (cell, face_number);
              neumann_bc.value_list (fe_face_values.get_quadrature_points(), nbc_values);

              for (unsigned int q=0; q<n_face_q_points; ++q)
                  for (unsigned int i=0; i<dofs_per_cell; ++i)
                    cell_rhs(i) += (fe_face_values.shape_value(i,q) *
                                    nbc_values[q] *
                                    fe_face_values.JxW(q));
            }

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (cell_matrix,
                                                cell_rhs,
                                                local_dof_indices,
                                                system_matrix,
                                                system_rhs);
      }
  }


  /**
   * The strategy to invert the stiffness matrix and solve the system is set here.
   * Because the problem is simple, we just use a Conjugate Gradient solver with the identity as
   * a preconditioner.
   */
  template <int dim>
  void ApShear<dim>::solve ()
  {
    SolverControl           solver_control (2000, 1e-12);
    SolverCG<>              solver (solver_control);

    solver.solve (system_matrix, solution, system_rhs,
                  PreconditionIdentity());

    constraints.distribute (solution);
  }


  /**
   * We need to calculate the size of each time step. We choose a size such that, if
   * the domain was three-dimensional, the distance that a particle can travel (in the direction
   * perpendicular to the modeled plane) within one time step is no larger than
   * the diameter of a single cell:
   * \f[ \Delta t = \frac{min(d)}{max(v_z)} \f]
   *
   * The function GridTools::minimal_cell_diameter computes the minimal diameter of all cells.
   * Since the cells are all squares, the minimal edge length is be the minimal
   * diameter divided by std::sqrt(dim+1.0).
   *
   * The maximal velocity is obtained by looping through every quadrature point.
   */
  template <int dim>
  void ApShear<dim>::get_time_step_and_displacement ()
  {
    Vector<double> inc_displacement;
    inc_displacement = solution;
    inc_displacement *= time_step;
    displacement += inc_displacement;

    const unsigned int   n_q_points
      = quadrature_formula.size();
    FEValues<dim> fe_values (*fe, quadrature_formula,
                             update_values);
    std::vector<double> solution_values(n_q_points);
    double max_velocity = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        fe_values.get_function_values (solution, solution_values);
        for (unsigned int q=0; q<n_q_points; ++q)
          {
            max_velocity = std::max (max_velocity,
                                     std::fabs(solution_values[q]));
          }
      }
    time_step = GridTools::minimal_cell_diameter(triangulation) /
                max_velocity /
                std::sqrt(dim + 1.0);
//    time_step = 1.0e3 * this->viscosity / this->elastic_modulus;
  }


  /**
   * This function will take care of computing and storing the difference between the numerical
   * solution and the analytical solution from Turcotte and Spence Model \cite turcotte_spence_74.
   * At the end we will output the results in several graphs. The comparison will stored in three
   * different variables:
   * - analytical_solution: Analytical solution computed at each quadrature point of the final
   * mesh.
   * - surface_numerical_solution: Numerical solution in each quadrature point at each refinement cycle.
   * - integrated_diff: Integrated difference between both solutions using the L1 norm.
   */
  template <int dim>
  void ApShear<dim>::compare_solutions_turcotte ()
  {

    QGauss<dim-1>   face_quadrature_formula(degree+1);
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_face_q_points    = face_quadrature_formula.size();

    const Solution_Turcotte<dim> exact_solution;

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                         update_values         | update_quadrature_points  |
                                         update_normal_vectors | update_JxW_values);
    /**
     * First we need to run a loop over every cell and every quadrature point to determine if
     * they are located at the surface, if so, we will compute the analytical solution and extract
     * the numerical solution to store them in analytical_solution and numerical_solution.
     * Once we have both solutions at the surface, we will compute the integrated
     * difference between them and we will store it in integrated_diff.
     */
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int face_number=0; face_number<n_faces; ++face_number)
        if (std::fabs(cell->face(face_number)->center()(1)) < 1e-12)
          {
            fe_face_values.reinit (cell, face_number);
            for (unsigned int q=0; q<n_face_q_points; ++q)
              {
                surface_q_points.push_back(fe_face_values.quadrature_point(q)[0]);
                surface_analytical_solution.push_back(exact_solution.value(fe_face_values.quadrature_point(q)));
                surface_numerical_solution.push_back(VectorTools::point_value (dof_handler, solution,
                                                                    fe_face_values.quadrature_point(q)));
                error += (std::fabs(surface_analytical_solution.back()-surface_numerical_solution.back())*
                                                    fe_face_values.JxW(q));

              }
          }
  }


  /**
   * This function will take care of computing and storing the difference between the numerical
   * solution and the analytical solution from Turcotte and Spence Model \cite turcotte_spence_74.
   * At the end we will output the results in several graphs. The comparison will stored in three
   * different variables:
   * - analytical_solution: Analytical solution computed at each quadrature point of the final
   * mesh.
   * - surface_numerical_solution: Numerical solution in each quadrature point at each refinement cycle.
   * - integrated_diff: Integrated difference between both solutions using the L1 norm.
   */
  template <int dim>
  void ApShear<dim>::compare_solutions_savage ()
  {

    QGauss<dim-1>   face_quadrature_formula(degree+1);
    const unsigned int n_faces = GeometryInfo<dim>::faces_per_cell;
    const unsigned int n_face_q_points    = face_quadrature_formula.size();

    const Solution_Savage<dim> exact_solution;

    FEFaceValues<dim> fe_face_values (*fe, face_quadrature_formula,
                                         update_values         | update_quadrature_points  |
                                         update_normal_vectors | update_JxW_values);
    /**
     * First we need to run a loop over every cell and every quadrature point to determine if
     * they are located at the surface, if so, we will compute the analytical solution and extract
     * the numerical solution to store them in analytical_solution and numerical_solution.
     * Once we have both solutions at the surface, we will compute the integrated
     * difference between them and we will store it in integrated_diff.
     */
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      for (unsigned int face_number=0; face_number<n_faces; ++face_number)
        if (std::fabs(cell->face(face_number)->center()(1)) < 1e-12)
          {
            fe_face_values.reinit (cell, face_number);
            for (unsigned int q=0; q<n_face_q_points; ++q)
              {
                surface_q_points.push_back(fe_face_values.quadrature_point(q)[0]);
                surface_analytical_solution.push_back(exact_solution.value(fe_face_values.quadrature_point(q)));
                surface_numerical_solution.push_back(VectorTools::point_value (dof_handler, solution,
                                                                    fe_face_values.quadrature_point(q)));
                error += (std::fabs(surface_analytical_solution.back()-surface_numerical_solution.back())*
                                                    fe_face_values.JxW(q));

              }
          }
  }


  /**
   * The results will be output in this step. Three different outputs are generated.
   * - The final mesh in eps format.
   * - The final displacement in the 2D domain in vtk format
   * - Information regarding the conparison between the numerical and the analytical solutions.
   */
  template <int dim>
  void ApShear<dim>::output_results ()
  {
    std::ostringstream filename;
    // The output file name depends on the refinement mode and the element type.
    switch (model)
      {
      case turcotte:
        filename << "solution-turcotte";
        break;
      case savage:
        filename << "solution-savage";
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    switch (refinement_mode)
      {
      case global_refinement:
        filename << "-global";
        break;
      case adaptive_refinement:
        filename << "-adaptive";
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    switch (fe->degree)
      {
      case 1:
        filename << "-q1";
        break;
      case 2:
        filename << "-q2";
        break;

      default:
        Assert (false, ExcNotImplemented());
      }

    filename << "(step" << step << ")";

    // Output final mesh
    std::ostringstream eps_filename;
    eps_filename << filename.str() << ".eps";
    std::ofstream output_mesh (eps_filename.str().c_str());

    GridOut grid_out;
    grid_out.write_eps (triangulation, output_mesh);

    // Output numerical solution
    std::ostringstream vtk_filename;
    vtk_filename << filename.str() << ".vtk";
    std::ofstream output_num_sol (vtk_filename.str().c_str());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "velocity");
    data_out.add_data_vector (displacement, "displacement");

    data_out.build_patches (fe->degree);
    data_out.write_vtk (output_num_sol);

    // Output point history
    std::ostringstream vtk_filename_history;
    vtk_filename_history << filename.str() << "_history.vtk";
    std::ofstream output_history (vtk_filename_history.str().c_str());

    DataOut<dim> data_out_history;
    data_out_history.attach_dof_handler (history_dof_handler);
    data_out_history.add_data_vector (elastic_stress[0], "elastic_stress_xz");
    data_out_history.add_data_vector (elastic_stress[1], "elastic_stress_yz");
    data_out_history.add_data_vector (history_field[0][2], "history_field_xz");
    data_out_history.add_data_vector (history_field[1][2], "history_field_yz");
//    data_out_history.add_data_vector (history_field[2][2], "history_field_zz");

    data_out_history.build_patches (history_fe->degree);
    data_out_history.write_vtk (output_history);


//    // Output comparison
//    std::ostringstream txt_filename_error;
//    txt_filename_error << filename.str() << "-error" << ".txt";
//    std::ofstream out_compare_error(txt_filename_error.str().c_str());
//
//    std::ostringstream txt_filename_analytical;
//    txt_filename_analytical << filename.str() << "-analytical" << ".txt";
//    std::ofstream out_compare_analytical (txt_filename_analytical.str().c_str());
//
//    std::ostringstream txt_filename_numerical;
//    txt_filename_numerical << filename.str() << "-numerical" << ".txt";
//    std::ofstream out_compare_numerical(txt_filename_numerical.str().c_str());
//
//    const int n_q_points = surface_q_points.size ();
//    for (unsigned int q = 0; q<n_q_points; ++q)
//      {
//        out_compare_analytical << surface_q_points [q] << "\t";
//        out_compare_analytical << surface_analytical_solution[q] << std::endl;
//
//        out_compare_numerical << surface_q_points [q] << "\t";
//        out_compare_numerical << surface_numerical_solution [q] << std::endl;
//      }
//
//    out_compare_error << error << std::endl;
  }
}





int main ()
{
  const unsigned int dim = 2;
  const unsigned int degree = 1;

  try
    {
      using namespace dealii;
      using namespace vsf;

      deallog.depth_console (0);
      /**
       * \note The equations are only valid for 2D geometries (anti-plane shear approximation). If the
       * user tries to create a 3D geometry, an exception is thrown and the program stops.
       */
      Assert (dim == 2, ExcNotImplemented ());

      /**
       * When an ApShear is created the user must specify the kind of elements and
       * refinement (adaptive_refinement or global_refinement) and Model (turcotte or savage)
       * are going oto be used.
       */
      {
        std::cout << "Solving with Q1 elements, adaptive refinement for Turcotte and Spence Model" << std::endl
                  << "==========================================================================" << std::endl
                  << std::endl;

        FE_Q<dim> fe(degree);
        FE_DGQ<dim> history_fe(degree);
        ApShear<dim>
        antiplane_shear_problem_2d (fe,
                                    history_fe,
                                    degree,
                                    ApShear<dim>::adaptive_refinement,
                                    ApShear<dim>::turcotte);

        antiplane_shear_problem_2d.run ();

        std::cout << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
