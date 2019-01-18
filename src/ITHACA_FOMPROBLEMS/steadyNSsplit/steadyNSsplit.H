/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------
License
    This file is part of ITHACA-FV
    ITHACA-FV is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    ITHACA-FV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.
    You should have received a copy of the GNU Lesser General Public License
    along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.
Class
    steadyNSsplit
Description
    Reduction class for a Stationary turbulent Navier-Stokes problem.
SourceFiles
    steadyNSsplit.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the steadyNSsplit class.
/// \dir
/// Directory containing the header and source files for the steadyNSsplit class.

#ifndef steadyNSsplit_H
#define steadyNSsplit_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "reductionProblem.H"
#include "ITHACAstream.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class steadyNSsplit Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> steady turbulent Navier Stokes problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a turbulent steady NS problem and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class steadyNSsplit: public steadyNS
{


    public:
        // Constructors
        steadyNSsplit();
        /// Construct with argc and argv
        steadyNSsplit(int argc, char* argv[]);

        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nuTmodes;

        // Create a SAMPLES for interpolation
        std::vector<SPLINTER::DataTable*> SAMPLES;
        // Create a RBF splines for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfsplines;

        /// List of vector fields for storing the x component of the velocity fields.
        PtrList<volVectorField> Ufieldx;

        /// List of vector fields for storing the y component of the velocity fields.
        PtrList<volVectorField> Ufieldy;

        /// List of vector fields for storing the z component of the velocity fields.
        PtrList<volVectorField> Ufieldz;


        /// List of vector fields for storing the modes corresponding to the x component of the velocity fields.
        PtrList<volVectorField> Umodesx;

        /// List of vector fields for storing the modes corresponding to the y component of the velocity fields.
        PtrList<volVectorField> Umodesy;

        /// List of vector fields for storing the modes corresponding to the z component of the velocity fields.
        PtrList<volVectorField> Umodesz;

        /// List of vector fields for storing the modes corresponding to the x component of the supremizer fields.
        PtrList<volVectorField> supmodesx;

        /// List of vector fields for storing the modes corresponding to the y component of the supremizer fields.
        PtrList<volVectorField> supmodesy;

        /// List of vector fields for storing the modes corresponding to the z component of the supremizer fields.
        PtrList<volVectorField> supmodesz;


        /// 3-D Vector to indicate the directions for which we have to consider doing the POD, 1 means consider and 0 means ignore
        Vector<double> split_directions;

        /// Number of velocity modes along the x direction used for the projection
        label NUmodesx;

        /// Number of velocity modes along the y direction used for the projection
        label NUmodesy;

        /// Number of velocity modes along the z direction used for the projection
        label NUmodesz;

        /// Number of supremizer modes along the x direction used for the projection
        label NSUPmodesx;

        /// Number of supremizer modes along the y direction used for the projection
        label NSUPmodesy;

        /// Number of supremizer modes along the z direction used for the projection
        label NSUPmodesz;

        /// Boundary term for penalty method - vector
        List <Eigen::MatrixXd> Q_vector;

        /// Boundary term for penalty method - matrix
        List <Eigen::MatrixXd> E_matrix;



        /** @name Reduced Matrices
        *
        */
        ///@{
        //
        /// Turbulent viscosity matrix
        Eigen::MatrixXd BT_matrix;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> CT2_matrix;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> CT1_matrix;

        /// Total C Matrix
        List <Eigen::MatrixXd> C_total_matrix;

        /// Total B Matrix
        Eigen::MatrixXd B_total_matrix;
        ///@}

        /// Number of viscoisty modes used for the projection
        label Nnutmodes;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;
        autoPtr<volScalarField> _nuTilda;

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);


        /// @brief      Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder used to save the reduced matrices.
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        void projectSUP(fileName folder, label NUx, label NUy, label NUz, label NP, label NSUP, label NSUPx,
        label NSUPy, label NSUPz, label Nnut);

        //--------------------------------------------------------------------------
        //  Projection Methods Momentum Equation
        /// Diffusive Term
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term.
        ///
        Eigen::MatrixXd diffusive_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        //--------------------------------------------------------------------------
        /// Gradient of pressure
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the Gradient of pressure term.
        ///
        Eigen::MatrixXd pressure_gradient_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        //--------------------------------------------------------------------------
        /// Convective Term
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced third order tensor in List <Eigen::MatrixXd> format for the convective term.
        ///
        List < Eigen::MatrixXd >  convective_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        //--------------------------------------------------------------------------
        /// Mass Term
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the mass matrix.
        ///
        Eigen::MatrixXd  mass_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        // Projection Methods Continuity Equation

        //--------------------------------------------------------------------------
        /// Divergence Term (supremizer approach)
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     reduced matrix in Eigen::MatrixXd format for the divergence term.
        ///
        Eigen::MatrixXd divergence_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        
        
        /// @brief      BT added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        ///
        /// @return     BT matrix for turbulence treatment
        ///
        Eigen::MatrixXd BT_turbulence(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);


        ///
        /// @brief      CT1 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     CT1 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulence_term1(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes,
    label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes);

        ///
        /// @brief      CT2 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
        ///
        /// @return     CT2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulence_term2(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes,
    label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes);

        
        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        ///
        /// @return     { description_of_the_return_value }
        ///
        List< Eigen::MatrixXd > Boundary_vector_list(label NUmodesx, label NUmodesy, label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

        ///
        /// @brief      { function_description }
        ///
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        /// @param[in]  <unnamed>  { parameter_description }
        ///
        /// @return     { description_of_the_return_value }
        ///
        List< Eigen::MatrixXd > Boundary_matrix_list(label NUmodesx, label NUmodesy, label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);

};

#endif









