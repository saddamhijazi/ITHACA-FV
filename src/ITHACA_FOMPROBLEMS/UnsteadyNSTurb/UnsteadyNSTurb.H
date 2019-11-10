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
    UnsteadyNSTurb
Description
    Reduction class for a NON-Stationary turbulent NS problem
SourceFiles
    UnsteadyNSTurb.C

\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the UnsteadyNSTurb class.
/// \dir
/// Directory containing the header and source files for the UnsteadyNSTurb class.

#ifndef UnsteadyNSTurb_H
#define UnsteadyNSTurb_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "IOporosityModelList.H"
#include "IOMRFZoneList.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "steadyNS.H"
#include "unsteadyNS.H"
#include <iostream>
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include <rbfspline.h>
#include <spline.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyNS Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady NS problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady NS problem
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the steadyNS class */
class UnsteadyNSTurb: public unsteadyNS
{
    public:
        // Constructors
        /// Construct Null
        UnsteadyNSTurb();

        /// Construct with argc and argv
        UnsteadyNSTurb(int argc, char* argv[]);

        // Specific variable for the unstationary case
        // Member Functions
        /// List of snapshots for the solution for eddy viscosity
        PtrList<volScalarField> nutFields;

        /// List of POD modes for eddy viscosity
        PtrList<volScalarField> nutModes;

        /// Create a Rbf splines for interpolation
        std::vector<SPLINTER::DataTable*> samplesInit;

        /// Create a Rbf splines for interpolation
        std::vector<SPLINTER::DataTable*> samples;

        /// Create a Rbf splines for interpolation
        std::vector<SPLINTER::DataTable*> samplesTimePar;

        /// Create a samples for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplines;

        /// Create a samples for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplinesTimePar;

        /// Create a samples for interpolation
        std::vector<SPLINTER::RBFSpline*> rbfSplinesInit;

        /// Turbulent viscosity term
        Eigen::MatrixXd btMatrix;

        /// Time-parameter combined matrix
        Eigen::MatrixXd z;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct2Matrix;
        Eigen::Tensor<double, 3 > ct2Tensor;

        /// Turbulent viscosity tensor
        List <Eigen::MatrixXd> ct1Matrix;
        Eigen::Tensor<double, 3 > ct1Tensor;

        /// Total Turbulent tensor
        List <Eigen::MatrixXd> cTotalMatrix;
        Eigen::Tensor<double, 3 > cTotalTensor;

        /// Velocity coefficients for RBF
        Eigen::MatrixXd velRBF;

        /// Initial velocity coefficients for RBF
        Eigen::MatrixXd velRBF0;

        /// Total B Matrix
        Eigen::MatrixXd bTotalMatrix;

        /// The matrix of L2 projection coefficients for the eddy viscosity
        Eigen::MatrixXd coeffL2;

        /// The matrix initial values of L2 projection coefficients for the eddy viscosity
        Eigen::MatrixXd coeff0L2;

        /// The vector of L2 projection coefficients for the eddy viscosity snapshot
        Eigen::VectorXd nutCoeff;

        /// Number of viscoisty modes used for the projection
        label nNutModes;

        /// Eddy viscosity field
        autoPtr<volScalarField> _nut;

        //--------------------------------------------------------------------------
        /// @brief      bt added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The N of velocity modes.s
        /// @param[in]  NSUPmodes  The N of supremizer modes.
        ///
        /// @return     bt matrix for turbulence treatment
        ///
        Eigen::MatrixXd btTurbulence(label NUmodes, label NSUPmodes);

        //--------------------------------------------------------------------------
        /// @brief      ct1 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm1(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      { function_description }
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct1 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor1(label NUmodes,
                label NSUPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      ct2 added matrix for the turbulence treatement
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 matrix for turbulence treatment
        ///
        List < Eigen::MatrixXd > turbulenceTerm2(label NUmodes, label NSUPmodes,
                label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      { function_description }
        ///
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        /// @return     ct2 tensor for the turbulence treatement
        ///
        Eigen::Tensor<double, 3 > turbulenceTensor2(label NUmodes,
                label NSUPmodes, label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      Perform a truthsolve
        /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
        /// to construct mu_interp matrix which is written out in the Offline folder, also for par file in
        /// the Parameters folder.
        ///
        void truthSolve(List<scalar> mu_now);

        //--------------------------------------------------------------------------
        /// Project using a supremizer approach
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        void projectSUP(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes);

        //--------------------------------------------------------------------------
        /// Project using the Poisson Equation for pressure
        ///
        /// @param[in]  folder     The folder name
        /// @param[in]  NUmodes    The number of velocity modes.
        /// @param[in]  NPmodes    The number of pressure modes.
        /// @param[in]  NSUPmodes  The number of supremizer modes which will be set to zero.
        /// @param[in]  nNutModes  The number of eddy viscosity modes.
        ///
        void projectPPE(fileName folder, label NUmodes, label NPmodes, label NSUPmodes,
                        label nNutModes);

        //--------------------------------------------------------------------------
        /// @brief      { function_description }
        ///
        /// @param[in]  A            { parameter_description }
        /// @param[in]  initSnapInd  The initialize snap ind
        ///
        /// @return     { description_of_the_return_value }
        ///
        List < Eigen::MatrixXd > rbfModifiedCoeffs(Eigen::MatrixXd A, Eigen::MatrixXd G,
                Eigen::VectorXd initSnapInd);

        ///
        List < Eigen::MatrixXd > rbfModifiedCoeffsZ(Eigen::MatrixXd A,
                Eigen::MatrixXd G,
                Eigen::VectorXd initSnapInd);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif










