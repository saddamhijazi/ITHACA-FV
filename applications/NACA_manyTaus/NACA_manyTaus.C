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
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    lift_and_drag

Author
    Giovanni Stabile, SISSA MathLab (International School for Advanced Studies) gstabile@sissa.it

Description
    Application to recover the lift and the drag after the simulation is performed

\*---------------------------------------------------------------------------*/

/// \file
/// \brief Application to recover the lift and the drag after the simulation is performed
/// \details In order to use this file one needs to prepare a FORCESdict file, in order to 
/// check the syntax one needs to check the \ref FORCESdict file.

/// \file FORCESdict
/// \brief Example of a FORCESdict file

#include "fvCFD.H"
#include "IOmanip.H"
#include "IFstream.H"
#include "primitiveFields.H"
#include "FieldFields.H"
#include "scalarMatrices.H"
#include "SortableList.H"
#include "volFieldsFwd.H"
#include "forces.H"
#include "forceCoeffs.H"
#include "volFields.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <stdio.h>
#include "ITHACAstream.H"
#include "ITHACAutilities.H"

#include <Eigen/Dense>

#include "Foam2Eigen.H"
#include <chrono>

#include "EigenFunctions.H"
#include <chrono>
#include <Eigen/SVD>
#include <GenEigsSolver.h>
#include <Eigen/SparseLU>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    std::string app = argv[1];
        word tau_var("./TAUS"); // the offline samples.
        Eigen::MatrixXd tau_setting = ITHACAstream::readMatrix(tau_var);
        int N_BC = tau_setting(0);
        int N_online = tau_setting(1);
        int N_tau = tau_setting.rows() - 2;
        Eigen::MatrixXd taus = tau_setting.bottomRows(N_tau);

        for (label k = 0; k < N_tau; k++)
        {
            Eigen::MatrixXd tau_now;
            tau_now.setOnes(N_online, N_BC);
            tau_now = taus(k,0) * tau_now;
            ITHACAstream::exportMatrix(tau_now, "tau_penalty", "eigen", "./");
            int status1 = system(app);
            int status2 = system("compute_error");
            word U_command = "cp ./ITHACAoutput/postProcessing/errorUSUP_mat.m ./ITHACAoutput/postProcessing/errorUSUP_mat"
            +name(k+1) + ".m";
            word p_command = "cp ./ITHACAoutput/postProcessing/errorPSUP_mat.m ./ITHACAoutput/postProcessing/errorPSUP_mat"
            +name(k+1) + ".m";
            word bc_online = "cp ./ITHACAoutput/Reconstruction/bc_online_mat.m ./ITHACAoutput/Reconstruction/bc_online_mat"
            +name(k+1) + ".m";
            int status3 = system(U_command);
            int status4 = system(p_command);
            int status5 = system(bc_online);




        }
        return 0;
    }

