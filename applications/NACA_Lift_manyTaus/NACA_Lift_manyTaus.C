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

// #include "setRootCase.H"
// #include "createTime.H"
// #include "createMesh.H"

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
            int status0 = system(app);

            word bc_online = "cp ./ITHACAoutput/Reconstruction/bc_online_mat.m ./ITHACAoutput/Reconstruction/bc_online_mat"
            +name(k+1) + ".m";
            int status1 = system(bc_online);

            int status2 = system("pwd");
            char s[100]; 

            printf("%s\n", getcwd(s, 100)); 
            chdir("ITHACAoutput/Reconstruction");
            printf("%s\n", getcwd(s, 100)); 
            //int status2 = system("cd ITHACAoutput/Reconstruction");
            int status3 = system("lift_and_drag");
            std::string s1 = std::to_string(N_online);
            char const *pchar1 = s1.c_str();
            char const *pchar2 = "postProcessing/FC/";

            char result[100];   // array to hold the result.

            strcpy(result,pchar2); // copy string one into the result.
            strcat(result,pchar1); // append string two to the result.[]
            chdir(result);

            word lift_command2 = "cp forces.dat forces_ddrom" + name(k+1) + ".dat";
            int status4 = system(lift_command2);
            int status5 = system("rm forces.dat");

            chdir("../../../../..");
        }
        return 0;
    }

