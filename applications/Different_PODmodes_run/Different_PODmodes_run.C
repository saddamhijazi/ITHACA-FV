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
        word toRun("./modes"); // the offline samples.
        Eigen::MatrixXd modesToRun = ITHACAstream::readMatrix(toRun);
        int N_run = modesToRun.rows();

        for (label k = 0; k < N_run; k++)
        {
        	word Ucommand = "sed -i 's/^NmodesU .*$/NmodesU " + name(modesToRun(k,0)) + ";/' ./system/ITHACAdict";
        	word Pcommand = "sed -i 's/^NmodesP .*$/NmodesP " + name(modesToRun(k,1)) + ";/' ./system/ITHACAdict";
        	word SUPcommand = "sed -i 's/^NmodesSUP .*$/NmodesSUP " + name(modesToRun(k,2)) + ";/' ./system/ITHACAdict";
        	word NUTcommand = "sed -i 's/^NmodesNUT .*$/NmodesNUT " + name(modesToRun(k,3)) + ";/' ./system/ITHACAdict";
        	int status0 = system(Ucommand);
        	int status1 = system(Pcommand);
        	int status2 = system(SUPcommand);
        	int status3 = system(NUTcommand);
            word tau_command = "cp tau_penalty_mat"+name(modesToRun(k,0))+".txt tau_penalty_mat.txt";
            int status33 = system(tau_command);

        	int status4 = system(app);
            int status5 = system("compute_error");
            word U_command = "cp ./ITHACAoutput/postProcessing/errorUSUP_mat.m ./ITHACAoutput/postProcessing/errorUSUP_mat"
            +name(modesToRun(k,0)) + ".m";
            word p_command = "cp ./ITHACAoutput/postProcessing/errorPSUP_mat.m ./ITHACAoutput/postProcessing/errorPSUP_mat"
            +name(modesToRun(k,0)) + ".m";
            word bc_online = "cp ./ITHACAoutput/Reconstruction/bc_online_mat.m ./ITHACAoutput/Reconstruction/bc_online_mat"
            +name(modesToRun(k,0)) + ".m";
            int status6 = system(U_command);
            int status7 = system(p_command);
            int status8 = system(bc_online);

            // word folderCommand1 = "cp -r ./ITHACAoutput/LiftandDragMatrices/ ./ITHACAoutput/LiftandDragMatrices"
            // + name(modesToRun(k,0)) + "_" + name(modesToRun(k,1)) + "_" + name(modesToRun(k,2)) + "_" 
            // + name(modesToRun(k,3));
            // word folderCommand2 = "cp -r ./ITHACAoutput/LiftandDragMatricesLam/ ./ITHACAoutput/LiftandDragMatricesLam"
            // + name(modesToRun(k,0)) + "_" + name(modesToRun(k,1)) + "_" + name(modesToRun(k,2)) + "_" 
            // + name(modesToRun(k,3));
            // int status5 = system(folderCommand1);
            // int status6 = system(folderCommand2);
        }
        return 0;
    }

