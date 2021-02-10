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
#include <Eigen/SparseLU>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

// #include "setRootCase.H"
// #include "createTime.H"
// #include "createMesh.H"

        word toRun("./modes"); // the offline samples.
        Eigen::MatrixXd modesToRun = ITHACAstream::readMatrix(toRun);
        int N_run = modesToRun.rows();
        word modesNames;
        for (label k = 0; k < N_run; k++)
        {
            if(modesToRun.cols()==1)
            {
                word Ucommand = "sed -i 's/^NmodesU .*$/NmodesU " + name(modesToRun(k,0)) + ";/' ./system/ITHACAdict";
                modesNames = name(modesToRun(k,0));
                int status0 = system(Ucommand);
            }
            else if(modesToRun.cols()==2)
            {
                word Ucommand = "sed -i 's/^NmodesU .*$/NmodesU " + name(modesToRun(k,0)) + ";/' ./system/ITHACAdict";
                word Pcommand = "sed -i 's/^NmodesP .*$/NmodesP " + name(modesToRun(k,1)) + ";/' ./system/ITHACAdict";
                modesNames = name(modesToRun(k,0)) + "_" + name(modesToRun(k,1));
                int status0 = system(Ucommand);
                int status1 = system(Pcommand);
            }
            else if(modesToRun.cols()==3)
            {
                word Ucommand = "sed -i 's/^NmodesU .*$/NmodesU " + name(modesToRun(k,0)) + ";/' ./system/ITHACAdict";
                word Pcommand = "sed -i 's/^NmodesP .*$/NmodesP " + name(modesToRun(k,1)) + ";/' ./system/ITHACAdict";
                word NUTcommand = "sed -i 's/^NmodesNUT .*$/NmodesNUT " + name(modesToRun(k,2)) + ";/' ./system/ITHACAdict";
                modesNames = name(modesToRun(k,0)) + "_" + name(modesToRun(k,1)) + "_" + name(modesToRun(k,2));
                int status0 = system(Ucommand);
                int status1 = system(Pcommand);
                int status3 = system(NUTcommand);
            }
            else if(modesToRun.cols()==4)
            {
              word Ucommand = "sed -i 's/^NmodesU .*$/NmodesU " + name(modesToRun(k,0)) + ";/' ./system/ITHACAdict";
              word Pcommand = "sed -i 's/^NmodesP .*$/NmodesP " + name(modesToRun(k,1)) + ";/' ./system/ITHACAdict";
              word SUPcommand = "sed -i 's/^NmodesSUP .*$/NmodesSUP " + name(modesToRun(k,2)) + ";/' ./system/ITHACAdict";
              word NUTcommand = "sed -i 's/^NmodesNUT .*$/NmodesNUT " + name(modesToRun(k,3)) + ";/' ./system/ITHACAdict";
              modesNames = name(modesToRun(k,0)) + "_" + name(modesToRun(k,1)) + "_" + name(modesToRun(k,2)) + "_" 
              + name(modesToRun(k,3));
              int status0 = system(Ucommand);
              int status1 = system(Pcommand);
              int status2 = system(SUPcommand);
              int status3 = system(NUTcommand);
          }
            
          
          word run_command = "mpirun -np 30 BoxLES -parallel > LES_Box.log";
          int status4 = system(run_command);
          


          word folders_command = "rm -r  ./ITHACAoutput/POD ./ITHACAoutput/supremizer";
          word log_command = "rm LES_Box.log";

          int status5 = system(folders_command);
          int status6 = system(log_command);


      }
      return 0;
  }
