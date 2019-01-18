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

\*---------------------------------------------------------------------------*/

#include "steadyNS.H"
#include "steadyNSsplit.H"
#include "viscosityModel.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
steadyNSsplit::steadyNSsplit() {}
steadyNSsplit::steadyNSsplit(int argc, char* argv[])
{
// #include "postProcess.H"
#include "setRootCase.H"
#include "createTime.H"
#include "createMesh.H"
    _simple = autoPtr<simpleControl>
    (
        new simpleControl
        (
          mesh
          )
        );
    simpleControl& simple = _simple();
#include "createFields.H"
#include "createFvOptions.H"
//
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "initContinuityErrs.H"
#pragma GCC diagnostic pop
//
    supex = ITHACAutilities::check_sup();
    turbulence->validate();
    ITHACAdict = new IOdictionary
    (
        IOobject
        (
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
            )
        );
    tolerance = ITHACAdict->lookupOrDefault<scalar>("tolerance", 1e-5);
    maxIter = ITHACAdict->lookupOrDefault<scalar>("maxIter", 1000);
    BCmethod = ITHACAdict->lookupOrDefault<word>("BCMethod", "lift");
    split_directions = ITHACAdict->lookupOrDefault<vector>("split_directions",Vector<double> (1, 1, 0));
    M_Assert(BCmethod == "lift" || BCmethod == "penalty",
        "The BC method must be set to lift or penalty in ITHACAdict");
    para = new ITHACAparameters;
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //

// Method to performa a truthSolve
void steadyNSsplit::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    volScalarField& p = _p();
    volVectorField& U = _U();
    surfaceScalarField& phi = _phi();
    fv::options& fvOptions = _fvOptions();
    simpleControl& simple = _simple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
#include "NLsolve.H"
    exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    volScalarField _nut(turbulence->nut());
//volScalarField nuTilda = mesh.lookupObject<volScalarField>("nuTilda");
    exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
//exportSolution(nuTilda, name(counter), "./ITHACAoutput/Offline/");
    
    if(split_directions[0]==1)
    {
        volScalarField Ux("Ux",U.component(0));
        volVectorField temp
        (
            IOobject
            (
                "temp",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(1,0,0))
            );

        volVectorField tempU
        (
            IOobject
            (
                "Ux",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
            );
        tempU = temp * Ux;
        Ufieldx.append(tempU);
    }
    if(split_directions[1]==1)
    {
        volScalarField Uy("Uy",U.component(1));
        volVectorField temp
        (
            IOobject
            (
                "temp",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(0,1,0))
            );

        volVectorField tempU
        (
            IOobject
            (
                "Uy",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
            );
        tempU = temp * Uy;
        Ufieldy.append(tempU);
    } 
    if(split_directions[2]==1)
    {
        volScalarField Uz("Uz",U.component(2));
        volVectorField temp
        (
            IOobject
            (
                "temp",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(0,0,1))
            );

        volVectorField tempU
        (
            IOobject
            (
                "Uz",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
                ),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero)
            );
        tempU = temp * Uz;
        Ufieldz.append(tempU);
    }

    Ufield.append(U);
    Pfield.append(p);
    nutFields.append(_nut);
    counter++;
    writeMu(mu_now);
// --- Fill in the mu_samples with parameters (mu) to be used for the PODI sample points
    mu_samples.conservativeResize(mu_samples.rows() + 1, mu_now.size());

    for (int i = 0; i < mu_now.size(); i++)
    {
        mu_samples(mu_samples.rows() - 1, i) = mu_now[i];
    }

// Resize to Unitary if not initialized by user (i.e. non-parametric problem)
    if (mu.cols() == 0)
    {
        mu.resize(1, 1);
    }

    if (mu_samples.rows() == mu.cols())
    {
        ITHACAstream::exportMatrix(mu_samples, "mu_samples", "eigen",
           "./ITHACAoutput/Offline");
    }
}

List < Eigen::MatrixXd > steadyNSsplit::turbulence_term1(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes,
    label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes)
{

    label Csize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    List < Eigen::MatrixXd > CT1_matrix;
    CT1_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT1_matrix[j].resize(Nnutmodes, Csize);
        CT1_matrix[j] = CT1_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

// Create PTRLIST with lift, velocities and supremizers

    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT1_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT1_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::laplacian(
                  nuTmodes[j], Together[k])).value();
            }
        }
    }

// Export the matrix
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "eigen",
        "./ITHACAoutput/Matrices/CT1");
    return CT1_matrix;
}







List < Eigen::MatrixXd > steadyNSsplit::turbulence_term2(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes,
    label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes)
{
    label Csize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    List < Eigen::MatrixXd > CT2_matrix;
    CT2_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT2_matrix[j].resize(Nnutmodes, Csize);
        CT2_matrix[j] = CT2_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

// Create PTRLIST with lift, velocities and supremizers

    

    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT2_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT2_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & (fvc::div(
                  nuTmodes[j] * dev((fvc::grad(Together[k]))().T())))).value();
            }
        }
    }

// Export the matrix
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "eigen",
        "./ITHACAoutput/Matrices/CT2");
    return CT2_matrix;
}

Eigen::MatrixXd steadyNSsplit::BT_turbulence(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{

    label BTsize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz ;
    Eigen::MatrixXd BT_matrix(BTsize, BTsize);

    BT_matrix = BT_matrix * 0;
// Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);

    

    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

// Project everything
    for (label i = 0; i < BTsize; i++)
    {
        for (label j = 0; j < BTsize; j++)
        {
            BT_matrix(i, j) = fvc::domainIntegrate(Together[i] & (fvc::div(dev((T(fvc::grad(
                Together[j]))))))).value();
        }
    }

// Export the matrix
    ITHACAstream::SaveDenseMatrix(BT_matrix, "./ITHACAoutput/Matrices/",
        "BT_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz));
    return BT_matrix;
}

// void steadyNSsplit::projectSUP(fileName folder, label NU, label NP, label NSUP,
//   label Nnut)
// {
//     if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
//     {
//         B_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/B_mat.txt");
//         C_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/C", "C");
//         K_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/K_mat.txt");
//         P_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/P_mat.txt");
//         M_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/M_mat.txt");
//         BT_matrix =
//         ITHACAstream::readMatrix("./ITHACAoutput/Matrices/BT_matrix_mat.txt");
//         CT1_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT1",
//           "CT1_matrix");
//         CT2_matrix = ITHACAstream::readMatrix("./ITHACAoutput/Matrices/CT2",
//           "CT2_matrix");
//     }
//     else
//     {
//         NUmodes = NU;
//         NPmodes = NP;
//         NSUPmodes = NSUP;
//         Nnutmodes = Nnut;
//         B_matrix = diffusive_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         C_matrix = convective_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         K_matrix = pressure_gradient_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         P_matrix = divergence_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         M_matrix = mass_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         BT_matrix = BT_turbulence(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         CT1_matrix = turbulence_term1(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes);
//         CT2_matrix = turbulence_term2(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz, label Nnutmodes);
//         if (BCmethod == "penalty")
//         {
//             BC_vel_vector = BC_velocity_vec(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//             BC_vel_matrix = BC_velocity_mat(label NUmodesx, label NUmodesy,label NUmodesz, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz);
//         }
//     }

//     B_total_matrix = B_matrix + BT_matrix;
//     C_total_matrix.setSize(C_matrix.size());

//     for (label i = 0; i < C_matrix.size(); i++)
//     {
//         C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
//     }

//     NUmodes = NU;
//     NPmodes = NP;
//     NSUPmodes = NSUP;
//     Nnutmodes = Nnut;
//     // Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
//     Eigen::MatrixXd Ncoeff = ITHACAutilities::get_coeffs_ortho(nutFields, nuTmodes);
//     ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python",
//      "./ITHACAoutput/Matrices/");
//     SAMPLES.resize(Nnutmodes);
//     rbfsplines.resize(Nnutmodes);

//     for (int i = 0; i < Nnutmodes; i++)
//     {
//         SAMPLES[i] = new SPLINTER::DataTable(1, 1);

//         for (int j = 0; j < Ncoeff.cols(); j++)
//         {
//             SAMPLES[i]->addSample(mu.row(j), Ncoeff(i, j));
//         }

//         rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i],
//             SPLINTER::RadialBasisFunctionType::GAUSSIAN);
//         std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
//     }
// }


void steadyNSsplit::projectSUP(fileName folder, label NUx, label NUy, label NUz, label NP, label NSUP, label NSUPx,
    label NSUPy, label NSUPz, label Nnut)
{
    NUmodesx = NUx;
    NUmodesy = NUy;
    NUmodesz = NUz;
    NPmodes = NP;
    NSUPmodes = NSUP;
    NSUPmodesx = NSUPx;
    NSUPmodesy = NSUPy;
    NSUPmodesz = NSUPz;
    Nnutmodes = Nnut;

    // if NSUPmodes is zero means you have to divide the supremizer modes into components
    if(NSUPmodes ==0)
    {
        for(int i = 0; i<supmodes.size(); i++)
        {
            if(split_directions[0]==1)
            {
                volScalarField Usupx("Usupx",supmodes[i].component(0));
                volVectorField temp
                (
                    IOobject
                    (
                        "temp",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(1,0,0))
                    );

                volVectorField tempU
                (
                    IOobject
                    (
                        "Usupx",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", supmodes[0].dimensions(), vector::zero)
                    );
                tempU = temp * Usupx;
                supmodesx.append(tempU);
            }
            if(split_directions[1]==1)
            {
                volScalarField Usupy("Usupy",supmodes[i].component(1));
                volVectorField temp
                (
                    IOobject
                    (
                        "temp",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(0,1,0))
                    );

                volVectorField tempU
                (
                    IOobject
                    (
                        "Usupy",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", supmodes[0].dimensions(), vector::zero)
                    );
                tempU = temp * Usupy;
                supmodesy.append(tempU);
            }
            if(split_directions[2]==1)
            {
                volScalarField Usupz("Usupz",supmodes[i].component(2));
                volVectorField temp
                (
                    IOobject
                    (
                        "temp",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", dimensionSet(0, 0, 0, 0, 0, 0, 0), vector(0,0,1))
                    );

                volVectorField tempU
                (
                    IOobject
                    (
                        "Usupz",
                        supmodes[0].time().timeName(),
                        supmodes[0].mesh(),
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                        ),
                    supmodes[0].mesh(),
                    dimensionedVector("zero", supmodes[0].dimensions(), vector::zero)
                    );
                tempU = temp * Usupz;
                supmodesz.append(tempU);
            }
        }
    }


    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {

        word B_str = "B_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        }
        else
        {
            B_matrix = diffusive_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }

        word BT_str = "BT_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + BT_str))
        {
            ITHACAstream::ReadDenseMatrix(BT_matrix, "./ITHACAoutput/Matrices/", BT_str);
        }
        else
        {
            BT_matrix = BT_turbulence( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);

        }

        word K_str = "K_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        }
        else
        {
            K_matrix = pressure_gradient_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }

        word P_str = "P_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz) + "_" + name(NPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        }
        else
        {
            P_matrix = divergence_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }

        word M_str = "M_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        }
        else
        {
            M_matrix = mass_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }

        C_matrix = convective_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        CT1_matrix = turbulence_term1( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz,  Nnutmodes);
        CT2_matrix = turbulence_term2( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes, NSUPmodesx,  NSUPmodesy,  NSUPmodesz,  Nnutmodes);

        if (BCmethod == "penalty")
        {
            Q_vector = Boundary_vector_list( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
            E_matrix = Boundary_matrix_list( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }
    }
    else
    {   
        B_matrix = diffusive_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        C_matrix = convective_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        K_matrix = pressure_gradient_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        P_matrix = divergence_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        M_matrix = mass_term( NUmodesx,  NUmodesy, NUmodesz,  NPmodes, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        BT_matrix = BT_turbulence( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        CT1_matrix = turbulence_term1( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz,  Nnutmodes);
        CT2_matrix = turbulence_term2( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz,  Nnutmodes);
        if (BCmethod == "penalty")
        {
            Q_vector = Boundary_vector_list( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
            E_matrix = Boundary_matrix_list( NUmodesx,  NUmodesy, NUmodesz, NSUPmodes,  NSUPmodesx,  NSUPmodesy,  NSUPmodesz);
        }
    }

// Export the matrices
    if (para->exportPython)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "python", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "python", "./ITHACAoutput/Matrices/");
    }

    if (para->exportMatlab)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "matlab", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "matlab", "./ITHACAoutput/Matrices/");
    }

    if (para->exportTxt)
    {
        ITHACAstream::exportMatrix(B_matrix, "B", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(K_matrix, "K", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(P_matrix, "P", "eigen", "./ITHACAoutput/Matrices/");
        ITHACAstream::exportMatrix(M_matrix, "M", "eigen", "./ITHACAoutput/Matrices/");
    }
    B_total_matrix = B_matrix + BT_matrix;
    C_total_matrix.setSize(C_matrix.size());

    for (label i = 0; i < C_matrix.size(); i++)
    {
        C_total_matrix[i] =  CT2_matrix[i] + CT1_matrix[i];
    }


// Get the coeffs for interpolation (the orthonormal one is used because basis are orthogonal)
    Eigen::MatrixXd Ncoeff = ITHACAutilities::get_coeffs_ortho(nutFields, nuTmodes);
    ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "matlab",
        "./ITHACAoutput/Matrices/");

    SAMPLES.resize(Nnutmodes);
    rbfsplines.resize(Nnutmodes);

    for (int i = 0; i < Nnutmodes; i++)
    {
        SAMPLES[i] = new SPLINTER::DataTable(1, 1);

        for (int j = 0; j < Ncoeff.cols(); j++)
        {
            SAMPLES[i]->addSample(mu.row(j), Ncoeff(i, j));
        }

        rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i],
            SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        std::cout << "Constructing RadialBasisFunction for mode " << i + 1 << std::endl;
    }
}

// * * * * * * * * * * * * * * Momentum Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNSsplit::diffusive_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    label Bsize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    Eigen::MatrixXd B_matrix;
    B_matrix.resize(Bsize, Bsize);
    PtrList<volVectorField> Together(0);



    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }



    // Project everything
    for (label i = 0; i < Bsize; i++)
    {
        for (label j = 0; j < Bsize; j++)
        {
            B_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::laplacian(
                dimensionedScalar("1", dimless, 1), Together[j])).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/",
      "B_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz));
    return B_matrix;
}

Eigen::MatrixXd steadyNSsplit::pressure_gradient_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    label K1size = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    label K2size = NPmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);



    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K_matrix(i, j) = fvc::domainIntegrate(Together[i] & fvc::grad(
                Pmodes[j])).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/",
      "K_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz) + "_" + name(NPmodes));
    return K_matrix;
}

List < Eigen::MatrixXd > steadyNSsplit::convective_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    label Csize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    List < Eigen::MatrixXd > C_matrix;
    C_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        C_matrix[j].resize(Csize, Csize);
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers


    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        for (label j = 0; j < Csize; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                C_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::div(
                    linearInterpolate(Together[j]) & Together[j].mesh().Sf(), Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(C_matrix, "C", "python", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(C_matrix, "C", "matlab", "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(C_matrix, "C", "eigen", "./ITHACAoutput/Matrices/C");
    return C_matrix;
}

Eigen::MatrixXd steadyNSsplit::mass_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    label Msize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    Eigen::MatrixXd M_matrix(Msize, Msize);
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);



    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    // Project everything
    for (label i = 0; i < Msize; i++)
    {
        for (label j = 0; j < Msize; j++)
        {
            M_matrix(i, j) = fvc::domainIntegrate(Together[i] & Together[j]).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/",
      "M_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz));
    return M_matrix;
}

// * * * * * * * * * * * * * * Continuity Eq. Methods * * * * * * * * * * * * * //

Eigen::MatrixXd steadyNSsplit::divergence_term(label NUmodesx, label NUmodesy,label NUmodesz, label NPmodes, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    label P1size = NPmodes;
    label P2size = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    Eigen::MatrixXd P_matrix(P1size, P2size);
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);



    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
            P_matrix(i, j) = fvc::domainIntegrate(Pmodes[i] * fvc::div (
                Together[j])).value();
        }
    }

    ITHACAstream::SaveDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/",
      "P_" + name(NUmodesx) + "_" + name(NUmodesy) + "_" + name(NUmodesz) + "_" + name(NSUPmodesx) + "_" + name(NSUPmodesy) + "_" + name(NSUPmodesz) + "_" + name(NPmodes));
    return P_matrix;
}

List< Eigen::MatrixXd > steadyNSsplit::Boundary_vector_list(label NUmodesx, label NUmodesy, label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    PtrList<volVectorField> Together(0);

    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }


    label Qsize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    List < Eigen::MatrixXd > Q_vector(inletIndex.rows());
    //List < Eigen::MatrixXd > BC_vel_vector(patches_penalty.size());

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        Q_vector[j].resize(Qsize, 1);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k,0);
        label BCcomp = inletIndex(k,1);

        for (label i = 0; i < Qsize; i++)
        {
            Q_vector[k](i, 0) = gSum(Together[i].boundaryField()[BCind]).component(
                BCcomp);   
        }

    }

    ITHACAstream::exportMatrix(Q_vector, "Q_vector", "eigen",
     "./ITHACAoutput/Matrices/Q_vector");
    return Q_vector;
}

List< Eigen::MatrixXd > steadyNSsplit::Boundary_matrix_list(label NUmodesx, label NUmodesy, label NUmodesz, label NSUPmodes, label NSUPmodesx, label NSUPmodesy, label NSUPmodesz)
{
    PtrList<volVectorField> Together(0);

    if (NUmodesx != 0)
    {
        for (label k = 0; k < NUmodesx; k++)
        {
            Together.append(Umodesx[k]);
        }
    }

    if (NUmodesy != 0)
    {
        for (label k = 0; k < NUmodesy; k++)
        {
            Together.append(Umodesy[k]);
        }
    }

    if (NUmodesz != 0)
    {
        for (label k = 0; k < NUmodesz; k++)
        {
            Together.append(Umodesz[k]);
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k]);
        }
    }
    else
    {
        if (NSUPmodesx != 0)
        {
            for (label k = 0; k < NSUPmodesx; k++)
            {
                Together.append(supmodesx[k]);
            }
        }

        if (NSUPmodesy != 0)
        {
            for (label k = 0; k < NSUPmodesy; k++)
            {
                Together.append(supmodesy[k]);
            }
        }

        if (NSUPmodesz != 0)
        {
            for (label k = 0; k < NSUPmodesz; k++)
            {
                Together.append(supmodesz[k]);
            }
        }
    }

    label Esize = NUmodesx + NUmodesy + NUmodesz + NSUPmodes + NSUPmodesx + NSUPmodesy + NSUPmodesz;
    label EUsize = inletIndex.rows();
    List < Eigen::MatrixXd > E_matrix(EUsize);

    for (label j = 0; j < inletIndex.rows(); j++)
    {
        E_matrix[j].resize(Esize, Esize);
    }

    for (label k = 0; k < inletIndex.rows(); k++)
    {
        label BCind = inletIndex(k,0);
        label BCcomp = inletIndex(k,1);

        for (label i = 0; i < Esize; i++)
        {
            for (label j = 0; j < Esize; j++)
            {
                E_matrix[k](i, j) = gSum(Together[i].boundaryField()[BCind].component(BCcomp) *
                  Together[j].boundaryField()[BCind].component(BCcomp));
            }
        }
    }

    ITHACAstream::exportMatrix(E_matrix, "E_matrix", "eigen",
     "./ITHACAoutput/Matrices/E_matrix");
    return E_matrix;
}