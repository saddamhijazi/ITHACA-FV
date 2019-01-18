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

#include "reducedSteadyNSturb.H"
#include "EigenFunctions.H"
#include "Foam2Eigen.H"


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSteadyNSturb::reducedSteadyNSturb()
{
}

reducedSteadyNSturb::reducedSteadyNSturb(steadyNSturb& FOMproblem)
:
problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
    Nphi_nut = problem->CT2_matrix[0].rows();

    for (label k = 0; k < problem->liftfield.size(); k++)
    {
        Umodes.append(problem->liftfield[k]);
    }

    for (label k = 0; k < problem->NUmodes; k++)
    {
        Umodes.append(problem->Umodes[k]);
    }

    for (label k = 0; k < problem->NSUPmodes; k++)
    {
        Umodes.append(problem->supmodes[k]);
    }

    newton_object = newton_steadyNSturb(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
        FOMproblem);
}

// int newton_steadyNSturb::operator()(const Eigen::VectorXd& x,
//     Eigen::VectorXd& fvec) const
// {
//     Eigen::VectorXd a_tmp(Nphi_u);
//     Eigen::VectorXd b_tmp(Nphi_p);
//     a_tmp = x.head(Nphi_u);
//     b_tmp = x.tail(Nphi_p);
//     // Convective term
//     Eigen::MatrixXd cc(1, 1);
//     // Mom Term
//     Eigen::VectorXd M1 = problem->B_total_matrix * a_tmp * nu;
//     // Gradient of pressure
//     Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
//     // Pressure Term
//     Eigen::VectorXd M3 = problem->P_matrix * a_tmp;

//     for (label i = 0; i < Nphi_u; i++)
//     {
//         cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp - nu_c.transpose() *
//         problem->C_total_matrix[i] * a_tmp;
//         // Info << "a_tmp.transpose() * problem->C_matrix[i] * a_tmp is " << a_tmp.transpose() * problem->C_matrix[i] * a_tmp << endl;
//         // Info << "nu_c.transpose() * problem->C_total_matrix[i] * a_tmp is " << nu_c.transpose() * problem->C_total_matrix[i] * a_tmp 
//         // << endl;

//         fvec(i) = M1(i) - cc(0, 0) - M2(i);
//     }
//     for (label j = 0; j < Nphi_p; j++)
//     {
//         label k = j + Nphi_u;
//         fvec(k) = M3(j);
//     }

//     for (label j = 0; j < N_BC; j++)
//     {
//         fvec(j) = x(j) - BC(j);
//     }

//     return 0;
// }


int newton_steadyNSturb::operator()(const Eigen::VectorXd& x,
    Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_p);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.tail(Nphi_p);
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Mom Term
    Eigen::VectorXd M1 = problem->B_total_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Pressure Term
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Penalty term
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    //Eigen::MatrixXd Temp(1, 1);
    // Term for penalty method
    if (problem->BCmethod == "penalty")
    {
        for (label l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = BC(l) * problem->BC_vel_vector[l] - problem->BC_vel_matrix[l] *
            a_tmp;
        }
    }
    for (label i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * problem->C_matrix[i] * a_tmp - nu_c.transpose() *
        problem->C_total_matrix[i] * a_tmp;
        // Info << "a_tmp.transpose() * problem->C_matrix[i] * a_tmp is " << a_tmp.transpose() * problem->C_matrix[i] * a_tmp << endl;
        // Info << "nu_c.transpose() * problem->C_total_matrix[i] * a_tmp is " << nu_c.transpose() * problem->C_total_matrix[i] * a_tmp 
        // << endl;
        fvec(i) = M1(i) - cc(0, 0) - M2(i);
        // std::cout << "fvec is " << fvec << std::endl;
        // std::cout << "cc(0, 0) is " << cc(0, 0) << std::endl;
        // std::cout << "M1 is " << M1 << std::endl;
        // std::cout << "M2 is " << M2 << std::endl;
        // std::cout << "problem->C_matrix[i] is " << problem->C_matrix[i] << std::endl;
        // std::cout << "problem->C_total_matrix[i] is " << problem->C_total_matrix[i] << std::endl;
        // std::cout << "nu_c.transpose() is " << nu_c.transpose() << std::endl;
        // exit(0);
        //Temp = penaltyU.row(i) * tauU ;
        //std::cout << "penaltyU rows is " << penaltyU.rows() << std::endl;
        //std::cout << "penaltyU cols is " << penaltyU.cols() << std::endl;
        //std::cout << "tauU size is " << tauU.size() << std::endl;
        //std::cout << "tauU is " << tauU << std::endl;
        //std::cout << "penaltyU.row(i) is " << penaltyU.row(i) << std::endl;
        //std::cout << "Temp is " << Temp << std::endl;
        //scalar Temp1 = Temp(0,0);
        //std::cout << typeid(Temp).name() <<  std::endl;
        //penaltyU = penaltyU.cwiseAbs();
        //fvec = fvec.cwiseAbs();

        if (problem->BCmethod == "penalty")
        {
            fvec(i) += ((penaltyU * tauU)(i,0));
        }

    }
    for (label j = 0; j < Nphi_p; j++)
    {
        label k = j + Nphi_u;
        fvec(k) = M3(j);
    }
    if (problem->BCmethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }
    }

    return 0;
}
int newton_steadyNSturb::df(const Eigen::VectorXd& x,
    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_steadyNSturb> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

void reducedSteadyNSturb::fixed_point_matrix_and_vector(Eigen::VectorXd& x, Eigen::VectorXd& BC, Eigen::VectorXd& nu_rbf, 
    Eigen::VectorXd& b, Eigen::MatrixXd& Aa)
{
    Aa.resize(Nphi_u+Nphi_p,Nphi_u+Nphi_p);
    b.resize(Nphi_u+Nphi_p,1);
    Eigen::MatrixXd cc(Nphi_u, Nphi_u);
    

    for (label i = 0; i < Nphi_u; i++)
    {
        cc.row(i) = x.transpose().head(Nphi_u) * problem->C_matrix[i] - nu_rbf.transpose() * problem->C_total_matrix[i];
    }
    Eigen::MatrixXd M1 = cc + problem->B_total_matrix * nu;
    Eigen::MatrixXd penaltyU1 = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    Eigen::MatrixXd penaltyU2 = Eigen::MatrixXd::Zero(Nphi_u, Nphi_u);

    // Term for penalty method
    if (problem->BCmethod == "penalty")
    {
        for (label l = 0; l < N_BC; l++)
        {

            penaltyU1.col(l) = BC(l) * problem->BC_vel_vector[l];
            penaltyU2 += - tauU(l,0) * problem->BC_vel_matrix[l];
        }
    }

    Aa.topLeftCorner(Nphi_u, Nphi_u) = M1 + penaltyU2;
    Aa.topRightCorner(Nphi_u,Nphi_p) = problem->K_matrix;
    Aa.bottomLeftCorner(Nphi_p, Nphi_u) = problem->P_matrix;
    Aa.bottomRightCorner(Nphi_p, Nphi_p) = Eigen::MatrixXd::Zero(Nphi_p, Nphi_p);
    b.head(Nphi_u) = - penaltyU1 * tauU;
    b.tail(Nphi_p) = Eigen::MatrixXd::Zero(Nphi_p, 1);
    //Info << "End of the method" << endl;


}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //


void reducedSteadyNSturb::solveOnline_sup(Eigen::MatrixXd vel_now)
{

    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    if (problem->BCmethod == "lift")
    {
        for (label j = 0; j < N_BC; j++)
        {
            y(j) = vel_now(j, 0);
        }
    }

    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newton_steadyNSturb> hnls(newton_object);
    newton_object.BC.resize(N_BC);
    newton_object.tauU = tauU;
    for (label j = 0; j < N_BC; j++)
    {
        newton_object.BC(j) = vel_now(j, 0);
    }


    if(problem->viscCoeff=="L2")
    {
        for (label i = 0; i < Nphi_nut; i++)
        {
            newton_object.nu_c = problem->nutCoeff;
        }
    }
    else if(problem->viscCoeff=="RBF")
    {
        for (label i = 0; i < Nphi_nut; i++)
        {
            newton_object.nu_c(i) = problem->rbfsplines[i]->eval(vel_now);
        }
    }
    else
    {
        Info << "The way to compute the eddy viscosity coefficients has to be either L2 or RBF" << endl;
        exit(0);
    }
    
    // for (label i = 0; i < Nphi_nut; i++)
    // {
    //     newton_object.nu_c(i) = problem->rbfsplines[i]->eval(vel_now);
    // }

    Ncoeff_DDROM = newton_object.nu_c;


    volScalarField nut_rec("nut_rec", problem->nuTmodes[0] * 0);

    for (label j = 0; j < Nphi_nut; j++)
    {
        nut_rec += problem->nuTmodes[j] * newton_object.nu_c(j);
    }

    nutREC.append(nut_rec);
    newton_object.nu = nu;


    Eigen::VectorXd nu_rbf = newton_object.nu_c;
    Eigen::VectorXd y1;
    Eigen::VectorXd res_fixedPoint;
    Eigen::VectorXd b;
    Eigen::VectorXd BC;
    Eigen::MatrixXd Aa;
    BC = newton_object.BC;


    //hnls.solve(y);
    int ite_max = 1000;
    int j = 0;

    while(j < ite_max)
    {
        fixed_point_matrix_and_vector(y,BC,nu_rbf,b,Aa);
        y1 = Aa.fullPivLu().solve(b);
        res_fixedPoint = y1 - y;
        if(res_fixedPoint.norm() < 1e-10)
        {
            break;
        }
        y = y1;
        j = j+1;
    }
    y = y1;
    // std::cout << "j is " << j << std::endl;
    // std::cout << "y is " << y << std::endl;
    // std::cout << "y1 is " << y1 << std::endl;
    // std::cout << "res_fixedPoint is " << res_fixedPoint << std::endl;
    // std::cout << "res_fixedPoint norm is " << res_fixedPoint.norm() << std::endl;

    hnls.solve(y);
    //std::cout << "y is " << y << std::endl;

    // exit(0);

    Eigen::VectorXd res(y);
    newton_object.operator()(y, res);
    std::cout << "################## Online solve N° " << count_online_solve <<
    " ##################" << std::endl;
    std::cout << "Solving for the parameter: " << vel_now << std::endl;

    if (res.norm() < 1e-5)
    {
        std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
        hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    else
    {
        std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " <<
        hnls.iter << " iterations " << def << std::endl << std::endl;
    }
    //exit(0);
    count_online_solve += 1;
}


void reducedSteadyNSturb::reconstruct_sup(fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;
    Eigen::MatrixXd bc_online(online_solution.size(),3);
    label BCind = problem->inletIndex(0,0);

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }
            // const fvPatchVectorField& BC = U_rec.boundaryField()[BCind];
            // List<Eigen::VectorXd> bcvalue = Foam2Eigen::field2EigenBC(U_rec);
            
            bc_online(counter,0) = U_rec.boundaryField()[BCind][0][0];
            bc_online(counter,1) = U_rec.boundaryField()[BCind][0][1];
            bc_online(counter,2) = U_rec.boundaryField()[BCind][0][2];
            

            problem->exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
            volScalarField P_rec("P_rec", problem->Pmodes[0] * 0);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            problem->exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
            nextwrite += printevery;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }
    ITHACAstream::exportMatrix(bc_online, "bc_online", "matlab", folder);
}

void reducedSteadyNSturb::reconstruct_sup_pLift(volScalarField PLiftfield,fileName folder, int printevery)
{
    mkDir(folder);
    ITHACAutilities::createSymLink(folder);
    int counter = 0;
    int nextwrite = 0;

    for (label i = 0; i < online_solution.size(); i++)
    {
        if (counter == nextwrite)
        {
            volVectorField U_rec("U_rec", Umodes[0] * 0);

            for (label j = 0; j < Nphi_u; j++)
            {
                U_rec += Umodes[j] * online_solution[i](j + 1, 0);
            }

            problem->exportSolution(U_rec, name(online_solution[i](0, 0)), folder);
            volScalarField P_rec("P_rec", PLiftfield);

            for (label j = 0; j < Nphi_p; j++)
            {
                P_rec += problem->Pmodes[j] * online_solution[i](j + Nphi_u + 1, 0);
            }

            problem->exportSolution(P_rec, name(online_solution[i](0, 0)), folder);
            nextwrite += printevery;
            UREC.append(U_rec);
            PREC.append(P_rec);
        }

        counter++;
    }
}
// ************************************************************************* //

