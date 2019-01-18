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

#include "reducedSteadyNSsplit.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
reducedSteadyNSsplit::reducedSteadyNSsplit()
{
}

reducedSteadyNSsplit::reducedSteadyNSsplit(steadyNSsplit& FOMproblem)
:
problem(&FOMproblem)
{

    N_BC = problem->inletIndex.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_p = problem->K_matrix.cols();
    Nphi_nut = problem->CT2_matrix[0].rows();
    //tauU = ITHACAdict->lookupOrDefault<vector>("penalty factor",Vector<double> (1, 1, 0));


    if (problem->NUmodesx != 0)
    {
        for (label k = 0; k < problem->NUmodesx; k++)
        {
            Umodes.append(problem->Umodesx[k]);
        }
    }

    if (problem->NUmodesy != 0)
    {
        for (label k = 0; k < problem->NUmodesy; k++)
        {
            Umodes.append(problem->Umodesy[k]);
        }
    }

    if (problem->NUmodesz != 0)
    {
        for (label k = 0; k < problem->NUmodesz; k++)
        {
            Umodes.append(problem->Umodesz[k]);
        }
    }

    if (problem->NSUPmodes != 0)
    {
        for (label k = 0; k < problem->NSUPmodes; k++)
        {
            Umodes.append(problem->supmodes[k]);
        }
    }
    else
    {

        if (problem->NSUPmodesx != 0)
        {
            for (label k = 0; k < problem->NSUPmodesx; k++)
            {
                Umodes.append(problem->supmodesx[k]);
            }
        }

        if (problem->NSUPmodesy != 0)
        {
            for (label k = 0; k < problem->NSUPmodesy; k++)
            {
                Umodes.append(problem->supmodesy[k]);
            }
        }

        if (problem->NSUPmodesz != 0)
        {
            for (label k = 0; k < problem->NSUPmodesz; k++)
            {
                Umodes.append(problem->supmodesz[k]);
            }
        }
    }

    
    newton_object = newton_steadyNSsplit(Nphi_u + Nphi_p, Nphi_u + Nphi_p,
        FOMproblem);
}

// int newton_steadyNSsplit::operator()(const Eigen::VectorXd& x,
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


int newton_steadyNSsplit::operator()(const Eigen::VectorXd& x,
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
            penaltyU.col(l) = BC(l) * problem->Q_vector[l] - problem->E_matrix[l] *
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

        if (problem->BCmethod == "penalty")
        {
            fvec(i) += (penaltyU * tauU)(i,0);
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
int newton_steadyNSsplit::df(const Eigen::VectorXd& x,
    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_steadyNSsplit> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //


void reducedSteadyNSsplit::solveOnline_sup(Eigen::MatrixXd vel_now)
{

    y.resize(Nphi_u + Nphi_p, 1);
    y.setZero();

    //if (problem->BCmethod == "lift")
    //{
    //    for (label j = 0; j < N_BC; j++)
    //    {
    //        y(j) = vel_now(j, 0);
    //    }
    //}
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    Eigen::HybridNonLinearSolver<newton_steadyNSsplit> hnls(newton_object);
    newton_object.BC.resize(N_BC);
    newton_object.tauU = tauU;

    for (label j = 0; j < N_BC; j++)
    {
        newton_object.BC(j) = vel_now(j, 0);
    }

    for (label i = 0; i < Nphi_nut; i++)
    {
        newton_object.nu_c(i) = problem->rbfsplines[i]->eval(vel_now);
    }

    Ncoeff_DDROM = newton_object.nu_c;


    volScalarField nut_rec("nut_rec", problem->nuTmodes[0] * 0);

    for (label j = 0; j < Nphi_nut; j++)
    {
        nut_rec += problem->nuTmodes[j] * newton_object.nu_c(j);
    }

    nutREC.append(nut_rec);
    newton_object.nu = nu;

    hnls.solve(y);

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

    count_online_solve += 1;
}


void reducedSteadyNSsplit::reconstruct_sup(fileName folder, int printevery)
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
}

void reducedSteadyNSsplit::reconstruct_sup_pLift(volScalarField PLiftfield,fileName folder, int printevery)
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

