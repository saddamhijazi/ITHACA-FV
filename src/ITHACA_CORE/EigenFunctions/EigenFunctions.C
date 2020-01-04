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

#include "EigenFunctions.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void EigenFunctions::sortEigenvalues(Eigen::VectorXd& eigenvalues,
                                     Eigen::MatrixXd& eigenvectors)
{
    labelList order;
    scalarField eigenValues(eigenvalues.size());

    for (int i = 0; i < eigenvalues.size(); i++)
    {
        eigenValues[i] = eigenvalues(i);
    }

    sortedOrder(eigenValues, order);
    scalarField eigenValues2(eigenValues);

    for (int i = 0; i < order.size(); i++)
    {
        eigenvalues(i) = eigenValues[order[order.size() - i - 1]];
    }

    Eigen::MatrixXd eigenvectors2 = eigenvectors;

    for (label i = 0; i < eigenvalues.size(); i++)
    {
        for (label k = 0; k < eigenvalues.size(); k++)
        {
            eigenvectors2(i, k) = eigenvectors(k, order[order.size() - i - 1]);
        }
    }

    eigenvectors = eigenvectors2;
}

template<typename T> Eigen::Matrix < T, -1,
         -1 > EigenFunctions::rowWiseNormalization(Eigen::Matrix < T, -1,
         -1 > matrix, Eigen::Matrix < T, -1, -1 > & rowWiseMax
         , Eigen::Matrix < T, -1, -1 > & rowWiseMin)
{
    Eigen::Matrix < T, -1, -1 > normalized;
    normalized = matrix;
    rowWiseMax.resize(matrix.rows(), 1);
    rowWiseMin.resize(matrix.rows(), 1);
    Eigen::Matrix < T, -1, -1 > temp;
    temp.setOnes(1, matrix.cols());

    for (int i = 0; i < matrix.rows(); i++)
    {
        rowWiseMin(i, 0) = matrix.row(i).minCoeff();
        rowWiseMax(i, 0) = matrix.row(i).maxCoeff();
        normalized.row(i) = normalized.row(i) - rowWiseMin(i,
                            0) * temp;
        normalized.row(i) /= rowWiseMax(i, 0) - rowWiseMin(i, 0);
    }

    return normalized;
}

template<typename T> Eigen::Matrix < T, -1,
         -1 > EigenFunctions::columnWiseNormalization(Eigen::Matrix < T, -1,
         -1 > matrix, Eigen::Matrix < T, -1, -1 > & columnWiseMax
         , Eigen::Matrix < T, -1, -1 > & columnWiseMin)
{
    Eigen::Matrix < T, -1, -1 > normalized;
    normalized = matrix;
    columnWiseMax.resize(matrix.cols(), 1);
    columnWiseMin.resize(matrix.cols(), 1);
    Eigen::Matrix < T, -1, -1 > temp;
    temp.setOnes(matrix.rows(), 1);

    for (int i = 0; i < matrix.cols(); i++)
    {
        columnWiseMin(i, 0) = matrix.col(i).minCoeff();
        columnWiseMax(i, 0) = matrix.col(i).maxCoeff();
        normalized.col(i) = normalized.col(i) - columnWiseMin(i, 0) * temp;
        normalized.col(i) /= columnWiseMax(i, 0) - columnWiseMin(i, 0);
    }

    return normalized;
}

template<typename T> Eigen::Vector < T,
         -1 > EigenFunctions::columnWiseNormalization(Eigen::Vector < T, -1 > vector,
         Eigen::Matrix < T, -1, -1 > columnWiseMax
         , Eigen::Matrix < T, -1, -1 > columnWiseMin)
{
    Eigen::Vector < T, -1 > normalized;
    normalized = vector;

    for (int i = 0; i < vector.size(); i++)
    {
        normalized(i) = normalized(i) - columnWiseMin(i, 0);
        normalized(i) /= columnWiseMax(i, 0) - columnWiseMin(i, 0);
    }

    return normalized;
}

template<typename T> Eigen::Matrix < T, -1,
         -1 > EigenFunctions::rowWiseRecovery(Eigen::Matrix < T, -1,
         -1 > normalized, Eigen::Matrix < T, -1, -1 > rowWiseMax
         , Eigen::Matrix < T, -1, -1 > rowWiseMin)
{
    Eigen::Matrix < T, -1, -1 > recovery;
    recovery = normalized;
    Eigen::Matrix < T, -1, -1 > temp;
    temp.setOnes(1, normalized.cols());

    for (int i = 0; i < normalized.rows(); i++)
    {
        recovery.row(i) *= rowWiseMax(i, 0) - rowWiseMin(i, 0);
        recovery.row(i) += rowWiseMin(i, 0) * temp;
    }

    return recovery;
}

template<typename T> Eigen::Matrix < T, -1,
         -1 > EigenFunctions::columnWiseRecovery(Eigen::Matrix < T, -1,
         -1 > normalized, Eigen::Matrix < T, -1, -1 > columnWiseMax
         , Eigen::Matrix < T, -1, -1 > columnWiseMin)
{
    Eigen::Matrix < T, -1, -1 > recovery;
    recovery = normalized;
    Eigen::Matrix < T, -1, -1 > temp;
    temp.setOnes(normalized.rows(), 1);

    for (int i = 0; i < normalized.cols(); i++)
    {
        recovery.col(i) *= columnWiseMax(i, 0) - columnWiseMin(i, 0);
        recovery.col(i) += columnWiseMin(i, 0) * temp;
    }

    return recovery;
}

template Eigen::Matrix < double, -1,
                         -1 > EigenFunctions::rowWiseNormalization(Eigen::Matrix < double, -1,
-1 > matrix, Eigen::Matrix < double, -1,
    -1 > & rowWiseMax
    , Eigen::Matrix < double, -1,
    -1 > & rowWiseMin);

template Eigen::Matrix < int, -1,
                         -1 > EigenFunctions::rowWiseNormalization(Eigen::Matrix < int, -1,
-1 > matrix, Eigen::Matrix < int, -1,
    -1 > & rowWiseMax
    , Eigen::Matrix < int, -1,
    -1 > & rowWiseMin);

template Eigen::Matrix < float, -1,
                         -1 > EigenFunctions::rowWiseNormalization(Eigen::Matrix < float, -1,
-1 > matrix, Eigen::Matrix < float, -1,
    -1 > & rowWiseMax
    , Eigen::Matrix < float, -1,
    -1 > & rowWiseMin);

template Eigen::Matrix < double, -1,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Matrix < double, -1,
-1 > matrix, Eigen::Matrix < double, -1, -1 > & columnWiseMax
    , Eigen::Matrix < double, -1, -1 > & columnWiseMin);

template Eigen::Matrix < int, -1,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Matrix < int, -1,
-1 > matrix, Eigen::Matrix < int, -1, -1 > & columnWiseMax
    , Eigen::Matrix < int, -1, -1 > & columnWiseMin);

template Eigen::Matrix < float, -1,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Matrix < float, -1,
-1 > matrix, Eigen::Matrix < float, -1, -1 > & columnWiseMax
    , Eigen::Matrix < float, -1, -1 > & columnWiseMin);

template Eigen::Matrix < double, -1,
                         -1 > EigenFunctions::rowWiseRecovery(Eigen::Matrix < double, -1,
-1 > normalized, Eigen::Matrix < double, -1,
    -1 > rowWiseMax
    , Eigen::Matrix < double, -1,
    -1 > rowWiseMin);

template Eigen::Matrix < int, -1,
                         -1 > EigenFunctions::rowWiseRecovery(Eigen::Matrix < int, -1,
-1 > normalized, Eigen::Matrix < int, -1,
    -1 > rowWiseMax
    , Eigen::Matrix < int, -1,
    -1 > rowWiseMin);

template Eigen::Matrix < float, -1,
                         -1 > EigenFunctions::rowWiseRecovery(Eigen::Matrix < float, -1,
-1 > normalized, Eigen::Matrix < float, -1,
    -1 > rowWiseMax
    , Eigen::Matrix < float, -1,
    -1 > rowWiseMin);

template Eigen::Matrix < double, -1,
                         -1 > EigenFunctions::columnWiseRecovery(Eigen::Matrix < double, -1,
-1 > normalized, Eigen::Matrix < double, -1, -1 > columnWiseMax
    , Eigen::Matrix < double, -1, -1 > columnWiseMin);

template Eigen::Matrix < int, -1,
                         -1 > EigenFunctions::columnWiseRecovery(Eigen::Matrix < int, -1,
-1 > normalized, Eigen::Matrix < int, -1, -1 > columnWiseMax
    , Eigen::Matrix < int, -1, -1 > columnWiseMin);

template Eigen::Matrix < float, -1,
                         -1 > EigenFunctions::columnWiseRecovery(Eigen::Matrix < float, -1,
-1 > normalized, Eigen::Matrix < float, -1, -1 > columnWiseMax
    , Eigen::Matrix < float, -1, -1 > columnWiseMin);

template Eigen::Vector < double,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Vector < double,
-1 > vector, Eigen::Matrix < double, -1, -1 > columnWiseMax
    , Eigen::Matrix < double, -1, -1 > columnWiseMin);

template Eigen::Vector < int,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Vector < int,
-1 > vector, Eigen::Matrix < int, -1, -1 > columnWiseMax
    , Eigen::Matrix < int, -1, -1 > columnWiseMin);

template Eigen::Vector < float,
                         -1 > EigenFunctions::columnWiseNormalization(Eigen::Vector < float,
-1 > vector, Eigen::Matrix < float, -1, -1 > columnWiseMax
    , Eigen::Matrix < float, -1, -1 > columnWiseMin);