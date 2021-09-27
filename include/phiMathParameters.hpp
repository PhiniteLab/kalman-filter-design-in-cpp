/* Copyright (C) 2021 Phinite Lab. All rights reserved.
 Contact information
 ------------------------------------------------------
 Developers: 
 Mehmet İşcan    : mehmetmkt@gmail.com
 Ali Burak Özden : ozdenaliburak@gmail.com 
 Company: 
 Phinite Engineering, Istanbul
 Web      :  http://phinitelab.com
 e-mail   :  info@phinitelab.com
 
 */
#ifndef __PHI_MATH_PARAMETERS_HPP__
#define __PHI_MATH_PARAMETERS_HPP__

// in order to generalize our code, we need to use
#ifdef __cplusplus
extern "C"
{
#endif

    ////////////////////////////////////////////////////////////
    // including some libraries for using input/output functions

#include "stdio.h"
#include "string.h"
#include "stdlib.h"

    // including some libraries for using input/output functions
    ////////////////////////////////////////////////////////////

#ifdef __cplusplus
}
#endif

////////////////////////////////////////////////////////////
// define error types

#define ALLOCATION_ERROR 1
#define INCONSISTENT_ROW_COLUMN 2
#define FILE_OPEN_ERROR 3
#define SAMPLING_RATE_ERROR 4

// define error types
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiMathParameters
{
private:
protected:
public:
    ////////////////////////////////////////////////////////////
    // general variables

    // general variables
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    phiMathParameters()
    {
    }

    ~phiMathParameters()
    {
    }

    // CONSTRUCTORS
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // function prototype decleration

    double **phiMathCreatingEmptyMatrices(int rows, int cols)
    {

        /*
         * creating dynamically empty matrices for general usage
         * output -> address of empty matrices
         * input  -> rows and columns values of matrices
         * */

        double **pd = (double **)malloc(rows * sizeof(double *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }
        for (int i = 0; i < rows; i++)
            pd[i] = (double *)malloc(cols * sizeof(double));

        return pd;
    }

    //////////////////////////////////////////////////////////////////////
    // error Handler

    void phiMathPhiErrorHandler(int errorType)
    {
        /*
         * notifying the error result
         * output -> return nothing
         * input  -> type of error
         * */

        switch (errorType)
        {
        case FILE_OPEN_ERROR:
            printf("System Dynamic Parameter files cannot be created!\n");
            break;
        case INCONSISTENT_ROW_COLUMN:
            printf("The rows and columns are not consistent!\n");
            break;
        case ALLOCATION_ERROR:
            printf("Memory allocation cannot be done!\n");
            break;
        case SAMPLING_RATE_ERROR:
            printf("Sampling period cannot be assigned to either negative or zero value!\n");
            break;

        default:
            break;
        }
    }

    // error Handler
    //////////////////////////////////////////////////////////////////////

    double **phiMathPhiVectorMatrixMultiplication(double **firstTerm,
                                                  double **SecondTerm,
                                                  int row1,
                                                  int col1,
                                                  int row2,
                                                  int col2)
    {
        /*
         * creating dynamically matrices produced by multiplication
         * output -> address of multiplied matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of first Matrices
         *           col of first Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        double **pd = this->phiMathCreatingEmptyMatrices(row1, col2);
        double sum = 0;

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row1; i++)
        {
            for (int j = 0; j < col2; j++)
            {
                pd[i][j] = 0;
            }
        }

        for (int i = 0; i < row1; i++) //row of first matrix
        {
            for (int j = 0; j < col2; j++) //column of second matrix
            {
                sum = 0;
                for (int k = 0; k < col1; k++)
                {
                    sum = sum + firstTerm[i][k] * SecondTerm[k][j];
                }
                pd[i][j] = sum;
            }
        }
        return pd;
    }

    double **phiMathPhiSkalarMatrixMultiplication(double skalarTerm,
                                                  double **SecondTerm,
                                                  int row,
                                                  int col)
    {
        /*
         * creating dynamically matrices produced by skalar multiplication
         * output -> address of multiplied matrices
         * input  -> value of skalar term
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        double **pd = this->phiMathCreatingEmptyMatrices(row, col);

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                pd[i][j] = skalarTerm * SecondTerm[i][j];
            }
        }

        return pd;
    }

    void phiMathPhiFree(double **pd,
                        int row,
                        int col)
    {
        /*
         * free the whole memory allocation
         * output -> return nothing
         * input  -> address of memory
         *           row value of matrices
         *           column value of matrices
         */

        for (int i = 0; i < row; i++)
            free(pd[i]);

        free(pd);
    }

    void phiMathPhiTranspose(double **phiMatrices, int row, int col, double **phiTransMatrices)
    {
        /*
         * creating inverse matrices of phiMatrices
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           value of row
         *           value of col
         *           address of transpose matrices
         */

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                phiTransMatrices[j][i] = phiMatrices[i][j];
            }
        }
    }

    double **phiMathPhiMatrixSummation(double **firstTerm,
                                       double **SecondTerm,
                                       int row,
                                       int col)
    {
        /*
         * creating dynamically matrices produced by summation
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        double **pd = this->phiMathCreatingEmptyMatrices(row, col);

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                pd[i][j] = firstTerm[i][j] + SecondTerm[i][j];
            }
        }

        return pd;
    }

    double **phiMathPhiMatrixNegation(double **firstTerm,
                                      double **SecondTerm,
                                      int row,
                                      int col)
    {
        /*
         * creating dynamically matrices produced by summation
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        double **pd = this->phiMathCreatingEmptyMatrices(row, col);

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                pd[i][j] = firstTerm[i][j] - SecondTerm[i][j];
            }
        }

        return pd;
    }

    void phiMathPhiMatrixAssignment(double **assignedTerm,
                                    double **SecondTerm,
                                    int row,
                                    int col)
    {
        /*
         * creating dynamically matrices produced by summation
         * output -> address of summed matrices
         * input  -> address of first Matrices
         *           address of second Matrices
         *           row of second Matrices
         *           col of second Matrices
         * */

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                assignedTerm[i][j] = SecondTerm[i][j];
            }
        }
    }

    double **phiMathEyeMatricesCreation(int row)
    {
        /*
         * creating dynamically eye matrices to produces discrete to continuos from conversion
         * output -> address of eye matrices
         * input  -> address of systemDynamicsParameter
         */

        double **pd = (double **)malloc(row * sizeof(double *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < row; i++)
            pd[i] = (double *)malloc(row * sizeof(double));

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < row; j++)
            {
                if (i == j)
                {
                    pd[i][j] = 1;
                }
                else
                {
                    pd[i][j] = 0;
                }
            }
        }

        return pd;
    }

    ////////////////////////////////////////////////////////////////////
    // inverse operation
    // Function to get cofactor of A[p][q] in temp[][]. n is current
    // dimension of A[][]
    void phiMathGetCofactor(double **A, double **temp, int p, int q, int n)
    {
        int i = 0, j = 0;

        // Looping for each element of the matrix
        for (int row = 0; row < n; row++)
        {
            for (int col = 0; col < n; col++)
            {
                // Copying into temporary matrix only those element
                // which are not in given row and column
                if (row != p && col != q)
                {
                    temp[i][j++] = A[row][col];

                    // Row is filled, so increase row index and
                    // reset col index
                    if (j == n - 1)
                    {
                        j = 0;
                        i++;
                    }
                }
            }
        }
    }

    /* Recursive function for finding determinant of matrix.
    n is current dimension of A[][]. */
    double phiMathDeterminant(double **A, int n)
    {
        double D = 0; // Initialize result

        // Base case : if matrix contains single element
        if (n == 1)
            return A[0][0];

        double **temp = this->phiMathCreatingEmptyMatrices(n, n); // To store cofactors

        double sign = 1; // To store sign multiplier

        // Iterate for each element of first row
        for (int f = 0; f < n; f++)
        {
            // Getting Cofactor of A[0][f]
            this->phiMathGetCofactor(A, temp, 0, f, n);
            D += sign * A[0][f] * this->phiMathDeterminant(temp, n - 1);

            // terms are to be added with alternate sign
            sign = -sign;
        }

        this->phiMathPhiFree(temp, n, n);

        return D;
    }

    // Function to get adjoint of A[N][N] in adj[N][N].
    void phiMathAdjoint(double **A, double **adj, int row)
    {
        if (row == 1)
        {
            adj[0][0] = 1;
            return;
        }

        // temp is used to store cofactors of A[][]
        double sign = 1;
        double **temp = this->phiMathCreatingEmptyMatrices(row, row);

        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < row; j++)
            {
                // Get cofactor of A[i][j]
                this->phiMathGetCofactor(A, temp, i, j, row);

                // sign of adj[j][i] positive if sum of row
                // and column indexes is even.
                sign = ((i + j) % 2 == 0) ? 1 : -1;

                // Interchanging rows and columns to get the
                // transpose of the cofactor matrix
                adj[j][i] = (sign) * (this->phiMathDeterminant(temp, row - 1));
            }
        }

        this->phiMathPhiFree(temp, row, row);
    }

    // Function to calculate and store inverse, returns false if
    // matrix is singular
    bool phiMathInverse(double **A, double **inverse, int row)
    {
        // Find determinant of A[][]
        double det = this->phiMathDeterminant(A, row);
        if (det == 0)
        {
            printf("Singular matrix, can't find its inverse\n");
            return false;
        }

        // Find adjoint
        double **adj = this->phiMathCreatingEmptyMatrices(row, row);
        this->phiMathAdjoint(A, adj, row);

        // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
        for (int i = 0; i < row; i++)
            for (int j = 0; j < row; j++)
                inverse[i][j] = adj[i][j] / double(det);

        this->phiMathPhiFree(adj, row, row);

        return true;
    }

    // inverse operation
    ////////////////////////////////////////////////////////////////////

    void phiMathDisplayMatrices(double **dp, int row, int col)
    {
        printf("\n\n");
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                printf("%lf ", dp[i][j]);
            }
            printf("\n");
        }
    }

    double phiMathPhiRand(double min, double max)
    {
        /*
         * generating random number in double number for a given min/max
         * output -> value of double number
         * input  -> value of min value
         *           value of max value
         */

        double scale = rand() / (double)RAND_MAX; /* [0, 1.0] */
        return min + scale * (max - min);         /* [min, max] */
    }

    // function prototype decleration
    ////////////////////////////////////////////////////////////
};

typedef phiMathParameters *phiMathParametersPtr;

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#endif