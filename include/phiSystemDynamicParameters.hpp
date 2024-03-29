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

#ifndef __PHI_SYSTEM_DYNAMIC_PARAMETERS_HPP__
#define __PHI_SYSTEM_DYNAMIC_PARAMETERS_HPP__

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
// including some libraries for using input/output functions
#include "phiMathParameters.hpp"

// including some libraries for using input/output functions
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiSystemDynamicParameters : public phiMathParameters
{
private:
protected:
public:
    ////////////////////////////////////////////////////////////
    // general variables

    /*
    * This structure represents the basic variables in terms of
    * system dynamics parameters.
    * 
    * There are only pointer types of object, the rest of them is
    * created in the solver blocks!
    * 
    */

    double **phiSystemDynamicStateMatrices; /*
                            * storing the state matrices of equation of motion
                            */

    double **phiSystemDynamicInputMatrices; /*
                            * storing the input matrices of equation of motion
                            */

    double **phiSystemDynamicDiscreteStateMatrices; /*
                            * storing discrete state matrices of equation of motion
                            */

    double **phiSystemDynamicDiscreteInputMatrices; /*
                            * storing discrete input matrices of equation of motion
                            */

    int phiSystemDynamicRows; /*
               * storing row value of state space
               */

    int phiSystemDynamicCols; /*
               * storing column value of state space
               * */

    double phiSystemDynamicInputValue; /*
                       * storing instanteous value of input
                       */

    int phiSystemDynamicInputNumber; /*
                      * storing the number of input value
                    */

    double **phiSystemDynamicStateNow; /*
                       * storing the state vectors in the next value
                       */

    double **phiSystemDynamicStatePre; /*
                       * storing the state vectors in the now value
                       */

    double **phiSystemDynamicInputNow; /*
                       * storing the input vectors in the now value
                       */

    double phiSystemDynamicDt; /*
               * storing the sampling period of the state space equation
               */

    FILE *phiSystemDynamicFp; /*
               * storing the text file name
               */

    /////////////////////////////////////////
    // time parameters

    double phiSystemDynamicTStartTime;

    double phiSystemDynamicTFinalTime;

    int phiSystemDynamicNumberOfIteration;

    // time parameters
    /////////////////////////////////////////

    // general variables
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    phiSystemDynamicParameters()
    {
    }

    phiSystemDynamicParameters(const char *fileName)
    {
        this->phiSystemDynamicFp = fopen(fileName, "w");

        if (this->phiSystemDynamicFp == NULL)
        {
            this->phiMathPhiErrorHandler(FILE_OPEN_ERROR);
            exit(EXIT_FAILURE);
        }
    }

    ~phiSystemDynamicParameters()
    {
        fclose(this->phiSystemDynamicFp);
    }

    // CONSTRUCTORS
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // function prototype decleration

    double **phiSystemDynamicCreatingEmptyStateMatrices()
    {
        /*
         * creating dynamically empty state matrices for state space equation 
         * output -> address of matrices (double **)
         * input  -> taking nothing
         */

        double **pd = (double **)malloc(this->phiSystemDynamicRows * sizeof(double *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        if (this->phiSystemDynamicRows != this->phiSystemDynamicCols)
        {
            this->phiMathPhiErrorHandler(INCONSISTENT_ROW_COLUMN);
            exit(EXIT_FAILURE);
        }

        if (this->phiSystemDynamicDt <= 0)
        {
            this->phiMathPhiErrorHandler(SAMPLING_RATE_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
            pd[i] = (double *)malloc(this->phiSystemDynamicCols * sizeof(double));

        return pd;
    }

    double **phiSystemDynamicCreatingEmptyInputMatrices()
    {
        /*
         * creating dynamically empty input matrices for state space equation
         * output -> address of matrices (double **)
         * input  -> taking nothing
         */

        double **pd = (double **)malloc(this->phiSystemDynamicRows * sizeof(double *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
            pd[i] = (double *)malloc(this->phiSystemDynamicCols * sizeof(double));

        return pd;
    }

    void phiSystemDynamicWriteTheMatrices()
    {
        /*
         * printing state and input matrices
         * output -> return nothing
         * input  -> taking nothing
         */

        printf("Printing state matrices...\n");
        for (int i = 0; i < this->phiSystemDynamicRows; i++)
        {
            for (int j = 0; j < this->phiSystemDynamicCols; j++)
            {
                printf("%f ", this->phiSystemDynamicStateMatrices[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        printf("Printing input matrices...\n");
        for (int i = 0; i < this->phiSystemDynamicRows; i++)
        {
            for (int j = 0; j < this->phiSystemDynamicInputNumber; j++)
            {
                printf("%f ", this->phiSystemDynamicInputMatrices[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    double **phiSystemDynamicEyeMatricesCreation()
    {
        /*
         * creating dynamically eye matrices to produces discrete to continuos from conversion
         * output -> address of eye matrices
         * input  -> address of systemDynamicsParameter
         */

        double **pd = (double **)malloc(this->phiSystemDynamicRows * sizeof(double *));

        if (pd == NULL)
        {
            this->phiMathPhiErrorHandler(ALLOCATION_ERROR);
            exit(EXIT_FAILURE);
        }

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
            pd[i] = (double *)malloc(this->phiSystemDynamicCols * sizeof(double));

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
        {
            for (int j = 0; j < this->phiSystemDynamicCols; j++)
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

    void phiSystemDynamicInitializeSystemParameters(double startTime, double finalTime)
    {
        ////////////////////////////////////////////////
        // creating dynamic matrices

        this->phiSystemDynamicStateNow = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicStatePre = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicInputNow = this->phiMathCreatingEmptyMatrices(1, 1);

        this->phiSystemDynamicDiscreteStateMatrices = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows,
                                                                                         this->phiSystemDynamicCols);

        this->phiSystemDynamicDiscreteInputMatrices = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows,
                                                                                         this->phiSystemDynamicInputNumber);

        // creating dynamic matrices
        ////////////////////////////////////////////////
        this->phiSystemDynamicTStartTime = startTime;
        this->phiSystemDynamicTFinalTime = finalTime;

        // calculate the length of model
        this->phiSystemDynamicNumberOfIteration = static_cast<int>((this->phiSystemDynamicTFinalTime - this->phiSystemDynamicTStartTime) / this->phiSystemDynamicDt) + 1;
    }

    void phiSystemDynamicNiteDynamicSolver(int iterator, double MNstd, double PNstd, double noiseRange)
    {
        /*
         * creating a solver for dynamic patterns!
         * output -> return nothing
         * input  -> iterator index
         *           measurement noise variance
         *           process noise variance
         *           noise range for phiRand function!
         * 
         * Important Note: This code creates a text file to keep the data in matrices form.
         *                 The format parameter of this file is
         *                 "%f %f %f"  -> x1 state, x2 state, elapsedTime 
         */

        ////////////////////////////////////////////////
        // internal terms
        double **DtA;
        double **eyeDtA;
        double **DtB;
        double **stateMultiplicationMatrices;
        double **inputMultiplicationMatrices;
        double **eye = this->phiSystemDynamicEyeMatricesCreation();

        // internal terms
        ////////////////////////////////////////////////

        ////////////////////////////////////////////////
        // Solving discrete form of state space

        // the whole operations
        this->phiSystemDynamicInputNow[0][0] = this->phiSystemDynamicInputValue;

        DtA = this->phiMathPhiSkalarMatrixMultiplication(this->phiSystemDynamicDt,
                                                         this->phiSystemDynamicStateMatrices,
                                                         this->phiSystemDynamicRows,
                                                         this->phiSystemDynamicCols);

        DtB = this->phiMathPhiSkalarMatrixMultiplication(this->phiSystemDynamicDt,
                                                         this->phiSystemDynamicInputMatrices,
                                                         this->phiSystemDynamicRows,
                                                         this->phiSystemDynamicInputNumber);
        eyeDtA = this->phiMathPhiMatrixSummation(eye, DtA,
                                                 this->phiSystemDynamicRows,
                                                 this->phiSystemDynamicCols);

        // updated kalman discrete dynamic
        this->phiMathPhiMatrixAssignment(this->phiSystemDynamicDiscreteStateMatrices,
                                         eyeDtA,
                                         this->phiSystemDynamicRows,
                                         this->phiSystemDynamicCols);

        this->phiMathPhiMatrixAssignment(this->phiSystemDynamicDiscreteInputMatrices,
                                         DtB,
                                         this->phiSystemDynamicRows,
                                         this->phiSystemDynamicInputNumber);

        stateMultiplicationMatrices = this->phiMathPhiVectorMatrixMultiplication(eyeDtA, this->phiSystemDynamicStatePre,
                                                                                 this->phiSystemDynamicRows, this->phiSystemDynamicCols, this->phiSystemDynamicRows, 1);
        inputMultiplicationMatrices = this->phiMathPhiVectorMatrixMultiplication(DtB, this->phiSystemDynamicInputNow,
                                                                                 this->phiSystemDynamicRows, 1, 1, 1);

        double currentTime = iterator * this->phiSystemDynamicDt;

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
        {
            this->phiSystemDynamicStateNow[i][0] = stateMultiplicationMatrices[i][0] + inputMultiplicationMatrices[i][0] + MNstd * this->phiMathPhiRand(-noiseRange, noiseRange);
        }

        ///////////////////////////////////////////////////////////
        // printing results

        printf("State values : ");

        for (int i = 0; i < this->phiSystemDynamicRows; i++)
        {
            printf("x[%d] : %f ", i, this->phiSystemDynamicStatePre[i][0]);
            fprintf(this->phiSystemDynamicFp, "%f ", this->phiSystemDynamicStatePre[i][0]);
        }
        printf("elapsed Time : %f seconds\n", currentTime);
        fprintf(this->phiSystemDynamicFp, "%f \n", currentTime);

        // printing results
        ///////////////////////////////////////////////////////////

        this->phiMathPhiMatrixAssignment(this->phiSystemDynamicStatePre,
                                         this->phiSystemDynamicStateNow,
                                         this->phiSystemDynamicRows, 1);

        // Solving discrete form of state space
        ////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // free the whole memory
        this->phiSystemDynamicPhiFree(DtA, this->phiSystemDynamicRows, this->phiSystemDynamicCols);
        this->phiSystemDynamicPhiFree(eyeDtA, this->phiSystemDynamicRows, this->phiSystemDynamicCols);
        this->phiSystemDynamicPhiFree(DtB, this->phiSystemDynamicRows, this->phiSystemDynamicInputNumber);
        this->phiSystemDynamicPhiFree(stateMultiplicationMatrices, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(inputMultiplicationMatrices, this->phiSystemDynamicRows, 1);
        // free the whole memory
        ///////////////////////////////////////////////////////
    }

    void phiSystemDynamicNiteStaticSolver(double finalTime, const char *fileName)
    {
        /*
         * creating a solver for static patterns!
         * output -> return nothing
         * input  -> finalTime in seconds
         *           file name to be stored data
         * 
         * Important Note: This code creates a text file to keep the data in matrices form.
         *                 The format parameter of this file is
         *                 "%f %f %f"  -> x1 state, x2 state, elapsedTime 
         */

        ////////////////////////////////////////////////
        // internal terms
        double **DtA;
        double **eyeDtA;
        double **DtB;
        double **stateMultiplicationMatrices;
        double **inputMultiplicationMatrices;
        double **eye = this->phiSystemDynamicEyeMatricesCreation();
        int numberOfLength = (int)(finalTime / this->phiSystemDynamicDt);

        // internal terms
        ////////////////////////////////////////////////

        ////////////////////////////////////////////////
        // file creation to write txt file
        FILE *fp;

        fp = fopen(fileName, "w");

        if (fp == NULL)
        {
            this->phiMathPhiErrorHandler(FILE_OPEN_ERROR);
            exit(EXIT_FAILURE);
        }

        // file creation to write txt file
        ////////////////////////////////////////////////

        ////////////////////////////////////////////////
        // creating dynamic matrices
        this->phiSystemDynamicStateNow = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicStatePre = this->phiMathCreatingEmptyMatrices(this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicInputNow = this->phiMathCreatingEmptyMatrices(1, 1);

        // creating dynamic matrices
        ////////////////////////////////////////////////

        ////////////////////////////////////////////////
        // Solving discrete form of state space

        // the whole operations
        this->phiSystemDynamicInputNow[0][0] = this->phiSystemDynamicInputValue;

        DtA = this->phiMathPhiSkalarMatrixMultiplication(this->phiSystemDynamicDt, this->phiSystemDynamicStateMatrices, this->phiSystemDynamicRows, this->phiSystemDynamicCols);
        DtB = this->phiMathPhiSkalarMatrixMultiplication(this->phiSystemDynamicDt, this->phiSystemDynamicInputMatrices, this->phiSystemDynamicRows, this->phiSystemDynamicInputNumber);
        eyeDtA = this->phiMathPhiMatrixSummation(eye, DtA, this->phiSystemDynamicRows, this->phiSystemDynamicCols);

        for (int iter = 0; iter < numberOfLength; iter++)
        {
            stateMultiplicationMatrices = this->phiMathPhiVectorMatrixMultiplication(eyeDtA, this->phiSystemDynamicStatePre,
                                                                                     this->phiSystemDynamicRows, this->phiSystemDynamicCols, this->phiSystemDynamicRows, 1);
            inputMultiplicationMatrices = this->phiMathPhiVectorMatrixMultiplication(DtB, this->phiSystemDynamicInputNow,
                                                                                     this->phiSystemDynamicRows, 1, 1, 1);

            double currentTime = iter * this->phiSystemDynamicDt;
            for (int i = 0; i < this->phiSystemDynamicRows; i++)
            {
                this->phiSystemDynamicStateNow[i][0] = stateMultiplicationMatrices[i][0] + inputMultiplicationMatrices[i][0];
            }

            ///////////////////////////////////////////////////////////
            // printing results

            printf("State values : ");

            for (int i = 0; i < this->phiSystemDynamicRows; i++)
            {
                printf("x[%d] : %f ", i, this->phiSystemDynamicStatePre[i][0]);
                fprintf(fp, "%f ", this->phiSystemDynamicStatePre[i][0]);
            }
            printf("elapsed Time : %f seconds\n", currentTime);
            fprintf(fp, "%f \n", currentTime);

            // printing results
            ///////////////////////////////////////////////////////////

            this->phiMathPhiMatrixAssignment(this->phiSystemDynamicStatePre, this->phiSystemDynamicStateNow, this->phiSystemDynamicRows, 1);
        }

        // Solving discrete form of state space
        ////////////////////////////////////////////////

        ///////////////////////////////////////////////////////
        // free the whole memory
        this->phiSystemDynamicPhiFree(DtA, this->phiSystemDynamicRows, this->phiSystemDynamicCols);
        this->phiSystemDynamicPhiFree(eyeDtA, this->phiSystemDynamicRows, this->phiSystemDynamicCols);
        this->phiSystemDynamicPhiFree(DtB, this->phiSystemDynamicRows, this->phiSystemDynamicInputNumber);
        this->phiSystemDynamicPhiFree(stateMultiplicationMatrices, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(inputMultiplicationMatrices, this->phiSystemDynamicRows, 1);

        this->phiSystemDynamicPhiFree(this->phiSystemDynamicStateNow, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(this->phiSystemDynamicStatePre, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(this->phiSystemDynamicInputNow, this->phiSystemDynamicInputNumber, 1);

        fclose(fp);

        // free the whole memory
        ///////////////////////////////////////////////////////
    }

    void phiSystemDynamicPhiFree(double **pd,
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

    void phiSystemDynamicPhiExit()
    {
        /*
         * free the systemDynamics parameter
         * output -> return nothing
         * input  -> taking nothing
         */

        this->phiSystemDynamicPhiFree(this->phiSystemDynamicStateNow, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(this->phiSystemDynamicStatePre, this->phiSystemDynamicRows, 1);
        this->phiSystemDynamicPhiFree(this->phiSystemDynamicInputNow, this->phiSystemDynamicInputNumber, 1);

        this->phiSystemDynamicPhiFree(this->phiSystemDynamicDiscreteStateMatrices,
                                      this->phiSystemDynamicRows,
                                      this->phiSystemDynamicCols);

        this->phiSystemDynamicPhiFree(this->phiSystemDynamicDiscreteInputMatrices,
                                      this->phiSystemDynamicRows,
                                      this->phiSystemDynamicInputNumber);

        free(this->phiSystemDynamicStateMatrices);
        free(this->phiSystemDynamicInputMatrices);
    }

    // function prototype decleration
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // demos

    void phiSystemDynamicEx1Demo()
    {
        /*
         * First order system example
         */

        this->phiSystemDynamicRows = 1;        // number of row in state space
        this->phiSystemDynamicCols = 1;        // number of col in state space
        this->phiSystemDynamicInputNumber = 1; // number of input in state space
        this->phiSystemDynamicDt = 0.001;      // sampling period

        this->phiSystemDynamicStateMatrices = this->phiSystemDynamicCreatingEmptyStateMatrices();
        // state matrices creation
        this->phiSystemDynamicInputMatrices = this->phiSystemDynamicCreatingEmptyInputMatrices();
        // input matrices creation

        // matrices value assignment
        this->phiSystemDynamicStateMatrices[0][0] = -1;
        this->phiSystemDynamicInputMatrices[0][0] = 1;

        printf("Writing state and input matrices \n\n");

        // write the state space matrices
        this->phiSystemDynamicWriteTheMatrices();

        printf("\n\n");

        /// solver solution
        // select the inputValue to be given to the system
        this->phiSystemDynamicInputValue = 10.0;

        // static time parameter solver
        this->phiSystemDynamicNiteStaticSolver(10, "ver1MCK.txt");

        // should be called to free memory
        this->phiSystemDynamicPhiExit();
    }

    void phiSystemDynamicEx2Demo()
    {
        /*
         * Second order system example
         */

        this->phiSystemDynamicRows = 2;        // number of row in state space
        this->phiSystemDynamicCols = 2;        // number of col in state space
        this->phiSystemDynamicInputNumber = 1; // number of input in state space
        this->phiSystemDynamicDt = 0.001;      // sampling period

        this->phiSystemDynamicStateMatrices = this->phiSystemDynamicCreatingEmptyStateMatrices();
        // state matrices creation
        this->phiSystemDynamicInputMatrices = this->phiSystemDynamicCreatingEmptyInputMatrices();
        // input matrices creation

        // matrices value assignment
        this->phiSystemDynamicStateMatrices[0][0] = 0;
        this->phiSystemDynamicStateMatrices[0][1] = 1;
        this->phiSystemDynamicStateMatrices[1][0] = -0.1;
        this->phiSystemDynamicStateMatrices[1][1] = -1;

        this->phiSystemDynamicInputMatrices[0][0] = 0;
        this->phiSystemDynamicInputMatrices[1][0] = 1;

        printf("Writing state and input matrices \n\n");

        // write the state space matrices
        this->phiSystemDynamicWriteTheMatrices();

        printf("\n\n");

        /// solver solution
        // select the inputValue to be given to the system
        this->phiSystemDynamicInputValue = 10.0;

        // static time parameter solver
        this->phiSystemDynamicNiteStaticSolver(10, "ver1MCK.txt");

        // should be called to free memory
        this->phiSystemDynamicPhiExit();
    }

    void phiSystemDynamicEx3Demo()
    {
        /*
         * nth order system example
         */

        this->phiSystemDynamicRows = 3;        // number of row in state space
        this->phiSystemDynamicCols = 3;        // number of col in state space
        this->phiSystemDynamicInputNumber = 1; // number of input in state space
        this->phiSystemDynamicDt = 0.001;      // sampling period

        this->phiSystemDynamicStateMatrices = this->phiSystemDynamicCreatingEmptyStateMatrices();
        // state matrices creation
        this->phiSystemDynamicInputMatrices = this->phiSystemDynamicCreatingEmptyInputMatrices();
        // input matrices creation

        // matrices value assignment
        this->phiSystemDynamicStateMatrices[0][0] = 0;
        this->phiSystemDynamicStateMatrices[0][1] = 1;
        this->phiSystemDynamicStateMatrices[0][2] = 0;

        this->phiSystemDynamicStateMatrices[1][0] = 0;
        this->phiSystemDynamicStateMatrices[1][1] = 0;
        this->phiSystemDynamicStateMatrices[1][2] = 1;

        this->phiSystemDynamicStateMatrices[2][0] = -2;
        this->phiSystemDynamicStateMatrices[2][1] = -3;
        this->phiSystemDynamicStateMatrices[2][2] = -4;

        this->phiSystemDynamicInputMatrices[0][0] = 0;
        this->phiSystemDynamicInputMatrices[1][0] = 0;
        this->phiSystemDynamicInputMatrices[2][0] = 1;

        printf("Writing state and input matrices \n\n");

        // write the state space matrices
        this->phiSystemDynamicWriteTheMatrices();

        printf("\n\n");

        /// solver solution
        // select the inputValue to be given to the system
        this->phiSystemDynamicInputValue = 2.0;

        // static time parameter solver
        this->phiSystemDynamicNiteStaticSolver(10, "ver1MCK.txt");

        // should be called to free memory
        this->phiSystemDynamicPhiExit();
    }

    // demos
    ////////////////////////////////////////////////////////////
};

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// type decleration for usage of system dynamic solver

#endif