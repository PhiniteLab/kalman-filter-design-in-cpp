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

#ifndef __PHI_KALMAN_DYNAMIC_PARAMETERS_HPP__
#define __PHI_KALMAN_DYNAMIC_PARAMETERS_HPP__

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
#include "phiSystemDynamicParameters.hpp"

// including some libraries for using input/output functions

////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
// CLASS DEFINITION

class phiKalmanDynamicParameters : public phiMathParameters
{
private:
protected:
public:
    ////////////////////////////////////////////////////////////
    // GENERAL VARIABLES

    double **phiKalmanDynamicStateMatrices; /*
                            * storing the state matrices of equation of motion
                            */

    double **phiKalmanDynamicInputMatrices; /*
                            * storing the input matrices of equation of motion
                            */

    int phiKalmanDynamicRows; /*
               * storing row value of state space
               */

    int phiKalmanDynamicCols; /*
               * storing column value of state space
               * */

    double phiKalmanDynamicInputValue; /*
                       * storing instanteous value of input
                       */

    int phiKalmanDynamicInputNumber; /*
                      * storing the number of input value
                    */

    double phiKalmanDynamicDt; /*
               * storing the sampling period of the state space equation
               */

    FILE *phiKalmanDynamicFp; /*
               * storing the text file name
               */

    /////////////////////////////////////////
    // time parameters

    double phiKalmanDynamicTStartTime;

    double phiKalmanDynamicTFinalTime;

    int phiKalmanDynamicNumberOfIteration;

    // time parameters
    /////////////////////////////////////////

    typedef struct kalman_filter
    {
        double **Amatrices; /*
                            * storing the state matrices of equation of motion
                            */

        double **Bmatrices; /*
                            * storing the input matrices of equation of motion
                            */

        int rows; /*
               * storing row value of state space
               */

        int cols; /*
               * storing column value of state space
               * */

        double inputValue; /*
                       * storing instanteous value of input
                       */

        int inputNumber; /*
                      * storing the number of input value
                    */

        double **stateNow; /*
                       * storing the state vectors in the next value
                       */

        double **statePre; /*
                       * storing the state vectors in the now value
                       */

        double **inputNow; /*
                       * storing the input vectors in the now value
                       */

        // kalman algorithm parameters

        double **Q;

        double **H;

        double **R;

        double **K;

        double **x;

        double **P;

        double detP;

        double **z;

    } phiKalmanDynamicKalmanFilter;

    phiKalmanDynamicKalmanFilter kalmanFilterParameters;

    // GENERAL VARIABLES
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // CONSTRUCTORS

    phiKalmanDynamicParameters(const char *fileName)
    {
        this->phiKalmanDynamicFp = fopen(fileName, "w");

        if (this->phiKalmanDynamicFp == NULL)
        {
            this->phiMathPhiErrorHandler(FILE_OPEN_ERROR);
            exit(EXIT_FAILURE);
        }
    }

    ~phiKalmanDynamicParameters()
    {
        fclose(this->phiKalmanDynamicFp);
    }

    // CONSTRUCTORS
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // FUNCTION PROTOTYPE DECLERATION

    void phiKalmanDynamicFirstXCalculation(phiSystemDynamicParameters *ptrSys)
    {
        /////////////////////
        // x calculation

        double **internalCalcA = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.Amatrices,
                                                                            this->kalmanFilterParameters.x,
                                                                            this->kalmanFilterParameters.rows,
                                                                            this->kalmanFilterParameters.cols,
                                                                            this->kalmanFilterParameters.rows, 1);

        double **internalCalcB = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.Bmatrices,
                                                                            this->kalmanFilterParameters.inputNow,
                                                                            this->kalmanFilterParameters.rows,
                                                                            this->kalmanFilterParameters.inputNumber,
                                                                            this->kalmanFilterParameters.inputNumber, 1);

        double **internalCalcAplusB = this->phiMathPhiMatrixSummation(internalCalcA, internalCalcB,
                                                                      this->kalmanFilterParameters.rows,
                                                                      this->kalmanFilterParameters.inputNumber);

        // obtaining the updated "x" value
        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.x,
                                         internalCalcAplusB,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.inputNumber);

        // x calculation
        /////////////////////

        // phi free
        this->phiKalmanDynamicPhiFree(internalCalcA, this->kalmanFilterParameters.rows, 1);

        this->phiKalmanDynamicPhiFree(internalCalcB, this->kalmanFilterParameters.rows, 1);

        this->phiKalmanDynamicPhiFree(internalCalcAplusB, this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.inputNumber);

        // phi free
    }

    void phiKalmanDynamicSecondPCalculation(phiSystemDynamicParameters *ptrSys)
    {
        ///////////////////////////////////
        // P calculation
        double **internalCalcP1 = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.Amatrices,
                                                                             this->kalmanFilterParameters.P,
                                                                             this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.cols,
                                                                             this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.rows);

        double **internalCalcAtranspose = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.cols);

        this->phiMathPhiTranspose(this->kalmanFilterParameters.Amatrices, this->kalmanFilterParameters.rows,
                                  this->kalmanFilterParameters.cols, internalCalcAtranspose);

        double **internalCalcP2 = this->phiMathPhiVectorMatrixMultiplication(internalCalcP1,
                                                                             internalCalcAtranspose,
                                                                             this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.cols,
                                                                             this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.rows);

        double **internalCalcPtotal = this->phiMathPhiMatrixSummation(internalCalcP2,
                                                                      this->kalmanFilterParameters.Q,
                                                                      this->kalmanFilterParameters.rows,
                                                                      this->kalmanFilterParameters.cols);

        // obtaining updated P values
        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.P,
                                         internalCalcPtotal,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.cols);

        // P calculation
        /////////////////////////////////////

        /////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(internalCalcP1,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);
        this->phiKalmanDynamicPhiFree(internalCalcAtranspose,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(internalCalcP2,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(internalCalcPtotal,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        // phi free
        /////////////////////////////////////
    }

    void phiKalmanDynamicThirdKCalculation(phiSystemDynamicParameters *ptrSys)
    {

        // obtain htranspose
        double **hTranspose = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                 this->kalmanFilterParameters.rows);

        this->phiMathPhiTranspose(this->kalmanFilterParameters.H,
                                  this->kalmanFilterParameters.rows,
                                  this->kalmanFilterParameters.rows,
                                  hTranspose);

        // multiplying process
        double **firstMultiplierFirstTerm = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.P,
                                                                                       hTranspose,
                                                                                       this->kalmanFilterParameters.rows,
                                                                                       this->kalmanFilterParameters.rows,
                                                                                       this->kalmanFilterParameters.rows,
                                                                                       this->kalmanFilterParameters.rows);

        double **secondMultiplierFirstTerm = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.H,
                                                                                        this->kalmanFilterParameters.P,
                                                                                        this->kalmanFilterParameters.rows,
                                                                                        this->kalmanFilterParameters.rows,
                                                                                        this->kalmanFilterParameters.rows,
                                                                                        this->kalmanFilterParameters.rows);

        double **secondMultiplierSecondTerm = this->phiMathPhiVectorMatrixMultiplication(secondMultiplierFirstTerm,
                                                                                         hTranspose,
                                                                                         this->kalmanFilterParameters.rows,
                                                                                         this->kalmanFilterParameters.rows,
                                                                                         this->kalmanFilterParameters.rows,
                                                                                         this->kalmanFilterParameters.rows);

        double **secondMultiplierThirdTerm = this->phiMathPhiMatrixSummation(secondMultiplierSecondTerm,
                                                                             this->kalmanFilterParameters.R,
                                                                             this->kalmanFilterParameters.rows,
                                                                             this->kalmanFilterParameters.rows);

        double **secondMultiplierTotalTerm = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                                this->kalmanFilterParameters.rows);

        this->phiMathInverse(secondMultiplierThirdTerm,
                             secondMultiplierTotalTerm,
                             this->kalmanFilterParameters.rows);

        double **totalTerm = this->phiMathPhiVectorMatrixMultiplication(firstMultiplierFirstTerm,
                                                                        secondMultiplierTotalTerm,
                                                                        this->kalmanFilterParameters.rows,
                                                                        this->kalmanFilterParameters.rows,
                                                                        this->kalmanFilterParameters.rows,
                                                                        this->kalmanFilterParameters.rows);

        // assignment to K value
        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.K,
                                         totalTerm,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.rows);

        ////////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(hTranspose,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(firstMultiplierFirstTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondMultiplierFirstTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondMultiplierSecondTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondMultiplierThirdTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondMultiplierTotalTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(totalTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        // phi free
        ////////////////////////////////////////
    }

    void phiKalmanDynamicFourthXCalculation(phiSystemDynamicParameters *ptrSys)
    {
        // assigning to x value
        double **firstTermFirstPart = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                         1);

        this->phiMathPhiMatrixAssignment(firstTermFirstPart,
                                         this->kalmanFilterParameters.x,
                                         this->kalmanFilterParameters.rows,
                                         1);

        // second part calculation

        double **secondTermFirstPart = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.H,
                                                                                  this->kalmanFilterParameters.x,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  1);

        double **secondTermSecondPart = this->phiMathPhiMatrixNegation(this->kalmanFilterParameters.z,
                                                                       secondTermFirstPart,
                                                                       this->kalmanFilterParameters.rows,
                                                                       1);

        double **secondTermTotalPart = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.K,
                                                                                  secondTermSecondPart,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  1);
        // assignment to updated x value
        double **totalTerm = this->phiMathPhiMatrixSummation(firstTermFirstPart,
                                                             secondTermTotalPart,
                                                             this->kalmanFilterParameters.rows,
                                                             1);

        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.x,
                                         totalTerm,
                                         this->kalmanFilterParameters.rows,
                                         1);
        /////////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(firstTermFirstPart,
                                      this->kalmanFilterParameters.rows,
                                      1);

        this->phiKalmanDynamicPhiFree(secondTermFirstPart,
                                      this->kalmanFilterParameters.rows,
                                      1);

        this->phiKalmanDynamicPhiFree(secondTermSecondPart,
                                      this->kalmanFilterParameters.rows,
                                      1);

        this->phiKalmanDynamicPhiFree(secondTermTotalPart,
                                      this->kalmanFilterParameters.rows,
                                      1);

        this->phiKalmanDynamicPhiFree(totalTerm,
                                      this->kalmanFilterParameters.rows,
                                      1);

        // phi free
        ////////////////////////////////////////
    }

    void phiKalmanDynamicFifthPCalculation(phiSystemDynamicParameters *ptrSys)
    {
        double **firstTermFirstPart = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                         this->kalmanFilterParameters.rows);

        this->phiMathPhiMatrixAssignment(firstTermFirstPart,
                                         this->kalmanFilterParameters.P,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.rows);

        double **secondTermFirstPart = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.K,
                                                                                  this->kalmanFilterParameters.H,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows);

        double **secondTermTotalPart = this->phiMathPhiVectorMatrixMultiplication(secondTermFirstPart,
                                                                                  this->kalmanFilterParameters.P,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows,
                                                                                  this->kalmanFilterParameters.rows);

        double **totalTerm = this->phiMathPhiMatrixNegation(firstTermFirstPart,
                                                            secondTermTotalPart,
                                                            this->kalmanFilterParameters.rows,
                                                            this->kalmanFilterParameters.rows);

        /// assignment to updated P value

        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.P,
                                         totalTerm,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.rows);

        ////////////////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(firstTermFirstPart,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondTermFirstPart,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(secondTermTotalPart,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        this->phiKalmanDynamicPhiFree(totalTerm,
                                      this->kalmanFilterParameters.rows,
                                      this->kalmanFilterParameters.rows);

        // phi free
        ////////////////////////////////////////////////
    }

    void phiKalmanDynamicCalculateFilter(int iterator, phiSystemDynamicParameters *ptrSys, double MNstd, double noiseRange)
    {

        // assigning discrete state/input matrices

        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.Amatrices,
                                         ptrSys->phiSystemDynamicDiscreteStateMatrices,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.cols);

        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.Bmatrices,
                                         ptrSys->phiSystemDynamicDiscreteInputMatrices,
                                         this->kalmanFilterParameters.rows,
                                         this->kalmanFilterParameters.inputNumber);

        // input values
        this->kalmanFilterParameters.inputValue = ptrSys->phiSystemDynamicInputValue;
        this->kalmanFilterParameters.inputNow[0][0] = ptrSys->phiSystemDynamicInputValue;

        // calculating the measurement form z
        double **measurementTerm = this->phiMathPhiVectorMatrixMultiplication(this->kalmanFilterParameters.H,
                                                                              ptrSys->phiSystemDynamicStateNow,
                                                                              this->kalmanFilterParameters.rows,
                                                                              this->kalmanFilterParameters.cols,
                                                                              ptrSys->phiSystemDynamicRows,
                                                                              1);

        this->phiMathPhiMatrixAssignment(this->kalmanFilterParameters.z,
                                         measurementTerm,
                                         ptrSys->phiSystemDynamicRows,
                                         1);

        // adding noise
        for (int i = 0; i < this->kalmanFilterParameters.rows; i++)
        {
            for (int j = 0; j < this->kalmanFilterParameters.inputNumber; j++)
            {
                this->kalmanFilterParameters.z[i][0] += MNstd * this->phiMathPhiRand(-noiseRange, noiseRange);
            }
        }

        ///////////////////////////////////////////
        // KALMAN UPDATE

        this->phiKalmanDynamicFirstXCalculation(ptrSys);

        this->phiKalmanDynamicSecondPCalculation(ptrSys);

        this->phiKalmanDynamicThirdKCalculation(ptrSys);

        this->phiKalmanDynamicFourthXCalculation(ptrSys);

        this->phiKalmanDynamicFifthPCalculation(ptrSys);

        // KALMAN UPDATE
        ///////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        // printing results

        double currentTime = iterator * ptrSys->phiSystemDynamicDt;

        for (int i = 0; i < this->kalmanFilterParameters.rows; i++)
        {
            fprintf(this->phiKalmanDynamicFp, "%lf ", this->kalmanFilterParameters.x[i][0]);
        }
        fprintf(this->phiKalmanDynamicFp, "%lf \n", currentTime);

        // printing results
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(measurementTerm,
                                      ptrSys->phiSystemDynamicRows,
                                      1);

        // phi free
        ///////////////////////////////////////////
    }

    void phiKalmanDynamicInitializeParameters(phiSystemDynamicParameters *ptrSys, double MNstd, double PNstd, double Qcoef, double Rcoef)
    {
        // firstly we need to assign the A and B matrices to the kalman parameters
        this->kalmanFilterParameters.rows = ptrSys->phiSystemDynamicRows;

        this->kalmanFilterParameters.cols = ptrSys->phiSystemDynamicCols;

        this->kalmanFilterParameters.inputNumber = ptrSys->phiSystemDynamicInputNumber;

        this->kalmanFilterParameters.Amatrices = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                                    this->kalmanFilterParameters.cols);

        this->kalmanFilterParameters.Bmatrices = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                                    this->kalmanFilterParameters.inputNumber);

        // state matrices
        this->kalmanFilterParameters.stateNow = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                                   1);

        this->kalmanFilterParameters.statePre = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                                   1);

        this->kalmanFilterParameters.inputNow = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.inputNumber,
                                                                                   1);

        // kalman parameters

        this->kalmanFilterParameters.x = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                            1);

        this->kalmanFilterParameters.K = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                            this->kalmanFilterParameters.rows);

        double detP = 1; // initial value of kalman determinant to update coefficient of this algorithm

        this->kalmanFilterParameters.z = this->phiMathCreatingEmptyMatrices(this->kalmanFilterParameters.rows,
                                                                            1);

        // kalman update parameters
        double **I = this->phiMathEyeMatricesCreation(this->kalmanFilterParameters.rows);
        this->kalmanFilterParameters.Q = this->phiMathPhiSkalarMatrixMultiplication(Qcoef * PNstd * PNstd, I,
                                                                                    this->kalmanFilterParameters.rows,
                                                                                    this->kalmanFilterParameters.cols);

        this->kalmanFilterParameters.R = this->phiMathPhiSkalarMatrixMultiplication(Rcoef * MNstd * MNstd, I,
                                                                                    this->kalmanFilterParameters.rows,
                                                                                    this->kalmanFilterParameters.cols);

        this->kalmanFilterParameters.P = this->phiMathPhiSkalarMatrixMultiplication(MNstd * MNstd, I,
                                                                                    this->kalmanFilterParameters.rows,
                                                                                    this->kalmanFilterParameters.cols);

        // set the specific values on H,R,Q parameters
        this->kalmanFilterParameters.H = this->phiMathEyeMatricesCreation(this->kalmanFilterParameters.rows); // eye matrices 2x2

        ////////////////////////////////////////////////////////////////////
        // initial value assignment

        // for x

        for (int i = 0; i < this->kalmanFilterParameters.rows; i++)
        {
            this->kalmanFilterParameters.x[i][0] = 0.0;
        }
        // for x
        for (int i = 0; i < this->kalmanFilterParameters.rows; i++)
        {
            this->kalmanFilterParameters.z[i][0] = 0.0;
        }

        // initial value assignment
        /////////////////////////////////////////////////////////////////////

        /////////////////////////////////////////
        // phi free

        this->phiKalmanDynamicPhiFree(I, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        // phi free
        /////////////////////////////////////////
    }

    void phiKalmanDynamicPhiExit()
    {
        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.Amatrices, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.Bmatrices, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.inputNumber);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.stateNow, this->kalmanFilterParameters.rows, 1);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.statePre, this->kalmanFilterParameters.rows, 1);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.inputNow, this->kalmanFilterParameters.inputNumber, 1);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.x, this->kalmanFilterParameters.rows, 1);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.K, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        //////////////////////////////////////////////////////////
        // free dynamic memory space

        // kalman parameters
        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.Q, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.H, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.R, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.P, this->kalmanFilterParameters.rows, this->kalmanFilterParameters.cols);

        this->phiKalmanDynamicPhiFree(this->kalmanFilterParameters.z, this->kalmanFilterParameters.rows, 1);

        // free dynamic memory space
        //////////////////////////////////////////////////////////
    }

    void phiKalmanDynamicPhiFree(double **pd,
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

    // FUNCTION PROTOTYPE DECLERATION
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // aided function prototype

    // aided function prototype
    ////////////////////////////////////////////////////////////
};

// CLASS DEFINITION
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////

#endif