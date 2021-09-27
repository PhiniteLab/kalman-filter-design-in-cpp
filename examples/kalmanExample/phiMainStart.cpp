#include "..\..\include\phiSystemDynamicSettings.hpp"

int main()
{
    phiSystemDynamicParameters phiSys("MCKResult.txt");
    phiKalmanDynamicParameters phiKalman("MCKKalmanResult.txt");

    phiSys.phiSystemDynamicRows = 2;        // number of row in state space
    phiSys.phiSystemDynamicCols = 2;        // number of col in state space
    phiSys.phiSystemDynamicInputNumber = 1; // number of input in state space
    phiSys.phiSystemDynamicDt = 0.001;      // sampling period

    phiSys.phiSystemDynamicStateMatrices = phiSys.phiSystemDynamicCreatingEmptyStateMatrices();
    // state matrices creation
    phiSys.phiSystemDynamicInputMatrices = phiSys.phiSystemDynamicCreatingEmptyInputMatrices();
    // input matrices creation

    // matrices value assignment
    phiSys.phiSystemDynamicStateMatrices[0][0] = 0;
    phiSys.phiSystemDynamicStateMatrices[0][1] = 1;
    phiSys.phiSystemDynamicStateMatrices[1][0] = -1;
    phiSys.phiSystemDynamicStateMatrices[1][1] = -2;

    phiSys.phiSystemDynamicInputMatrices[0][0] = 0;
    phiSys.phiSystemDynamicInputMatrices[1][0] = 1;

    printf("Writing state and input matrices \n\n");

    // write the state space matrices
    phiSys.phiSystemDynamicWriteTheMatrices();

    printf("\n\n");

    /// solver solution
    // select the inputValue to be given to the system
    phiSys.phiSystemDynamicInputValue = 10.0;

    //////////////////////////////////////////////////////
    // initialize the model

    double MNstd = 0.004; // measurement noise variance
    double PNstd = 0.002; // process noise variance

    phiSys.phiSystemDynamicInitializeSystemParameters(0.0, 10.0);
    phiKalman.phiKalmanDynamicInitializeParameters(&phiSys, MNstd, PNstd, 0.00001, 0.01);

    // initialize the model
    //////////////////////////////////////////////////////

    // dynamic time parameter solver
    for (int i = 0; i < phiSys.phiSystemDynamicNumberOfIteration; i++)
    {
        phiSys.phiSystemDynamicNiteDynamicSolver(i, MNstd, PNstd, 5.0);
        phiKalman.phiKalmanDynamicCalculateFilter(i, &phiSys, MNstd, 5.0);
    }

    // should be called to free memory
    phiSys.phiSystemDynamicPhiExit();
    phiKalman.phiKalmanDynamicPhiExit();

    return 0;
}
/*
#include "..\..\include\phiSystemDynamicSettings.hpp"

int main()
{
    phiSystemDynamicParameters phiSys("MCKResult.txt");

    phiSys.phiSystemDynamicRows = 2;        // number of row in state space
    phiSys.phiSystemDynamicCols = 2;        // number of col in state space
    phiSys.phiSystemDynamicInputNumber = 1; // number of input in state space
    phiSys.phiSystemDynamicDt = 0.001;      // sampling period

    phiSys.phiSystemDynamicStateMatrices = phiSys.phiSystemDynamicCreatingEmptyStateMatrices();
    // state matrices creation
    phiSys.phiSystemDynamicInputMatrices = phiSys.phiSystemDynamicCreatingEmptyInputMatrices();
    // input matrices creation

    // matrices value assignment
    phiSys.phiSystemDynamicStateMatrices[0][0] = 0;
    phiSys.phiSystemDynamicStateMatrices[0][1] = 1;
    phiSys.phiSystemDynamicStateMatrices[1][0] = -1;
    phiSys.phiSystemDynamicStateMatrices[1][1] = -2;

    phiSys.phiSystemDynamicInputMatrices[0][0] = 0;
    phiSys.phiSystemDynamicInputMatrices[1][0] = 1;

    printf("Writing state and input matrices \n\n");

    // write the state space matrices
    phiSys.phiSystemDynamicWriteTheMatrices();

    double **inverseA = phiSys.phiMathCreatingEmptyMatrices(phiSys.phiSystemDynamicRows,
                                                            phiSys.phiSystemDynamicCols);

    phiSys.phiMathInverse(phiSys.phiSystemDynamicStateMatrices,
                          inverseA,
                          phiSys.phiSystemDynamicRows);

    printf("Inverse A\n");
    for (int i = 0; i < phiSys.phiSystemDynamicRows; i++)
    {
        for (int j = 0; j < phiSys.phiSystemDynamicCols; j++)
        {
            printf("%lf ", inverseA[i][j]);
        }
        printf("\n\n");
    }

    printf("\n\n");

    phiSys.phiSystemDynamicPhiFree(inverseA,
                                   phiSys.phiSystemDynamicRows,
                                   phiSys.phiSystemDynamicCols);

    return 0;
}
*/