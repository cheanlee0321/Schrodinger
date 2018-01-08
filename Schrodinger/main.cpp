#include <iostream>
#include "schrodinger.h"

using namespace std;

int main()
{

    /*
    //1D single well Schrodinger, single well and double well
    Schrodinger *test=new Schrodinger();
    test->Schro_Parameter1D();
    test->Schro_BlockMeshingMesh1D();
    test->Schro_InitialGuess1D();
    test->Schro_PrintMaterial1D("Struct.txt");
    test->Schro_DeclareHamiltonianArray1D();
    test->Schro_AssignHamiltonianArray1D();
    test->Schro_PrintH("H.txt");
    test->Schro_SchrodingerSolver();
    test->Schro_SortEigen_Merge();
    test->Schro_PrintEigenValues("EigenValues.txt");
    test->Schro_PrintEigenVectors1D("EigenVectors0.txt", 0);
    test->Schro_PrintEigenVectors1D("EigenVectors1.txt", 1);
    test->Schro_PrintEigenVectors1D("EigenVectors2.txt", 2);
    test->Schro_PrintEigenVectors1D("EigenVectors3.txt", 3);
    */

    //2D Squre well Schrodinger
    Schrodinger *test=new Schrodinger();
    test->Schro_Parameter2D();
    test->Schro_BlockMeshingMesh2D();
    test->Schro_InitialGuess2D();
    test->Schro_PrintMaterial2D("Struct.txt");
    test->Schro_DeclareHamiltonianArray2D();
    test->Schro_AssignHamiltonianArray2D();
    test->Schro_PrintH("H.txt");
    test->Schro_SchrodingerSolver();
    test->Schro_SortEigen_Merge();
    test->Schro_PrintEigenValues("EigenValues.txt");
    test->Schro_PrintEigenVectors2D("EigenVectors0.txt", 0);
    test->Schro_PrintEigenVectors2D("EigenVectors1.txt", 1);
    test->Schro_PrintEigenVectors2D("EigenVectors2.txt", 2);
    test->Schro_PrintEigenVectors2D("EigenVectors3.txt", 3);
    /*
    */
}

