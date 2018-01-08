#ifndef SCHRODINGER_H
#define SCHRODINGER_H

# include <eigen3/Eigen/Dense>
using namespace Eigen;


struct Material1 {
    double coordX; //coordinate
    double coordY;
    double coordZ;
    double k;   // permitivitty
    double dop; // doping, -Nai Ndi
    double phi; // potential, Ei
    double u;   // quasi-fermi level, exp(-Ef/Vt)
    double v;   //                  , exp(Ef/Vt)
    double r;   // recombination rate, SRH
    double rho; // charge density, C/cm3
    double Crho; // charge density, C/cm3
    double mun; // mobility, nm^2/v.s
    double mup;
    double tau;
    double Ex;  // electric field
    double Ey;
    double Ez;
    double flag;
};

struct Nq {
    double nq_1;
    double nq_2;
    double nq_3;
    double nq_4;
    double nq_5;
    double nq_6;
    double nq_7;
    double nq_8;
    double nq_9;
};

struct Psi {
    double psi_1;
    double psi_2;
    double psi_3;
    double psi_4;
    double psi_5;
    double psi_6;
    double psi_7;
    double psi_8;
    double psi_9;
};

class Schrodinger
{

public:
    Schrodinger();
    ~Schrodinger();

    void Schro_Parameter1D();
    void Schro_Parameter2D();
    void Schro_BlockMeshingMesh1D();
    void Schro_BlockMeshingMesh2D();
    void Schro_Initialize();
    void Schro_InitialGuess1D();
    void Schro_InitialGuess2D();

    void Schro_DeclareHamiltonianArray1D();
    void Schro_AssignHamiltonianArray1D();
    void Schro_DeclareHamiltonianArray2D();
    void Schro_AssignHamiltonianArray2D();

    void Schro_SchrodingerSolver();
    void Schro_SortEigen_Bubble();
    void Schro_SortEigen_Merge();

    double Schro_PoissonSolver2D();
    double Schro_PoissonGaussSeidel2D();
    double Schro_PoissonGaussSeidelInner2D(int i, int j);

    void Schro_ElectronConcentration();
    void Schro_PrintNq2D(const char *path);

    void Schro_PrintMaterial1D(const char *path);
    void Schro_PrintMaterial2D(const char *path);
    void Schro_PrintH(const char *path);
    void Schro_PrintEigenValues(const char *path);
    void Schro_PrintEigenVectors1D(const char *path, int nb);
    void Schro_PrintEigenVectors2D(const char *path, int nb);

protected:
    Material1 *sample1;

private:

    void Schro_Merge(int low, int mid, int high);
    void Schro_MergeSort(int low, int high);

    int Mx,My,Mz,px,py,pz,*xb,*yb,*zb,L1,loop;
    double lx,ly,lz,*xpin,*ypin,*zpin,*meshx,*meshy,*meshz;
    double NWRradiusx,NWRradiusz,NWRlength,NWRcenterx,NWRcenterz;
    int NWRleft,NWRright,NWRtop,NWRbottom,NWRcenterxP,NWRcenterzP;

    double SimTolPoisson;
    double ni_nm;

    //*********************************
    //Schrodinger Poisson
    //*********************************
    MatrixXd phi_m, EigenValues_m, EigenVectors_m;
    MatrixXd H_m;
    EigenSolver<MatrixXd> es;
    int gridx,gridy,gridL;
    Psi *Psi_C;
    Nq *Nq_C;

};

#endif // SCHRODINGER_H
