# include "schrodinger.h"
# include "Parameter.h"
# include <math.h>
# include <fstream>
# include <iostream>
# include <cstdlib>
# include <math.h>
# include <time.h>
# include <Eigen/Dense>

using namespace std;
using namespace Eigen;

Schrodinger::Schrodinger()
{
    //control how many cpu for parallel computing.

    omp_set_num_threads(8); //for C++ code
    Eigen::initParallel();
    Eigen::setNbThreads(8); //for Eigen Code
}

Schrodinger::~Schrodinger(){

}

void Schrodinger::Schro_Parameter1D(){

    //for square well

    SimTolPoisson=1e-6;

    // x in nm.
    lx=10;
    // x points iitial value
    px=1;
    // x block numbers.
    Mx=1;
    // set x pins
    xpin=new double [Mx+1];

    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }

    // set x mesh steps
    meshx=new double [Mx];
    // mesh unit = 1/nm = 1/step
    for(int i=0;i<Mx;i++){
        meshx[i]=2*(i+1);
    }

    // points calculation
    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }

    // set x point numbers till each block
    xb=new int [Mx+1];
    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }

    L1=px;
}

void Schrodinger::Schro_Parameter2D(){

    //for square well

    SimTolPoisson=1e-6;

    // xy in nm.
    lx=10;
    ly=10;
    // xy points iitial value
    px=py=1;
    // xy block numbers.
    Mx=1;
    My=1;
    // set xy pins
    xpin=new double [Mx+1];
    ypin=new double [My+1];
    for(int i=0;i<Mx+1;i++){
        xpin[i]=0+lx/Mx*i;
    }

    /*
    xpin[0]=0;
    xpin[1]=lx/2-NWR-(1000-fmod(NWR,1000));
    xpin[2]=lx/2+NWR+(1000-fmod(NWR,1000));;
    xpin[3]=lx;
    */

    for(int i=0;i<My+1;i++){
        ypin[i]=0+ly/My*i;
    }

    /*
    ypin[0]=0;
    ypin[1]=50;
    ypin[2]=1000;
    ypin[3]=200000;
    */

    // set xy mesh steps
    meshx=new double [Mx];
    meshy=new double [My];
    // mesh unit = 1/nm = 1/step
    for(int i=0;i<Mx;i++){
        meshx[i]=1*(i+1);
    }

    /*
    meshx[0]=1e-3;
    meshx[1]=1e-2;
    meshx[2]=1e-3;
    */
    for(int i=0;i<My;i++){
        meshy[i]=1*(i+1);
    }
    /*
    meshx[3]=1e-2;
    meshx[4]=1e-3;

    meshy[0]=1;
    meshy[1]=2e-2;
    meshy[2]=1e-3;
    */
    // points calculation
    for(int i=0;i<Mx;i++){
        px=px+meshx[i]*(xpin[i+1]-xpin[i]);
    }
    for(int i=0;i<My;i++){
        py=py+meshy[i]*(ypin[i+1]-ypin[i]);
    }
    // set xy  point numbers till each block
    xb=new int [Mx+1];
    yb=new int [My+1];
    for(int i=1;i<Mx+1;i++){
        xb[0]=0;
        xb[i]=xb[i-1]+(xpin[i]-xpin[i-1])*meshx[i-1];
    }
    for(int i=1;i<My+1;i++){
        yb[0]=0;
        yb[i]=yb[i-1]+(ypin[i]-ypin[i-1])*meshy[i-1];
    }
    L1=px*py;
}

void Schrodinger::Schro_BlockMeshingMesh1D(){

    sample1 = new Material1[L1];

    //assign all coordinate
    Schro_Initialize();

    for(int m=0;m<Mx;m++){

        double a=xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            int pointer = (i);
            sample1[pointer].coordX=a+(i-xb[m])/meshx[m];
        }
    }
}

void Schrodinger::Schro_BlockMeshingMesh2D(){

    sample1 = new Material1[L1];

    //assign all coordinate
    Schro_Initialize();

    for(int m=0;m<Mx;m++){

        double a= xpin[m];

        for(int i=xb[m];i<xb[m+1]+1;i++){
            for (int j=0;j<py;j++){
                int pointer = (px)*(j) + (i);
                sample1[pointer].coordX=a+(i-xb[m])/meshx[m];
            }
        }
    }

    for(int m=0;m<My;m++){

        double a= ypin[m];

        for (int i=0;i<px;i++){
            for(int j=yb[m];j<yb[m+1]+1;j++){
                int pointer = (px)*(j) + (i);
                sample1[pointer].coordY=a+(j-yb[m])/meshy[m];
            }
        }
    }
}

void Schrodinger::Schro_Initialize(){

    #pragma omp parallel for
    for(int i=0;i<L1;i++){
        sample1[i].coordX=0;
        sample1[i].coordY=0;
        sample1[i].coordZ=0;
        sample1[i].k=80;
        sample1[i].dop=0;
        sample1[i].phi=0;
        sample1[i].u=0;
        sample1[i].v=0;
        sample1[i].r=0;
        sample1[i].rho=0;
        sample1[i].Crho=0;
        sample1[i].mun=0;
        sample1[i].mup=0;
        sample1[i].tau=0;
        sample1[i].Ex=0;
        sample1[i].Ey=0;
        sample1[i].Ez=0;
        sample1[i].flag=0;
    }

}

void Schrodinger::Schro_InitialGuess1D(){

    // no use

    /*

#pragma omp parallel for
    for (int i=0;i<px;i++){

        int pointer = i;

        //Oxide
        sample1[pointer].k=SiO2_permi;
        sample1[pointer].flag=2;

        //N Wire
        if(sample1[pointer].coordX>=(NWRcenterx-NWRradiusx) && sample1[pointer].coordX<=(NWRcenterx+NWRradiusx)){
            if(sample1[pointer].coordY>=(NWRcenterz-NWRradiusz) && sample1[pointer].coordY<=(NWRcenterz+NWRradiusz)){
                sample1[pointer].k=Si_permi;
                sample1[pointer].flag=1;
                //sample1[pointer].dop=Ndi;
                //sample1[pointer].phi=0;
                //sample1[pointer].u=exp((-1)*0/VT);
                //sample1[pointer].v=exp(0/VT);
                //sample1[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                //sample1[pointer].mup=mupCal(Tamb, 0, 1);
                //sample1[pointer].tau=tauPCal(0);
                //sample1[pointer].r=SRHrecomb2D(i,j);
            }
        }
    }
    */
}

void Schrodinger::Schro_InitialGuess2D(){

    // no use

    /*
#pragma omp parallel for
    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){

            int pointer = (px)*(j) + (i);

            //Oxide
            sample1[pointer].k=SiO2_permi;
            sample1[pointer].flag=2;

            //N Wire
            if(sample1[pointer].coordX>=(NWRcenterx-NWRradiusx) && sample1[pointer].coordX<=(NWRcenterx+NWRradiusx)){
                if(sample1[pointer].coordY>=(NWRcenterz-NWRradiusz) && sample1[pointer].coordY<=(NWRcenterz+NWRradiusz)){
                    sample1[pointer].k=Si_permi;
                    sample1[pointer].flag=1;
                    //sample1[pointer].dop=Ndi;
                    //sample1[pointer].phi=0;
                    //sample1[pointer].u=exp((-1)*0/VT);
                    //sample1[pointer].v=exp(0/VT);
                    //sample1[pointer].mun=munCal(Tamb, 0, 1); // max Na Nd
                    //sample1[pointer].mup=mupCal(Tamb, 0, 1);
                    //sample1[pointer].tau=tauPCal(0);
                    //sample1[pointer].r=SRHrecomb2D(i,j);
                }
            }
        }
    }
    */
}

double Schrodinger::Schro_PoissonSolver2D(){

    loop=0;

    double errPhi(0),errPhi_max(0);

    do{
        loop++;

        errPhi=Schro_PoissonGaussSeidel2D();

        if(errPhi_max < errPhi) {errPhi_max=errPhi;}

        if(loop<20 || loop%100==0)
        cerr <<"PS:"<< loop <<"\t" <<errPhi<<"\t"<<errPhi_max<<endl;


    }while(errPhi>SimTolPoisson);

    return errPhi_max;

}

double Schrodinger::Schro_PoissonGaussSeidel2D(){

    double max_val=0;

#pragma omp parallel for reduction(max:max_val)
    for (int i=1; i<px-1; i++) {
        for (int j=1; j<py-1; j++) {

            int pointer = (px)*(j) + (i);

            double phik=sample1[pointer].phi;

            sample1[pointer].phi=Schro_PoissonGaussSeidelInner2D(i,j);

            double error=abs(sample1[pointer].phi-phik);

            error=error/(abs(phik)+1);

            if(error>max_val)
                max_val=error;
        }
    }

    return max_val;

}

double Schrodinger::Schro_PoissonGaussSeidelInner2D(int i, int j){

    int pointer = (px)*(j) + (i);
    int pointer_ip =   (px)*(j) + (i+1);
    int pointer_in =   (px)*(j) + (i-1);
    int pointer_jp =   (px)*(j+1) + (i);
    int pointer_jn =   (px)*(j-1) + (i);

    double permitivity_ip=sample1[pointer].k*sample1[pointer_ip].k / (0.5*sample1[pointer_ip].k+0.5*sample1[pointer].k);
    double permitivity_in=sample1[pointer].k*sample1[pointer_in].k / (0.5*sample1[pointer_in].k+0.5*sample1[pointer].k);
    double permitivity_jp=sample1[pointer].k*sample1[pointer_jp].k / (0.5*sample1[pointer_jp].k+0.5*sample1[pointer].k);
    double permitivity_jn=sample1[pointer].k*sample1[pointer_jn].k / (0.5*sample1[pointer_jn].k+0.5*sample1[pointer].k);

    double deltax=abs(sample1[pointer_ip].coordX-sample1[pointer_in].coordX)/2;
    double deltay=abs(sample1[pointer_jp].coordY-sample1[pointer_jn].coordY)/2;
    double xstep_p=abs(sample1[pointer_ip].coordX-sample1[pointer].coordX);
    double xstep_n=abs(sample1[pointer_in].coordX-sample1[pointer].coordX);
    double ystep_p=abs(sample1[pointer_jp].coordY-sample1[pointer].coordY);
    double ystep_n=abs(sample1[pointer_jn].coordY-sample1[pointer].coordY);

    double f,df,phik;
    double volume=deltax*deltay;

    phik=sample1[pointer].phi;

    if(sample1[pointer].flag==0){ // initialized state
        f=volume*sample1[pointer].rho/e0;
        df=0;
    }
    else if(sample1[pointer].flag==1||sample1[pointer].flag==4){ // semiconductor
        f=volume*ni_nm*(sample1[pointer].u*exp(phik/VT)-sample1[pointer].v*exp((-1)*phik/VT)-sample1[pointer].dop)*(-1)*q0/e0;
       df=volume*ni_nm*(sample1[pointer].v*exp((-1)*phik/VT))*(-1)*q0/e0/VT;
    }
    else if(sample1[pointer].flag==2||sample1[pointer].flag==5){ // insulator, Glass substrate
        f=volume*sample1[pointer].rho/e0;
        df=0;
    }
    else{
        cerr << "error, undetermined flag for poisson behavior."<<endl;
    }

    return (((permitivity_ip*sample1[pointer_ip].phi/xstep_p+permitivity_in*sample1[pointer_in].phi/xstep_n)*deltay
            +(permitivity_jp*sample1[pointer_jp].phi/ystep_p+permitivity_jn*sample1[pointer_jn].phi/ystep_n)*deltax
            + f - df*phik )
            /
            ((permitivity_ip/xstep_p+permitivity_in/xstep_n)*deltay
            +(permitivity_jp/ystep_p+permitivity_jn/ystep_n)*deltax - df));
}

void Schrodinger::Schro_ElectronConcentration(){

    Nq_C = new Nq[gridL];
    double C1(0),C2(0),C3(0),C4(0),C5(0),C6(0),C7(0),C8(0),C9(0);

    for(int i=0;i<gridx;i++){
        for(int j=0;j<gridy;j++){
            int pointer =(gridy)*(i) + (j);
            Nq_C[pointer].nq_1=pow(Psi_C[pointer].psi_1,2);
            C1=C1+pow(Psi_C[pointer].psi_1,2);
            Nq_C[pointer].nq_2=pow(Psi_C[pointer].psi_2,2);
            C2=C2+pow(Psi_C[pointer].psi_2,2);
            Nq_C[pointer].nq_3=pow(Psi_C[pointer].psi_3,2);
            C3=C3+pow(Psi_C[pointer].psi_3,2);
            Nq_C[pointer].nq_4=pow(Psi_C[pointer].psi_4,2);
            C4=C4+pow(Psi_C[pointer].psi_4,2);
            Nq_C[pointer].nq_5=pow(Psi_C[pointer].psi_5,2);
            C5=C5+pow(Psi_C[pointer].psi_5,2);
            Nq_C[pointer].nq_6=pow(Psi_C[pointer].psi_6,2);
            C6=C6+pow(Psi_C[pointer].psi_6,2);
            Nq_C[pointer].nq_7=pow(Psi_C[pointer].psi_7,2);
            C7=C7+pow(Psi_C[pointer].psi_7,2);
            Nq_C[pointer].nq_8=pow(Psi_C[pointer].psi_8,2);
            C8=C8+pow(Psi_C[pointer].psi_8,2);
            Nq_C[pointer].nq_9=pow(Psi_C[pointer].psi_9,2);
            C9=C9+pow(Psi_C[pointer].psi_9,2);
        }
    }
}

void Schrodinger::Schro_PrintNq2D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tPsi_1(3)\tPsi_2(4)\tPsi_3(5)\tPsi_4(6)\tPsi_5(7)\tPsi_6(8)\tPsi_7(9)\tPsi_8(10)\tPsi_9(11)\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<gridx;i++){
        for (int j=0;j<gridy;j++){
            int pointer =(gridy)*(i) + (j);
            int pointer2 =(py)*(i+NWRleft) + (j+NWRbottom);
            output << sample1[pointer2].coordX << '\t' << sample1[pointer2].coordY << '\t'
                   << Nq_C[pointer].nq_1 << '\t' << Nq_C[pointer].nq_2 << '\t' << Nq_C[pointer].nq_3 << '\t'
                   << Nq_C[pointer].nq_4 << '\t' << Nq_C[pointer].nq_5 << '\t' << Nq_C[pointer].nq_6 << '\t'
                   << Nq_C[pointer].nq_7 << '\t'<< Nq_C[pointer].nq_8 << '\t'<< Nq_C[pointer].nq_9 << '\t' << endl;
        }
    }

    output.close();
}

void Schrodinger::Schro_PrintMaterial2D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tY(2)\tK(3)\tdop(4)\tphi(5)\tu(6)\tv(7)\tr(8)\trho(9)\tmun(10)\tmup(11)\ttau(12)\tEx(13)\tEy(14)\tflag(19)\tCrho(20)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        for (int j=0;j<py;j++){
            int pointer =(px)*(j) + (i);
            output << sample1[pointer].coordX << '\t' << sample1[pointer].coordY << '\t'
                   << sample1[pointer].k << '\t' <<sample1[pointer].dop << '\t' <<sample1[pointer].phi << '\t'
                   << sample1[pointer].u << '\t' << sample1[pointer].v << '\t' << sample1[pointer].r << '\t'
                   << sample1[pointer].rho << '\t'<< sample1[pointer].mun << '\t'<< sample1[pointer].mup << '\t'
                   << sample1[pointer].tau << '\t'<< sample1[pointer].Ex << '\t'<< sample1[pointer].Ey<< '\t'
                   << sample1[pointer].flag << '\t'<< sample1[pointer].Crho <<endl;

        }
    }

    output.close();
}

void Schrodinger::Schro_PrintMaterial1D(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(6);


    output << "X(1)\tK(2)\tdop(3)\tphi(4)\tu(5)\tv(6)\tr(7)\trho(8)\tmun(9)\tmup(10)\ttau(11)\tEx(12)\tEy(13)\tflag(14)\tCrho(15)#"<<endl;
    output << "[nm]\t[nm]\t[nm]\t#"<<endl;
    output <<"--------------------------------------------------------------------------------------------------------------------------------#" << endl;

    for (int i=0;i<px;i++){
        int pointer = (i);
        output << sample1[pointer].coordX << '\t'
               << sample1[pointer].k << '\t' <<sample1[pointer].dop << '\t' <<sample1[pointer].phi << '\t'
               << sample1[pointer].u << '\t' << sample1[pointer].v << '\t' << sample1[pointer].r << '\t'
               << sample1[pointer].rho << '\t'<< sample1[pointer].mun << '\t'<< sample1[pointer].mup << '\t'
               << sample1[pointer].tau << '\t'<< sample1[pointer].Ex << '\t'<< sample1[pointer].Ey<< '\t'
               << sample1[pointer].flag << '\t'<< sample1[pointer].Crho <<endl;

    }

    output.close();
}

void Schrodinger::Schro_DeclareHamiltonianArray1D(){

    gridx=px;

    gridL=gridx;

    H_m = MatrixXd::Zero(gridL,gridL);

}

void Schrodinger::Schro_AssignHamiltonianArray1D(){

    // uniform mesh is assumed, 0.5nm
    // the potential unit is eV.
    // Assuming a square well having 0 potential.
    // Outside of the well having infinite potential barrier.

    double delta = 1/meshx[0];

    double a = (-1)*hbar*hbar_eV/(2*m0*delta*1e-9*delta*1e-9);

    int row = 0;

    //H_m is L*L matrix, in the following loop row = i*gridy+j.

    for (int i=0;i<gridL;i++){

        //H_m(row,i*gridy+j) = ((-1)*phi_m(i,j)+0.55) - 4*a; // (-1)*V + Eg/2 -4a

        H_m(row,i) = (-1)*0 - 2*a; // -2a, V=0

        //set potential barrier at boundary point
        if (i==0 || i==gridL-1 ){
            H_m(row,i) = (-1)*(-999) - 2*a; // -2a, V=0
        }

        /*
        //set potential barrier for double well
        if (i==0 || i== gridL-1 || ((i > 2*gridL/5) && (i < 3*gridL/5)) ){
            H_m(row,i) = (-1)*(-999) - 2*a; // -2a, V=0
        }
        */

        //(i-1)
        if (0 <= (i-1) && (i-1) <= gridL-1){
            H_m(row,i-1) = a;
        }

        //(i+1)
        if (0 <= (i+1) && (i+1) <= gridL-1){
            H_m(row,i+1) = a;
        }

        row = row + 1;
    }
}

void Schrodinger::Schro_DeclareHamiltonianArray2D(){

    gridx=px;
    gridy=py;

    gridL=gridx*gridy;

    H_m = MatrixXd::Zero(gridL,gridL);

}

void Schrodinger::Schro_AssignHamiltonianArray2D(){

    // uniform mesh is assumed, 1nm
    // the potential unit is eV.
    // Assuming a square well having 0 potential.
    // Outside of the well having infinite potential barrier.

    double delta = 1/meshx[0];

    double a = (-1)*hbar*hbar_eV/(2*m0*delta*1e-9*delta*1e-9);

    int row = 0;

    //H_m is L*L matrix, in the following loop row = i*gridy+j.

    for (int i=0;i<gridx;i++){
        for (int j=0;j<gridy;j++){

            // (-1)*V + Eg/2 -4a

            H_m(row,i*gridy+j) = (0) - 4*a; // -4a, V=0


            //set potential barrier at boundary point
            if (i==0 || i==gridx-1 || j==0 || j==gridy-1 ){
                H_m(row,i*gridy+j) = (-1)*(-999) - 4*a; // -2a, V=0
            }

            //(i-1)
            if (0 <= (i-1)*gridy+j && (i-1)*gridy+j <= gridL-1){
                H_m(row,(i-1)*gridy+j) = a;
            }

            //(i+1)
            if (0 <= (i+1)*gridy+j && (i+1)*gridy+j <= gridL-1){
                H_m(row,(i+1)*gridy+j) = a;
            }

            //(j-1)
            if (0 <= (i)*gridy+j-1 && (i)*gridy+j-1 <= gridL-1){
                H_m(row,(i)*gridy+j-1) = a;
            }

            //(j+1)
            if (0 <= (i)*gridy+j+1 && (i)*gridy+j+1 <= gridL-1){
                H_m(row,(i)*gridy+j+1) = a;
            }

            //the assignment above will let the j direction boundary have excess coefficient.
            //j-1 at j=0 should be 0, and j+j at j=gridy-1 should be 0.

            if ((j%gridy) == 0 && (i)*gridy+j != 0){
                H_m(row, (i)*gridy+j-1) = 0;
            }

            if ((j%gridy) == (gridy-1) && (i)*gridy+j != gridL-1){
                H_m(row, (i)*gridy+j+1) = 0;
            }

            /*x direction do not have excess coefficient is because at (i-1)
             *
             *if (0 <= (i-1)*gridy+j) already exclude i=0
             *
             *at (i+1)
             *if((i+1)*gridy+j <= gridL-1) already exclude i=gridx
             *
             */

            row = row + 1;
        }
    }
}

void Schrodinger::Schro_PrintH(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    output<< H_m;

    output.close();
}

void Schrodinger::Schro_SchrodingerSolver(){

    cerr << "Solving Eigen, "<< nbThreads() << " parallel threads." <<endl;

    es.compute(H_m);

    EigenValues_m.resize(gridL,1);
    EigenVectors_m.resize(gridL,gridL);


    #pragma omp parallel for
    for(int i=0;i<gridL;i++){
        EigenValues_m(i,0)=real(es.eigenvalues()(i));
    }

    #pragma omp parallel for
    for(int i=0;i<gridL;i++){
        for(int j=0;j<gridL;j++){
            EigenVectors_m(i,j)=real(es.eigenvectors()(i,j));
        }
    }
}

void Schrodinger::Schro_SortEigen_Bubble(){

    // bubble sort, not efficient.
    for(int i=0;i<es.eigenvalues().rows();i++){
        for(int j=0;j<es.eigenvalues().rows()-i-1;j++){
            if(EigenValues_m(j,0) > EigenValues_m(j+1,0)){

                double temp=EigenValues_m(j,0);
                EigenValues_m(j,0)=EigenValues_m(j+1,0);
                EigenValues_m(j+1,0)=temp;

                MatrixXd temp_m(gridL,1);
                temp_m=EigenVectors_m.col(j);
                EigenVectors_m.col(j)=EigenVectors_m.col(j+1);
                EigenVectors_m.col(j+1)=temp_m;

            }
        }
    }
}

void Schrodinger::Schro_SortEigen_Merge(){

    //Merge Sort
    cerr << "Sorting Eigenvalues and Eigenvectors."<<endl;
    Schro_MergeSort(0,gridL-1);

}

void Schrodinger::Schro_PrintEigenValues(const char *path){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    output<< EigenValues_m;

    output.close();
}

void Schrodinger::Schro_PrintEigenVectors2D(const char *path, int nb){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    //output<< EigenVectors_m;

    for (int i=0;i<gridx;i++){
        for (int j=0;j<gridy;j++){
            int pointer =(gridy)*(i) + (j);
            output<<i<<'\t'<<j<<'\t'<< EigenVectors_m(pointer,nb)<<endl;
        }
    }

    output.close();
}

void Schrodinger::Schro_PrintEigenVectors1D(const char *path, int nb){

    fstream output;

    output.open(path, fstream::out | fstream::trunc);

    output.precision(4);

    //output<< EigenVectors_m;

    for (int i=0;i<gridx;i++){
            int pointer = (i) ;
            output<<i<<'\t'<< EigenVectors_m(pointer,nb)<<endl;
    }

    output.close();
}

void Schrodinger::Schro_Merge(int low,int mid,int high){

    int h,i,j,k;
    MatrixXd b(gridL,1);
    MatrixXd temp_b(gridL,gridL);
    h=low;
    i=low;
    j=mid+1;

    while((h<=mid)&&(j<=high)){
        if(EigenValues_m(h,0)<=EigenValues_m(j,0)){
            b(i,0)=EigenValues_m(h,0);
            temp_b.col(i)=EigenVectors_m.col(h);
            h++;
        }
        else{
            b(i,0)=EigenValues_m(j,0);
            temp_b.col(i)=EigenVectors_m.col(j);
            j++;
        }
        i++;
    }
    if(h>mid){
        for(k=j;k<=high;k++){
            b(i,0)=EigenValues_m(k,0);
            temp_b.col(i)=EigenVectors_m.col(k);
            i++;
        }
    }
    else{
        for(k=h;k<=mid;k++){
            b(i,0)=EigenValues_m(k,0);
            temp_b.col(i)=EigenVectors_m.col(k);
            i++;
        }
    }
    for(k=low;k<=high;k++){
        EigenValues_m(k,0)=b(k,0);
        EigenVectors_m.col(k)=temp_b.col(k);
    }
}

void Schrodinger::Schro_MergeSort(int low, int high){

    int mid;
    if(low<high){
        mid=(low+high)/2;
        Schro_MergeSort(low,mid);
        Schro_MergeSort(mid+1,high);
        Schro_Merge(low,mid,high);
    }
}
