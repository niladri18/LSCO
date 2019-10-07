#define ARMA_DONT_USE_WRAPPER
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <complex>
#include <cmath>
#include <math.h>
#include "lattice.h"
#include "matrix.h"
#include "gauleg_double.cpp"
#include <armadillo>
#include <time.h>
#include <omp.h>
//#include </opt/armadillo/9.200.7/gcc-4.8.5/include/armadillo>
using namespace std;
using namespace arma;

	/* LATTICE 
	 * COPYRIGHT: NILADRI GOMES, 
	 * EMIAL: niladri.gomes@gmail.com
	
	   		       y
	   		       ^
	   		       |
   	O---Cu---O---Cu---O---Cu-->x
	     |	      |        |
	     O	      O        O
	     |	      |        |     
   	O---Cu---O---Cu---O---Cu--
	     |	      |        |
	     O	      O        O
	     |	      |        |     
   	O---Cu---O---Cu---O---Cu--
	     |	      |        |
	     O	      O        O
	     |	      |        |     
   	O---Cu---O---Cu---O---Cu--

	Orbital ordering :
	|0> = |px up>, |1> = |px dn>, |2> = |py up>, |3> = |py dn>, |4> = |pz up>, |5> = |pz dn>



	*/


int flag_omp = 1;
/* Define parameters */
double Haxis[3] = {0,0,1}; //Direction of the magnetic field

	//double ts = 1.0, tp = 0.2, U = 1.0;

double a = 1.0; // lattice constant
double **Hintra, **Hxy, **Hxx, **Hyy; 
Cell *data;
double **superL, **unitL;
double rho = 11.0/12.0;

cx_mat DENSITY_MATRIX, DENSITY_MATRIX_tmp, H_U, H_mag;
cx_cube GU;
cube Eigval;

int ConstructLattice(Cell *data);
int PrintLattice(Cell *data, int num_check);
int initHop(double ts, double tp);
int initLattice();
int initQuadrature(Quad *quadrature);
int compute_hopmatrix(cx_mat& hopmatrix, double kx, double ky,  Cell *data, int lines);
double HartreeFock(Cell *data, Quad *quadrature, int lines, double U, double J);
double HartreeFockParallel(Cell *data, Quad *quadrature, int lines, double U, double J);
double FindFermi(Cell *data, Quad *quadrature, int lines);
int DecoupleHF( cx_mat& H_U, double U, double J);
double errorDensityMatrix();
int ComputeMagMoment(mat& moment);
double maxDiff(mat& moment1, mat& moment2);
void printMoment(mat& moment);
int ComputeCharge(vec& charge);
//void saveBand(Quad *quadrature );
double saveBand(Cell *data, Quad *quadrature, int lines, double U, double J, double Hz);
void saveMoment(mat& moment, vec& charge, double U, double Hz);
void saveFullBand( Quad *quadrature, double U, double Hz);
void saveArpes(Cell *data, Quad *quadrature, int lines, double U, double J, double Hz);
int constructHmag(double Hz, double direction[3]);
void saveBandGap(Quad *quadrature, double U, double Hz);

/* Utilitites */
int SpinorGenerator(double J, cx_cube& G);
int makeSpinOperator();







int delta(double m1, double m2){
	/* Delta function */
	double x = 0;
	if (m1 == m2) x = 1;
	return x;
}

double findmax(double a, double b){
	if (a > b){
		return a;
	}
	else{
		return b;
	}
}

int main(int argc, char *argv[]){
	
	time_t start_t, end_t;
	start_t = time(NULL);
	//int omp_num = omp_get_max_threads();
	//int omp_num = omp_set_max_threads();
	//printf("Number of threads: %d\n",omp_num);
	//exit(1);
	double ts = 1.0, tp = 0.2, U, J, Hz;
	if (argc !=3){
		printf("USAGE: %s <U> <Hz> \n",argv[0]);
		exit(1);
	}
	U = atof(argv[1]);
	J = 0.0*U;
	Hz = atof(argv[2]);
	//J = 0;
	Hintra = mkmatrixd(NU,NU);
	Hxy = mkmatrixd(NU,NU);
	Hxx = mkmatrixd(NU,NU);
	Hyy = mkmatrixd(NU,NU);
	superL = mkmatrixd(6,2);
	unitL = mkmatrixd(2,2);

	int num = 4000,lines;
	data = (Cell*) malloc(sizeof(Cell)*num);
	initHop(ts,tp);
	printmatrixd("Hintra",Hintra,NU,NU);
	printmatrixd("Hxy",Hxy,NU,NU);
	printmatrixd("Hxx",Hxx,NU,NU);
	printmatrixd("Hyy",Hyy,NU,NU);

	initLattice();
	printmatrixd("superL",superL,6,2);
	printmatrixd("unitL",unitL,2,2);

	lines = ConstructLattice(data);

	PrintLattice(data,lines);

	Quad quadrature;
	initQuadrature(&quadrature);

	/* Create the spinor operator GU for each site */
	makeSpinOperator();
	//GU.print("Generators");
	
	/* Inititalize density matrix */
	//FILE *fden;
	//fden = fopen("density.rho","rw");
	//ifstream fb("density.rho");
	//if (fden != NULL){
	//	DENSITY_MATRIX.load("density.rho");
	//}	
	//else{
	//	DENSITY_MATRIX.eye(NU,NU);
	//}
	//DENSITY_MATRIX.randu(NU,NU);
	//DENSITY_MATRIX.zeros(NU,NU);
	DENSITY_MATRIX.eye(NU,NU);
	DENSITY_MATRIX_tmp.eye(NU,NU);
	Eigval.zeros(Nintx, Ninty, NU);
	//DENSITY_MATRIX_tmp.randu(NU,NU);
	//DENSITY_MATRIX.print("Density Matrix");
	
	mat moment(2,3), oldmoment(2,3);	
	vec charge(2);
	ComputeMagMoment(moment);
	oldmoment = moment;

	/* Construct Magnetic field contribution */
	constructHmag(Hz, Haxis);
	H_mag.print("Hmag");
	//ComputeMagMoment(moment);
	//moment.print("Moment");
	//HartreeFock(data, &quadrature, lines);
	double Etot, error;
	int iter;
	for (int i = 0; i < 1000; i++){
		if (flag_omp){
			printf("Using parallel thread\n");
			Etot = HartreeFockParallel(data, &quadrature, lines, U, J);
		}
		else{
			Etot = HartreeFock(data, &quadrature, lines, U, J);
		}
		ComputeMagMoment(moment);
		//error = maxDiff(moment, oldmoment);
		iter = i;
		//oldmoment = moment;
		error = errorDensityMatrix();
		printf("Energy = %lf  at step %d Error: %f\n",Etot, i, error);
		if(error < 10e-10) break;
		oldmoment = moment;
		//cout<<"Error"<<endl;
		//cout<<error<<endl;
	}
	printf("Converged in %d steps. \n",iter);

	//if(fden == NULL){	
	//	DENSITY_MATRIX.save("density.rho");
	//}
	//saveDensityMatrix();

	ComputeMagMoment(moment);
	//printMoment(moment);
	ComputeCharge(charge);
	saveMoment(moment, charge, U, Hz);
	charge.print("Charge:");

	saveBand(data, &quadrature, lines, U, J, Hz);
	saveArpes(data, &quadrature, lines, U, J, Hz);
	saveFullBand(&quadrature, U, Hz);
	saveBandGap (&quadrature, U, Hz);

	end_t = time(NULL);
	printf("Finished in %lds.\n",end_t-start_t);
	//cx_mat hopmatrix(NU,NU);
	//hopmatrix.zeros();
	//for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) hopmatrix(i,j) = {0.0,-0.5};
	//hopmatrix.print("hopmatrix:");


	//compute_hopmatrix( hopmatrix, 2, 2, data, lines);
	//vec eval;
	//cx_mat evec;
	//eig_sym(eval,evec,hopmatrix);

	//hopmatrix.print("hopmatrix:");
	//cout<<hopmatrix.is_hermitian()<<endl;

	//for (int i = 0; i < 10; i++){
	//	cout << "Hi" << endl;
	//}
	return 0;
}

int initLattice(){
	/* initialize the coordinates of the unit cell and the nearest neighbors */
	double stmp[6][2] = { {0,1},
		{1,1},
		{1,0},
		{0,-1},
		{-1,-1},
		{-1,0}
	};

	// coordinates for O1- and O2-
	double utmp[2][2] = { {-0.5,0.0},{0.0,0.5} }; 	
	for (int i = 0; i < 6; i++){
		for (int j = 0; j<2; j++){
			superL[i][j] = stmp[i][j];
		}
	}

	for (int i = 0; i < 2; i++){
		for (int j = 0; j<2; j++){
			unitL[i][j] = utmp[i][j];
		}
	}
	return 0;
}

int ConstructLattice(Cell *data){


	int numcheck = 0;
	double x,y, dist;

	/*	
	for (int mu = 0; mu < NU; mu++) {
		for (int nu = 0; nu < NU; nu++){
			data[numcheck].x = 0;
			data[numcheck].y = 0;
			data[numcheck].mu = mu;	
			data[numcheck].nu = nu;
			data[numcheck].tmunu = Hxy[mu][nu];
			data[numcheck].dist = 0;	
			//numcheck += 1;

			if (mu/NC!=nu/NC){
				x0 = 0.0;
				y0 = 0.0;
				data[numcheck].x = 0;
				data[numcheck].y = 0;
				data[numcheck].mu = mu;
				data[numcheck].nu = nu;
				data[numcheck].tmunu = Hxy[mu][nu];
				data[numcheck].dist = sqrt(0.5*sqrt(2));
				numcheck += 1;

			}
		}
	}
	*/
	
	/* Intra-unit cell hop */
	for (int mu = 0; mu < NU; mu++) {
		for (int nu = 0; nu < NU; nu++){
			//x = unitL[site1][0] - unitL[site0][0];
			//y = unitL[site1][1] - unitL[site0][1];
			dist = sqrt(0.5);
			data[numcheck].x = 0;
			data[numcheck].y = 0;
			data[numcheck].mu = mu;
			data[numcheck].nu = nu;
			data[numcheck].tmunu = Hxy[mu][nu];
			data[numcheck].dist = dist;
			numcheck += 1;
		}
	}


	/* Nearest neighbor hop */
	int site0, site1;

	/* Loop over six nearest neighbor */
	for (int i = 0; i < 6; i++){

		/* For the orbitals at the origin  unit*/
		for (int mu = 0; mu < NU; mu++) {

			/* For the orbitals at the n.n. unit */
			for (int nu = 0; nu < NU; nu++){
				site0 = mu/NC;
				site1 = nu/NC;
				x = superL[i][0] + unitL[site1][0] - unitL[site0][0];
				y = superL[i][1] + unitL[site1][1] - unitL[site0][1];
				dist = sqrt(x*x + y*y);
				if(dist > 1.0 ) continue;
				data[numcheck].x = superL[i][0];
				data[numcheck].y = superL[i][1];
				data[numcheck].mu = mu;
				data[numcheck].nu = nu;
				data[numcheck].tmunu = 0.0;
				/* Three cases */
				if (superL[i][0] == 0 && superL[i][1] != 0){
					/* Hopping along y axis : only O1->O1 and O1->O2 allowed */
					data[numcheck].tmunu = Hyy[mu][nu];
				}
				else if (superL[i][0] != 0 && superL[i][1] == 0 ){
					/* Hopping along x axis : only O2->O2 and O2->O1 allowed */
					data[numcheck].tmunu = Hxx[mu][nu];
				}
				else if (superL[i][0] != 0 && superL[i][1] != 0 ){
					data[numcheck].tmunu = Hxy[mu][nu];
				}
				data[numcheck].dist = dist;
				numcheck += 1;
			}
		}
	}

	


	return numcheck;

}


int PrintLattice(Cell *data, int num_check){
	for(int i = 0; i < num_check; i++){
		//if (data[i].mu/NC == data[i].nu/NC){
			printf("%d\t%d\t%d\t%d\t%d\t%f\t%f\n",i,data[i].x,data[i].y,data[i].mu, data[i].nu, data[i].tmunu,data[i].dist);
		//}
	}
	return 0;
}

int initHop(double ts, double tp){

	/* Initialize the hopping matrices 
	
	   		       y
	   		       ^
	   		       |
   	O---Cu---O---Cu---O---Cu-->x
	     |	      |        |
	     O	      O        O
	     |	      |        |     
   	O---Cu---O---Cu---O---Cu--

	Orbital ordering :
	|0> = |px up>, |1> = |px dn>, |2> = |py up>, |3> = |py dn>, |4> = |pz up>, |5> = |pz dn>

	*/
	double tmpxy[6][6] = {{0,0,ts,0,0,0},
			   {0,0,0,ts,0,0},	
			   {ts,0,0,0,0,0},
			   {0,ts,0,0,0,0},
			   {0,0,0,0,tp,0},
			   {0,0,0,0,0,tp}
			};
	double tmpxx[6][6] ={{ts,0,0,0,0,0},
			   {0,ts,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,0,0},
			   {0,0,0,0,tp,0},
			   {0,0,0,0,0,tp}
			};
        double tmpyy[6][6] ={{0,0,0,0,0,0},
                           {0,0,0,0,0,0},
                           {0,0,ts,0,0,0},
                           {0,0,0,ts,0,0},
                           {0,0,0,0,tp,0},
                           {0,0,0,0,0,tp}
                        };
	
	/* Onsite energies of two oxygens */
	double eps[2] = {0,0};


	for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) Hintra[i][j] = 0.0;
	for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) Hxy[i][j] = 0.0;
	for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) Hxx[i][j] = 0.0;
	for (int i = 0; i < NU; i++) for (int j = 0; j < NU; j++) Hyy[i][j] = 0.0;
	

	/* Define O-Cu-O x-y (intra unit cell) hopping */
	for (int i = 0; i < NU/2; i++) for (int j = NU/2; j < NU; j++){
		Hintra[i][j] = tmpxy[i][j%(NU/2)];
	}
	for (int i = NU/2; i < NU; i++) for (int j = 0; j < NU/2; j++){
		Hintra[i][j] = tmpxy[i%(NU/2)][j];
	}
	for (int i = 0; i < NU; i ++){
		Hintra[i][i] = eps[i/NC];
	}


	/* Define x-x unit cell hopping */
        for (int i = 0; i < NC; i++) for (int j = NC; j < NU; j++){
                Hxx[i][j] = tmpxy[i][j%NC];
        }
        for (int i = NC; i < NU; i++) for (int j = 0; j < NC; j++){
                Hxx[i][j] = tmpxy[i%NC][j];
        }
        for (int i = 0; i < NC; i++) for (int j = 0; j < NC; j++){
                Hxx[i][j] = tmpxx[i%NC][j];
        }

	/* Define y-y unit cell hopping */
        for (int i = 0; i < NC; i++) for (int j = NC; j < NU; j++){
                Hyy[i][j] = tmpxy[i][j%NC];
        }
        for (int i = NC; i < NU; i++) for (int j = 0; j < NC; j++){
                Hyy[i][j] = tmpxy[i%NC][j];
        }
	for (int i = NC; i < NU; i++) for (int j = NC; j < NU; j++){
		Hyy[i][j] = tmpyy[i%NC][j%NC];
	}

	/* Define x-y unit cell hopping */
	for (int i = 0; i < NU/2; i++) for (int j = NU/2; j < NU; j++){
		Hxy[i][j] = tmpxy[i][j%(NU/2)];
	}
	for (int i = NU/2; i < NU; i++) for (int j = 0; j < NU/2; j++){
		Hxy[i][j] = tmpxy[i%(NU/2)][j];
	}

	return 0;
}



int initQuadrature(Quad *quadrature){
       
	/* Generate Gaussian-Legendre quadrature to integrate over the discretized k-space 
	 *
	 * position: stores Nint(x/y) points between [-PI,PI ]
	 *
	 * weight: stores corresponding Gaussian-Legendre weight
	 *
	 * */
	int Dim = 2;
        quadrature->Nint[0] = Nintx;
        quadrature->Nint[1] = Ninty;
        //quadrature->Nint[2] = Nintx;

        for(int k=0; k<Dim; k++){
                quadrature->position[k] = mkvectord( quadrature->Nint[k] );
                quadrature->weight[k] = mkvectord( quadrature->Nint[k] );
                gauleg( -PI, PI, (quadrature->position[k])-1, (quadrature->weight[k])-1, quadrature->Nint[k] );
                //gauleg( -PI, PI, (quadrature->position[k]), (quadrature->weight[k]), quadrature->Nint[k] );
                printf("axis%d: %d points generated\n", k, quadrature->Nint[k]);
        }
        return 0;
}


int compute_hopmatrix(cx_mat& hopmatrix, double kx, double ky, Cell *data, int lines){
        
	
        /* Computes the fourier transform of the Hopping matrix H_0(k) */
	
	cx_double ex;	
	hopmatrix.zeros();


	int mu, nu;
        for(int i=0; i<lines; i++) {
                mu = data[i].mu;
                nu = data[i].nu;
                ex = cx_double( cos( kx*data[i].x + ky*data[i].y ), -sin( kx*data[i].x + ky*data[i].y ) );
		hopmatrix(mu,nu) +=  cx_double( data[i].tmunu  * real(ex) , data[i].tmunu * imag(ex) );
        }
        return hopmatrix.is_hermitian();

}
double FindFermi(Cell *data, Quad *quadrature, int lines, double U, double J){
        cx_mat hamil(NU,NU), H_U(NU,NU);
        double Efermi = 0;
        int hermitian, idx, Nk, Nbin = 1<<24;
	Nk = Nintx*Ninty; // No of k-points
	vec etmp(NU); // Save all eigenvalues
	//mat eval(Nk,NU);
	cube eval(Nintx, Ninty, NU);
	cx_mat evec;
	vec dos(Nbin);//density of states
	dos.zeros();
	eval.zeros();

	double emin = 1e6, emax = -1e6;
	//double emax = 1e6;
	DecoupleHF(H_U, U, J);


        for (int kx = 0; kx < Nintx; kx++) {
                for (int ky = 0; ky < Ninty; ky++){
			//idx = kx*quadrature->Nint[1]*quadrature->Nint[2] + ky*quadrature->Nint[2] + kz;
			idx = kx*quadrature->Nint[1] + ky;
			hamil.zeros();
                        hermitian = compute_hopmatrix(hamil,quadrature->position[0][kx], quadrature->position[1][ky],data,lines);
			if (hermitian != 1) {
				cout<<"hopmatrix is non-hermitian!"<<endl;
				exit(1);
			}
			hamil += H_U;
			hamil += H_mag;
			//H_mag.print("Hag");
			//hamil.print("hamil");
			eig_sym(etmp, evec, hamil);
			//for (int i = 0; i < NU; i++) eval(idx,i) = etmp(i);
			for (int i = 0; i < NU; i++) eval(kx,ky,i) = etmp(i);
		
			//cout<<"---"<<endl;
			//cout<<etmp[0]<<endl;	
			//cout<<etmp[NU-1]<<endl;	
			if (etmp(0) < emin ) emin = etmp(0);
			if (etmp(NU-1) > emax) emax = etmp(NU-1);
			//emax = etmp[0];
			//eig_gen(eval, evec, hamil);
			//eig_sym(eval, evec, hamil);
			//etmp.print("eval");
                }
        }
	emax += 1e-6;
	emin -= 1e-6;
	double dE = (emax-emin)/Nbin;
	printf("------------------<find fermi energy>::------------------\nFermi level in (%lf:%lf), dE = %f with %d grids\n", emin, emax, dE, Nbin);
	//cout<<emax<<endl;
	//cout<<emin<<endl;
	//cout<<dE<<endl;
	//for (int k = 0; k < Nk; k++){
	//	for(int mu = 0; mu < NU; mu++ ){
	//		idx = (int)(eval(k,mu) - emin)/dE;
	//		dos(idx) += 1;
			//cout<<idx<<endl;
	//	}
	//}

	for (int kx = 0; kx < Nintx; kx++){
		for (int ky = 0; ky < Ninty; ky++){
			for(int mu = 0; mu < NU; mu++ ){
				idx = (int)(eval(kx,ky,mu) - emin)/dE;
				dos(idx) += 1;
				//cout<<idx<<endl;
			}
		}
	}
	double occupancy = 0;//, norm = 1.0/(NU*NK);
	//cout<<"sum"<<endl;
	//cout<<sum(dos)<<endl;
	//dos *= (1.0/(NU*Nk));
	//dos /= (NU*Nk);
	idx = 0;
	//cout<<rho*Nk*NU<<endl;
	//cout<<"occupancy"<<endl;
	//double pop;
	for (int i = 0; i < Nbin; i++){
		occupancy += dos(i);
		idx = i;
		if(occupancy > rho*Nk*NU) break;
		//idx = i;
		//pop = occupancy;
		//occupancy += dos(i);
		//idx = i;
		//if(occupancy > rho*Nk*NU) break;
	}
	//while (occupancy < rho*Nk*NU){
	//	cout<<"----"<<endl;
	//	cout<<idx<<endl;
	//	cout<<dos(idx)<<endl;
	//	occupancy += dos(idx);
	//	cout<<occupancy<<endl;
	//	idx += 1;
	//}
	//printf("Occupancy:%f, %f\n",dos(idx), occupancy);
	//cout<<"idx"<<endl;
	//cout<<idx<<endl;
	//cout<<"Nbin"<<endl;
	//cout<<Nbin<<endl;i
	Efermi = emin + idx*dE;
	//cout<<"Efermi"<<endl;
	//cout<<( emin + idx*dE )<<endl;
	//dos.print();
	return Efermi;
}

double HartreeFock(Cell *data, Quad *quadrature, int lines, double U, double J){


	cx_mat hamil(NU,NU), H_U(NU,NU);
	vec eval(NU);
	cx_mat evec(NU,NU),my_rho, psi;
	double weight, energy = 0;
	int hermitian, ifermi;
	double Efermi;
	//double deltaMu;
	DecoupleHF(H_U, U, J);
        Efermi = FindFermi(data, quadrature, lines, U, J);
	//DecoupleHF(H_U, U, J);
	//H_U.print("HU");
	//cout<<"Fermi Energy"<<endl;
	//cout<<Efermi<<endl;
	printf("Fermi Energy: %f\n",Efermi);
	my_rho.zeros(NU,NU);
	for (int kx = 0; kx < Nintx; kx++) {
		for (int ky = 0; ky < Ninty; ky++){
			weight = quadrature->weight[0][kx] * quadrature->weight[1][ky]/Factor ;
			hermitian = compute_hopmatrix(hamil,quadrature->position[0][kx], quadrature->position[1][ky],data,lines);
			if (hermitian != 1) {
				cout<<"hopmatrix is non-hermitian!"<<endl;
				exit(1);
			}
			/* Add electron correlation */
			hamil += H_U;
			hamil += H_mag;
			eig_sym(eval,evec,hamil);
			for (int i = 0; i < NU; i++) Eigval(kx,ky,i) = eval(i)-Efermi;
			//eval.print();
			for (int i = 0; i < NU; i++){
				ifermi = i;
				//energy += eval(i)*weight;
				if(eval(i)>=Efermi) break;
				//ifermi = i;
				energy += eval(i)*weight;
			}
			//for (int i = 0; i < NU; i++) Eigval(kx,ky,i) = eval(i);
			//psi.zeros(ifermi,NU);	
			//psi.zeros(NU,ifermi-1);	
			//cout<<"idx"<<endl;
			//deltaMu = 10 - ifermi + 1;
			//cout<<"delta:"<<deltaMu<<endl;
			//ifermi = 10;
			//cout<<ifermi<<endl;
			// compute density matrix
			//evec.print("evec");
			
			//psi = evec.submat(0,0,ifermi,NU-1);
			psi = evec.submat(0,0,NU-1,ifermi);
			//cout<<size(psi)<<"\t"<<ifermi<<endl;
			
			//my_rho += psi.t()*psi*weight;
			//my_rho += psi*conj(psi.t())*weight;
			my_rho += psi*psi.t()*weight;
			
			//cout<<size(my_rho)<<"\t"<<ifermi<<endl;
			//my_rho += evec.submat(0,0,NU-1,ifermi) * evec.submat(0,0,NU-1,ifermi).t() * weight;
			//cout<<my_rho.n_rows<<endl;
			//cout<<my_rho.n_cols<<endl;
			//my_rho.print("myrho");	
			//cout<<hermitian<<endl;
		}
	}
	DENSITY_MATRIX = my_rho;
	printf("Trace: %f\n",real(trace(DENSITY_MATRIX)));
	
	return energy;
}


double HartreeFockParallel(Cell *data, Quad *quadrature, int lines, double U, double J){

	/* Declare OMP variable */
	int omp_num = omp_get_max_threads();
	//omp_set_num_threads(1);
	//int omp_num = omp_get_num_threads();
	printf("Number of threads: %d\n",omp_num);

	//exit(1);
	//cx_mat hamil(NU,NU), H_U(NU,NU);
	cx_mat total_rho;
	total_rho.zeros(NU,NU);
	///cx_cube hamil(NU,NU,omp_num);
	//vec eval(NU);
	//mat eval(omp_num, NU);
	//cx_mat evec(NU,NU),my_rho, psi;
	double Efermi, energy;
	//double deltaMu;
	DecoupleHF(H_U, U, J);
        Efermi = FindFermi(data, quadrature, lines, U, J);
	printf("Fermi Energy: %f\n",Efermi);
	//total_rho.print("rho");
	//DecoupleHF(H_U, U, J);
	//H_U.print("HU");
	//cout<<"Fermi Energy"<<endl;
	//cout<<Efermi<<endl;
	//
	//
//#pragma omp parallel default(none) shared ( H_U, H_mag, Efermi, quadrature, data, lines, Eigval, DENSITY_MATRIX, energy, total_rho, hamil)
#pragma omp parallel default(none) shared ( H_U, H_mag, Efermi, quadrature, data, lines, Eigval, DENSITY_MATRIX, energy, total_rho)
	{
		//int mythread = omp_get_thread_num();
		cx_mat my_rho, psi, evec(NU,NU);
		cx_mat hamil(NU,NU);
		vec eval(NU);
		double weight, my_energy = 0;
		int hermitian, ifermi;
		my_rho.zeros(NU,NU);
#pragma omp for collapse(2)
//#pragma omp for
		for (int kx = 0; kx < Nintx; kx++) for (int ky = 0; ky < Nintx; ky++){
		//for (int kx = 0; kx < Nintx; kx++) for (int ky = 0; ky < Nintx; ky++)printf("%d,%d, thread: %d\n",kx,ky,mythread);
			//printf("%d,%d, thread: %d\n",kx,ky,mythread);
			//hamil.zeros(NU,NU);
			weight = quadrature->weight[0][kx] * quadrature->weight[1][ky]/Factor ;
			//weight = quadrature->weight[0][kx] * quadrature->weight[1][ky];
			//hamil.slice(mythread).zeros();
			//hermitian = compute_hopmatrix(hamil.slice(mythread),quadrature->position[0][kx], quadrature->position[1][ky],data,lines);
			hermitian = compute_hopmatrix(hamil,quadrature->position[0][kx], quadrature->position[1][ky],data,lines);
			//hermitian = compute_hopmatrix(hamil,quadrature->position[0][kx], quadrature->position[1][ky],data,lines);
			if (hermitian != 1) {
				//cout<<"hopmatrix is non-hermitian!"<<endl;
				printf("hopmatrix is non-hermitian! \n");
				exit(1);
			}

			//hamil.slice(mythread) += H_U;
			//hamil.slice(mythread) += H_mag;
			hamil += H_U;
			hamil += H_mag;
			//printf("Diagonalizing..\n");
			//eig_sym(eval,evec);
			
			eig_sym(eval,evec,hamil);
			//eig_sym(eval,evec,hamil.slice(mythread));
			//printf("Done..\n");

			for (int i = 0; i < NU; i++) Eigval(kx,ky,i) = eval(i)-Efermi;
			//eval.print();
			for (int i = 0; i < NU; i++){
				ifermi = i;
				//energy += eval(i)*weight;
				if(eval(i)>=Efermi) break;
				//ifermi = i;
				my_energy += eval(i)*weight;
			}

			psi = evec.submat(0,0,NU-1,ifermi);
			//psi = evec.submat(0,0,ifermi-1,NU-1);
			//cout<<size(psi)<<"\t"<<ifermi<<endl;

			//my_rho += conj(psi.t())*psi*weight;
			my_rho = (psi)*(psi.t())*weight;
#pragma omp critical 
		{
			energy += my_energy;
			//DENSITY_MATRIX += my_rho;
			total_rho += my_rho;
		}




		}
	}
	//total_rho.print("new rho");
	DENSITY_MATRIX = total_rho;

	printf("Trace: %f\n",real(trace(DENSITY_MATRIX)));
	
	return energy;
	//return 0;
}


int DecoupleHF( cx_mat& H_U, double U, double J){

	//int i,j;
	//gsl_complex **HU = mkgscmatrixd(NU,NU);
	cx_mat HU(NU,NU);
	int map[6] = {4,5,0,1,2,3};

	double U_dd[15][3] = {
		{1, 2, U},
		{1, 3, U - 3*J},
		{1, 4, U - 2*J},
		{1, 5, U-3*J},
		{1, 6, U-2*J},
		{2, 3, U - 2*J},
		{2, 4, U - 3*J},
		{2, 5, U-2*J},
		{2, 6, U-3*J},
		{3, 4, U},
		{3, 5, U-3*J},
		{3, 6, U-2*J},
		{4, 5, U -2*J},
		{4, 6, U -3*J},
		{5, 6, U},
	};

	double U_hop[12][5] = {
		{1, 2, 4, 3, J},
		{1, 2, 6, 5, J},
		{1, 4, 2, 3,  J},
		{1, 6, 5, 2, -J},
		{3, 4, 6, 5, J},
		{3, 6, 5, 4, -J},
		{2, 3, 4, 1, -J},
		{2, 5, 6, 1, -J},
		{3, 4, 2, 1, J},
		{4, 5, 6, 3, -J},
		{5, 6, 2, 1, J},
		{5, 6, 4, 3, J},
	};

	HU.zeros();

	//for(i = 0; i < NU; i++) for(j = 0; j < NU; j++) HU[i][j] = zero;
	cx_double tmp;
	//int i1, i2, i3, i4, site;
	int i1, i2, i3, i4;

	/* Rotate the DENSITY MATRIX to the local octahedral coordinates */
	//iunitary(DENSITY_MATRIX,target,RUtotal,NU);
	
	for (int site = 0; site < 2; site++){
		for (int i = 0; i < 15; i++){
			i1 = 6*site + map[ (int)(U_dd[i][0] -1) ];
			i2 = 6*site + map[ (int)(U_dd[i][1] -1) ];
			tmp = DENSITY_MATRIX(i2,i2)*U_dd[i][2];
			HU(i1,i1) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i2][i2], U_dd[i][2]);
			//HU[i1][i1] =  gsl_complex_add(HU[i1][i1], tmp);

			
			tmp = DENSITY_MATRIX(i1,i1)*U_dd[i][2];
			HU(i2,i2) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i1][i1], U_dd[i][2]);
			//HU[i2][i2] =  gsl_complex_add(HU[i2][i2], tmp);

			tmp = DENSITY_MATRIX(i2,i1)*(-U_dd[i][2]);
			HU(i1,i2) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i2][i1], -U_dd[i][2]);
			//HU[i1][i2] =  gsl_complex_add(HU[i1][i2], tmp);

			tmp = DENSITY_MATRIX(i1,i2)*(-U_dd[i][2]);
			HU(i2,i1) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i1][i2], -U_dd[i][2]);
			//HU[i2][i1] = gsl_complex_add(HU[i2][i1], tmp);
		}
	}

	/* add pair-hopping contribution */
	for (int site = 0; site < 2; site++){
		for (int i = 0; i < 12; i++){
			i1 = 6*site + map[ (int)(U_hop[i][0] -1) ];
			i2 = 6*site + map[ (int)(U_hop[i][1] -1) ];
			i3 = 6*site + map[ (int)(U_hop[i][2] -1) ];
			i4 = 6*site + map[ (int)(U_hop[i][3] -1) ];


			/* positive terms */
			tmp = DENSITY_MATRIX(i2,i3)*U_hop[i][4];
			HU(i1,i4) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i2][i3], U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i2][i3], U_hop[i][4]);
			//HU[i1][i4] =  gsl_complex_add(HU[i1][i4], tmp);


			tmp = DENSITY_MATRIX(i1,i4)*U_hop[i][4];
			HU(i2,i3) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i1][i4], U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i1][i4], U_hop[i][4]);
			//HU[i2][i3] =  gsl_complex_add(HU[i2][i3], tmp);

			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i4][i1], U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i4][i1], U_hop[i][4]);
			//HU[i3][i2] =  gsl_complex_add(HU[i3][i2], tmp);


			/* negative terms */
			tmp = DENSITY_MATRIX(i2,i4)*(-U_hop[i][4]);
			HU(i1,i3) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i2][i4], -U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i2][i4], -U_hop[i][4]);
			//HU[i1][i3] =  gsl_complex_add(HU[i1][i3], tmp);

			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i4][i2], -U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i4][i2], -U_hop[i][4]);
			//HU[i3][i1] =  gsl_complex_add(HU[i3][i1], tmp);

			tmp = DENSITY_MATRIX(i1,i3)*(-U_hop[i][4]);
			HU(i2,i4) += tmp;
			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i1][i3], -U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i1][i3], -U_hop[i][4]);
			//HU[i2][i4] = gsl_complex_add(HU[i2][i4], tmp);

			//tmp = gsl_complex_mul_real(DENSITY_MATRIX[i3][i1], -U_hop[i][4]);
			//tmp = gsl_complex_mul_real(target[i3][i1], -U_hop[i][4]);
			//HU[i4][i2] = gsl_complex_add(HU[i4][i2], tmp);
		}
	}

	H_U = HU;
	//for (i=0; i<NU; i++) for(j=0;j<NU;j++) H_U[i][j] = HU[i][j];
	//freegscmatrixd(HU,NU);
	//freegscmatrixd(target,NU);
	return 0;
}	

double errorDensityMatrix(){

	/* Computes the maximum error from density matrix between two consecutive iterations */
	cx_mat tmp(NU,NU);
	mat rtmp(NU,NU), ctmp;
	//tmp = (DENSITY_MATRIX - DENSITY_MATRIX_tmp).t()*(DENSITY_MATRIX - DENSITY_MATRIX_tmp  );
	tmp = (DENSITY_MATRIX - DENSITY_MATRIX_tmp);
	rtmp = real(tmp);
	ctmp = imag(tmp);

	//rtmp.print("rtmp");
	//ctmp.print("ctmp");

	double e1, e2;
	e1 = max(max(rtmp));
	e2 = max(max(ctmp));
	//double error = findmax(max(max(rtmp)),max(max(ctmp))) ;
	double error = findmax(e1,e2) ;
	DENSITY_MATRIX_tmp = DENSITY_MATRIX;
	return error;
}

int SpinorGenerator(double J, cx_cube& G){
	
	/* Computes the Generator of spin-J rotation 
	 *
	 * Note: the ordering of the spin basis
	 * 	 considered here is {m,..,-m}
	 * 	
	 *
	 *******************************************/ 

	double term1,term2,term3, m1,m2;
	int span = 2*J + 1;
	//G = mkgsctritensord(3,span,span); // Generator
	G.zeros(span,span,3);
	/* Perpare the generator */
	m1 = J;
	for(int i = 0; i < span; i++){
		m2 = J;
		for(int j = 0; j < span; j++){
			term1 = 0.5*sqrt(J*(J+1) - m2*(m2+1))*delta(m1,m2+1) + 0.5*sqrt(J*(J+1) - m2*(m2-1))*delta(m1,m2-1);
			term2 = 0.5*sqrt(J*(J+1) - m2*(m2+1))*delta(m1,m2+1) - 0.5*sqrt(J*(J+1) - m2*(m2-1))*delta(m1,m2-1);
			term3 = m2*delta(m1,m2);
			G(i,j,0) = cx_double(term1,0);
			G(i,j,1) = cx_double(0,-term2);
			G(i,j,2) = cx_double(term3,0);
			m2 -= 1;
		}
		m1 -= 1;
	}

	return 0;
}



int makeSpinOperator(){

	/* Creates the spin operator matrices from Generators */
	//int i,j,k;
	//GU = mkgsctritensord(Dim,6,6);
	GU.zeros(NC,NC,3);
	cx_cube G2(2,2,3);
	//gsl_complex ***G2;
	//G2 = mkgsctritensord(Dim,2,2);


	SpinorGenerator(0.5,G2);

	//for(k = 0; k<3; k++)for(i = 0; i<NC; i++)for(j = 0; j<6; j++) GU[k][i][j] = zero;
	
	for(int k = 0; k<3; k++){
		for(int i = 0; i < 2; i++){
			for(int j = 0; j < 2; j++){
				//GU[0][i+2*k][j+2*k] = G2[0][i][j];
				GU(i+2*k,j+2*k,0) = G2(i,j,0);
				GU(i+2*k,j+2*k,1) = G2(i,j,1);
				GU(i+2*k,j+2*k,2) = G2(i,j,2);
				//GU[1][i+2*k][j+2*k] = G2[1][i][j];
				//GU[2][i+2*k][j+2*k] = G2[2][i][j];
			}
		}
	}



	//freegsctritensord(G2,Dim,2);
	return 0;

}

int ComputeMagMoment(mat& moment){
	/* Computes the magnetic moment */
	moment.zeros();
	cx_mat tmp;
	
	for (int site = 0; site < 2; site++){
		for (int x = 0; x < 3; x++){
			//cout << size( DENSITY_MATRIX.submat(NC*site,NC*site,NC*(site + 1)-1, NC*(site + 1)-1) ) << endl;
			tmp = DENSITY_MATRIX.submat(NC*site,NC*site,NC*(site + 1)-1, NC*(site + 1)-1).t() ;
			tmp *= GU.slice(x);
			//(tmp*GU.slice(x)).print();
			//tmp.print();
			//cout<<trace(tmp)<<endl;
			moment(site,x) =  real(trace(tmp));
		}
	}
	return 0;

}

int ComputeCharge(vec& charge){

	/* Computes charge at each site */
	charge.zeros();	
	cx_mat charge_operator(NC,NC);
	charge_operator.eye();
	cx_mat tmp;
	for (int site = 0; site < 2; site++){
			//cout << size( DENSITY_MATRIX.submat(NC*site,NC*site,NC*(site + 1)-1, NC*(site + 1)-1) ) << endl;
			tmp = DENSITY_MATRIX.submat(NC*site,NC*site,NC*(site + 1)-1, NC*(site + 1)-1) ;
			tmp *= charge_operator;
			charge(site) =  real(trace(tmp));
	}
	return 0;

}

double maxDiff(mat& moment1, mat& moment2){
	/* Computes the max difference between two moments */
	mat diff(2,3);
	double error;
	diff = abs(moment2 - moment1);
	error = max(max(diff));
	return error;
}

void printMoment(mat& moment){
	/*Prints the magnetic moment at each site */
	for (int i = 0; i < 2; i++){
		printf("Site: [%d]; total moment = %lf\n",i, norm(moment.row(i)));
	}
	for (int i = 0; i < 2; i++){
		moment.row(i).print();
	}

}

void saveMoment(mat& moment, vec& charge, double U, double Hz){
	/*Saves the magnetic moment at each site */
	char paramoment[1024];
	FILE *fmoment;
        sprintf(paramoment,"moment_H%.2f_%d%d%d_U%.2f.txt",Hz,int(Haxis[0]),int(Haxis[1]),int(Haxis[2]), U );
        fmoment = fopen(paramoment,"w");
	for (int i = 0; i < 2; i++){
		printf("Site: [%d]; total moment = %lf\n",i, norm(moment.row(i)));
		fprintf(fmoment,"Site: [%d]; total moment = %lf\n",i, norm(moment.row(i)));
	}
	for (int i = 0; i < 2; i++){
		moment.row(i).print();
		for (int x = 0; x < 3; x++) fprintf(fmoment,"%lf\t",moment(i,x));
		fprintf(fmoment,"\n");
	}
	for (int i = 0; i < 2; i++){
		printf("Site: [%d]; total charge = %lf\n",i, charge(i) );
		fprintf(fmoment,"Site: [%d]; charge = %lf\n",i, charge(i));
	}

	fclose(fmoment);

}
double saveBand(Cell *data, Quad *quadrature, int lines, double U, double J, double Hz){
	/* Save band structure along high symmetry lines 
	 *
	 *
	 * NOTE: High symmetry points in the HIYM file below is given in cartesian coordinates
	 *
	 * */
	int Naxis, Dim = 2;
	mat sympoint;
	char *indices = NULL;
	char paraband[1024];

	FILE *fb = fopen("inputs/HISYM", "r");
	
	if( fb != NULL ){
		if( fscanf(fb, "%d", &Naxis) == 0 ){
			printf("reading error\n");
			exit(1);
		}
		//sym_temp = mkmatrixd(Naxis, Dim);
		cout<<Naxis<<endl;
		sympoint.zeros(Naxis,Dim);
		indices = mkvectorc(Naxis);
		for(int i=0; i<Naxis; i++){
			for(int k=0; k<Dim; k++){
				//if( fscanf(fb, "%lf", &sym_temp[i][k]) == 1 )
				if( fscanf(fb, "%lf", &sympoint(i,k)) == 1 )
					//sym_temp[i][k] *= 2*PI;
					sympoint(i,k) *= 2*PI;
				//else exit(1);//	printf("reading error\n");
			}
			if( fscanf(fb, " %c", &indices[i]) == 1 )
				continue;
			//else exit(1);
		}
		fclose(fb);
	}
	else {
		printf("HISYM is required!\n");
		exit(1);
	}
	sympoint.print("Sympoint");

	vec dk(Dim); // Change in K between different high-symmetry points
	vec Nkdx(Naxis-1); // No of K points between two high-symmetry points
	Nkdx.fill(100);
	double kx,ky;

	cx_mat hamil(NU,NU), H_U(NU,NU), evec(NU,NU);
	vec eval(NU);
	double Efermi;
	int hermitian;
	
	DecoupleHF(H_U, U, J);
        Efermi = FindFermi(data, quadrature, lines, U, J);
	//DecoupleHF(H_U, U, J);
	int num = 0;
	FILE *fband;
	//ofstream fband;
	sprintf(paraband,"band_H%.2f_%d%d%d_U%.2f.txt",Hz,int(Haxis[0]),int(Haxis[1]),int(Haxis[2]), U);
	fband = fopen(paraband,"w");
	//fband.open(paraband,"w");
	//fband.open(paraband);
	//H_U.print("H_U");
	for(int axis=0; axis<Naxis-1; axis++) {
		for(int i=0; i<Dim; i++)	dk(i) = sympoint(axis+1,i)-sympoint(axis,i);
		//for(k=0; k<Nkdx[axis]; k++){
		for(int k=0; k<Nkdx(axis); k++){
			//kx = sympoint(axis,0) + k*dk(0)/Nkdx[axis];
			kx = sympoint(axis,0) + k*dk(0)/Nkdx(axis);
			ky = sympoint(axis,1) + k*dk(1)/Nkdx(axis);
			//printf("%f\t%f\n",kx,ky);
			//printf("%d\t",count);
			fprintf(fband,"%d\t",num);
			hermitian = compute_hopmatrix(hamil,kx, ky,data,lines);
			if (hermitian != 1) {
				cout<<"hopmatrix is non-hermitian!"<<endl;
				exit(1);
			}
			/* Add electron correlation */
			hamil += H_U;
			/* Add magnetic field contribution */
			//H_mag.print("Hmag");
			hamil += H_mag;
			eig_sym(eval,evec,hamil);
			//for (int mu = 0; mu < NU; mu++) printf("%lf\t",eval(mu)-Efermi);
			for (int mu = 0; mu < NU; mu++) fprintf(fband,"%lf\t",eval(mu)-Efermi);
			fprintf(fband,"\n");
			//ky = sympoint[axis][1] + k*dk[1]/Nkdx[axis];
			//kz = sympoint[axis][2] + k*dk[2]/Nkdx[axis];
			num += 1;
		}
	}
	fclose(fband);

	return Efermi;
	//fband.close();
	/*
	double x,y;
	for (int kx = 0; kx < Nintx; kx++){
	for (int ky = 0; ky < Ninty; ky++){
		x = quadrature->position[0][kx];
		y = quadrature->position[1][ky];
		printf("%f\t %f\t",x,y);
		for(int i = 0; i < NU; i++)printf("%f\t", Eigval(kx,ky,i));
		printf("\n");

	}
	}
	*/

}

void saveFullBand(Quad *quadrature, double U, double Hz){
	/* Save top two bands for entire B-zone */
	double x,y;
	//int Nband = 3;
	char paraband[1024];
	FILE *fband;
	sprintf(paraband,"Fullband_H%.2f_%d%d%d_U%.2f.txt",Hz,int(Haxis[0]),int(Haxis[1]),int(Haxis[2]),U);
	fband = fopen(paraband,"w");

	for (int kx = 0; kx < Nintx; kx++){
		for (int ky = 0; ky < Ninty; ky++){
			x = quadrature->position[0][kx];
			y = quadrature->position[1][ky];
			//printf("%f\t %f\t",x,y);
			fprintf(fband,"%f\t %f\t",x,y);
			for(int i = 0; i < NU; i++)fprintf(fband,"%f\t", Eigval(kx,ky,i)  );
			//fprintf(fband,"%f \t",  Eigval(kx,ky,NU-1)-Eigval(kx,ky,NU-2)  );
			//fprintf(fband,"%f \t %f\t %f", Eigval(kx,ky,NU-2), Eigval(kx,ky,NU-1), Eigval(kx,ky,NU-1)-Eigval(kx,ky,NU-2)  );
			//printf("\n");
			fprintf(fband,"\n");

		}
	}

	fclose(fband);
}

void saveBandGap(Quad *quadrature, double U, double Hz){
	/* Save band gap of the top two bands for entire B-zone */
	double x,y;
	//int Nband = 3;
	char paraband[1024];
	FILE *fband;
	sprintf(paraband,"Gapband_H%.2f_%d%d%d_U%.2f.txt",Hz,int(Haxis[0]),int(Haxis[1]),int(Haxis[2]),U);
	fband = fopen(paraband,"w");

	for (int kx = 0; kx < Nintx; kx++){
		for (int ky = 0; ky < Ninty; ky++){
			x = quadrature->position[0][kx];
			y = quadrature->position[1][ky];
			//printf("%f\t %f\t",x,y);
			fprintf(fband,"%f\t %f\t",x,y);
			//for(int i = 0; i < NU; i++)fprintf(fband,"%f\t", Eigval(kx,ky,i)  );
			fprintf(fband,"%f \t",  Eigval(kx,ky,NU-1)-Eigval(kx,ky,NU-3)  );// spin up
			fprintf(fband,"%f \t",  Eigval(kx,ky,NU-2)-Eigval(kx,ky,NU-4)  );// spin down
			fprintf(fband,"%f \t",  Eigval(kx,ky,NU-1)-Eigval(kx,ky,NU-2)  );// spin up - spin down
			//fprintf(fband,"%f \t %f\t %f", Eigval(kx,ky,NU-2), Eigval(kx,ky,NU-1), Eigval(kx,ky,NU-1)-Eigval(kx,ky,NU-2)  );
			//printf("\n");
			fprintf(fband,"\n");

		}
	}

	fclose(fband);
}
void saveArpes(Cell *data, Quad *quadrature, int lines, double U, double J, double Hz){
	/* plot energy gap to trak ARPES data */
	//FILE *fb = fopen("inputs/ARPES", "r");

	cx_mat hamil(NU,NU), H_U(NU,NU), evec(NU,NU);
	vec eval(NU);
	double Efermi;
	int hermitian;
	
	DecoupleHF(H_U, U, J);
        Efermi = FindFermi(data, quadrature, lines, U, J);

	double kx,ky, theta, frac = 0.8;
	char paraband[1024];
	//, buffer[1024];
	sprintf(paraband,"arpes_H%.2f_%d%d%d_U%.2lf.txt",Hz,int(Haxis[0]),int(Haxis[1]),int(Haxis[2]), U);
	ifstream fb("inputs/ARPES");
	//ifstream fb("inputs/contourE1.txt");
	//ifstream fb("inputs/contourEup1U0.00.txt");
	ofstream fo(paraband);
	while (fb >> kx >> ky){
		//cout<<kx<<"\t"<<ky<<endl;
		//printf("%f\t%f\n",kx,ky);
		theta = (180/PI)*atan( (PI - kx)/(PI-ky));
		//printf("%f\t",theta);
		//sprintf(buffer,"%f\t",theta);
		fo << theta <<"\t"<<kx<<"\t"<<ky<<"\t";
		hermitian = compute_hopmatrix(hamil,kx*frac, ky*frac,data,lines);
		if (hermitian != 1) {
			cout<<"hopmatrix is non-hermitian!"<<endl;
			exit(1);
		}
		/* Add electron correlation */
		hamil += H_U;
		/* Add magnetic field contribution */
		hamil += H_mag;
		eig_sym(eval,evec,hamil);

		//for (int mu = 0; mu < NU; mu++) printf("%lf\t",eval(mu)-Efermi);
		//fprintf(fband,"%f \t %f\t %f", eval(NU-2), eval(NU-1), eval(NU-1)-eval(NU-2)  );
		//printf("%f \t %f\t %f", eval(NU-2), eval(NU-1), eval(NU-1)-eval(NU-2)  );
		//sprintf(buffer,"%f \t %f\t ", eval(NU-1) - eval(NU-3), eval(NU-2) - eval(NU-4)  );
		//for (int mu = 0; mu < NU; mu++) sprintf(buffer,"%lf\t",eval(mu)-Efermi);
		for (int mu = 0; mu < NU; mu++) fo <<(eval(mu)-Efermi)<<"\t";
		
		//fo << eval(NU-2)<<"\t"<< eval(NU-1)<<"\t"<< eval(NU-1)-eval(NU-2)  <<endl;
		fo << "" << endl;
		//printf("\n");
		//for (int mu = 0; mu < NU; mu++) fprintf(fband,"%lf\t",eval(mu)-Efermi);
	}
	fb.close();
	fo.close();

}


int constructHmag(double Hz, double direction[3]){
	
	/* Computes the Magnetic field contribution 
		Hz : magnitude of the magnetic field
		direction: Direction of the magnetic field
	*/
	
	vec V(3); 
	//V.zeros(3);
	//double normV;	
	//double spinProj[NC] = {0.5,-0.5,0.5,-0.5,0.5,-0.5};
	
	for (int i = 0; i < 3; i++) {
		V(i) = direction[i];
	}
	
	V = normalise(V);
	//V /= norm;
	vec H(3);// = zeros(3);
	for (int i = 0; i < 3; i++){
		 H(i) = V(i)*Hz;
	}
	
	H.print("H="); 
	
	// S*H
	//H_mag = zeros(NU,NU);
	cx_mat tmpmag(NC,NC);
	tmpmag.zeros();
	for (int i = 0; i < 3; i++){
		//tmp.slice(i) = GU.slice(i)*H(i);
		tmpmag += GU.slice(i)*H(i);
	}
	//tmpmag.print("tmpmag");
	H_mag.zeros(NU,NU);
	for (int site = 0; site < 2; site ++){
		for (int i = 0; i < NC; i++) for(int j = 0; j < NC; j++) H_mag(site*NC + i,site*NC + j) = tmpmag(i,j);
	}
	//H_mag.print("Hmag");
	//for (int i = 0; i < NU; i++) H_mag(i,i) = cx_double(spinProj[i],0);
	
	return 0;
}

