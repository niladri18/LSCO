#include <complex>
#include <cmath>
#include <math.h>
#include <omp.h>

#define NU      12 
#define NC      6 
#define Nintx	128
#define Ninty	128
#define PI      3.14159265358979323846264338
#define Factor	(4*PI*PI)

using namespace std;

class Cell {
	public:
		int x; 
		int y;
		int mu;
		int nu;
		//complex<long double> tmunu;
		double tmunu;
		double dist;
};

class Quad{
	public:
		double *position[2];
		double *weight[2];
		int Nint[2];

};

/*
class Latt{
	public:
		Cell *data;
		int num = 4000;
		int ConstructLattice();
		int PrintLattice();
		data = (Cell*) malloc(sizeof(Cell)*num);
		int Construct();
		int PrintLattice(int num_check);

};
*/
