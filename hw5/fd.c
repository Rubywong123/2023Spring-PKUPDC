#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// TO run with tolerance 1e-4 and 4x4 loop currents
//    ./fd 4 1e-4

#define PI 3.14159265359
#define MAX(a,b) (((a)>(b))?(a):(b))

// solve for solution vector u
int solve(const int N, const double tol, double * u, double * f){

  double *unew = (double*)calloc((N+2)*(N+2),sizeof(double));
  
  double w = 1.0;
  double invD = 1./4.;  // factor of h cancels out

  double res2 = 1.0;
  unsigned int iter = 0;
  while(res2>tol*tol){

    res2 = 0.0;

    // update interior nodes using Jacobi
    for(int i=1; i<=N; ++i){
      for(int j=1; j<=N; ++j){
	
	const int id = i + j*(N+2); // x-index first
	const double Ru = -u[id-(N+2)]-u[id+(N+2)]-u[id-1]-u[id+1];
	const double rhs = invD*(f[id]-Ru);
	const double oldu = u[id];
	const double newu = w*rhs + (1.0-w)*oldu;
	res2 += (newu-oldu)*(newu-oldu);
	unew[id] = newu;
      }
    }

    for (int i = 0; i < (N+2)*(N+2); ++i){
      u[i] = unew[i];
    }
    
    ++iter;
    if(!(iter%500)){
      printf("at iter %d: residual = %g\n", iter, sqrt(res2));
    }
  }

  return iter;
}

int main(int argc, char **argv){
  
  if(argc!=3){
    printf("Usage: ./main N tol\n");
    exit(-1);
  }
  
  int N = atoi(argv[1]);
  double tol = atof(argv[2]);

  double *u = (double*) calloc((N+2)*(N+2), sizeof(double));
  double *f = (double*) calloc((N+2)*(N+2), sizeof(double));
  double h = 2.0/(N+1);
  for (int i = 0; i < N+2; ++i){
    for (int j = 0; j < N+2; ++j){
      const double x = -1.0 + i*h;
      const double y = -1.0 + j*h;
      f[i + j*(N+2)] = sin(PI*x)*sin(PI*y) * h*h;
    }
  }

   
  int iter = solve(N, tol, u, f);

  double err = 0.0;
  for (int i = 0; i < (N+2)*(N+2); ++i){
    err = MAX(err,fabs(u[i] - f[i]/(h*h*2.0*PI*PI)));
  }
  
  printf("Iters: %d\n", iter);
  printf("Max error: %lg\n", err);
  printf("Memory usage: %lg GB\n", (N+2)*(N+2)*sizeof(double)/1.e9);
  
  free(u);
  free(f);  

}
  
