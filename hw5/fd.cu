#include <math.h>
#include <stdio.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdlib.h>
#include "CycleTimer.h"

// TO run with tolerance 1e-4 and 4x4 loop currents
//    ./fd 4 1e-4

#define PI 3.14159265359
#define MAX(a,b) (((a)>(b))?(a):(b))

//kernel function
__global__ void compute_unew(double* u, double* f, double* unew, int N, double w, double invD){
  int xid = blockIdx.x * blockDim.x + threadIdx.x;
  int yid = blockIdx.y * blockDim.y + threadIdx.y;
  if(xid < 1 || xid >N || yid < 1 || yid > N)return;
  int id = xid + (yid) * (N+2);
  
  const double Ru = -u[id-(N+2)]-u[id+(N+2)]-u[id-1]-u[id+1];
	const double rhs = invD*(f[id]-Ru);
	const double oldu = u[id];
  // similar to gradient descent?
	const double newu = w*rhs + (1.0-w)*oldu;
  unew[id] = newu;
  __syncthreads();
}

// solve for solution vector u
// host function
int solve(const int N, const double tol, double * u, double * f){

  double start = CycleTimer::currentSeconds();
  double *unew = (double*)calloc((N+2)*(N+2),sizeof(double));
  size_t size = (N+2)*(N+2) * sizeof(double);

  //device memory
  double* d_u = NULL;
  cudaMalloc((void**)&d_u, size);
  double* d_f = NULL;
  cudaMalloc((void**)&d_f, size);
  double* d_unew = NULL;
  cudaMalloc((void**)&d_unew, size);

  double malloc_time = CycleTimer::currentSeconds() - start;
  printf("Malloc time: %.4f\n", malloc_time);

  cudaMemcpy(d_f, f, size, cudaMemcpyHostToDevice);
  double res2 = 1.0;
  unsigned int iter = 0;
  double w = 1.0;
  double invD = 1./4.;  // factor of h cancels out
  while(res2>tol*tol){

    res2 = 0.0;

    // copy from host to device
    // start = CycleTimer::currentSeconds();
    cudaMemcpy(d_u, u, size, cudaMemcpyHostToDevice);
    // double copy_time = CycleTimer::currentSeconds() - start;
    // printf("Copy_time: %.4f\n", copy_time);

    //setup device blocks & threads, then call kernel function
    // start = CycleTimer::currentSeconds();
    dim3 threadPerBlock(16, 16);
    dim3 blocksPerGrid((N + threadPerBlock.x - 1) / threadPerBlock.x,( N + threadPerBlock.y - 1) / threadPerBlock.y);
    compute_unew <<< blocksPerGrid, threadPerBlock >>> (d_u, d_f, d_unew, N, w, invD);
    
    cudaGetLastError();
    // double compute_time = CycleTimer::currentSeconds() - start;
    // printf("Computation time: %.4f\n", compute_time);

    //copy from device to host
    cudaMemcpy(unew, d_unew, size, cudaMemcpyDeviceToHost);

    // start = CycleTimer::currentSeconds();
    for(int i=1; i<=N; ++i){
      for(int j=1; j<=N; ++j){
        int id = i + j * (N+2);
        double newu = unew[id];
        double oldu = u[id];
	      res2 += (newu-oldu)*(newu-oldu);
      }
    }
    // double res_time = CycleTimer::currentSeconds() - start;
    // printf("Compute res time: %.4f\n", res_time);

    // start = CycleTimer::currentSeconds();
    for (int i = 0; i < (N+2)*(N+2); ++i){
      u[i] = unew[i];
    }
    // double update_time =  CycleTimer::currentSeconds() - start;
    // printf("Update time: %.4f\n", update_time);

    ++iter;
    if(!(iter%500)){
      printf("at iter %d: residual = %g\n", iter, sqrt(res2));
    }
  }

  cudaFree(d_f);
  cudaFree(d_u);
  cudaFree(d_unew);

  // free(unew);

  return iter;
}

int main(int argc, char **argv){
  
  if(argc!=3){
    printf("Usage: ./main N tol\n");
    exit(-1);
  }
  
  int N = atoi(argv[1]);
  double tol = atof(argv[2]);

  // flatten array
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
  double err = 0.0;
  for (int i = 0; i < (N+2)*(N+2); ++i){
    // solution: 1/2PI^2 * sin(PI*x) * sin(PI*y)
    // if (fabs(u[i] - f[i]/(h*h*2.0*PI*PI)) > err){
    //   printf("%d, %.4f, %.4f\n", i, fabs(u[i] - f[i]/(h*h*2.0*PI*PI)), err);
    // }
    err = MAX(err,fabs(u[i] - f[i]/(h*h*2.0*PI*PI)));
  }
  printf("Max error: %lg\n", err);

  //u: random array at initial state
  double start = CycleTimer::currentSeconds();
  int iter = solve(N, tol, u, f);
  double solve_time = CycleTimer::currentSeconds() - start;
  printf("Solve time: %.4f\n", solve_time);

  err = 0.0;
  for (int i = 0; i < (N+2)*(N+2); ++i){
    // if (fabs(u[i] - f[i]/(h*h*2.0*PI*PI)) > err){
    //   printf("%d, %.4f, %.4f\n", i, fabs(u[i] - f[i]/(h*h*2.0*PI*PI)), err);
    // }
    // solution: 1/2PI^2 * sin(PI*x) * sin(PI*y)
    err = MAX(err,fabs(u[i] - f[i]/(h*h*2.0*PI*PI)));
  }
  
  printf("Iters: %d\n", iter);
  printf("Max error: %lg\n", err);
  printf("Memory usage: %lg GB\n", (N+2)*(N+2)*sizeof(double)/1.e9);
  
  free(u);
  free(f);  

}
  
