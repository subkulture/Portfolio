// C++ / CUDA Program to perform matrix maths using GPUs via CUDA


// Includes
#include <iostream>  // cout
#include <iomanip>  // setprecision
#include <stdlib.h>  // atoi
using namespace std;

//--------------------------------------------------------------------
// CUDA Kernel function to add the elements of two arrays on the GPU
//--------------------------------------------------------------------
__global__ // all kernels are preceded by __global__ keyword
void add(int n, // number of elements in an array
         float *A, // device pointer to array
         float *B, // device pointer to another array
         float *C) // device pointer to another array
{
    // determine thread ID within block
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    // determine stride of loop (more elements in array than threads)
    int stride = blockDim.x * gridDim.x;

    // each thread loops over array elements in steps of grid size
    for (int i = index; i < n; i += stride){
        A[i] = B[i] + C[i];
    }
}

int main(int argc, char** argv)
{
    // Get N from command line argument
    int N = 1000000;
    if (argc > 1) N = atoi(argv[1]); // specify array size
    std::cout << "Running with N = " << N << std::endl;
    float *A, *B, *C; // pointers for arrays

    // Allocate Unified Memory – accessible from CPU or GPU
    cudaMallocManaged(&A, N*sizeof(float));
    cudaMallocManaged(&B, N*sizeof(float));
    cudaMallocManaged(&C, N*sizeof(float));

    // initialize arrays on the host
    for(int i = 0; i < N; i++) A[i] = 0.0f;
    for(int i = 0; i < N; i++) B[i] = 0.1f;
    for(int i = 0; i < N; i++) C[i] = 0.2f;

    // Check the arithmetic on the CPU.
    double a;
    float b, c;
    b = 0.1f;
    c = 0.2f;
    a = b + c;
    std::cout << "a = " << std::setprecision(16) << a << std::endl;
    std::cout << "delta = " << std::setprecision(16) << 0.3 - a << std::endl;

    int num_threads_in_block = 256; // set number of threads (multiple of 32)
    int num_blocks = 32*2; // for 2 SMs, set a multiple blocks for each one

    // Run kernel on N elements on the GPU
    add<<<num_blocks, num_threads_in_block>>>(N, A, B, C);

    // cudaDeviceSynchronize call added here as the CPU continues through the rest of the
    // program so the printing to std::cout and cudaFree calls will operate on the arrays before the
    // GPU has finished doing its work.
    cudaDeviceSynchronize();

    std::cout << "A[0] = " << std::setprecision(16) << double(A[0]) << std::endl;
    std::cout << "delta = " << std::setprecision(16) <<  0.3 - double(A[0]) << std::endl;

    cudaFree(A);
    cudaFree(B);
    cudaFree(C);
}
