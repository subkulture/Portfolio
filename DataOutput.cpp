//------------------------------------------------------------------------------
//                       Data output using MPI-IO
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <vector>

// mpi header
#include <mpi.h>

// openmp header
#include <omp.h>

// simple wrapper class for rng
class rng{

   // std::random variables (internal to class)
   std::mt19937 mt; // mersenne twister
   std::uniform_real_distribution<double> dist;

public:

   // seed rng
   void seed(unsigned int random_seed){
      dist = std::uniform_real_distribution<double>(0.0,1.0); // uniform distribution [0:1)
      std::mt19937::result_type mt_seed = random_seed;
      mt.seed(mt_seed);
   }

   // generate a uniform random number between 0 and 1
   double grnd() {
      return dist(mt);
   }

};

//------------------------------------------------------------------------------
//  Program to estimate the value of PI using Monte Carlo
//
//  The code is parallelised with hybrid MPI and OpenMP allowing thread and
//  process level parallelism.
//
//------------------------------------------------------------------------------
int main(int argument_count, char *argument_values[]){

   // initialise mpi
   MPI_Init(&argument_count, &argument_values);

   // mpi variables
   int rank = 0; // my rank in all processors
   int num_processors = 1; // total processors in set
   bool root = false; // flag identifying root process

   // omp variables
   int thread_rank = 0; // my thread rank in all processors
   int num_threads = 1; // total processors in set

   // set number of omp threads from command line
   if(argument_count == 2 ){

      // set number of threads from argument (converting to stringstream)
      std::string sizestr = argument_values[1];
      std::stringstream sizess(sizestr.c_str());
      sizess >> num_threads;

      // check for valid number of threads, and print an error message if not
      if(num_threads <1 || num_threads>32){
         std::cerr << "Error number of threads specified on the command line must be between 1 and 32" << std::endl;
         exit(1); // exit program
      }

   }

   // Get number of processors and rank
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
   if(rank == 0) root = true; // set root process as rank 0

   // Set OpenMP runtime parameters
   omp_set_dynamic(false); // disable dynamic alteration of number of threads
   omp_set_num_threads(num_threads); // explicitly set a number of threads per mpi proc

   // Set parameters for calculating packing fraction
   double box_size = 1000.0; // size of box
   double radius   = 10.0;   // radius of circles

   // Initialise array of random number generators (one per thread)
   std::vector<rng> random_generators(num_threads);

   // Seed thread parallel random number generators
   for(int t = 0; t< num_threads; t++){
      // define seed (using standard compliant wraparound behaviour for unsigned int)
      unsigned int seed = 37852+(1+rank)*(num_processors+num_threads) + (1+t)*num_threads;
      random_generators[t].seed(seed);
   }

   uint64_t in_count = 0; // total number of darts landing inside the circle

   // calculate number of darts per mpi process per thread (use unsigned 64 bit ints to avoid overflow)
   const uint64_t num_darts = 10000000000;
   const uint64_t num_darts_per_thread = num_darts/(num_processors*num_threads);

   // start the timer
   double start_time = MPI_Wtime();

   // omp parallel region
   #pragma omp parallel
   {

      const int thread = omp_get_thread_num(); // thread ID
      uint64_t thread_in_count = 0; // number of darts landing inside the circle (thread private)

      // for loop, but omp pragma not needed as using separate number of iterations for each
      for( int d=0; d<num_darts_per_thread; ++d){

         // generate random x and y positions
         const double x = random_generators[thread].grnd();
         const double y = random_generators[thread].grnd();

         // calculate radius squared
         const double radius_sq = x*x + y*y;

         // check if coordinates are inside circle of radius 1 and increment counter
         if(radius_sq <= 1.0) thread_in_count++;

      }

      // now add all counts together
      #pragma omp critical
      in_count += thread_in_count;

   } // end omp parallel

   // reduce the counts from all mpi processors on root
   uint64_t total_count = 0;
   MPI_Reduce(&in_count, &total_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

   // stop the timer
   double end_time = MPI_Wtime();

   // calculate actual number of darts thrown (to account for rounding errors)
   uint64_t total_darts = 0;
   MPI_Reduce(&num_darts_per_thread, &total_darts, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

   // print calculated value of pi and error to screen
   const double calculated_pi = 4.0*double(total_count)/double(total_darts*num_threads);
   if(root){
      std::cout << std::setprecision(15) << "Estimated value of pi is " << calculated_pi << std::setprecision(6) << " with a fractional error of " << fabs(calculated_pi-M_PI)*100.0/M_PI << " %" << std::endl;
      std::cout << "   MPI processes: " << num_processors << std::endl;
      std::cout << "   OMP threads  : " << num_threads << std::endl;
      std::cout << "   Run time     : " << end_time - start_time << " s" << std::endl;
   }

   // Wait for all processes
   MPI_Barrier(MPI_COMM_WORLD);

   // Finalize MPI
   MPI_Finalize();

}
