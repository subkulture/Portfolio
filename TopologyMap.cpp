//------------------------------------------------------------------------------
//                          Mapping topology in MPI
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <iostream>
#include <vector>

// mpi library
#include "mpi.h"

//==============================================================================
// Allocate processess (threads) to cores
//--------------------------------------------------
//
// Global thread id       0 1 2 3 4 5 6 7     8 9 10 11 12 13 14 15
//                            socket 0     |         socket 1
//                           cores  0-7    |        cores 0-7
//
// The following allocations assumes MPI process allocation by core and
// requires the following optons to mpirun:
//
//   --map-by core
//   --bind-to-core
//
//------------------------------------------------------------------------------
// Calculate mpi bandwidth and latency between different processes in topology
//
//==============================================================================
//
int main(int argc, char *argv[]){

   // mpi variables
   int rank = 0; // my rank in all processors
   int num_processors = 1; // total processors in set
   bool root = false; // flag identifying root process

   // initialise mpi
   MPI_Init(&argc, &argv);

   // Get number of processors and rank
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
   if(rank == 0) root = true; // set root process as rank 0

   // constant variables for number of cores per socket
   const unsigned int threads_per_core = 1;
   const unsigned int cores_per_socket = 8;
   const unsigned int sockets_per_node = 2;

   // Calculate node id
   const unsigned int threads_per_node = threads_per_core * cores_per_socket * sockets_per_node;
   const unsigned int node_id = rank / threads_per_node;

   // Calculate socket id
   const unsigned int threads_per_socket = threads_per_core * cores_per_socket;
   const unsigned int socket_id = (rank - threads_per_node * node_id) / threads_per_socket;

   // Calculate core id
   const unsigned int core_id = (rank - threads_per_node * node_id - threads_per_socket * socket_id)/threads_per_core;

   // Calculate thread id
   const unsigned int thread_id = (rank - threads_per_node * node_id - threads_per_socket * socket_id - core_id);

   // MPI Barrier
   MPI_Barrier(MPI_COMM_WORLD);

   //----------------------------------------------------
   // Determine local topology mappings for each process
   //----------------------------------------------------

   // determine topology matrix from my process
   std::vector<unsigned int> nodes(num_processors,0);
   std::vector<unsigned int> sockets(num_processors,0);
   std::vector<unsigned int> cores(num_processors,0);
   std::vector<unsigned int> threads(num_processors,0);

   // set and reduce node ids
   nodes[rank] = node_id;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &nodes[0], num_processors, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

   // set and reduce socket ids
   sockets[rank] = socket_id;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &sockets[0], num_processors, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

   // set and reduce core ids
   cores[rank] = core_id;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &cores[0], num_processors, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

   // set and reduce thread ids
   threads[rank] = thread_id;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &threads[0], num_processors, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

   // print to screen from root to check node, socket core allocations
   if(root){
      for(int p=0; p<num_processors; p++){
         std::cout << p << "\tnode id: " << nodes[p] << "\tsocket id: " << sockets[p] << "\tcore id: " << cores[p] << "\tthread id: " << threads[p] << std::endl;
      }
   }

   //--------------------------------------------------
   // Computing inter core, socket and node bandwidth
   //--------------------------------------------------

   const unsigned int buffer_size = 10000000; // buffer for sending and receiving data
   const double bytes = double(buffer_size)*sizeof(double); // calclulate number of bytes for buffer
   const unsigned int messages = 10; // number of messages to send

   // Print informative message on buffer size
   if(root) std::cout << "Buffer size: " << bytes/1e6 << " MB" << std::endl;

   // set up buffer for data on all processes
   std::vector<double> buffer;
   buffer.resize(buffer_size,1.23456);

   // data for storing timings for transfers and latency
   std::vector<double> timings(num_processors,0.0);
   std::vector<double> latency(num_processors,0.0);

   // determine topology from my rank to others
   // (0 = core-core, 1 = socket-socket, 2 = node-node)
   std::vector<unsigned int> topology(num_processors,0);
   for(int p=0; p<num_processors; p++){
      // determine inter communications types
      bool inode = (node_id == nodes[p]);
      bool isocket = (socket_id == sockets[p]);
      bool icore = (core_id == cores[p]);
      // inter core (same node, same socket)
      if(inode && isocket) topology[p] = 0;
      // inter socket (same node, different socket)
      else if(inode && !isocket) topology[p] = 1;
      // inter node
      else if(!inode) topology[p] = 2;
      // error
      else std::cerr << "Error: unknown inter process communication category" << std::endl;
   }

   // Declare status variable for MPI_Recv calls
   MPI_Status status;

   // Loop over all processes to determine mpi bandwidth and latency from process 0
   for(int p=1; p<num_processors; p++){

      // Print informative message to screen
      if(root) std::cout << "Sending data to and recieving data from process " << p << std::endl;

      //---------------------------------------------------------
      // MPI bandwidth calculation
      //---------------------------------------------------------

      // start the timer
      double start_time = MPI_Wtime();

      // Repeat calculations for better statistics
      for(unsigned int m = 0; m < messages; m++){
         // send the buffer on root
         if(root) MPI_Send(&buffer[0], buffer_size, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
         // receive the buffer on rank p
         else if(rank == p) MPI_Recv(&buffer[0], buffer_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      }

      // end the timer
      double end_time = MPI_Wtime();

      // Calculate timing for bandwidth calculation
      timings[p] = end_time - start_time;

      //---------------------------------------------------------
      // MPI latency calculation (1 value from buffer)
      //---------------------------------------------------------

      // start the timer
      start_time = MPI_Wtime();

      // Repeat calculations 1M times for better statistics
      for(unsigned int m = 0; m < 1000000; m++){
         // send the buffer on root
         if(root) MPI_Send(&buffer[0], 1, MPI_DOUBLE, p, 0, MPI_COMM_WORLD);
         // receive the buffer on rank p
         else if(rank == p) MPI_Recv(&buffer[0], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
      }

      // end the timer
      end_time = MPI_Wtime();

      // Calculate mean latency
      latency[p] = (end_time - start_time)/1000000.0;

      // Wait for all processes to ensure only one process p posts receives at one time
      MPI_Barrier(MPI_COMM_WORLD);

   }

   // Wait for all processes to finish send/receive operations
   MPI_Barrier(MPI_COMM_WORLD);

   // On root output beautified topology map, bandwidth and latency
   if(root){
      std::vector<std::string> key(3);
      key[0] = "c-c";
      key[1] = "s-s";
      key[2] = "n-n";
      std::cout << "------------------------------------------------" << std::endl;
      for(int p=1; p<num_processors; p++){
         std::cout << p << "\t" << key[topology[p]] << "\t" << timings[p] << "\t" << (double(messages)*bytes/timings[p])/1.e9 << " GB/s" << "\t" << latency[p]*1.e9 << " ns" << std::endl;
      }
      std::cout << "------------------------------------------------" << std::endl;
   }

   // Wait for all processes
   MPI_Barrier(MPI_COMM_WORLD);

   // Finalize MPI
   MPI_Finalize();

   return 0;

}
