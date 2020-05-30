//------------------------------------------------------------------------------
//               MPI-IO
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// mpi library
#include "mpi.h"

//------------------------------------------------------------------------------
// MPI Parallel program to read in a plain text file input.txt on the root
// process and broadcast its contents to all other processes. A copy of the file
// is then output from each process into input<rank>.txt to verify the correct
// broadcast of the data. The code then distributes a large array among all
// processes and outputs it to disk using MPI-IO routines.
//
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){

   //---------------------------------------------------------------------------
   // Initialise MPI
   //---------------------------------------------------------------------------

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

   //---------------------------------------------------------------------------
   // Open plain text file and broadcast to all processors
   //---------------------------------------------------------------------------

   // Print out informative message to user
   if(root) std::cout << "Reading input file on root process and broadcasting to all" << std::endl;

   // Define vector of character array to store characters
   std::vector<char> character_array(0);

   // Open input file on root process
   if(root){

      std::ifstream ifile("input.txt");
      std::string line; // declare a string to hold line of text
      std::string file_contents; // declare a string to hold file contents

      // constant string for manually appending newlines \n to file_contents
      const std::string newline = "\n";

      // Loop over all lines until the end of the file
      while( getline(ifile,line) ){
         // append line to file contents preserving new lines
         file_contents += line+newline;
      }

      // close input file
      ifile.close();

      // For debugging print file to screen
      // std::cout << file_contents << std::flush;

      // resize character array on root process
      character_array.resize(file_contents.size());

      // copy string to character array on root process only
      for(unsigned int idx = 0; idx < file_contents.size(); idx++){
         character_array[idx] = file_contents[idx];
      }

   }

   // Broadcast size of array to all processors
   unsigned int character_array_size = character_array.size();
   MPI_Bcast(&character_array_size, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

   // resize character array on all processors to be the correct size
   if(!root) character_array.resize(character_array_size);

   // Broadcast actual character array to all processes
   MPI_Bcast(&character_array[0], character_array.size(), MPI_CHAR, 0, MPI_COMM_WORLD);

   // set filename on each process
   std::stringstream filename_ss; // stringstream
   filename_ss << "input" << rank << ".txt"; // output text to string stream
   std::string filename = filename_ss.str(); // convert stringstream to text

   // open output file on each process
   std::ofstream ofile; // define output file stream
   ofile.open(filename.c_str()); // open file (converting filename to c-style string)

   // output all characters in character array to file
   for(unsigned int idx = 0; idx < character_array.size(); idx++){
      ofile << character_array[idx];
   }

   // close file
   ofile.close();

   // Wait for all processes
   MPI_Barrier(MPI_COMM_WORLD);

   //-----------------------------------------------------------------------------
   //  Distribute array of size N over all processors and output to shared file
   //-----------------------------------------------------------------------------

   // Print out informative message to user
   if(root) std::cout << "Writing data file data.bin to disk with MPI-IO" << std::endl;

   // Specify number of data points in array
   const unsigned int N = 103; // number of data across all processes

   // determine number of data points on my rank
   int num_data = N/num_processors;
   // allocate stray data points to last process
   if(rank == num_processors-1) num_data = N-(num_data*(num_processors-1));

   // Print out informative message showing all online processes
   std::cout << "Processor " << rank+1 << " of " << num_processors << " online with " << num_data << " data points\n" << std::flush;

   // wait here for all to finish
   MPI_Barrier(MPI_COMM_WORLD);

   // check total data points are correct
   int total_data = 0;
   MPI_Reduce(&num_data, &total_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   if(root) std::cout << "Total data points on all processors: " << total_data << std::endl;

   // initialise data with rank to see source of data
   std::vector<int> data_array(num_data, rank);

   // define mpi file variables
   MPI_File handle; // file handle (pointer) for shared MPI-IO
   MPI_Status status;

   // Open MPI file handle on all processes in MPI_COMM_WORLD (collective operation)
   MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_WRONLY|MPI_MODE_CREATE, MPI_INFO_NULL, &handle);

   // Calculate offsets (in bytes) for data so that each process writes to the
   // correct part of the file
   MPI_Offset offset = sizeof(int)*rank*(N/num_processors);

   // All processors now engage in collective write
   MPI_File_write_at_all(handle, offset, &data_array[0], num_data, MPI_INT, &status);

   // Now close MPI file handle
   MPI_File_close(&handle);

   // Wait here for all to finish
   MPI_Barrier(MPI_COMM_WORLD);

   //-----------------------------------------------------------------------------
   //  Verify correctness of shared file
   //-----------------------------------------------------------------------------

   // Print out informative message to user
   if(root) std::cout << "Reading in data.bin to verify correctness" << std::endl;

   // specify a buffer to lead data
   std::vector<int> buffer(0);

   // Collectively open file for reading (note we are repurposing handle for a different file)
   MPI_File_open(MPI_COMM_WORLD, "data.bin", MPI_MODE_RDONLY, MPI_INFO_NULL, &handle);

   // Only read file on root process
   if(root){

      // resize buffer to correct size for data and initialise to -1
      buffer.resize(N,-1);

      // read data from file
      MPI_File_read(handle, &buffer[0], N, MPI_INT, &status);

   }

   // Close input file
   MPI_File_close(&handle);

   // Verify correctness of data for this problem on root process
   if(root){

      // boolean variable to define if all data is correct
      bool correct = true; // assume all data is correct initially

      // loop over all data points
      for(int i=0; i<buffer.size(); i++){

         // calculated expected value of rank
         int expected = i/num_data;
         if(expected > num_processors-1) expected = num_processors-1; // remember surplus points

         // optionally print to screen
         // std::cout << expected << "\t" << buffer[i] << std::endl;

         // ensure read data is what we expect
         if(expected != buffer[i]) correct = false;

      }

      // print informative message to screen
      if(correct) std::cout << "Data file verified as correct" << std::endl;
      else std::cout << "Error in data file" << std::endl;

   }

   // Wait for all processes
   MPI_Barrier(MPI_COMM_WORLD);

   // Finalize MPI
   MPI_Finalize();

   return 0;

}
