//------------------------------------------------------------------------------
//                   Hierarchical decomposition
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

//------------------------------------------------------------------------------
// Program to generate a square grid of 2D points within a circle of radius 64.
// The points are decomposed into a blocks of size 4 units, where blocks with no
// points are ignored in the decomposition. The blocks are then assigned to
// processors in linear order so that each processor has approximately the same
// number of blocks.
//
// The probram outputs a list of generated points (points.txt) and a list of
// populated blocks including the block coordinates and rank (blocks.txt).
//
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){

   //---------------------------------------------------------------------------
   // Initialise number of processors
   //---------------------------------------------------------------------------
   int num_processors = 8; // total processors in set

   //---------------------------------------------------------------------------
   // Generate 1D list of points 128 x 128 within a circle of radius 64
   //---------------------------------------------------------------------------

   // pair of 1D arrays for storing x-y positions of point within a circle
   std::vector<double> coordinate_x(0);
   std::vector<double> coordinate_y(0);

   // Area fraction of a circle pi*r^2/l^2 with 2% margin for error
   const unsigned int estimated_points = (unsigned int) (3.14159*64.0*64.0*1.02);

   // Print informative message for user
   std::cout << "Estimated " << estimated_points << " data points within circle" << std::endl;

   // reserve memory for 1D arrays
   coordinate_x.reserve(estimated_points);
   coordinate_y.reserve(estimated_points);

   // calculate square of radius of circle
   const double r2 = 64.0*64.0;

   // Allocate points within circle looping over all potential grid points
   for(int i=0; i<128; i++){
      for(int j=0; j<128; j++){
         // calculate offset from origin of circle (assume 64,64)
         const double di = double(i)-64.0;
         const double dj = double(j)-64.0;
         // check if point is within radius
         if(di*di + dj*dj < r2){
            // save coordinates to 1D arrays
            coordinate_x.push_back(double(i));
            coordinate_y.push_back(double(j));
         }
      }
   }

   // Print information for the user
   std::cout << "Actual number for points in circle " << coordinate_x.size() << std::endl;

   // save actual points in circle to disk
   std::ofstream ofile;
   ofile.open("points.txt");
   for(unsigned int i=0; i<coordinate_x.size(); i++){
      ofile << coordinate_x[i] << "\t" << coordinate_y[i] << std::endl;
   }
   ofile.close();

   //---------------------------------------------------------------------------
   // Discretise points into blocks generating a list of points in each block
   //---------------------------------------------------------------------------

   // Print information for the user
   std::cout << "Calculating Hierarchical decomposition" << std::endl;

   // 2D storage of block_id and list of points in block
   std::vector<std::vector<unsigned int> > point_list;

   // 1D storage of block coordinates
   std::vector<double> block_coordinate_x(0);
   std::vector<double> block_coordinate_y(0);

   // set number of blocks in x and y based on resolution
   const double block_resolution = 4.0;
   const unsigned int blocks_in_xy = ceil(128.0/block_resolution); // rounded up
   const unsigned int max_blocks = blocks_in_xy*blocks_in_xy;

   // resize storage for all blocks
   point_list.resize(max_blocks);
   block_coordinate_x.resize(max_blocks);
   block_coordinate_y.resize(max_blocks);

   // temprary variable to count block number in 1D
   unsigned int block_id = 0;

   // loop over each block and points to determine if points are in block
   for(int i = 0; i < blocks_in_xy; i++){
      for(int j = 0; j < blocks_in_xy; j++){

         // calculate block coordinates
         const double bx = i*block_resolution + block_resolution*0.5;
         const double by = j*block_resolution + block_resolution*0.5;

         // calculate block minima and maxima
         const double min_x = bx - block_resolution*0.5;
         const double min_y = by - block_resolution*0.5;
         const double max_x = bx + block_resolution*0.5;
         const double max_y = by + block_resolution*0.5;

         // save coordinates in block
         block_coordinate_x[block_id] = bx;
         block_coordinate_y[block_id] = by;

         // loop over all points and allocate to blocks
         for(unsigned int p=0; p<coordinate_x.size(); p++){

            // save point coordinates into temporary variables
            const double x = coordinate_x[p];
            const double y = coordinate_y[p];

            // check if point is in this block
            if(x >= min_x && x < max_x && y >= min_y && y < max_y){
               // push back point list with point number
               point_list[block_id].push_back(p);
            }

         }

         // increment block id
         block_id++;

      }
   }

   // calculate number of non-empty blocks
   unsigned int num_full_blocks = 0;
   for(unsigned int b=0; b<point_list.size(); b++){
      if(point_list[b].size()>0){
         num_full_blocks++;
      }
   }

   // Print information for the user
   std::cout << "Number of blocks containing points " << num_full_blocks << std::endl;

   // Allocate an even share of blocks to processors
   const unsigned int num_blocks_per_processor = ceil(num_full_blocks/num_processors);
   const unsigned int num_blocks_on_last_processor = num_full_blocks-((num_processors-1)*num_blocks_per_processor);
   const unsigned int total_blocks = num_blocks_per_processor*(num_processors-1) + num_blocks_on_last_processor;

   // Print information for the user
   std::cout << "Number of blocks per processor " << num_blocks_per_processor << std::endl;
   std::cout << "Number of blocks on last processor " << num_blocks_on_last_processor << std::endl;
   std::cout << "Total number of blocks in decomposition " << total_blocks << std::endl;

   // array to store which block is on each processor
   std::vector<unsigned int> block_rank(point_list.size(),-1);

   // repurpose block id again to calculate which block is on which rank
   block_id = 0;

   // keep track of how many blocks are assigned to each processor
   unsigned int cpu_count = 0;
   unsigned int cpu_rank = 0;

   // loop over all blocks and allocate to processors in order
   for(unsigned int b=0; b<block_coordinate_x.size(); b++){

      // get num points in block
      const unsigned int num_points_in_block = point_list[b].size();

      // for non-zero blocks, assign to current rank
      if(num_points_in_block > 0){

         block_rank[b] = cpu_rank;

         // increment cpu_count
         cpu_count++;

         // increment cpu_rank after num_blocks are assigned as long as rank is
         // valid (assumes last points are allocated to last processor)
         if(cpu_count > num_blocks_per_processor-1 && cpu_rank < num_processors-1){
            cpu_rank++;
            cpu_count=0;
         }
      }

   }

   // Save list of blocks and points into file
   ofile.open("blocks.txt");

   // loop over all blocks
   for(unsigned int b=0; b<block_coordinate_x.size(); b++){

      // get num points in block
      const unsigned int num_points_in_block = point_list[b].size();

      // Only output blocks with at least 1 data point
      if(num_points_in_block > 0){

         // x, y, num_points, rank, point list
         ofile << block_coordinate_x[b] << "\t" << block_coordinate_y[b] << "\t";
         ofile << num_points_in_block << "\t" << block_rank[b] << "\t";

         // loop over all points in block and output to file
         for(unsigned int p = 0; p < num_points_in_block; p++) ofile << point_list[b][p] << "\t";

         // output end of line after each list of points
         ofile << std::endl;

      }
   }

   // close output file
   ofile.close();

   return 0;

}
