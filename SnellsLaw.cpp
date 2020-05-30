//------------------------------------------------------------------------------
//                 Snell's law with automated testing suite
//------------------------------------------------------------------------------
//
// C++ standard library headers
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

//------------------------------------------------------------------------------
// Function to calculate refracted angle using snells law
//------------------------------------------------------------------------------
double snells_law(double theta, double n1, double n2){
   return asin(n1*sin(theta)/n2);
}

//------------------------------------------------------------------------------
// Function to calculate ray equation with gradient m and intercept c
//------------------------------------------------------------------------------
void update_ray_equation(double x, double y, double theta, double& m, double& c){

   // update gradient based on angle of ray
   m = tan(theta);

   // update intercept by solving for c given current ray position x and y
   c = y - m*x;

   return;

}

// forward function declaration for test suite
void test_suite();

//------------------------------------------------------------------------------
// Function to calculate propagation of a light beam through a series of lenses
// Propagation is determined by snells law given by:
//
//                n_1 sin theta_1 = n_2 sin theta_2
//
// where n1 is the refractive index of the medium for the incoming light ray,
// theta_1 is the angle of incidence from the surface normal of the incoming
// light ray, n_2 is the refractive index of the medium for the outgoing ray,
// and theta_2 is the angle of refraction from the surface normal of the
// outgoing light ray.
// The code calculates the trajectory of an incoming light ray with angle
// theta_0 through a series of lenses. The code then outputs the intersection
// points of the light ray with the lenses.
//
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){

   // capture command line arguments for testing
   if(argc>1){
      std::string arg = argv[1];
      std::string teststr = "--test";
      if(arg==teststr){
         // call test suite
         test_suite();
         return 0;
      }
   }

   // initialise beam
   double x = 0.0; // x-position of light ray
   double y = 0.0; // y position of light ray
   double theta = 20*M_PI/180.0; // initial angle of light ray
   double n1 = 1.0; // initial refractive index

   // set up lense interface positions and refractive indices
   std::vector<double> lense_x(5); // x-position of lenses
   std::vector<double> lense_n(5); // refractive index of lense

   lense_x[0] = 5.0;
   lense_n[0] = 1.5; // glass

   lense_x[1] = 6.0;
   lense_n[1] = 1.0; // air

   lense_x[2] = 10.0;
   lense_n[2] = 1.5; // glass

   lense_x[3] = 11.0;
   lense_n[3] = 1.0; //air

   lense_x[4] = 20.0;
   lense_n[4] = 1.0; //air

   // open data file for output
   std::ofstream ofile;
   ofile.open ("data.txt");

   // write initial point to file
   ofile << x << "\t" << y << std::endl;

   // temporary variable to stire equation of light ray
   double m = 0.0;
   double c = 0.0;

   // determine initial equation of light ray
   update_ray_equation(x, y, theta, m, c);

   // loop over all interfaces to determine path of light ray
   for(int i=0; i< lense_x.size(); i++){

      // get refractive index of material beyond interface
      double n2 = lense_n[i];

      // update x and y coordinates of beam
      x = lense_x[i];
      y = m*x + c;

      // calculate new angle for light ray
      theta = snells_law(theta, n1, n2);

      // write refraction point to file
      ofile << x << "\t" << y << "\t" << theta << std::endl;

      // calculate new equation for light ray
      update_ray_equation(x, y, theta, m, c);

      // update n1 to new medium value
      n1 = lense_n[i];

   }

   // close output file
   ofile.close();

   return 0;

}

//------------------------------------------------------------------------------
// Test suite
//------------------------------------------------------------------------------
void test_suite(){

   // temporary variables used for test
   double theta = 0.0;

   // numerical epsilon for numerical comparisons
   double epsilon = 1.0e-14;

   //----------------------------------------------------------
   // test snells law
   //----------------------------------------------------------
   // Test normal incidence
   //----------------------------------------------------------
   std::cerr << "Testing snells law for normal incidence...                       " << std::flush;
   theta = snells_law(0.0,1.0,1.0);
   // pass
   if(theta > -epsilon && theta < epsilon) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test 30 degree incidence at imaginary interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for 30 degree incidence at imaginary interface" << std::flush;
   const double theta30 = 30.0*M_PI/180.0;
   theta = snells_law(theta30,1.0,1.0);
   // pass
   if(theta > theta30-epsilon && theta < theta30+epsilon) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test 30 degree incidence at n = 1/ n = 2 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for 30 degrees and n = 1 / n = 2 interface... " << std::flush;
   theta = snells_law(theta30,1.0,2.0);
   // pass
   if(theta > 0.252680254 && theta < 0.252680256) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test 30 degree incidence at n = 2/ n = 1 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for 30 degrees and n = 2 / n = 1 interface... " << std::flush;
   theta = snells_law(theta30,2.0,1.0);
   // pass
   if(theta > 1.57079630 && theta < 1.57079632) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test -30 degree incidence at n = 1/ n = 2 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for -30 degrees and n = 1 / n = 2 interface..." << std::flush;
   theta = snells_law(-theta30,1.0,2.0);
   // pass
   if(theta > -0.252680256 && theta < -0.252680254) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test -30 degree incidence at n = 1/ n = 2 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for -30 degrees and n = 2 / n = 1 interface..." << std::flush;
   theta = snells_law(-theta30,2.0,1.0);
   // pass
   if(theta > -1.57079632 && theta < -1.57079630) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test 0 degree incidence for negative n = 1/ n = 2 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for negative refractive index n1 zero degrees " << std::flush;
   theta = snells_law(0,-1.0,1.0);
   // pass
   if(theta > -epsilon && theta < epsilon) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;
   //----------------------------------------------------------
   // Test 30 degree incidence for negative n = 1/ n = 2 interface
   //----------------------------------------------------------
   std::cerr << "Testing snells law for negative refractive index n1 30 degrees   " << std::flush;
   theta = snells_law(theta30,-1.0,1.0);
   std::cout << theta << std::endl;
   // pass
   if(theta > theta30-epsilon && theta < theta30+epsilon) std::cerr << " [pass]" << std::endl;
   else std::cerr << " [fail]" << std::endl;

}
