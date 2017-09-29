#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

#include <iostream>         // standard library for reading inputs
#include <sstream>
#include <fstream>          // standard library for showing outputs


#include <boost/lexical_cast.hpp>


// constants
const double Pi = 3.14159265359;



// structure that holds all variables that are externally inputed
struct all_t {     double delta, mu, kT_min, kT_increment, kT_max, U_min, U_increment, U_max, mu_min, mu_increment, mu_max, x_max, disp_range, disp_mult;
                      int n, nn, disp_density;
              std::string model, print_mode, string_mode, output_mode, physics;
             };


// structure that holds only mesh export properties
struct exp_prop_t {int n; double xMax; double kT; int div; double range; double mult; std::string mode; };

// structure that holds only basic mesh properties
struct array_prop_t {int n; double xMax; double kT; };



// function to import initial parameters of the mesh and the physical model as well
void initialize_parameters(all_t &, std::string);




#endif // CONSTANTS_H
