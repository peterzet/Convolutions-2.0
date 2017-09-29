#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

#include <iostream>         // standard library for reading inputs
#include <sstream>
#include <fstream>          // standard library for showing outputs


#include <boost/lexical_cast.hpp>

const double Pi = 3.14159265359;

enum class exp_mode_t: int  {EXP_RAW = 0, EXP_AWT = 1, EXP_DFT = 2};

// structure that holds all variables that are externally inputed
struct all_t {     double delta, mu, kT_min, kT_increment, kT_max, U_min, U_increment, U_max, mu_min, mu_increment, mu_max, x_max, disp_range, disp_mult;
                      int n, nn, disp_density;
              std::string model, print_mode, string_mode, output_mode, physics;
             };

void initialize_parameters(all_t &, std::string);


// only variables with export properties
struct exp_prop_t {int n; double xMax; double kT; int div; double range; double mult; std::string mode; };

// only variables with array properties
struct array_prop_t {int n; double xMax; double kT; };


#endif // CONSTANTS_H
