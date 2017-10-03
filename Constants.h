#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

#include <iostream>         // standard library for reading inputs
#include <sstream>
#include <fstream>          // standard library for showing outputs


#include <boost/lexical_cast.hpp>

// defines the arrays output method depending on c++ version
#define CX11
// #define OLD

// this allow constructors and destructors to write on the console
#define CALL

#define _(x) std::cout << #x << std::endl; x

// constants
const double Pi = 3.14159265359;



// structure that holds all variables that are externally inputed
struct all_t {     double delta, kT, kT_min, kT_increment, kT_max, U, U_min, U_increment, U_max, mu, mu_min, mu_increment, mu_max, x_max, disp_range, disp_mult;
                      int n, disp_density;
              std::string model, print_mode, string_mode, output_mode, physics;
             };


// structure that holds only basic mesh properties
struct array_prop_t
{
                    int n; double x_max; double kT;

                    void set_prop_t(int _n, double _x_max, double _kT){n = _n; x_max = _x_max; kT = _kT; };

                    void set_prop_all(all_t & all){n = all.n; x_max = all.x_max; kT = all.kT; };
};

// structure that holds only mesh export properties
struct exp_prop_t {int div; double range; double mult; std::string mode; };


// function to import initial parameters of the mesh and the physical model as well
void initialize_parameters(all_t &, std::string);
void store_calc_parameters(all_t &, std::string);




#endif // CONSTANTS_H
