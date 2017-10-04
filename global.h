#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>

#include <iostream>         // standard library for reading inputs
#include <sstream>
#include <fstream>          // standard library for showing outputs

#include <boost/lexical_cast.hpp>

//////////////////////////////////////////////////////////////////////
////////                        MACROS                      /////////
////////////////////////////////////////////////////////////////////

// defines C++ standard that affects output procedures in the program
#define CX11                                // CX11 standard
// #define OLD                              // old standard

// this allow constructors and destructors to write on the console
#define CALL                                // commenting the macro disallows constructor/destructor calls



////////////////    SELECTING FFTW algorithm     /////////////////
// Fourier transforms are calcultated using the fftw3 library, the following macros specify which method of calculation is applied

#define FFTW_METHOD FFTW_ESTIMATE           // instead of actual measurements of different algorithms,
                                            // a simple heuristic is used to pick a (probably sub-optimal)
                                            // plan quickly.

// #define FFTW_METHOD FFTW_MEASURE
// #define FFTW_METHOD FFTW_PATIENT
// #define FFTW_METHOD FFTW_EXHAUSTIVE

//////////////////////////////////////////////////////////////////////
////////                 GLOBAL CONSTANTS                   /////////
////////////////////////////////////////////////////////////////////
const double Pi = 3.14159265359;


//////////////////////////////////////////////////////////////////////
////////              GLOBALLY USED STRUCTURES              /////////
////////////////////////////////////////////////////////////////////

// structure that holds all variables that are externally inputed by user
struct all_t
{
    // members of the struct
    double delta, kT, kT_min, kT_increment, kT_max, U, U_min, U_increment, U_max, mu, mu_min, mu_increment, mu_max, x_max, disp_range, disp_mult;
    int n, disp_density;
    std::string model, print_mode, string_mode, output_mode, physics;
};


// structure that holds only basic mesh properties
struct array_prop_t
{
    // members of the struct
    int n; double x_max; double kT;
    // function to set manually all members of the struct
    void set_prop_t(int _n, double _x_max, double _kT){n = _n; x_max = _x_max; kT = _kT; };
    // function to set all members of the struct from struct all_t
    void set_prop_all(all_t & all){n = all.n; x_max = all.x_max; kT = all.kT; };
};

// structure that holds only mesh export properties
struct exp_prop_t
{
    // members of the struct
    int div; double range; double mult; std::string mode;
};


//////////////////////////////////////////////////////////////////////
////////              GLOBALLY USED FUNCTIONS               /////////
////////////////////////////////////////////////////////////////////


// function to import initial parameters of the mesh and the physical model as well
void initialize_parameters(all_t &, std::string);

// function to export file with parameters of the calculation
void store_calc_parameters(all_t &, std::string);




#endif // GLOBAL_H
