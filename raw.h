/////////////////////////////////////////////////////////////////////////////////
//
//      raw is an array class to store function from -|x_max| to +|x_max|
//      0th element contains value of the function at -|x_max|
//      last element contains value of the function at +|x_max|
//      elements between contain function values in the corresponding sequence
//
////////////////////////////////////////////////////////////////////////////////


#ifndef RAWFUNCTION_H
#define RAWFUNCTION_H

#include "Constants.h"
#include "base.h"

#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <vector>
#include <complex>          // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting



class raw: public base_array
{
    protected:




    public:
        raw();
        ~raw();


        void init(int, double, double);

        void output(std::string, all_t);

        // set operations
        void set_zero();
        void set_real(double);
        void set_imag(double);

        void set_FD();
        void set_BE();
        void set_K3();
        void set_copy(raw &);



        // operations changing existing values


};

#endif // RAWFUNCTION_H
