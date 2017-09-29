#ifndef BASE_H
#define BASE_H

#include "Constants.h"

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



class base_array
{
    protected:
        base_array();
        ~base_array();


    public:

        virtual void init(int, double, double) = 0;
        virtual void output(std::string, all_t) = 0;
        virtual void conjugate() = 0;

        // set operations
        virtual void set_zero() = 0;
        virtual void set_real(double) = 0;
        virtual void set_imag(double) = 0;

        virtual void set_FD() = 0;
        virtual void set_BE() = 0;
        virtual void set_K3() = 0;


        // parameters of the mesh
        int n;
        int nn;
        double xMax;
        double kT;
        bool dftKnown;

        std::complex<double> * y;      // array of function values
        std::complex<double> * yDFT;   // DFT of y, DFT means discrete Fourier transform



};

#endif // BASE_H

