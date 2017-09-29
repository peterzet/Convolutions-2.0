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
        // parameters of the mesh
        int n;
        int nn;
        double xMax;
        double kT;
        bool dftKnown;

        std::complex<double> * y;      // array of function values
        std::complex<double> * yDFT;   // DFT of y, DFT means discrete Fourier transform

    public:
        base_array();
        ~base_array();

        // pure virtual function
        virtual void init(int, double, double) = 0;
        virtual void set_copy();

        virtual void output(std::string, all_t);


        virtual void set_zero();
        virtual void set_real(double);
        virtual void set_imag(double);

        virtual void set_FD();
        virtual void set_BE();
        virtual void set_K3();

        // operations on arrays
        virtual void conjugate_y();
        virtual void conjugate_dft();


        // get information functions
        void print_stat();


};

#endif // BASE_H

