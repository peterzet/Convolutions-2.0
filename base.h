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
        double xMax;
        double kT;

    public:
        base_array();
        ~base_array();

        int nn;
        bool dftKnown;

        std::complex<double> * y;      // array of function values
        std::complex<double> * yDFT;   // DFT of y, DFT means discrete Fourier transform

        // pure virtual function
        virtual void init(int, double, double)=0;

        // copy functions
        virtual void  copy_all();
        virtual void copy_real();
        virtual void copy_imag();


        virtual void output(std::string, all_t &);


        virtual void set_zero() =0;
        virtual void set_real(double) =0;
        virtual void set_imag(double) =0;

        virtual void set_FD() =0;
        virtual void set_BE() =0;
        virtual void set_K3()=0;

        // operations on arrays that change these
        virtual void conjugate_y();
        virtual void conjugate_dft();

        // get information functions
        void print_stat();


};

#endif // BASE_H

