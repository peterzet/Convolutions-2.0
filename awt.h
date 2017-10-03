#ifndef AWT_H
#define AWT_H

#include "Constants.h"
#include "raw.h"

#include <vector>
#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <string>           // standard library for manipultaing strings
#include <cassert>          // error handling library, function assert to terminate the program

#include <complex> // implements the complex class to contain complex numbers in cartesian form
#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting
#include <stdlib.h>



//// Acronym AWT stands for ArrayWithTails

class AWT: public base_array
{
    public:
        AWT(std::string, const array_prop_t);
        void init(array_prop_t stat);     // function to actually construct the object class

        ~AWT();                                    // deconstructor of the class object


        // n, nn, xMax, kT are inherited from base


        ////////////////  FOURIER TRANSFORMs  ////////////////////////////////////////////
        fftw_complex * yINTER;    // DFT of y
        fftw_plan forwardFFT;
        fftw_plan backwardFFT;

        void forwardDFT();        // constructor of an empty function forwardDFT
        void backwardDFT();       // constructor of an empty function backwardDFT


        //////////////  OUTPUT AWTs  ///////////////////////////
        void output(std::string, all_t& );

        /////////////// COPY AWTs   //////////////////////////
        void  copy_all(AWT &);
        void copy_real(AWT &);
        void copy_imag(AWT &);

        /////////////// SET AWTs    //////////////////////////
        void set_zero();
        void set_real(double);
        void set_imag(double);

        void set_FD();
        void set_BE();
        void set_K3();

        /////////////// OPERATIONS on AWTs   /////////////////
        void conjugate_y(AWT &);
        void conjugate_dft(AWT &);
};

#endif // AWT_H
