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
        AWT();
        void init(int,double,double);     // function to actually construct the object class

        ~AWT();                                    // deconstructor of the class object


        // n, nn, xMax, kT are inherited from base

        /////////////   MEMBERS   //////////////////////////////////////////////////////
        fftw_complex * yINTER;    // DFT of y

        ////////////////  FOURIER TRANSFORMs  ////////////////////////////////////////////
        fftw_plan forwardFFT;
        fftw_plan backwardFFT;

        void forwardDFT();        // constructor of an empty function forwardDFT
        void backwardDFT();       // constructor of an empty function backwardDFT



        ////////////////  SET AWTs  ///////////////////////////
        void output(std::string, all_t);
        void set_zero();
        void set_real(double);
        void set_imag(double);

        void set_FD();
        void set_BE();
        void set_K3();

        void set_values(int, int, double, double);


        void conjugate();





};

#endif // AWT_H
