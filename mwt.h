/////////////////////////////////////////////////////////////////////////////////
//
//      mwt is an array class to store matrices periodically
//
//      0th element contains value of the function at 0
//      sequentially positive domain of the function is stored
//      nth element stores therefore f( +|x_max| ) value
//
//      next 2*n + 3 elements contain just zeros
//
//      (3*n + 4)th element stores f( -|x_max| ) value
//      sequentially negative domain of the function is stored
//      (4*n + 3)th element stores therefore f( -|x_max|/n ) value
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MWT_H
#define MWT_H

#include "global.h"
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



//// Acronym MWT stands for ArrayWithTails

class MWT: public base_array
{
    protected:
        void init();     // function to actually construct the object class

    public:
        MWT(std::string, const array_prop_t);     // constructor of the class
        ~MWT();                                   // deconstructor of the class object


        ////////////////  FOURIER TRANSFORMs  ////////////////////////////////////////////
        fftw_complex * yINTER;    // DFT of y
        fftw_plan forwardFFT;
        fftw_plan backwardFFT;

        void forwardDFT();        // constructor of an empty function forwardDFT
        void backwardDFT();       // constructor of an empty function backwardDFT


        //////////////  OUTPUT MWTs  ///////////////////////////
        void output(std::string, all_t& );

        /////////////// COPY MWTs   //////////////////////////
        void  copy_all(MWT &);
        void copy_real(MWT &);
        void copy_imag(MWT &);

        /////////////// SET MWTs    //////////////////////////
        void set_zero();
        void set_real(double);
        void set_imag(double);

        void set_FD();
        void set_BE();
        void set_K3();

        /////////////// OPERATIONS on MWTs   /////////////////
        void conjugate_y(MWT &);
        void conjugate_dft(MWT &);
};

#endif // MWT_H

