#include "raw.h"

#include <string>           // standard library for manipultaing strings
#include <vector>
#include <complex> // implements the std::complex class to contain std::complex numbers in cartesian form


#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <cassert>          // error handling library, function assert to terminate the program


#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting




base_array::base_array()
{
    init(0, 0, 0);
}


void base_array::init(int _n, double _xMax, double _kT)
{
    // these are the basic para,eters of the function
    n    = _n;                  // 2*n + 1 is the number of the mesh points
    nn   = 2 * n + 1;
    xMax = _xMax;               // <-xMax, xMax> is the interval on which the function is defined
    kT = _kT;

}

base_array::~base_array()
{
    //dtor
}

// the array is filled with std::complex zero
void base_array::set_zero()
{

}

void base_array::set_real(double numb)
{


}

void base_array::set_imag(double numb)
{

}

// The initializer of FERMI-DIRAC and BOSE-EINSTEIN distribution
void base_array::set_FD()
{

}

void base_array::set_BE()
{

}

// Kernel3 for KRAMMERS KRONIG
void base_array::set_K3()
{


}





void base_array::output(std::string name, all_t all)
{

}

void base_array::conjugate()
{

}







