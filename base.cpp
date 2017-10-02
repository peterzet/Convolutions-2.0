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
    #ifdef CALL
        std::cout << "Constructor of virtual base" << std::endl;
    #endif // OLD
}

void base_array::init(int _n, double _xMax, double _kT)
{
    std::cout << "Call from virtual base class base_array: no implementation of initialize" << std::endl;
}

base_array::~base_array()
{
    #ifdef CALL
        std::cout << "Destructor of virtual base" << std::endl;
    #endif // OLD
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                    OUTPUT FUNCTIONS                                             ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void base_array::output(std::string name, all_t & all)
{
    std::cout << "Call from virtual base class base_array: no implementation of output" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                     SET FUNCTIONS                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

// the array is filled with std::complex zero
void base_array::set_zero()
{
    std::cout << "Call from virtual base class base_array: no implementation of set_zero" << std::endl;
}

void base_array::set_real(double numb)
{
    std::cout << "Call from virtual base class base_array: no implementation of set_real" << std::endl;
}

void base_array::set_imag(double numb)
{
    std::cout << "Call from virtual base class base_array: no implementation of set_imag" << std::endl;
}

// The initializer of FERMI-DIRAC and BOSE-EINSTEIN distribution
void base_array::set_FD()
{
    std::cout << "Call from virtual base class base_array: no implementation of set_FD" << std::endl;
}

void base_array::set_BE()
{
    std::cout << "Call from virtual base class base_array: no implementation of set_BE" << std::endl;
}

// Kernel3 for KRAMMERS KRONIG
void base_array::set_K3()
{
    std::cout << "Call from virtual base class base_array: no implementation of set_K3" << std::endl;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                               FUNCTIONS CHANGING ARRAYS                                         ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void base_array::conjugate_y()
{
    std::cout << "Call from virtual base class base_array: no implementation of conjugate_y" << std::endl;
}

void base_array::conjugate_dft()
{
    std::cout << "Call from virtual base class base_array: no implementation of conjugate_y" << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                                     COPY FUNCTIONS                                              ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void base_array::copy_all()
{
    std::cout << "Call from virtual base class base_array: no implementation of set_copy" << std::endl;
}

void base_array::copy_real()
{
    std::cout << "Call from virtual base class base_array: no implementation of delete_real" << std::endl;
}

void base_array::copy_imag()
{
    std::cout << "Call from virtual base class base_array: no implementation of delete_imag" << std::endl;
}


void base_array::print_stat()
{
    std::cout << "n is: " << n << " nn is: " << nn <<  " xMax is: " << xMax <<  " kT is: " << kT << std::endl;
}




