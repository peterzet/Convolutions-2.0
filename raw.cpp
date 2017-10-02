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




raw::raw()
{
    #ifdef CALL
        std::cout << "constructing raw" << std::endl;
    #endif // OLD

    init(0, 0, 0);
}


void raw::init(int _n, double _xMax, double _kT)
{
    // these are the basic para,eters of the function
    n    = _n;                  // 2*n + 1 is the number of the mesh points
    nn   = 2 * n + 1;
    xMax = _xMax;               // <-xMax, xMax> is the interval on which the function is defined
    kT = _kT;

    y    = (std::complex<double> *) malloc(nn * sizeof(std::complex<double>));
    yDFT = (std::complex<double> *) malloc(nn * sizeof(std::complex<double>));


    // we always initialize the function as constant zero
    std::complex<double> u(0,0);
    for(int i=0; i < nn; i++)    y[i] = u;

    #ifdef CALL
        std::cout << "initializing raw" << std::endl;
    #endif // OLD
}

raw::~raw()
{
    #ifdef CALL
        std::cout << "destructing raw" << std::endl;
    #endif // OLD
}


//////////////////////////////////////////////////////////////////////
/////////////            OUTPUT raws            /////////////////////
////////////////////////////////////////////////////////////////////

void raw::output(std::string name, all_t & all)
{
    // output stream file is open according to old standard
    #ifdef OLD
        std::ofstream output(name.c_str());
    #endif // OLD
    // output stream file is open according to cx11 standard
    #ifdef CX11
        std::ofstream output;
        output.open(name);
    #endif // CX11


    // counter to ensure that writting was successful
    int control = 0;

    // x coordinate of maximum range is calculated
    int limit = all.n * all.disp_range / all.x_max;


    // output of y values in the form of a function is selected
    if(all.output_mode == "raw")
    {
        for(int i = 0; i < nn; i++)
        {
              output << (i-n)* all.x_max/n << "  "
                   << all.disp_mult * real(yDFT[i]) << "  "
                   << all.disp_mult * imag(yDFT[i]) << "  "
                   << std::endl;
        }
        control += 1;
    }

    // output of y values in the raw form is selected
    if(all.output_mode == "raw")
    {
        for(int i = 0; i <= n; i++)
        {
              output << i << "  "
                   << all.disp_mult * real(yDFT[i+all.n]) << "  "
                   << all.disp_mult * imag(yDFT[i+all.n]) << "  "
                   << std::endl;
        }
        for(int i = 0; n+1 < 3*n+4; i++)
        {
              output << i << "  "
                   << 0 << "  "
                   << 0 << "  "
                   << std::endl;
        }
        for(int i = 3*n+4; i < 4*n+4; i++)
        {
              output << i << "  "
                   << all.disp_mult * real(yDFT[i-3*n-4]) << "  "
                   << all.disp_mult * imag(yDFT[i-3*n-4]) << "  "
                   << std::endl;
        }
        control += 1;
    }

    // output of DFT values is selected
    // output of DFT values is selected
    if(all.output_mode == "dft")
    {
        // first the negative values are printed
        for(int i=0;   i<nn;     i += all.disp_range)
        {
            output << (i-n)* all.x_max/n << "  "
                   << all.disp_mult * real(yDFT[i]) << "  "
                   << all.disp_mult * imag(yDFT[i]) << "  "
                   << std::endl;
        }
        control += 1;
    }

    // if output was successfull control equals to one
    if(control != 1)
    {
        std::cout << "output of " << name << " was not successful" << std::endl;
        output << "Did you define the output_mode in the input file correctly?" << std::endl;
    }
}



//////////////////////////////////////////////////////////////////////
/////////////              COPY raws            /////////////////////
////////////////////////////////////////////////////////////////////

// the array is filled with std::complex zero

void raw::copy_all(raw & inX)
{
    // checking if mesh properties are the same
    int compatibility = 1;
    if(n != inX.n)
    {
        std::cerr << "cannot copy raw to another raw: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(xMax != inX.xMax)
    {
        std::cerr << "cannot copy raw to another raw: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the raw!" << std::endl;
        exit(1);
    }

    // actual copying of the data
    for(int i = 0;   i<inX.nn;     i += 1)      y[i] = inX.y[i];

}

// real part is copied
void raw::copy_real(raw & inX)
{
    // checking if mesh properties are the same
    int compatibility = 1;
    if(n != inX.n)
    {
        std::cerr << "cannot copy raw to another raw: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(xMax != inX.xMax)
    {
        std::cerr << "cannot copy raw to another raw: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the raw!" << std::endl;
        exit(1);
    }

    // actual copying of the data
    for(int i=0; i < nn; i++)
    {
        std::complex<double> u(0, imag(inX.y[i]) );
        y[i]= u;
    }
    dftKnown = false;
}

// imag part is copied
void raw::copy_imag(raw & inX)
{
        // checking if mesh properties are the same
    int compatibility = 1;
    if(n != inX.n)
    {
        std::cerr << "cannot copy raw to another raw: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(xMax != inX.xMax)
    {
        std::cerr << "cannot copy raw to another raw: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the raw!" << std::endl;
        exit(1);
    }

    // actual copying of the data
    for(int i=0; i < nn; i++)
    {
        std::complex<double> u(real(inX.y[i]), 0);
        y[i]= u;
    }
    dftKnown = false;
}


//////////////////////////////////////////////////////////////////////
/////////////              SET raws            //////////////////////
////////////////////////////////////////////////////////////////////

// the array is filled with std::complex zero
void raw::set_zero()
{
    std::complex<double> u(0,0);
    for(int i=0;  i<nn;    i++)    y[i] = u;

}

void raw::set_real(double numb)
{
    std::complex<double> u(numb,0);
    for(int i=0;  i<nn;    i++)    y[i] = u;
}

void raw::set_imag(double numb)
{
    std::complex<double> u(0,numb);
    for(int i=0;  i<nn;    i++)    y[i] = u;
}

// The initializer of FERMI-DIRAC distribution
void raw::set_FD()
{
    if( kT == 0 )
    {
        for(int i=0;     i<n;   i++)     y[i] = 1.0;
                                         y[n] = 0.5;
        for(int i=n+1;  i<nn;   i++)     y[i] = 0;
    }
    else
    {
        for(int i=0;     i<n;    i++)    y[i] = 1.0/( exp( (i - n) * xMax/ (n * kT)  ) + 1 );
                                         y[n] = 0.5;
        for(int i=n+1;  i<nn;    i++)    y[i] = 1.0/( exp( (i - n) * xMax / (n *kT)  ) + 1 );
    }
}


// The initializer of BOSE-EINSTEIN distribution
void raw::set_BE()
{
    if( kT == 0 )
    {
        for(int i=0;     i<n;   i++)     y[i] = 1.0;
                                         y[n] = 0.5;
        for(int i=n+1;  i<nn;   i++)     y[i] = 0;
    }
    else
    {
        for(int i=0;     i<n;    i++)    y[i] = 1.0/( exp( (i - n) * xMax/ (n * kT)  ) - 1 );
                                         y[n] = 0.5;
        for(int i=n+1;  i<nn;    i++)    y[i] = 1.0/( exp( (i - n) * xMax / (n *kT)  ) - 1 );
    }
}

// Kernel3 for KRAMMERS KRONIG
void raw::set_K3()
{

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                               FUNCTIONS CHANGING ARRAYS                                         ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void raw::conjugate_y(raw & inX)
{
    for(int i = 0;      i<nn;    i++)      y[i]  =  conj(inX.y[i]);
}

void raw::conjugate_dft(raw & inX)
{
    for(int i = 0;      i<nn;    i++)      yDFT[i]  =  conj(inX.yDFT[i]);
}

