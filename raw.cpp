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
}

raw::~raw()
{
    //dtor
}


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





void raw::output(std::string name, all_t all)
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

    all.nn           = 4*all.n + 4;

    // counter to ensure that writting was successful
    int control = 0;

    // x coordinate of maximum range is calculated
    int limit = all.n * all.disp_range / all.x_max;


    // output of y values in the form of a function is selected
    if(all.output_mode == "raw")
    {

        control += 1;
    }

    // output of y values in the AWT form is selected
    if(all.output_mode == "awt")
    {

        control += 1;
    }

    // output of DFT values is selected
    if(all.output_mode == "dft")
    {

        control += 1;
    }

    // if output was successfull control equals to one
    if(control != 1)        std::cout << "output of " << name << " was not successful" << std::endl;

}

