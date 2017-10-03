#include "awt.h"

#include "Constants.h"

#include <string>           // standard library for manipultaing std::strings
#include <vector>
#include <complex>          // implements the std::complex class to contain std::complex numbers in cartesian form


#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <cassert>          // error handling library, function assert to terminate the program

#include <fftw3.h>          // FFTW library
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting
#include <stdlib.h>


////////////////    SELECTING FFTW algorithm     /////////////////



#define FFTW_METHOD FFTW_ESTIMATE           // instead of actual measurements of different algorithms,
                                            // a simple heuristic is used to pick a (probably sub-optimal)
                                            // plan quickly.

// #define FFTW_METHOD FFTW_MEASURE
// #define FFTW_METHOD FFTW_PATIENT
// #define FFTW_METHOD FFTW_EXHAUSTIVE

// Formal constructor of AWT
AWT::AWT(std::string initName, const array_prop_t initStat): base_array(initName, initStat)
{
    init(stat);

    #ifdef CALL
        std::cout << "constructing AWT " << initName << std::endl;
    #endif // OLD
}

// The initializer of AWT
void AWT::init(array_prop_t _stat)
{
    stat.n    = _stat.n;              // 2_n + 1 is the number of the mesh points
    stat.x_max = _stat.x_max;              // <-xMax, xMax> is the interval on which the function is defined
    nn   = 4 * stat.n  + 4;      // nn is the size of the complete AWT array
    stat.kT = _stat.kT;               // kT is the temperature associated with the AWT

    // data is allocated for y in the heap, the argument of fftw_malloc defines the size of the input array as nn = 4*n +4
    y = (std::complex<double> *) fftw_malloc(nn * sizeof(std::complex<double>));

    // this informs you whenever there is no space to store y
    if(y == NULL)
    {
        std::cerr << "insufficient memory!" << std::endl;  // std::cerr
        exit(1);
    }

    // data is allocated, the argument of fftw_malloc defines the size of the input array as nnDFT = 4*n +4
    yINTER = (fftw_complex *) fftw_malloc(nn * sizeof(fftw_complex));

    // this informs you whenever there is no space to store yDFT
    if(yINTER == NULL)
    {
        std::cerr << "insufficient memory!" << std::endl;
        exit(1);
    }

    //required because of variable type mismatch
    yDFT = (std::complex<double> *) yINTER;

    // A plan forwardFFT is created for forward fast Fourier transform
    forwardFFT  = fftw_plan_dft_1d(nn, (fftw_complex *) y, (fftw_complex *) yDFT, FFTW_FORWARD, FFTW_METHOD | FFTW_PRESERVE_INPUT);

    // A plan fbackwardFFT is created forbackward fast Fourier transform
    backwardFFT = fftw_plan_dft_1d(nn, (fftw_complex *) yDFT, (fftw_complex *) y, FFTW_BACKWARD, FFTW_METHOD | FFTW_PRESERVE_INPUT);

    dftKnown = false;

}


// The destructor of AWT
AWT::~AWT()
{
    fftw_free(y);
    fftw_free(yINTER);
    fftw_destroy_plan(forwardFFT);
    fftw_destroy_plan(backwardFFT);

    #ifdef CALL
        std::cout << "destructing AWT " << name << std::endl;
    #endif // OLD
}

//////////////////////////////////////////////////////////////////////////
/////////////          FOURIER TRANSFORMs          //////////////////////
////////////////////////////////////////////////////////////////////////

// Forward and Backward discrete Fourier transform
void AWT::forwardDFT()
{
    if(dftKnown)
    return;
    fftw_execute(forwardFFT);
    dftKnown = true;
}

void AWT::backwardDFT()
{
    int n = stat.n;

    assert(dftKnown);
    fftw_execute(backwardFFT);

    // normalize with the size of the transformed array and erase garbage at the place of the zero padding
    for(int i=0;     i<=n; i++)      y[i] /= nn;
    for(int i=n+1; i<nn-n; i++)      y[i] = 0;
    for(int i=nn-n;  i<nn; i++)      y[i] /= nn;
}

//////////////////////////////////////////////////////////////////////
/////////////            OUTPUT AWTs           //////////////////////
////////////////////////////////////////////////////////////////////
void AWT::output(std::string name, all_t & all)
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

    int nn           = 4*all.n + 4;
    int n = all.n;
    double x_max = all.x_max;
    double disp_mult = all.disp_mult;


    // counter to ensure that writting was successful
    int control = 0;

    // x coordinate of maximum range is calculated
    int limit = all.n * all.disp_range / all.x_max;


    // output of y values in the form of a function is selected
    if(all.output_mode == "raw")
    {
        // first the negative values are printed
        for(int i=nn-limit;   i<nn;     i += all.disp_range)
        {
            output << (i-nn)* x_max/n << "  "
                   << disp_mult * real(y[i]) << "  "
                   << disp_mult * imag(y[i]) << "  "
                   << std::endl;
        }
        for(        int i=0;      i<limit;      i += all.disp_range)
        {
            output << i * x_max/n << "  "
                   << disp_mult * real(y[i]) << "  "
                   << disp_mult * imag(y[i]) << "  "
                   << std::endl;
        }
        control += 1;
    }

    // output of y values in the AWT form is selected
    if(all.output_mode == "awt")
    {
        // the whole AWT with zero padding is outputed
        for(int i = 0;   i<nn;     i += all.disp_range)
        {
            output << i << "  "
                   << all.disp_mult * real(y[i]) << "  "
                   << all.disp_mult * imag(y[i]) << "  "
                   << std::endl;
        }
        control += 1;
    }

    // output of DFT values is selected
    if(all.output_mode == "dft")
    {
        // first the negative values are printed
        for(int i=nn-limit;   i<nn;     i += all.disp_range)
        {
            output << (i-nn)* all.x_max/n << "  "
                   << all.disp_mult * real(yDFT[i]) << "  "
                   << all.disp_mult * imag(yDFT[i]) << "  "
                   << std::endl;
        }
        for(        int i=0;      i<limit;      i += all.disp_range)
        {
            output << i * all.x_max/all.n << "  "
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
/////////////              COPY AWTs            /////////////////////
////////////////////////////////////////////////////////////////////

// the array is filled with std::complex zero

void AWT::copy_all(AWT & inX)
{
    // checking if mesh properties are the same
    int compatibility = 1;
    if(stat.n != inX.stat.n)
    {
        std::cerr << "cannot set AWT to another AWT: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(stat.x_max != inX.stat.x_max)
    {
        std::cerr << "cannot set AWT to another AWT: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the AWT " << inX.name << " to " << name << "!"<< std::endl;
        exit(1);
    }

    // actual copying of the data
    for(int i = 0;   i<inX.nn;     i += 1)      y[i] = inX.y[i];

}

// real part is copied
void AWT::copy_real(AWT & inX)
{
        // checking if mesh properties are the same
    int compatibility = 1;
    if(stat.n != inX.stat.n)
    {
        std::cerr << "cannot copy AWT to another AWT: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(stat.x_max != inX.stat.x_max)
    {
        std::cerr << "cannot copy AWT to another AWT: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the AWT " << inX.name << " to " << name << "!"<< std::endl;
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
void AWT::copy_imag(AWT & inX)
{
        // checking if mesh properties are the same
    int compatibility = 1;
    if(stat.n != inX.stat.n)
    {
        std::cerr << "cannot copy AWT to another AWT: incompatible number of mesh points!" << std::endl;
        compatibility -= 1;
    }

    if(stat.x_max != inX.stat.x_max)
    {
        std::cerr << "cannot copy AWT to another AWT: incompatible range of the functions!" << std::endl;
        compatibility -= 1;
    }

    if(compatibility != 1)
    {
        std::cout << "Calculation corrupted, could not copy the AWT " << inX.name << " to " << name << "!"<< std::endl;
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
/////////////              SET AWTs            //////////////////////
////////////////////////////////////////////////////////////////////

// the array is filled with std::complex zero
void AWT::set_zero()
{
    std::complex<double> u(0,0);
    for(int i=0; i<nn; i++)         y[i] = u;
}

void AWT::set_real(double numb)
{
    int n =stat.n;

    std::complex<double> u(numb, 0);
    for(int i=0;     i<=n; i++)      y[i] = u;
    for(int i=n+1; i<nn-n; i++)      y[i] = 0;
    for(int i=nn-n;  i<nn; i++)      y[i] = u;

}

void AWT::set_imag(double numb)
{
    int n =stat.n;

    std::complex<double> u(0, numb);
    for(int i=0;     i<=n; i++)      y[i] = u;
    for(int i=n+1; i<nn-n; i++)      y[i] = 0;
    for(int i=nn-n;  i<nn; i++)      y[i] = u;
}

// FERMI-DIRAC distribution
void AWT::set_FD()
{
    int n =stat.n;
    double kT = stat.kT;
    double xMax = stat.x_max;

    if( kT == 0 )
    {
                                           y[0] = 0.5;
        for(int i=1;     i<nn-n;  i++)     y[i] = 0;
        for(int i=nn-n;  i<nn;    i++)     y[i] = 1.0;
    }
    else
    {
                                            y[0] = 1.0/( exp( +1e-12) + 1 );
        for(int i=0;     i<=n;    i++)      y[i] = 1.0/( exp(i*xMax/ (n * kT)  ) + 1 );
        for(int i=n+1; i<nn-n;    i++)      y[i] = 0;
        for(int i=nn-n;  i<nn;    i++)      y[i] = 1.0/( exp(  (i - 4*n -4)*xMax / (n *kT)  ) + 1 );
    }
}

// BOSE-EINSTEIN distribution
void AWT::set_BE()
{
    int n =stat.n;
    double kT = stat.kT;
    double xMax = stat.x_max;


    if( kT == 0 )
    {
                                            y[0] = 0.5;
        for(int i=1;     i<nn-n;  i++)      y[i] = 0;
        for(int i=nn-n;  i<nn;    i++)      y[i] = -1.0;
    }
    else
    {
                                            y[0] = 1.0/( exp( +1e-12/kT ) - 1 );
        for(int i=0;     i<=n;    i++)      y[i] = 1.0/( exp(i*xMax/ (n * kT)  ) - 1 );
        for(int i=n+1; i<nn-n;    i++)      y[i] = 0;
        for(int i=nn-n;  i<nn;    i++)      y[i] = 1.0/( exp( (i - 4*n -4)* xMax / (n*kT) ) - 1 );
    }
}

// Kernel3 for KRAMMERS KRONIG
void AWT::set_K3()
{
    int n =stat.n;
    double kT = stat.kT;
    double xMax = stat.x_max;


    y[1] =  0;
    y[1] =  4*log(3.0/2.0);
    y[2] = 10*log(4.0/3.0) -6*log(3.0/2.0);

    int i;
    for(i=3; i<=n; i++)
    {
        if(i<1000)
        {
            y[i] =     (1-i*i)*(i-2)         *atanh(1.0/(1 - 2*i))
                   +   (1-i*i)*(i+2)         *atanh(1.0/(1 + 2*i))
                   +   (i*i*i-6*i*i+11*i-6)  *atanh(1.0/(3 - 2*i))/3.0
                   +   (i+3)*(i*i+3*i+2)     *atanh(1.0/(3 + 2*i))/3.0;
        }
        else    y[i] = 1.0/i;

    }
    for(int i=n+1;   i<nn-n;   i++)         y[i] =  0;
    for(int i=nn-n;  i<nn;     i++)         y[i] =  -y[nn-i];
    for(int i=0;     i<nn;     i++)      yDFT[i] =  0;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////                               FUNCTIONS CHANGING ARRAYS                                         ////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
void AWT::conjugate_y(AWT & inX)
{
    int n =stat.n;


    for(int i = 0;      i<n;    i++)      y[i]  =  conj(inX.y[i]);
    for(int i = 3*n+4;  i<nn;   i++)      y[i]  =  conj(inX.y[i]);
}

void AWT::conjugate_dft(AWT & inX)
{
    int n =stat.n;


    for(int i = 0;      i<n;    i++)      yDFT[i]  =  conj(inX.yDFT[i]);
    for(int i = 3*n+4;  i<nn;   i++)      yDFT[i]  =  conj(inX.yDFT[i]);
}


