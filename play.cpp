#include "play.h"

#include "global.h"

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





// Formal constructor of play
play::play(const char *value)
{
    #ifdef CALL
        std::cout << "constructing play " << this << std::endl;
    #endif // OLD

    if(value)
    {
        data = new char[strlen(value) + 1];
        strcpy(data, value);
    }
    else
    {
        data = new char[1];
        *data = '\0';
    }

}

// The destructor of play
play::~play()
{
    #ifdef CALL
        std::cout << "destructing play " << this << std::endl;
    #endif // OLD

    delete [] data;
}

std::string play::get_data()
{
    return data;
}

int play::get_length()
{
    return strlen(data);
}

char** play::get_ptr()
{
    return &data;
}

int play::compare(play & Y)
{
    std::cout << "length: " << this->get_length() << std::endl;
    return this->get_length() > Y.get_length();

}
