#ifndef play_H
#define play_H

#include "global.h"
#include "play.h"

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



class play
{
    private:
        char *data;



    public:

        play(const char *value = 0);
        ~play();

        std::string get_data();

        int get_length();

        char** get_ptr();

        int compare(play &);

};

#endif // play_H


