#ifndef game_H
#define game_H

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



class game
{

    public:
        game(array_prop_t);
        ~game();

        int nn;
        bool dftKnown;

        int get_n();
        double get_xMax();
        double get_kT();
        std::string get_name();

        void init(array_prop_t&);


};

#endif // game_H


