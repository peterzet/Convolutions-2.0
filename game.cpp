#include "game.h"

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





// Formal constructor of game
game::game(array_prop_t initStat)
{
    #ifdef CALL
        std::cout << "constructing game" << std::endl;
    #endif // OLD

    init(initStat);
}

// The destructor of game
game::~game()
{

    #ifdef CALL
        std::cout << "destructing game" << std::endl;
    #endif // OLD
}

void game::init(array_prop_t& _stat)
{


}
