#include "Constants.h"

#include "base.h"
#include "raw.h"
#include "awt.h"
#include "aux.h"

/*
#include "Calculator.h"
#include "physics.h"
#include "test.h"


#include "Methods.h"
#include "EffectiveT0.h"
#include "Effective_freq.h"
*/

#include <boost/math/tools/minima.hpp>
#include <cassert>          // error handling library, function assert to terminate the program
#include <cmath>            // declares some common mathematical operations and transformation
#include <cstdlib>          // several general purpose functions, including dynamic memory management,
                            // random number generation, communication with the environment,
                            // integer arithmetics, searching, sorting and converting
#include <fftw3.h>          // FFTW library
#include <fstream>          // standard library for showing outputs

#include <iostream>         // standard library for reading inputs

#include <string>           // standard library for manipultaing std::strings
#include <sstream>
#include <vector>
#include <iomanip> // setprecision

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

//#define BOOST_MATH_INSTRUMENT





int main()
{
    // the name of the file to import desired definitions and values
    std::string name = "input";

    // structure that holds all variables that need to be inputed
    all_t all;

    // function that imports all the desired initial definitions and values
    //import_initial(name, all);




    //AWT test;
    //test.init(1000,10,0);
    //test.set_FD();


    std::string callme;
    callme = "FD";


    //all.mult = 1.0;
    //test.output(callme, all);





    return 0;
}



