#include "Constants.h"

#include "base.h"
#include "raw.h"
#include "awt.h"
#include "aux.h"


//#include "play.h"
#include "game.h"

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


#include <iostream>         // standard library for reading inputs
#include <sstream>
#include <fstream>          // standard library for showing outputs

#include <string>           // standard library for manipultaing std::strings
#include <vector>
#include <complex>

#include <iomanip>          // setprecision

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>

//#define BOOST_MATH_INSTRUMENT





int main()
{

    // the name of the file to import desired definitions and values
    all_t params;
    initialize_parameters(params, "input");
    store_calc_parameters(params, "calc");

    // properties array
    array_prop_t stat;
    stat.set_prop_t(10, 1, 0);
    stat.set_prop_all(params);


    // actual working area
    AWT A("A", stat);

    std::cout << A.get_name() << std::endl;
    std::cout << A.get_n() << std::endl;
    std::cout << A.get_x_max() << std::endl;
    std::cout << A.get_kT() << std::endl;

    std::cout << "first copy attempt!" << std::endl;
    AWT B("B", stat);

    std::cout << B.get_name() << std::endl;
    std::cout << B.get_n() << std::endl;
    std::cout << B.get_x_max() << std::endl;
    std::cout << B.get_kT() << std::endl;

    B.copy_all(A);

    raw C("C", stat);



    return 0;
}



