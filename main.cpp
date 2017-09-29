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
    std::string name = "input";
    all_t params;
    initialize_parameters(params, name);

    std::cout << params.n << std::endl;

    AWT A;
    A.init(params.n,params.x_max,params.kT_min);
    A.conjugate_y(A);
    A.set_imag(1.0);
    A.conjugate_y(A);


    AWT AA;
    AA.init(2,2,params.kT_min);
    AA.set_copy(A);



    raw B;
    B.init(10,2,0);
    B.conjugate_y();

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



