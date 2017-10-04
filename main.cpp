#include "global.h"

#include "base.h"

#include "raw.h"
#include "awt.h"
#include "mwt.h"

#include "aux.h"


#include "play.h"
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

/*
play* get_pt_play(play in)
{
    return &in;
};

play** get_ptpt_play(play in)
{
    play *a;
    a = get_pt_play(in);
    return &a;
}

play*** get_ptptpt_play(play in)
{
    play** a;
    a = get_ptpt_play(in);
    return &a;
}
*/

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

/*
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


    char aa;
    std::cout << sizeof(aa)  << std::endl;

    int b;
    std::cout << sizeof(b)  << std::endl;

    double c;
    std::cout << sizeof(c)  << std::endl;


    play s  = "dadaism";
    std::cout << s.get_data() << std::endl;
    std::cout << s.get_ptr()  << std::endl;
    std::cout << sizeof(s.get_ptr())  << std::endl;
    std::cout << s.get_length()  << std::endl;


    play d;
    std::cout << d.get_data() << std::endl;
    std::cout << d.get_ptr()  << std::endl;
    std::cout << sizeof(d.get_ptr())  << std::endl;
    std::cout << d.get_length()  << std::endl;


    std::cout << s.compare(d) << std::endl;

    play* s_pt = &s;
    play* d_pt = &d;

    play** r;

    play a;
    a = *s_pt;



    std::cout << a.get_data() << std::endl;
    std::cout << a.get_ptr()  << std::endl;
    std::cout << sizeof(a.get_ptr())  << std::endl;
    std::cout << a.get_length()  << std::endl;

    std::cout << &s << std::endl;
    std::cout << &d << std::endl;


    int i = 5;
    int* i_pt = &i;

    std::cout <<  i << std::endl;
    std::cout <<  *i_pt << std::endl;
    */

    AWT B("B", stat);

    std::cout << B.get_name() << std::endl;
    std::cout << B.get_n() << std::endl;
    std::cout << B.get_x_max() << std::endl;
    std::cout << B.get_kT() << std::endl;


    return 0;
}



