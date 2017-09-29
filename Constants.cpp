#include "Constants.h"

#include <string>           // standard library for manipultaing strings
#include <vector>
#include <complex> // implements the std::complex class to contain std::complex numbers in cartesian form


#include <iostream>         // standard library for reading inputs
#include <fstream>          // standard library for showing outputs
#include <cassert>          // error handling library, function assert to terminate the program



void initialize_parameters(all_t & all, std::string name)
{

    std::ifstream ListFile;
    ListFile.open (name.c_str());   // list contains all the file names to be processed
    std::string item;                    // used for temporary store the data

    if (ListFile.fail()) {std::cerr<< "Error opening the file" << std::endl; exit(1);}

    for(int i=0; i<58; i++)
    {
        ListFile >> item;

        if(i == 3)   all.delta        = boost::lexical_cast<double>(item);
        if(i == 5)   all.mu           = boost::lexical_cast<double>(item);
        if(i == 7)   all.model        = boost::lexical_cast<std::string>(item);

        if(i == 12)  all.kT_min       = boost::lexical_cast<double>(item);
        if(i == 15)  all.kT_increment = boost::lexical_cast<double>(item);
        if(i == 18)  all.kT_max       = boost::lexical_cast<double>(item);

        if(i == 21)  all.U_min        = boost::lexical_cast<double>(item);
        if(i == 24)  all.U_increment  = boost::lexical_cast<double>(item);
        if(i == 27)  all.U_max        = boost::lexical_cast<double>(item);

        if(i == 30)  all.mu_min       = boost::lexical_cast<double>(item);
        if(i == 33)  all.mu_increment = boost::lexical_cast<double>(item);
        if(i == 36)  all.mu_max       = boost::lexical_cast<double>(item);


        if(i == 40)  all.x_max        = boost::lexical_cast<double>(item);
        if(i == 42)  all.n            = boost::lexical_cast<int>(item);


        if(i == 47)  all.disp_density = boost::lexical_cast<int>(item);
        if(i == 49)  all.disp_range   = boost::lexical_cast<int>(item);
        if(i == 52)  all.print_mode   = boost::lexical_cast<std::string>(item);
        if(i == 55)  all.string_mode  = boost::lexical_cast<std::string>(item);
        if(i == 58)  all.output_mode  = boost::lexical_cast<std::string>(item);
        if(i == 60)  all.physics      = boost::lexical_cast<std::string>(item);

        all.disp_mult = 1.0;

    }

    ListFile.close();
};

