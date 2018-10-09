#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include "smartastro.h"

#include <string>

#include <cspice/SpiceUsr.h>

#include <stdio.h>

int main()
{

    // String date
    std::string strDate = "Thu Mar 20 12:53:29 PST 1997";
    char * date = new char[strDate.size() + 1];
    std::copy(strDate.begin(), strDate.end(), date);
    date[strDate.size()] = '\0'; // don't forget the terminating 0

    // Kernel
    std::string strLeap = "AstroData/SpiceKernels/lsk/naif0012.tls";
    char * leap = new char[strLeap.size() + 1];
    std::copy(strLeap.begin(), strLeap.end(), leap);
    leap[strLeap.size()] = '\0'; // don't forget the terminating 0

    SpiceDouble et;

    furnsh_c ( leap      );
    str2et_c ( date, &et );

    printf ( "%f\n", et );

    std::cout << "\nNothing coded here. See examples.\n" << std::endl;

    // don't forget to free the string after finished using it
    delete[] date;
    delete[] leap;

    return 0;
}
