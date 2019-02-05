#include <iostream>

#include "Astro-Core/conversion_frames.h"

using namespace std;


int main()
{
    // try conversion frame
    vector<double> scCar(6), vRth(3), vCar(3);

    scCar[0] = -0.2215;
    scCar[1] = +0.0469;
    scCar[2] = +0.4575;
    scCar[3] = +0.4649;
    scCar[4] = -0.3424;
    scCar[5] = +0.4706;

    vRth[0] = +0.4572;
    vRth[1] = -0.0146;
    vRth[2] = +0.3003;

    smartastro::astrocore::conversion_frames::rth2car(scCar,vRth,vCar);

    cout << "vCar = ";
    for (auto i = 0 ; i < 3; i++)
        cout << vCar[i] << " " ;
    cout << endl;

    // Convert back
    smartastro::astrocore::conversion_frames::car2rth(scCar,vCar,vRth);

    cout << "vRth = ";
    for (auto i = 0 ; i < 3; i++)
        cout << vRth[i] << " " ;
    cout << endl;

    cout << endl << "Goodbye from SMART-ASTRO main" << endl;

    return 0;
}
