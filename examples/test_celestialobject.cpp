#define CATCH_CONFIG_MAIN

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>

#include "smartastro.h"
#include "AstroData/SpiceKernels/spiceKernelNames.h"
#include "../src/AstroBodies/celestial_object.h"
#include "../src/Astro-Core/spice_general_functions.h"
#include "catch.hpp"

using namespace std;
using namespace smartastro;
using namespace spiceKernels;

string SPK (smartastro::spiceKernels::planetsEph);
	
std::vector<std::string> params{"J2000", SPK, "moon", "NONE"};
vector<double> positn{1, 2, 3};

smartastro::astrobodies::Celestial_Object celestialObject("Earth", params, positn, 1);

TEST_CASE( "npedln tested", "[npedln]" ) {
	vector<double> pnear(3);
	vector<double> linept{ 1.0e6,  2.0e6,  3.0e6 };
	vector<double> linedr{ -4.472091234e-1, -8.944182469e-1, -4.472091234e-3 }; 
	double dist;
	double a = 7.0e5; 
        double b = 7.0e5; 
        double c = 6.0e5;

	celestialObject.npedln(a, b, c, linept, linedr, pnear, dist);

	REQUIRE(std::fabs(pnear[0] - -.16333110792340931E+04) < 1);
	REQUIRE(std::fabs(pnear[1] - -.32666222157820771E+04) < 1);
	REQUIRE(std::fabs(pnear[2] - .59999183350006724E+06) < 1);
	REQUIRE(std::fabs(dist - .23899679338299707E+06) < 2);
}

TEST_CASE( "npjelpl tested", "[npjelpl]" ) {
	vector<double> normalNVC{ 0.,  0.,  1. , 0.};
	vector<double> ellipseVector{ 1., 1., 1., 2., 0., 0., 0., 1., 1.};
	vector<double> elout(9);

	celestialObject.pjelpl(ellipseVector, normalNVC, 1, elout);

	REQUIRE(elout[0] == 1.);
	REQUIRE(elout[1] == 1.);
	REQUIRE(elout[2] == 0.);
	REQUIRE(elout[3] == 2.);
	REQUIRE(elout[4] == 0.);
	REQUIRE(elout[5] == 0.);
	REQUIRE(elout[6] == 0.);
	REQUIRE(elout[7] == 1.);
	REQUIRE(elout[8] == 0.);
	
}

TEST_CASE( "saelgv tested", "[saelgv]" ) {
	vector<double> smajor(3);
	vector<double> sminor(3);
	vector<double> vec1 { 1.,  1.,  1. };
        vector<double> vec2 { 1., -1.,  1. };

	celestialObject.saelgv(vec1, vec2, smajor, sminor);

	REQUIRE(std::fabs(smajor[0] - -1.414213562373095) < 3);
	REQUIRE(std::fabs(smajor[1] - 0.0) < 3);
	REQUIRE(std::fabs(smajor[2] - -1.414213562373095) < 3);
	REQUIRE(std::fabs(sminor[0] - -2.4037033579794549D-17) < 15);
	REQUIRE(std::fabs(sminor[1] - 1.414213562373095) < 1);
	REQUIRE(std::fabs(sminor[2] - -2.4037033579794549D-17) < 15);
}

TEST_CASE( "nplnpt tested", "[nplnpt]" ) {
	double dist;
	vector<double> pnear(3);
	vector<double> LINPT{1.0, 2.0, 3.0}; 
        vector<double> LINDIR{0.0, 1.0, 1.0};  
        vector<double> POINT{6.0, 9.0, 10.0};

	celestialObject.nplnpt(LINPT, LINDIR, POINT, pnear, dist);

	REQUIRE(pnear[0] == 1.);
	REQUIRE(pnear[1] == 9.);
	REQUIRE(pnear[2] == 10.);
	REQUIRE(dist == 5.0);
}
