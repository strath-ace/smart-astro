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

void test_npedln(){
	vector<double> pnear(3);
	vector<double> linept{ 1.0e6,  2.0e6,  3.0e6 };
	vector<double> linedr{ -4.472091234e-1, -8.944182469e-1, -4.472091234e-3 }; 
	double dist;
	double a = 7.0e5; 
        double b = 7.0e5; 
        double c = 6.0e5;

	celestialObject.npedln(a, b, c, linept, linedr, pnear, dist);
	
	printf ("npedln: Find nearest point on a triaxial ellipsoid to a specified line, and the distance from the ellipsoid to the line. \n"
		"Input \n"
		"a: %1f\n"
		"b: %1f\n"
		"c: %1f\n"
		"linept[0]: %1f\n"
             	"linept[1]: %1f\n"
	     	"linept[2]: %1f\n"
		"linedr[0]: %1f\n"
             	"linedr[1]: %1f\n"
	     	"linedr[2]: %1f\n"
	     	"Output\n"
             	"pnear[0]: %1f\n"
             	"pnear[1]: %1f\n"
	     	"pnear[2]: %1f\n"
	     	"dist: %1f\n"
		"\n",
             	a, b, c, linept[0], linept[1], linept[2], linedr[0], linedr[1], linedr[2], pnear[0], pnear[1], pnear[2], dist);
}

void test_pjelpl(){
	vector<double> normalNVC{ 0.,  0.,  1. , 0.};
	vector<double> ellipseVector{ 1., 1., 1., 2., 0., 0., 0., 1., 1.};
	vector<double> elout(9);

	celestialObject.pjelpl(ellipseVector, normalNVC, 1, elout);

	printf ("pjelpl: Project an ellipse onto a plane, orthogonally. \n"
		"Input \n"
		"ellipseVector[0]: %1f\n"
		"ellipseVector[1]: %1f\n"
		"ellipseVector[2]: %1f\n"
		"ellipseVector[3]: %1f\n"
             	"ellipseVector[4]: %1f\n"
	     	"ellipseVector[5]: %1f\n"
		"ellipseVector[6]: %1f\n"
             	"ellipseVector[7]: %1f\n"
	     	"ellipseVector[8]: %1f\n"
		"normalNVC[0]: %1f\n"
		"normalNVC[1]: %1f\n"
		"normalNVC[2]: %1f\n"
		"normalNVC[3]: %1f\n"
		"planetype: 1\n"
	     	"Output\n"
             	"elout[0]: %1f\n"
             	"elout[1]: %1f\n"
	     	"elout[2]: %1f\n"
		"elout[3]: %1f\n"
		"elout[4]: %1f\n"
		"elout[5]: %1f\n"
		"elout[6]: %1f\n"
		"elout[7]: %1f\n"
	     	"elout[8]: %1f\n"
		"\n",
             	ellipseVector[0], ellipseVector[1], ellipseVector[2], ellipseVector[3], ellipseVector[4], ellipseVector[5], ellipseVector[6], ellipseVector[7], ellipseVector[8], normalNVC[0], normalNVC[1], normalNVC[2], normalNVC[3], elout[0], elout[1], elout[2], elout[3], elout[4], elout[5], elout[6], elout[7], elout[8]);
}

void test_sealgv(){
	vector<double> smajor(3);
	vector<double> sminor(3);
	vector<double> vec1 { 1.,  1.,  1. };
        vector<double> vec2 { 1., -1.,  1. };

	celestialObject.saelgv(vec1, vec2, smajor, sminor);

	printf ("saelgv: Find semi-axis vectors of an ellipse generated by two arbitrary three-dimensional vectors. \n"
		"Input \n"
		"vec1[0]: %1f\n"
		"vec1[1]: %1f\n"
		"vec1[2]: %1f\n"
		"vec2[0]: %1f\n"
		"vec2[1]: %1f\n"
		"vec2[2]: %1f\n"
		"Output \n"
		"smajor[0]: %1f\n"
		"smajor[1]: %1f\n"
		"smajor[2]: %1f\n",
		vec1[0], vec1[1], vec1[2], vec2[0], vec2[1], vec2[2], smajor[0], smajor[1], smajor[2]);
		cout << "sminor[0]: " << sminor[0] << endl;
		cout << "sminor[1]: " << sminor[1] << endl;
		cout << "sminor[2]: " << sminor[2] << endl;
		cout << endl;
}

void test_nplnpt(){
	double dist;
	vector<double> pnear(3);
	vector<double> LINPT{1.0, 2.0, 3.0}; 
        vector<double> LINDIR{0.0, 1.0, 1.0};  
        vector<double> POINT{6.0, 9.0, 10.0};

	celestialObject.nplnpt(LINPT, LINDIR, POINT, pnear, dist);

	printf ("nplnpt: Find the nearest point on a line to a specified point, and find the distance between the two points.  \n"
		"Input \n"
		"LINPT[0]: %1f\n"
		"LINPT[1]: %1f\n"
		"LINPT[2]: %1f\n"
		"LINDIR[0]: %1f\n"
		"LINDIR[1]: %1f\n"
		"LINDIR[2]: %1f\n"
		"POINT[0]: %1f\n"
		"POINT[1]: %1f\n"
		"POINT[2]: %1f\n"
		"Output \n"
		"pnear[0]: %1f\n"
		"pnear[1]: %1f\n"
		"pnear[2]: %1f\n"
		"dist: %1f\n"
		"\n",
		LINPT[0], LINPT[1], LINPT[2], LINDIR[0], LINDIR[1], LINDIR[2], POINT[0], POINT[1], POINT[2], pnear[0], pnear[1], pnear[2], dist);
}

void test_surfnm(){
	double a = 1;
	double b = 2;
	double c = 3;
	vector<double> point {4, 5, 6};
	vector<double> normal(3);

	celestialObject.surfnm(a, b, c, point, normal);

	printf ("surfnm: This routine computes the outward-pointing, unit normal vector from a point on the surface of an ellipsoid. \n"
		"Input \n"
		"a: %1f\n"
		"b: %1f\n"
		"c: %1f\n"
		"point[0]: %1f\n"
		"point[1]: %1f\n"
		"point[2]: %1f\n"
		"Output \n"
		"normal[0]: %1f\n"
		"normal[1]: %1f\n"
		"normal[2]: %1f\n"
		"\n",
		a, b, c, point[0], point[1], point[2], normal[0], normal[1], normal[2]);
}

void test_surfpt(){
	vector<double> u {4, 5, 6};
	double a = 7;
	double b = 8;
	double c = 9;
	vector<double> point(3);
	int found;

	celestialObject.surfpt(u, a, b, c, point, found);

	printf ("surfpt: Determine the intersection of a line-of-sight vector with the surface of an ellipsoid. \n"
		"Input \n"
		"u[0]: %1f\n"
		"u[1]: %1f\n"
		"u[2]: %1f\n"
		"a: %1f\n"
		"b: %1f\n"
		"c: %1f\n"
		"Output \n"
		"point[0]: %1f\n"
		"point[1]: %1f\n"
		"point[2]: %1f\n"
		"found: %1f\n"
		"\n",
		u[0], u[1], u[2], a, b, c, point[0], point[1], point[2]);
}

void test_edlimb(){
	double a = 1;
	double b = 2;
	double c = 3;
	vector<double> viewpt {4, 5, 6};
	vector<double> limb (9);

	celestialObject.edlimb(a, b, c, viewpt, limb);

	printf ("edlimb: Find the limb of a triaxial ellipsoid, viewed from a specified point.  \n"
		"Input \n"
		"a: %1f\n"
		"b: %1f\n"
		"c: %1f\n"
		"viewpt[0]: %1f\n"
		"viewpt[1]: %1f\n"
		"viewpt[2]: %1f\n"
		"Output \n"
		"limb[0]: %1f\n"
		"limb[1]: %1f\n"
		"limb[2]: %1f\n"
		"limb[3]: %1f\n"
		"limb[4]: %1f\n"
		"limb[5]: %1f\n"
		"limb[6]: %1f\n"
		"limb[7]: %1f\n"
		"limb[8]: %1f\n"
		"\n",
		a, b, c, viewpt[0], viewpt[1], viewpt[2], limb[0], limb[1], limb[2], limb[3], limb[4], limb[5], limb[6], limb[7], limb[8]);
}

void test_inelpl(){
	vector<double> ellips {1, 2, 3, 4, 5, 6, 7, 8, 9};
	vector<double> plane {10, 11, 12, 13, 14, 15, 16, 17, 18};
	int nxpts = 19;
	vector<double> xpt1 {19, 20, 21};
	vector<double> xpt2 {22, 23, 24};

	celestialObject.inelpl(ellips, plane, 3, nxpts, xpt1, xpt2);

	printf ("inelpl: Find the intersection of an ellipse and a plane. \n"
		"Input \n"
		"ellips[0]: %1f\n"
		"ellips[1]: %1f\n"
		"ellips[2]: %1f\n"
		"ellips[3]: %1f\n"
		"ellips[4]: %1f\n"
		"ellips[5]: %1f\n"
		"ellips[6]: %1f\n"
		"ellips[7]: %1f\n"
		"ellips[8]: %1f\n"
		"plane[0]: %1f\n"
		"plane[1]: %1f\n"
		"plane[2]: %1f\n"
		"plane[3]: %1f\n"
		"plane[4]: %1f\n"
		"plane[5]: %1f\n"
		"plane[6]: %1f\n"
		"plane[7]: %1f\n"
		"plane[8]: %1f\n"
		"planeType: 3"
		"Output \n"
		"nxpts: %1f\n"
		"xpt1[0]: %1f\n"
		"xpt1[1]: %1f\n"
		"xpt1[2]: %1f\n"
		"xpt2[0]: %1f\n"
		"xpt2[1]: %1f\n"
		"xpt2[2]: %1f\n"
		"\n",
		ellips[0], ellips[1], ellips[2], ellips[3], ellips[4], ellips[5], ellips[6], ellips[7], ellips[8], plane[0], plane[1], plane[2], plane[3], plane[4], plane[5], plane[6], plane[7], plane[8], nxpts, xpt1[0], xpt1[1], xpt1[2], xpt2[0], xpt2[1], xpt2[2]);
}

void test_conics(){
	vector<double> elts {1, 2, 3, 4, 5, 6, 7, 8};
	double et = 9;
	vector<double> state {10, 11, 12, 13, 14, 15};

	celestialObject.conics(elts, et, state);

	printf ("conics: Determine the state (position, velocity) of an orbiting body from a set of elliptic, hyperbolic, or parabolic orbital elements.\n"
		"Input \n"
		"elts[0]: %1f\n"
		"elts[1]: %1f\n"
		"elts[2]: %1f\n"
		"elts[3]: %1f\n"
		"elts[4]: %1f\n"
		"elts[5]: %1f\n"
		"elts[6]: %1f\n"
		"elts[7]: %1f\n"
		"et: %1f\n"
		"Output \n"
		"state[0]: %1f\n"
		"state[1]: %1f\n"
		"state[2]: %1f\n"
		"state[3]: %1f\n"
		"state[4]: %1f\n"
		"state[5]: %1f\n"
		"\n",
		elts[0], elts[1], elts[2], elts[3], elts[4], elts[5], elts[6], elts[7], et, state[0], state[1], state[2], state[3], state[4], state[5]);
}

void test_oscelt(){
	vector<double> state {1, 2, 3, 4, 5, 6};
	double et = 7;
	vector<double> elts {8, 9, 10, 11, 12, 13, 14, 15};

	celestialObject.oscelt(state, et, elts);

	printf ("oscelt: Determine the set of osculating conic orbital elements that corresponds to the state (position, velocity) of a body at some epoch.\n"
		"Input \n"
		"state[0]: %1f\n"
		"state[1]: %1f\n"
		"state[2]: %1f\n"
		"state[3]: %1f\n"
		"state[4]: %1f\n"
		"state[5]: %1f\n"
		"et: %1f\n"
		"Output \n"
		"elts[0]: %1f\n"
		"elts[1]: %1f\n"
		"elts[2]: %1f\n"
		"elts[3]: %1f\n"
		"elts[4]: %1f\n"
		"elts[5]: %1f\n"
		"elts[6]: %1f\n"
		"elts[7]: %1f\n"
		"\n",
		state[0], state[1], state[2], state[3], state[4], state[5], et, elts[0], elts[1], elts[2], elts[3], elts[4], elts[5], elts[6], elts[7]);
}

int main(){
	test_npedln();
	test_pjelpl();
	test_sealgv();
	test_nplnpt();
	test_surfnm();
	test_surfpt();
	test_edlimb();
	test_inelpl();
	test_conics();
	test_oscelt();
	return 0;
}


