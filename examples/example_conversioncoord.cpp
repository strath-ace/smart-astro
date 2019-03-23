#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <string>

#include "smartastro.h"
#include "../src/Astro-Core/conversion_coordinates.h"
#include "catch.hpp"

using namespace std;
using namespace smartastro;

void test_radrec(){
	double range = 1;
	double ra = 2;
	double dec = 3;
	vector<double> rectan(3);

	smartastro::astrocore::conversion_coordinates::radrec(range, ra, dec, rectan);

	printf ("radrec: Convert from range, right ascension, and declination to rectangular coordinates.\n"
		"Input \n"
		"range: %1f\n"
		"ra: %1f\n"
		"dec: %1f\n"
	     	"Output\n"
           	"rectan[0]: %1f\n"
             	"rectan[1]: %1f\n"
	     	"rectan[2]: %1f\n"
		"\n",
             	range, ra, dec, rectan[0], rectan[1], rectan[2]);
}

void test_recrad(){
	vector<double> rectan{ 1, 2, 3 };
	double range;
	double ra;
	double dec;

	smartastro::astrocore::conversion_coordinates::recrad(rectan, range, ra, dec);

	printf ("recrad: Convert rectangular coordinates to range, right ascension, and declination.\n"
		"Input \n"
		"rectan[0]: %1f\n"
             	"rectan[1]: %1f\n"
	     	"rectan[2]: %1f\n"
	     	"Output\n"
           	"range: %1f\n"
		"ra: %1f\n"
		"dec: %1f\n"
		"\n",
             	rectan[0], rectan[1], rectan[2], range, ra, dec);
}

void test_sphlat(){
	double r = 1;
	double colat = 2;
	double lons = 3;
	double radius;
	double lon;
	double lat;
	
	smartastro::astrocore::conversion_coordinates::sphlat(r, colat, lons, radius, lon, lat);

	printf ("sphlat: Convert from spherical coordinates to latitudinal coordinates.\n"
		"Input \n"
		"r: %1f\n"
             	"colat: %1f\n"
	     	"lons: %1f\n"
	     	"Output\n"
           	"radius: %1f\n"
		"lon: %1f\n"
		"lat: %1f\n"
		"\n",
             	r, colat,lons, radius, lon, lat);
}

void test_latsph(){
	double radius = 1;
	double lon = 2;
	double lat = 3;
	double rho;
	double colat;
	double lons;

	smartastro::astrocore::conversion_coordinates::latsph(radius, lon, lat, rho, colat, lons);

	printf ("latsph: Convert from latitudinal coordinates to spherical coordinates.\n"
		"Input \n"
		"radius: %1f\n"
		"lon: %1f\n"
		"lat: %1f\n"
		"Output\n"
		"rho: %1f\n"
             	"colat: %1f\n"
	     	"lons: %1f\n"	
		"\n",
             	radius, lon, lat, rho, colat, lons);
}

void test_sphcyl(){
	double radius = 1;
	double colat = 2;
	double slon = 3;
	double r;
	double lon;
	double z;

	smartastro::astrocore::conversion_coordinates::sphcyl(radius, colat, slon, r, lon, z);

	printf ("sphcyl: This routine converts from spherical coordinates to cylindrical coordinates.\n"
		"Input \n"
		"radius: %1f\n"
		"colat: %1f\n"
		"slon: %1f\n"
		"Output\n"
		"r: %1f\n"
		"lon: %1f\n"
		"z: %1f\n"	
		"\n",
             	radius, colat, slon, r, lon, z);
}

void test_latcyl(){
	double radius = 1;
	double lon = 2;
	double lat = 3;
	double r;
	double lonc;
	double z;

	smartastro::astrocore::conversion_coordinates::latcyl(radius, lon, lat, r, lonc, z);

	printf ("latcyl: Convert from latitudinal coordinates to cylindrical coordinates.\n"
		"Input \n"
		"radius: %1f\n"
		"lon: %1f\n"
		"lat: %1f\n"
		"Output\n"
		"r: %1f\n"
		"lonc: %1f\n"
		"z: %1f\n"	
		"\n",
             	radius, lon, lat, r, lonc, z);
}

int main(){
	test_radrec();
	test_recrad();
	test_sphlat();
	test_latsph();
	test_sphcyl();
	test_latcyl();
	return 0;
}
