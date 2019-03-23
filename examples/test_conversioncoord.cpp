#define CATCH_CONFIG_MAIN

#include <vector>
#include "catch.hpp"
#include "../src/Astro-Core/conversion_coordinates.h"

TEST_CASE( "cylrec tested", "[cylrec]" ) {
	double r = 1.4142;
	double lon = 0.785398;
	double z = 1.0000;
	std::vector<double> rectan(3);

	smartastro::astrocore::conversion_coordinates::cylrec(r, lon, z, rectan);

	REQUIRE(std::fabs(rectan[0] - 1.0000) < 0.01);
	REQUIRE(std::fabs(rectan[1] - 1.0000) < 0.01);
	REQUIRE(std::fabs(rectan[1] - 1.0000) < 0.01);
}

TEST_CASE( "reccyl tested", "[reccyl]" ) {
	double r;
	double lon;
	double z;
	std::vector<double> rectan{1.0000, 1.0000, 1.0000};

	smartastro::astrocore::conversion_coordinates::reccyl(rectan, r, lon, z);

	REQUIRE(std::fabs(r == 1.4142) < 0.01);
	REQUIRE(std::fabs(lon - 0.785398) < 0.01);
	REQUIRE(std::fabs(z - 1.0000) < 0.01);
}

TEST_CASE( "sphrec tested", "[sphrec]" ) {
	double r = 1.7320;
	double colat = 0.9553164381;
	double lon = 0.785398;
	std::vector<double> rectan(3);

	smartastro::astrocore::conversion_coordinates::sphrec(r, colat, lon, rectan);

	REQUIRE(std::fabs(rectan[0] - 1.0000) < 0.01);
	REQUIRE(std::fabs(rectan[1] - 1.0000) < 0.01);
	REQUIRE(std::fabs(rectan[1] - 1.0000) < 0.01);
}

TEST_CASE( "recsph tested", "[recsph]" ) {
	double r;
	double colat;
	double lon;
	std::vector<double> rectan {1, 1, 1};

	smartastro::astrocore::conversion_coordinates::recsph(rectan, r, colat, lon);

	REQUIRE(std::fabs(r - 1.7320) < 0.01);
	REQUIRE(std::fabs(colat - 0.9553164381) < 0.01);
	REQUIRE(std::fabs(lon - 0.785398) < 0.01);
}

TEST_CASE( "cylsph tested", "[cylsph]" ) {
	double r = 1;
	double lonc = 3.14159;
	double z = -1.000;
	double radius;
	double colat;
	double lon;

	smartastro::astrocore::conversion_coordinates::cylsph(r, lonc, z, radius, colat, lon);

	REQUIRE(std::fabs(radius - 1.4142) < 0.01);
	REQUIRE(std::fabs(lon - 3.14159) < 0.01);
	REQUIRE(std::fabs(colat - 2.35619) < 0.01);
}

TEST_CASE( "cyllat tested", "[cyllat]" ) {
	double r = 1.0000;
	double lonc = 3.14159;
	double z = 1.000;
	double radius;
	double lon;
	double lat;

	smartastro::astrocore::conversion_coordinates::cyllat(r, lonc, z, radius, lon, lat);

	REQUIRE(std::fabs(radius - 1.4142) < 0.01);
	REQUIRE(std::fabs(lon - 3.14159) < 0.01);
	REQUIRE(std::fabs(lat + -0.785398) < 0.01);
}

TEST_CASE( "reclat tested", "[reclat]" ) {
	std::vector<double> rectan {1, 1, 1};
	double r;
	double lon;
	double lat;

	smartastro::astrocore::conversion_coordinates::reclat(rectan, r, lon, lat);

	REQUIRE(std::fabs(r - 1.7320) < 0.01);
	REQUIRE(std::fabs(lon - 0.785397999997775) < 0.01);
	REQUIRE(std::fabs(lat - 0.6154781434) < 0.01);
}

TEST_CASE( "latrec tested", "[latrec]" ) {
	std::vector<double> rectan (3);
	double r = 1.7320;
	double lon = 0.785397999997775;
	double lat = 0.6154781434;

	smartastro::astrocore::conversion_coordinates::latrec(r, lon, lat, rectan);

	REQUIRE(std::fabs(rectan[0] == 1) < 0.01);
	REQUIRE(std::fabs(rectan[1] == 1) < 0.01);
	REQUIRE(std::fabs(rectan[2] == 1) < 0.01);
}
