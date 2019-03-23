#define CATCH_CONFIG_MAIN

#include <vector>
#include "catch.hpp"
#include "../src/Astro-Core/spice_general_functions.h"

TEST_CASE( "Vector converted.", "[vector2array]" ) {
	std::vector<double> vector {1, 2, 3, 4, 5, 6};
	double array [2][3];

	smartastro::astrocore::spice_general_functions::doubleVectorTo2dArray3(vector, array);

	REQUIRE(array[0][0] == 1);
	REQUIRE(array[0][1] == 2);
	REQUIRE(array[0][2] == 3);
	REQUIRE(array[1][0] == 4);
	REQUIRE(array[1][1] == 5);
	REQUIRE(array[1][2] == 6);
}

TEST_CASE( "Array converted with 1 row.", "[array2vector]" ) {
	std::vector<double> vector(3);
	double array [][3] = 
	{	
		{1,2,3}
	};
	
	smartastro::astrocore::spice_general_functions::double2dArray3ToVector(array, vector);

	REQUIRE(vector[0] == 1);
	REQUIRE(vector[1] == 2);
	REQUIRE(vector[2] == 3);
}

TEST_CASE( "Array converted with 2 rows.", "[array2vector]" ) {
	std::vector<double> vector(6);
	double array [][3] = 
	{	
		{1,2,3},
		{4,5,6}
	};
	
	smartastro::astrocore::spice_general_functions::double2dArray3ToVector(array, vector);

	REQUIRE(vector[0] == 1);
	REQUIRE(vector[1] == 2);
	REQUIRE(vector[2] == 3);
	REQUIRE(vector[3] == 4);
	REQUIRE(vector[4] == 5);
	REQUIRE(vector[5] == 6);
}

TEST_CASE( "DLADvector converted and back.", "[vector2DLAD2vector]" ) {
	std::vector<int> vector {1, 2, 3, 4, 5, 6, 7, 8};
	SpiceDLADescr descr;

	smartastro::astrocore::spice_general_functions::vectorToSpiceDLADescr(vector, descr);
	smartastro::astrocore::spice_general_functions::spiceDLADescrToVector(descr, vector);

	REQUIRE(vector[0] == 1);
	REQUIRE(vector[1] == 2);
	REQUIRE(vector[2] == 3);
	REQUIRE(vector[3] == 4);
	REQUIRE(vector[4] == 5);
	REQUIRE(vector[5] == 6);
	REQUIRE(vector[6] == 7);
	REQUIRE(vector[7] == 8);
}

TEST_CASE( "DSKvector converted and back.", "[vector2DSK2vector]" ) {
	std::vector<double> vector {1, 2, 3, 4, 5, 6, 7.1, 8.1, 9.1, 10.1, 11.1, 12.1,
				 13.1, 14.1, 15.1, 16.1, 17.1, 18.1, 19.1, 20.1, 21.1, 22.1, 23.1, 24.1};
	SpiceDSKDescr descr;

	smartastro::astrocore::spice_general_functions::vectorToSpiceDSKDescr(vector, descr);
	smartastro::astrocore::spice_general_functions::spiceDSKDescrToVector(descr, vector);

	REQUIRE(vector[0] == 1);
	REQUIRE(vector[1] == 2);
	REQUIRE(vector[2] == 3);
	REQUIRE(vector[3] == 4);
	REQUIRE(vector[4] == 5);
	REQUIRE(vector[5] == 6);
	REQUIRE(vector[6] == 7.1);
	REQUIRE(vector[7] == 8.1);
	REQUIRE(vector[8] == 9.1);
	REQUIRE(vector[9] == 10.1);
	REQUIRE(vector[10] == 11.1);
	REQUIRE(vector[11] == 12.1);
	REQUIRE(vector[12] == 13.1);
	REQUIRE(vector[13] == 14.1);
	REQUIRE(vector[14] == 15.1);
	REQUIRE(vector[15] == 16.1);
	REQUIRE(vector[16] == 17.1);
	REQUIRE(vector[18] == 19.1);
	REQUIRE(vector[20] == 21.1);
	REQUIRE(vector[21] == 22.1);
	REQUIRE(vector[22] == 23.1);
	REQUIRE(vector[23] == 24.1);
}

TEST_CASE( "bodc2n tested", "[bodc2n]" ) {
	std::string name;
	int found;
	int lenout = 100;
	int code = 399;

	smartastro::astrocore::spice_general_functions::bodc2n(code, lenout, name, found);

	REQUIRE(name.compare("EARTH") == 0);
	REQUIRE(found == 1);
}

TEST_CASE( "bodn2c tested", "[bodn2c]" ) {
	std::string name = "Earth";
	int found;
	int code;

	smartastro::astrocore::spice_general_functions::bodn2c(name, code, found);

	REQUIRE(code == 399);
	REQUIRE(found == 1);
}
