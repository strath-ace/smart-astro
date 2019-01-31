#include "Astro-Core/conversion_coordinates.h"
#include <stdio.h>

int main()
{
	doublereal r = 1;
	doublereal lon = 1;
	doublereal lat = 1;
	doublereal jacobi = 1;
	
	std::cout << "" << smartastro::astrocore::conversion_coordinates::drdlat_(r, lon, lat, jacobi) << endl;

	return 0;
}
