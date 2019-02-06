/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#include "Astro-Core/spice_general_functions.h"

using namespace smartastro;
using namespace smartastro::astrocore;

char* stringToCharArray(std::string string)
{
   char *charArray = string.c_string();
  return charArray;
}

std::vector<int> intArrayToVector(int array[])
{
  std::vector<int> vector(array, array + sizeof array / sizeof array[0]);
  return vector;
}

int* intVectorToArray(std::vector<int> vector)
{
  int* arrayPointer = &vector[0];
  return arrayPointer;
}

std::vector<char> charArrayToVector(char array[])
{
  std::vector<char> vector(array, array + sizeof array / sizeof array[0]);
  return vector;
}

char* charVectorToArray(std::vector<char> vector)
{
  char* arrayPointer = &vector[0];
  return arrayPointer;
}

std::vector<double> doubleArrayToVector(double array[])
{
  std::vector<double> vector(array, array + sizeof array / sizeof array[0]);
  return vector;
}

double* intVectorToArray(std::vector<double> vector)
{
  double* arrayPointer = &vector[0];
  return arrayPointer;
}

int spice_general_functions::subslr_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len)
{
	return subslr_(method, target, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

int spice_general_functions::subpnt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len)
{
	return subpnt_(method, target, et, fixref, abcorr, obsrvr, spoint, trgepc, srfvec, method_len, target_len, fixref_len, abcorr_len, obsrvr_len);
}

double spice_general_functions::lspcn_( std::string &body,  std::vector<double> &et,  std::string &abcorr, int body_len, int abcorr_len)
{
	return lspcn_(body, et, abcorr, body_len, abcorr_len);
}

int spice_general_functions::furnsh_( std::string &file, int file_len)
{
	return furnsh_(file, file_len);
}

int spice_general_functions::unload_( std::string &file, int file_len)
{
	return unload_(file, file_len);
}
