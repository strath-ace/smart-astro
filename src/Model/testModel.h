#ifndef TESTMODEL_H
#define TESTMODEL_H

#include <vector>
#include <string>

/**
* @brief Simple demonstration model.
*/

class TestModel
{
    public:

        /**
        * @brief Reverses a string.
        *
        * @param[in/out] inputString String to be reversed.
        * @return void
        *
        */

        static void reverseString(std::string &inputString);

        /**
        * @brief Doubles a double.
        *
        * @param[in/out] inputDouble Double to be doubled.
        * @return void
        *
        */

        static void doubleDouble(double &inputDouble);

        /**
        * @brief Squares an integer.
        *
        * @param[in/out] inputInt Integer to be squared.
        * @return void
        *
        */

        static void squareInt(int &inputInt);

        /**
        * @brief doubles every element of a vector of doubles.
        *
        * @param[in/out] doubleVector vector whose elements are to be doubled.
        * @return void
        *
        */

    	static void doubleDoubleVector(std::vector<double> &doubleVector);

        /**
        * @brief squares every element of a vector of integers.
        *
        * @param[in/out] intVector vector whose elements are to be squared.
        * @return void
        *
        */

        static void squareIntVector(std::vector<int> &intVector);
};

#endif /* TESTMODEL_H */
