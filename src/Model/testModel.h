#ifndef TESTMODEL_H
#define TESTMODEL_H

#include <vector>
#include <string>

class TestModel
{
    public:
        static void reverseString(std::string &inputString);

        static void doubleDouble(double &inputDouble);

        static void squareInt(int &inputInt);

    	static void doubleDoubleVector(std::vector<double> &doubleVector);

        static void squareIntVector(std::vector<int> &intVector);
};

#endif /* TESTMODEL_H */
