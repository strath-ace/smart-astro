#include "testModel.h"

void TestModel::reverseString(std::string &inputString){
        int n = inputString.length();

        for (int i = 0; i < n / 2; i++)
            std::swap(inputString[i], inputString[n - i - 1]);
}

void TestModel::doubleDouble(double &inputDouble){
	inputDouble = inputDouble * 2;
}

void TestModel::squareInt(int &inputInt){
	int temp = inputInt * inputInt;
	inputInt = temp;
}

void TestModel::doubleDoubleVector(std::vector<double> &doubleVector){
    for(int i = 0; i < doubleVector.size(); ++i){
       doubleVector[i]  = doubleVector[i] * 2;
    }
}

void TestModel::squareIntVector(std::vector<int> &intVector){
    for(int i = 0; i < intVector.size(); ++i){
            int temp = intVector[i] * intVector[i];
            intVector[i] = temp;
    }
}
