#include "parameterinputcontroller.h"
#include <iostream>
#include <string>

ParameterInputController::ParameterInputController()
{
}

ParameterInputController::~ParameterInputController()
{

}

void ParameterInputController::calculate(std::vector<QString> &args, int index, bool &conversionFlag)
{
    switch (index) {
        case 0:
        {
            double bubble = args[0].toDouble(&conversionFlag);
            TestModel::doubleDouble(bubble);
            args[0] = QString::number(bubble);
           }
            break;
        case 1:
            {
            std::vector<double> vector;
            qstringToDoubleVector(args[0], vector);
            if(vector.size() > 0){
                TestModel::doubleDoubleVector(vector);
                doubleVectorToQstring(vector, args[0]);
                conversionFlag = true;
            } else {
                conversionFlag = false;
            }
            }
            break;
        case 2:
        {
            std::string string = args[0].toStdString();
            TestModel::reverseString(string);
            args[0] = QString::fromStdString(string);
            conversionFlag = true;
        }
            break;
        case 3:
        {
            int integer = args[0].toInt(&conversionFlag);
            TestModel::squareInt(integer);
            args[0] = QString::number(integer);
        }
            break;
        case 4:
        {
            std::vector<int> vector;
            qstringToIntVector(args[0], vector);
            if(vector.size() > 0){
                TestModel::squareIntVector(vector);
                intVectorToQstring(vector, args[0]);
                conversionFlag = true;
            } else {
                conversionFlag = false;
            }
        }
            break;
    }
}

void ParameterInputController::qstringToIntVector(QString &string, std::vector<int> &vector)
{
    std::string stringy = string.toStdString();
    std::stringstream iss(stringy);

    int number;
    while ( iss >> number )
      vector.push_back( number );
}

void ParameterInputController::qstringToDoubleVector(QString &string, std::vector<double> &vector)
{
    std::string stringy = string.toStdString();
    std::stringstream iss(stringy);

    double number;
    while ( iss >> number )
      vector.push_back( number );
}

void ParameterInputController::intVectorToQstring(std::vector<int> &vector, QString &string)
{
    std::string stringer;
    stringer = stringer + std::to_string(vector[0]);

    for(int i = 1; i < vector.size(); ++i){
        stringer = stringer + ' ' + std::to_string(vector[i]);
    }

    string = QString::fromStdString(stringer);
}

void ParameterInputController::doubleVectorToQstring(std::vector<double> &vector, QString &string)
{
    std::string stringer;
    stringer = stringer + std::to_string(vector[0]);

    for(int i = 1; i < vector.size(); ++i){
        stringer = stringer + ' ' + std::to_string(vector[i]);
    }

    string = QString::fromStdString(stringer);
}




