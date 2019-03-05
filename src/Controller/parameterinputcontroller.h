#ifndef PARAMETERINPUTCONTROLLER_H
#define PARAMETERINPUTCONTROLLER_H

#include <QWidget>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <string>
#include "../Model/testModel.h"

class ParameterInputController
{
public:
    ParameterInputController();

    ~ParameterInputController();

    void calculate(std::vector<QString> &args, int index, bool &conversionFlag);

    void qstringToIntVector(QString &string, std::vector<int> &vector);

    void qstringToDoubleVector(QString &string, std::vector<double> &vector);

    void intVectorToQstring(std::vector<int> &vector, QString &string);

    void doubleVectorToQstring(std::vector<double> &vector, QString &string);
};

#endif // PARAMETERINPUTCONTROLLER_H
