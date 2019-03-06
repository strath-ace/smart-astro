#ifndef PARAMETERINPUTCONTROLLER_H
#define PARAMETERINPUTCONTROLLER_H

#include <QWidget>
#include <vector>
#include <stdlib.h>
#include <sstream>
#include <string>
#include "../Model/testModel.h"

/**
* @brief Acts as the interface between the model and the view.
*/

class ParameterInputController
{
public:
    ParameterInputController();

    ~ParameterInputController();

    /**
    * @brief Calls the function from the model the user has requested.
    *
    * @param[in/out] args arguments to be passed to the model's function.
    * @param[in] index index corresponding to the requested model function.
    * @param[in] inputString flag to determine if conversion from QString to relevant data type was successful.
    * @return void
    *
    */

    void calculate(std::vector<QString> &args, int index, bool &conversionFlag);

    /**
    * @brief Converts from a QString to a vector of integers.
    *
    * @param[in] string QString to be converted from.
    * @param[out] vector vector to be converted to.
    * @return void
    *
    */

    void qstringToIntVector(QString &string, std::vector<int> &vector);

    /**
    * @brief Converts from a QString to a vector of doubles.
    *
    * @param[in] string QString to be converted from.
    * @param[out] vector vector to be converted to.
    * @return void
    *
    */

    void qstringToDoubleVector(QString &string, std::vector<double> &vector);

    /**
    * @brief Converts from a vector of integers to a QString.
    *
    * @param[in] vector vector to be converted from.
    * @param[out] string QString to be converted to.
    * @return void
    *
    */

    void intVectorToQstring(std::vector<int> &vector, QString &string);

    /**
    * @brief Converts from a vector of doubles to a QString.
    *
    * @param[in] vector vector to be converted from.
    * @param[out] string QString to be converted to.
    * @return void
    *
    */

    void doubleVectorToQstring(std::vector<double> &vector, QString &string);
};

#endif // PARAMETERINPUTCONTROLLER_H
