#ifndef PARAMETERINPUTCONTROLLER_H
#define PARAMETERINPUTCONTROLLER_H

#include "../../Astro-Core/conversion_coordinates.h"
#include "../../Astro-Core/spice_general_functions.h"
#include <QObject>
#include <QWidget>
#include <vector>

#include "stdlib.h"

class ParameterInputController
{
public:
    ParameterInputController();

    ~ParameterInputController();

    QString calculate(std::vector<QString> args, int index);

private:

    double qstringToDouble(QString string);
};

#endif // PARAMETERINPUTCONTROLLER_H
