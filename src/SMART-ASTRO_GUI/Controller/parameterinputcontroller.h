#ifndef PARAMETERINPUTCONTROLLER_H
#define PARAMETERINPUTCONTROLLER_H

#include <QWidget>
#include <vector>

class ParameterInputController
{
public:
    ParameterInputController();

    ~ParameterInputController();

    void calculate(std::vector<QString> &args, int index);
};

#endif // PARAMETERINPUTCONTROLLER_H
