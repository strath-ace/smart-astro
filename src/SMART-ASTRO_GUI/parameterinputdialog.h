#ifndef PARAMETERINPUTDIALOG_H
#define PARAMETERINPUTDIALOG_H

#include <QDialog>
#include "Controller/parameterinputcontroller.h"
#include "ui_parameterinputdialog.h"

namespace Ui {
class ParameterInputDialog;
}

class ParameterInputDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ParameterInputDialog(QWidget *parent = nullptr, int index = 0);
    ~ParameterInputDialog();

private slots:
    void on_pushButton_clicked();

private:
    Ui::ParameterInputDialog *ui;
    int index;
    ParameterInputController parameterInputController;
    std::vector<QString> getArgs();
    void setInputFields();
    void setVisibility(int fieldNum);
    void setParams(std::vector<std::string> &params);
};

#endif // PARAMETERINPUTDIALOG_H
