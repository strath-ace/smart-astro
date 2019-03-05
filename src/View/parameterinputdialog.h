#ifndef PARAMETERINPUTDIALOG_H
#define PARAMETERINPUTDIALOG_H

#include <QDialog>
#include <QMessageBox>
#include "../Controller/parameterinputcontroller.h"
#include "answerdialog.h"

namespace Ui {
class ParameterInputDialog;
}

class ParameterInputDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ParameterInputDialog(QWidget *parent = nullptr, QWidget *genFunctDialog = nullptr, int index = 0);
    ~ParameterInputDialog();
    void setIndex(int index);

private slots:
    void on_calculateButton_clicked();

    void on_backButton_clicked();

private:
    Ui::ParameterInputDialog *ui;
    int index;
    bool conversionFlag = false;
    ParameterInputController *parameterInputController;
    AnswerDialog *answerDialog = nullptr;
    QWidget *parentWindow;
    QWidget *genFunctDialog;
    std::vector<QString> args;
    std::vector<QString> resultsLabels;
    void getArgs();
    void setInputFields();
    void setVisibility(int fieldNum);
    void setParams(std::vector<std::string> &params);
};

#endif // PARAMETERINPUTDIALOG_H
