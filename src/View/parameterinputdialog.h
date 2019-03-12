#ifndef PARAMETERINPUTDIALOG_H
#define PARAMETERINPUTDIALOG_H

#include <QDialog>
#include <QMessageBox>
#include "../Controller/parameterinputcontroller.h"
#include "answerdialog.h"

namespace Ui {
class ParameterInputDialog;
}

/**
* @brief Window for inputting parameters for calculation.
*/

class ParameterInputDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ParameterInputDialog(QWidget *parent = nullptr, QWidget *genFunctDialog = nullptr, int index = 0);
    ~ParameterInputDialog();

    /**
    * @brief Sets the function index.
    *
    * @param[in] index function index to be set.
    * @return void
    *
    */

    void setIndex(int index);

private slots:

    /**
    * @brief Performs the requested calculation, displays the results on the answer dialog and hides this window.
    *
    * @return void
    *
    */

    void on_calculateButton_clicked();

    /**
    * @brief Shows the window for selcting a function and hides this one.
    *
    * @return void
    *
    */

    void on_backButton_clicked();

private:
    Ui::ParameterInputDialog *ui;

    // determines the function selected
    int index;

    // determines if type conversion was successful
    bool conversionFlag = false;
    ParameterInputController *parameterInputController;
    AnswerDialog *answerDialog = nullptr;
    QWidget *parentWindow;
    QWidget *genFunctDialog;

    // stores arguments to be passed to model functions and the answer
    std::vector<QString> args;

    // labels for answers
    std::vector<QString> resultsLabels;

    /**
    * @brief Gets the QString values from the parameter input fields and places them in the args vector.
    *
    * @return void
    *
    */

    void getArgs();

    /**
    * @brief Sets labels for the parameters and results.
    *
    * @return void
    *
    */

    void setInputFields();

    /**
    * @brief Sets input fields and their respective labels visibilty depending on the number required.
    *
    * @param[in] fieldNum Number of input fields that are to be visible.
    * @return void
    *
    */

    void setVisibility(int fieldNum);

    /**
    * @brief Sets input parameters labels values.
    *
    * @param[in] params values to be assigned to labels.
    * @return void
    *
    */

    void setParams(std::vector<std::string> &params);
};

#endif // PARAMETERINPUTDIALOG_H
