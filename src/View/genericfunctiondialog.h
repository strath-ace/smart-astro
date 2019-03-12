#ifndef GENERICFUNCTIONDIALOG_H
#define GENERICFUNCTIONDIALOG_H

#include <QDialog>
#include "QButtonGroup"
#include "parameterinputdialog.h"

namespace Ui {
class GenericFunctionDialog;
}

/**
* @brief Window for selecting a function.
*/

class GenericFunctionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit GenericFunctionDialog(QWidget *parent = nullptr);
    ~GenericFunctionDialog();

private slots:

    /**
    * @brief If a button has been checked, the parameter input dialog is shown and this window is hidden.
    *
    * @return void
    *
    */

    void on_confirmButton_clicked();

    /**
    * @brief Shows the main window and hides this one.
    *
    * @return void
    *
    */

    void on_backButton_clicked();

private:
    Ui::GenericFunctionDialog *ui;
    ParameterInputDialog *parameterInputDialog = nullptr;
    QWidget *parentWindow;

    //container for buttons for easy access
    QButtonGroup *qButtonGroup;

    /**
    * @brief Displays the output of the calculation.
    *
    * @param[in] qButtonGroup ButtonGroup buttons are to be added to.
    * @return void
    *
    */

    void addButtonsToGroup(QButtonGroup *qButtonGroup);
};

#endif // GENERICFUNCTIONDIALOG_H

