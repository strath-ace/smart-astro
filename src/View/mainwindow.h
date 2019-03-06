#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "genericfunctiondialog.h"


namespace Ui {
class MainWindow;
}

/**
* @brief Window for selecting the type of calulation the user wants to carry out.
*/

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:

    /**
    * @brief Shows the function selection dialog corresponding to the type the user selected and hides this window when button is clicked.
    *
    * @return void
    *
    */

    void on_showFunctionsButton_clicked();

private:
    Ui::MainWindow *ui;
    GenericFunctionDialog *genericFunctionDialog;
};

#endif // MAINWINDOW_H
