#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "genericfunctiondialog.h"


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:

    void on_showFunctionsButton_clicked();

private:
    Ui::MainWindow *ui;
    GenericFunctionDialog *genericFunctionDialog;
};

#endif // MAINWINDOW_H
