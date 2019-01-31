#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <string>
#include <iostream>

using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_showFunctionsButton_clicked()
{
    int index = ui->functionTypeBox->currentIndex();

    switch(index) {
    case 0:
        break;
    case 1:
    {
        genericFunctionDialog = new  GenericFunctionDialog(this);
        genericFunctionDialog->show();
    }
        break;
    case 2:
        break;
    case 3:
        break;
    case 4:
        break;
    case 5:
        break;
    case 7:
        break;
    }
}
