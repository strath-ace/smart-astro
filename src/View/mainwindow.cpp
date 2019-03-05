#include "mainwindow.h"
#include "ui_mainwindow.h"

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
    {
        genericFunctionDialog = new GenericFunctionDialog(this);
        this->hide();
        genericFunctionDialog->show();
    }
        break;
    }
}
