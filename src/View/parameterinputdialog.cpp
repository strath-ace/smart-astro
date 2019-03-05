#include "parameterinputdialog.h"
#include "ui_parameterinputdialog.h"
#include <iostream>
#include <string>

ParameterInputDialog::ParameterInputDialog(QWidget *parent, QWidget *genFunctDialog, int index) :
    QDialog(parent),
    ui(new Ui::ParameterInputDialog)
{
    ui->setupUi(this);
    this->index = index;
    setInputFields();
    parentWindow = parent;
    this->genFunctDialog = genFunctDialog;
    parameterInputController = new ParameterInputController();
}

ParameterInputDialog::~ParameterInputDialog()
{
}

void ParameterInputDialog::setVisibility(int fieldNum)
{
    if (fieldNum <= 0){
        ui->label->setVisible(false);
        ui->lineEdit->setVisible(false);
    }
}


void ParameterInputDialog::setParams(std::vector<std::string> &params)
{
    if(params.size() >= 1){
        QString string;
        ui->label->setText(QString::fromStdString(params[0]));
    }
}

void ParameterInputDialog::setInputFields()
{
    switch(index) {
    case 4:
    {

        std::vector<std::string> params{"intVector"};
        resultsLabels.push_back(QString::fromStdString("intVector"));
        setParams(params);
        setVisibility(1);
    }
        break;
    case 3:
    {
        std::vector<std::string> params{"inputInt"};
        resultsLabels.push_back(QString::fromStdString("inputInt"));
        setParams(params);
        setVisibility(1);
    }
        break;
    case 2:
    {
        std::vector<std::string> params{"inputString"};
        resultsLabels.push_back(QString::fromStdString("inputString"));
        setParams(params);
        setVisibility(1);
    }
        break;
    case 1:
    {
        std::vector<std::string> params{"doubleVector"};
        resultsLabels.push_back(QString::fromStdString("doubleVector"));
        setParams(params);
        setVisibility(1);
    }
        break;
    case 0:
    {
        std::vector<std::string> params{"inputDouble"};
        resultsLabels.push_back(QString::fromStdString("inputDouble"));
        setParams(params);
        setVisibility(1);
    }
        break;
    }
}

void ParameterInputDialog::getArgs()
{
    if(ui->label->text() != NULL){
    args.push_back(ui->lineEdit->text());
    }
}

void ParameterInputDialog::setIndex(int index){
    this->index = index;
    setInputFields();
}

void ParameterInputDialog::on_calculateButton_clicked()
{
    getArgs();
    setInputFields();
    parameterInputController->calculate(args, index, conversionFlag);
    if(conversionFlag == true){
        if(answerDialog == nullptr){
            answerDialog = new AnswerDialog(parentWindow, this);
        }
        answerDialog->displayResults(args, resultsLabels);
        answerDialog->show();

        this->hide();
    } else {
        QMessageBox::warning(this, "Warning", "Incorrect parameters.");
    }
    args.clear();
    resultsLabels.clear();
}

void ParameterInputDialog::on_backButton_clicked()
{
    genFunctDialog->show();
    this->hide();
}
