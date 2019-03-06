#include "answerdialog.h"
#include "ui_answerdialog.h"
#include <iostream>
#include <string>

AnswerDialog::AnswerDialog(QWidget *parent, QWidget *paramDialog) :
    QDialog(parent),
    ui(new Ui::AnswerDialog)
{
    ui->setupUi(this);
    parentWindow = parent;
    this->paramDialog = paramDialog;
}

AnswerDialog::~AnswerDialog()
{
    delete ui;
}

void AnswerDialog::displayResults(std::vector<QString> args, std::vector<QString> resultsLabels){
    if(args.size() >= 1){
        ui->answerLabel->setText(resultsLabels[0]);
        ui->answerField->setText(args[0]);
    }
}

void AnswerDialog::on_okButton_clicked()
{
 this->paramDialog->show();
 this->hide();
}
