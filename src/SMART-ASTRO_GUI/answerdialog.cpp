#include "answerdialog.h"
#include "ui_answerdialog.h"

AnswerDialog::AnswerDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::AnswerDialog)
{
    ui->setupUi(this);
}

AnswerDialog::~AnswerDialog()
{
    delete ui;
}

void AnswerDialog::displayResults(std::vector<QString> args, std::vector<QString> resultsLabels){
    ui->answerLabel1->setText(resultsLabels[0]);
    ui->answerLabel2->setText(resultsLabels[1]);
    ui->answerLabel3->setText(resultsLabels[2]);

    ui->answerField1->setText(args[0]);
    ui->answerField2->setText(args[1]);
    ui->answerField3->setText(args[2]);
}

void AnswerDialog::on_okButton_clicked()
{
    this->close();
}
