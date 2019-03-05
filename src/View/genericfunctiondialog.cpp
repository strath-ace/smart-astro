#include "genericfunctiondialog.h"
#include "ui_genericfunctiondialog.h"



GenericFunctionDialog::GenericFunctionDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::GenericFunctionDialog)
{
    ui->setupUi(this);

    qButtonGroup = new QButtonGroup();
    addButtonsToGroup(qButtonGroup);
    parentWindow = parent;
}

GenericFunctionDialog::~GenericFunctionDialog()
{
}

void GenericFunctionDialog::addButtonsToGroup(QButtonGroup *qButtonGroup)
{
    qButtonGroup->addButton(ui->doubleDoubleButton, 0);
    qButtonGroup->addButton(ui->doubleDoubleVectorButton, 1);
    qButtonGroup->addButton(ui->reverseStringButton, 2);
    qButtonGroup->addButton(ui->squareIntButton, 3);
    qButtonGroup->addButton(ui->squareIntVectorButton, 4);
}

void GenericFunctionDialog::on_confirmButton_clicked()
{
    if(qButtonGroup->checkedButton() != NULL){
        if(parameterInputDialog == nullptr){
            parameterInputDialog = new ParameterInputDialog(parentWindow, this, qButtonGroup->checkedId());
        }
    parameterInputDialog->setIndex(qButtonGroup->checkedId());
    this->hide();
    parameterInputDialog->show();   
    }
}

void GenericFunctionDialog::on_backButton_clicked()
{
    this->hide();
    parentWindow->show();
}
