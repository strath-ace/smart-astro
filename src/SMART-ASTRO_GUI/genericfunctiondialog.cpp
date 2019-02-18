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
    delete ui;
}

void GenericFunctionDialog::addButtonsToGroup(QButtonGroup *qButtonGroup)
{
    qButtonGroup->addButton(ui->cyllat_Button, 0);
    qButtonGroup->addButton(ui->cylrec_Button, 1);
    qButtonGroup->addButton(ui->cylsph_Button, 2);
    qButtonGroup->addButton(ui->dcyldr_Button, 3);
    qButtonGroup->addButton(ui->dgeodr_Button, 4);
    qButtonGroup->addButton(ui->dlatdr_Button, 5);
    qButtonGroup->addButton(ui->dpgrdr_Button, 6);
    qButtonGroup->addButton(ui->drdcyl_Button, 7);
    qButtonGroup->addButton(ui->drdgeo_Button, 8);
    qButtonGroup->addButton(ui->drdlat_Button, 9);
    qButtonGroup->addButton(ui->drdpgr_Button, 10);
    qButtonGroup->addButton(ui->drdsph_Button, 11);
    qButtonGroup->addButton(ui->dsphdr_Button, 12);
    qButtonGroup->addButton(ui->furnsh_Button, 13);
    qButtonGroup->addButton(ui->georec_Button, 14);
    qButtonGroup->addButton(ui->latcyl_Button, 15);
    qButtonGroup->addButton(ui->latrec_Button, 16);
    qButtonGroup->addButton(ui->latsph_Button, 17);
    qButtonGroup->addButton(ui->pgrrec_Button, 18);
    qButtonGroup->addButton(ui->radrec_Button, 19);
    qButtonGroup->addButton(ui->reccyl_Button, 20);
    qButtonGroup->addButton(ui->recgeo_Button, 21);
    qButtonGroup->addButton(ui->reclat_Button, 22);
    qButtonGroup->addButton(ui->recpgr_Button, 23);
    qButtonGroup->addButton(ui->recrad_Button, 24);
    qButtonGroup->addButton(ui->recsph_Button, 25);
    qButtonGroup->addButton(ui->sphcyl_Button, 26);
    qButtonGroup->addButton(ui->sphlat_Button, 27);
    qButtonGroup->addButton(ui->sphrec_Button, 28);
    qButtonGroup->addButton(ui->subpnt_Button, 29);
    qButtonGroup->addButton(ui->subslr_Button, 30);
    qButtonGroup->addButton(ui->unload_Button, 31);
    qButtonGroup->addButton(ui->xfmsta_Button, 31);
}

void GenericFunctionDialog::on_confirmButton_clicked()
{
    if(qButtonGroup->checkedButton() != NULL){
    parameterInputDialog = new ParameterInputDialog(parentWindow, qButtonGroup->checkedId());
    parameterInputDialog->show();
    this->close();
    }
}
