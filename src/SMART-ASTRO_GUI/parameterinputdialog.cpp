#include "parameterinputdialog.h"
#include "ui_parameterinputdialog.h"

ParameterInputDialog::ParameterInputDialog(QWidget *parent, int index) :
    QDialog(parent),
    ui(new Ui::ParameterInputDialog)
{
    ui->setupUi(this);
    this->index = index;
}

ParameterInputDialog::~ParameterInputDialog()
{
    delete ui;
}

void ParameterInputDialog::setVisibility(int fieldNum)
{
    if(fieldNum <= 14){
        ui->label_24->setVisible(false);
        ui->label_23->setVisible(false);
        ui->label_22->setVisible(false);
        ui->label_21->setVisible(false);
        ui->label_20->setVisible(false);
        ui->label_19->setVisible(false);
        ui->label_18->setVisible(false);
        ui->label_17->setVisible(false);
        ui->label_16->setVisible(false);
        ui->label_15->setVisible(false);
        ui->label_14->setVisible(false);
        ui->lineEdit_24->setVisible(false);
        ui->lineEdit_23->setVisible(false);
        ui->lineEdit_22->setVisible(false);
        ui->lineEdit_21->setVisible(false);
        ui->lineEdit_20->setVisible(false);
        ui->lineEdit_19->setVisible(false);
        ui->lineEdit_18->setVisible(false);
        ui->lineEdit_17->setVisible(false);
        ui->lineEdit_16->setVisible(false);
        ui->lineEdit_15->setVisible(false);
        ui->lineEdit_14->setVisible(false);
    }
    if (fieldNum <= 7){
        ui->label_13->setVisible(false);
        ui->label_12->setVisible(false);
        ui->label_11->setVisible(false);
        ui->label_10->setVisible(false);
        ui->label_9->setVisible(false);
        ui->label_8->setVisible(false);
        ui->lineEdit_13->setVisible(false);
        ui->lineEdit_12->setVisible(false);
        ui->lineEdit_11->setVisible(false);
        ui->lineEdit_10->setVisible(false);
        ui->lineEdit_9->setVisible(false);
        ui->lineEdit_8->setVisible(false);
    }
    if (fieldNum <= 6){
        ui->label_7->setVisible(false);
        ui->lineEdit_7->setVisible(false);
    }
    if (fieldNum <= 5){
        ui->label_6->setVisible(false);
        ui->lineEdit_6->setVisible(false);
    }
    if (fieldNum <= 4){
        ui->label_5->setVisible(false);
        ui->lineEdit_5->setVisible(false);
    }
    if (fieldNum <= 3){
        ui->label_4->setVisible(false);
        ui->lineEdit_4->setVisible(false);
    }
}


void ParameterInputDialog::setParams(std::vector<std::string> &params)
{
    if(params.size() >= 3){
        QString string;
        ui->label->setText(QString::fromStdString(params[0]));
        ui->label_2->setText(QString::fromStdString(params[1]));
        ui->label_3->setText(QString::fromStdString(params[2]));

    }
    if(params.size() >= 4){
        ui->label_4->setText(QString::fromStdString(params[3]));

    }
    if(params.size() >= 5){
        ui->label_5->setText(QString::fromStdString(params[4]));

    }
    if(params.size() >= 6){
            ui->label_6->setText(QString::fromStdString(params[5]));

        }
    if(params.size() >= 7){
        ui->label_7->setText(QString::fromStdString(params[6]));

    }
    if(params.size() >= 14){
        ui->label_8->setText(QString::fromStdString(params[7]));
        ui->label_9->setText(QString::fromStdString(params[8]));
        ui->label_10->setText(QString::fromStdString(params[9]));
        ui->label_11->setText(QString::fromStdString(params[10]));
        ui->label_12->setText(QString::fromStdString(params[11]));
        ui->label_13->setText(QString::fromStdString(params[12]));
        ui->label_14->setText(QString::fromStdString(params[13]));
    }
}

void ParameterInputDialog::setInputFields()
{
    setVisibility(14);

    switch(index) {
    case 31:
    case 30:

        break;
    case 29:
        setVisibility(7);
        break;
    case 28:
    {
        std::vector<std::string> params{"r", "colat", "lons", "radius", "lon", "lat"};
        resultsLabels.push_back(QString::fromStdString("radius"));
        resultsLabels.push_back(QString::fromStdString("lon"));
        resultsLabels.push_back(QString::fromStdString("lat"));
        setParams(params);
        setVisibility(6);
    }
        break;
    case 27:
    {
        std::vector<std::string> params{"radius", "colat", "slon", "r", "lon", "z"};
        resultsLabels.push_back(QString::fromStdString("r"));
        resultsLabels.push_back(QString::fromStdString("lon"));
        resultsLabels.push_back(QString::fromStdString("z"));
        setParams(params);
        setVisibility(6);
    }
        break;
    case 26:
        setVisibility(7);
        break;
    case 25:
        setVisibility(5);
        break;
    case 24:
        setVisibility(5);
        break;
    case 23:
        setVisibility(5);
        break;
    case 22:
        setVisibility(5);
        break;
    case 21:
        setVisibility(5);
        break;
    case 20:
        setVisibility(5);
        break;
    case 19:
        setVisibility(5);
        break;
    case 18:
    {
        std::vector<std::string> params{"body", "et", "abcorr"};
        resultsLabels.push_back(QString::fromStdString("Answer"));
        setParams(params);
        setVisibility(3);
    }
        break;
    case 17:
    {
        std::vector<std::string> params{"radius", "lon", "lat", "r", "lonc", "z"};
        resultsLabels.push_back(QString::fromStdString("r"));
        resultsLabels.push_back(QString::fromStdString("lonc"));
        resultsLabels.push_back(QString::fromStdString("z"));
        setParams(params);
        setVisibility(6);
    }
        break;
    case 16:
        setVisibility(5);
        break;
    case 15:
    {
        std::vector<std::string> params{"radius", "lon", "lat", "r", "lonc", "z"};
        resultsLabels.push_back(QString::fromStdString("r"));
        resultsLabels.push_back(QString::fromStdString("lonc"));
        resultsLabels.push_back(QString::fromStdString("z"));
        setParams(params);
        setVisibility(6);
    }
        break;
    case 14:
        setVisibility(3);
        break;
    case 13:
        setVisibility(3);
        break;
    case 12:
        setVisibility(3);
        break;
    case 11:
        setVisibility(3);
        break;   
    case 10:
        setVisibility(3);
        break;
    case 9:
        setVisibility(3);
        break;
    case 8:
        setVisibility(3);
        break;
    case 7:
        setVisibility(3);
        break;
    case 6:
        setVisibility(3);
        break;
    case 5:
        setVisibility(3);
        break;
    case 4:
        setVisibility(3);
        break;
    case 3:
        setVisibility(3);
        break;
    case 2:
    {
        std::vector<std::string> params{"r", "lonc", "z", "radius", "colat", "lon"};
        resultsLabels.push_back(QString::fromStdString("radius"));
        resultsLabels.push_back(QString::fromStdString("colat"));
        resultsLabels.push_back(QString::fromStdString("lon"));
        setParams(params);
        setVisibility(6);
    }
        break;
    case 1:
        setVisibility(5);
        break;
    case 0:
    {
        std::vector<std::string> params{"r", "lonc", "z", "radius", "lon", "lat"};
        resultsLabels.push_back(QString::fromStdString("radius"));
        resultsLabels.push_back(QString::fromStdString("lon"));
        resultsLabels.push_back(QString::fromStdString("lat"));
        setParams(params);
        setVisibility(6);
    }
        break;
    }
}

void ParameterInputDialog::getArgs()
{

    if(ui->label->text() != NULL){
    args.push_back(ui->label->text());
    }
    if(ui->label_2->text() != NULL){
    args.push_back(ui->label_2->text());
    }
    if(ui->label_3->text() != NULL){
    args.push_back(ui->label_3->text());
    }
    if(ui->label_4->text() != NULL){
    args.push_back(ui->label_3->text());
    }
    if(ui->label_5->text() != NULL){
    args.push_back(ui->label_3->text());
    }
    if(ui->label_6->text() != NULL){
    args.push_back(ui->label_2->text());
    }
    if(ui->label_7->text() != NULL){
    args.push_back(ui->label_7->text());
    }
    if(ui->label_8->text() != NULL){
    args.push_back(ui->label_8->text());
    }
    if(ui->label_9->text() != NULL){
    args.push_back(ui->label_9->text());
    }
    if(ui->label_10->text() != NULL){
    args.push_back(ui->label_10->text());
    }
    if(ui->label_11->text() != NULL){
    args.push_back(ui->label_11->text());
    }
    if(ui->label_12->text() != NULL){
    args.push_back(ui->label_12->text());
    }
    if(ui->label_13->text() != NULL){
    args.push_back(ui->label_13->text());
    }
    if(ui->label_14->text() != NULL){
    args.push_back(ui->label_14->text());
    }
    if(ui->label_15->text() != NULL){
    args.push_back(ui->label_15->text());
    }
    if(ui->label_16->text() != NULL){
    args.push_back(ui->label_16->text());
    }
    if(ui->label_17->text() != NULL){
    args.push_back(ui->label_17->text());
    }
    if(ui->label_18->text() != NULL){
    args.push_back(ui->label_18->text());
    }
    if(ui->label_19->text() != NULL){
    args.push_back(ui->label_19->text());
    }
    if(ui->label_20->text() != NULL){
    args.push_back(ui->label_20->text());
    }
    if(ui->label_21->text() != NULL){
    args.push_back(ui->label_21->text());
    }
}

void ParameterInputDialog::on_pushButton_clicked()
{
    answerDialog.show();
    parameterInputController.calculate(args, index);
    answerDialog.displayResults(args, resultsLabels);
    args.clear();
    resultsLabels.clear();
    this->close();
}
