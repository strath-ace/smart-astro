#include "parameterinputdialog.h"

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
    if(index <= 14){
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
    if (index <= 7){
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
    if (index <= 5){
        ui->label_7->setVisible(false);
        ui->label_6->setVisible(false);
        ui->lineEdit_7->setVisible(false);
        ui->lineEdit_6->setVisible(false);
    }
    if (index <= 4){
        ui->label_5->setVisible(false);
        ui->lineEdit_5->setVisible(false);
    }
    if (index <= 3){
        ui->label_4->setVisible(false);
        ui->lineEdit_4->setVisible(false);
    }
}


void ParameterInputDialog::setParams(int &params)
{
    if(*params.size() >= 3){
        ui->label->setText(args[0]);
        ui->label_2->setText(args[1]);
        ui->label_3->setText(args[2]);

    }
    if(args.size() >= 4){
        ui->label_4->setText(args[3]);

    }
    if(args.size() >= 5){
        ui->label_5->setText(args[4]);

    }
    if(args.size() >= 7){
        ui->label_6->setText(args[5]);
        ui->label_7->setText(args[6]);

    }
    if(args.size() >= 14){
        ui->label_8->setText(args[7]);
        ui->label_9->setText(args[8]);
        ui->label_10->setText(args[9]);
        ui->label_11->setText(args[10]);
        ui->label_12->setText(args[11]);
        ui->label_13->setText(args[12]);
        ui->label_14->setText(args[13]);
    }
}

void ParameterInputDialog::setInputFields()
{
    std::vector<std::string> args;
    setVisibility(14);

    switch(index) {
    case 31:
    case 30:
        args = {"method:", "target:", "et:", "fixref:", "abcorr:", "obsrvr:", "spoint:", "trgepc:", "srfvec:", "method_len:", "target_len:", "fixref_len:", "abcorr_len:", "obsrvr_len:"};

        break;
    case 24:
        setVisibility(7);
        ui->label->setText("body:");
        ui->label_2->setText("target:");
        ui->label_3->setText("et:");
        ui->label_4->setText("fixref:");
        ui->label_5->setText("abcorr:");
        ui->label_6->setText("obsrvr:");
        ui->label_7->setText("spoint:");
        break;
    case 27:
        setVisibility(7);
        break;
    case 26:
        setVisibility(7);
        break;
    case 25:
        setVisibility(5);
        break;
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
        setVisibility(5);
        break;
    case 17:
        setVisibility(5);
        break;
    case 16:
        setVisibility(5);
        break;
    case 15:
        setVisibility(4);
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
        setVisibility(3);
        break;
    case 1:
        setVisibility(3);
        break;
    case 0:
        setVisibility(3);
        break;
    }
}

std::vector<QString> ParameterInputDialog::getArgs()
{
    std::vector<QString> vector;

    if(ui->label->text() != NULL){
    vector.push_back(ui->label->text());
    }
    if(ui->label_2->text() != NULL){
    vector.push_back(ui->label_2->text());
    }
    if(ui->label_3->text() != NULL){
    vector.push_back(ui->label_3->text());
    }
    if(ui->label_4->text() != NULL){
    vector.push_back(ui->label_3->text());
    }
    if(ui->label_5->text() != NULL){
    vector.push_back(ui->label_3->text());
    }
    if(ui->label_6->text() != NULL){
    vector.push_back(ui->label_2->text());
    }
    if(ui->label_7->text() != NULL){
    vector.push_back(ui->label_7->text());
    }
    if(ui->label_8->text() != NULL){
    vector.push_back(ui->label_8->text());
    }
    if(ui->label_9->text() != NULL){
    vector.push_back(ui->label_9->text());
    }
    if(ui->label_10->text() != NULL){
    vector.push_back(ui->label_10->text());
    }
    if(ui->label_11->text() != NULL){
    vector.push_back(ui->label_11->text());
    }
    if(ui->label_12->text() != NULL){
    vector.push_back(ui->label_12->text());
    }
    if(ui->label_13->text() != NULL){
    vector.push_back(ui->label_13->text());
    }
    if(ui->label_14->text() != NULL){
    vector.push_back(ui->label_14->text());
    }
    if(ui->label_15->text() != NULL){
    vector.push_back(ui->label_15->text());
    }
    if(ui->label_16->text() != NULL){
    vector.push_back(ui->label_16->text());
    }
    if(ui->label_17->text() != NULL){
    vector.push_back(ui->label_17->text());
    }
    if(ui->label_18->text() != NULL){
    vector.push_back(ui->label_18->text());
    }
    if(ui->label_19->text() != NULL){
    vector.push_back(ui->label_19->text());
    }
    if(ui->label_20->text() != NULL){
    vector.push_back(ui->label_20->text());
    }
    if(ui->label_21->text() != NULL){
    vector.push_back(ui->label_21->text());
    }

    return vector;
}

void ParameterInputDialog::on_pushButton_clicked()
{
    std::vector<QString> argsVector = getArgs();
    ui->answerLine->setText(parameterInputController.calculate(argsVector, index));
}
