#ifndef GENERICFUNCTIONDIALOG_H
#define GENERICFUNCTIONDIALOG_H

#include <QDialog>
#include "QButtonGroup"
#include "parameterinputdialog.h"

namespace Ui {
class GenericFunctionDialog;
}

class GenericFunctionDialog : public QDialog
{
    Q_OBJECT

public:
    explicit GenericFunctionDialog(QWidget *parent = nullptr);
    ~GenericFunctionDialog();

private slots:
    void on_confirmButton_clicked();

    void on_backButton_clicked();

private:
    Ui::GenericFunctionDialog *ui;
    ParameterInputDialog *parameterInputDialog = nullptr;
    QWidget *parentWindow;
    QButtonGroup *qButtonGroup;

    void addButtonsToGroup(QButtonGroup *qButtonGroup);
};

#endif // GENERICFUNCTIONDIALOG_H

