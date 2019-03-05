#ifndef ANSWERDIALOG_H
#define ANSWERDIALOG_H

#include <QDialog>

namespace Ui {
class AnswerDialog;
}

class AnswerDialog : public QDialog
{
    Q_OBJECT

public:
    explicit AnswerDialog(QWidget *parent = nullptr, QWidget *paramDialog = nullptr);
    ~AnswerDialog();
    void displayResults(std::vector<QString> args, std::vector<QString> resultsLabels);

private slots:
    void on_okButton_clicked();

private:
    Ui::AnswerDialog *ui;
    QWidget *parentWindow;
    QWidget *paramDialog;

};

#endif // ANSWERDIALOG_H
