#ifndef ANSWERDIALOG_H
#define ANSWERDIALOG_H

#include <QDialog>

namespace Ui {
class AnswerDialog;
}

/**
* @brief Window for Displaying the output of a calculation.
*/

class AnswerDialog : public QDialog
{
    Q_OBJECT

public:
    explicit AnswerDialog(QWidget *parent = nullptr, QWidget *paramDialog = nullptr);
    ~AnswerDialog();

    /**
    * @brief Displays the output of the calculation.
    *
    * @param[in] args Output from the model's function.
    * @param[in] resultsLabels Labels for the output.
    * @return void
    *
    */

    void displayResults(std::vector<QString> args, std::vector<QString> resultsLabels);

private slots:

    /**
    * @brief Shows the parameter input dialog and hides this window when button is clicked.
    *
    * @return void
    *
    */

    void on_okButton_clicked();

private:
    Ui::AnswerDialog *ui;
    QWidget *parentWindow;
    QWidget *paramDialog;

};

#endif // ANSWERDIALOG_H
