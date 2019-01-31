#include "parameterinputcontroller.h"

ParameterInputController::ParameterInputController()
{

}

ParameterInputController::~ParameterInputController()
{

}

doublereal *qstringToDoublerealP(QString string)
{
    std::string arg1 = string.toUtf8().constData();
    doublereal arg1d = atof(arg1.c_str());
    doublereal *arg1dp = &arg1d;
    return arg1dp;
}

char *qstringToCharP(QString string)
{
    std::string arg1 = string.toUtf8().constData();
    std::vector<char> writable(arg1.begin(), arg1.end());
    writable.push_back('\0');
    return &writable[0];
}

ftnlen qstringToFtnlen(QString string)
{
    int convertedInt = std::stoi(string);
    return convertedInt;
}

QString calculate(std::vector<QString> args, int index)
{
    QString s;
    int result;
    doublereal doublerealResult;

    switch(index) {
    case 0:
        result = smartastro::astrocore::conversion_coordinates::cyllat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 1:
        result = smartastro::astrocore::conversion_coordinates::cylrec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 2:
        result = smartastro::astrocore::conversion_coordinates::cylsph_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 3:
        result = smartastro::astrocore::conversion_coordinates::dcyldr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 4:
        result = smartastro::astrocore::conversion_coordinates::dgeodr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 5:
        result = smartastro::astrocore::conversion_coordinates::dlatdr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]),stringToDoublerealP(args[3]));
        break;
    case 7:
        result = smartastro::astrocore::conversion_coordinates::dpgrdr_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), stringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));
        break;
    case 8:
        result = smartastro::astrocore::conversion_coordinates::drdcyl_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 9:
        result = smartastro::astrocore::conversion_coordinates::drdgeo_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]),qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 10:
        result = smartastro::astrocore::conversion_coordinates::drdlat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 11:
        result = smartastro::astrocore::conversion_coordinates::drdpgr_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), stringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));
        break;
    case 12:
        result = smartastro::astrocore::conversion_coordinates::drdsph_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 13:
        result = smartastro::astrocore::conversion_coordinates::dsphdr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), stringToDoublerealP(args[3]));
        break;
    case 14:
        result = smartastro::astrocore::conversion_coordinates::georec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 15:
        result = smartastro::astrocore::conversion_coordinates::latcyl_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 16:
        result = smartastro::astrocore::conversion_coordinates::latrec_(stringToDoublerealP(args[0]),qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 17:
        result = smartastro::astrocore::conversion_coordinates::latsph_(stringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 18:
        doublerealResult = smartastro::astrocore::spice_general_functions::lspcn_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToCharP(args[2]), qstringToFtnlen(args[3]), qstringToFtnlen(args[4]));
        break;
    case 19:
        result = smartastro::astrocore::conversion_coordinates::pgrrec_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));
        break;
    case 20:
        result = smartastro::astrocore::conversion_coordinates::radrec_(stringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 21:
        result = smartastro::astrocore::conversion_coordinates::reccyl_(stringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 22:
        result = smartastro::astrocore::conversion_coordinates::recgeo_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 23:
        result = smartastro::astrocore::conversion_coordinates::reclat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 24:
        result = smartastro::astrocore::conversion_coordinates::recpgr_(qstringToCharP(args[0]),  qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));
        break;
    case 25:
        result = smartastro::astrocore::conversion_coordinates::recrad_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 26:
        result = smartastro::astrocore::conversion_coordinates::recsph_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 27:
        result = smartastro::astrocore::conversion_coordinates::sphcyl_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 28:
        result = smartastro::astrocore::conversion_coordinates::sphlat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));
        break;
    case 29:
        result = smartastro::astrocore::conversion_coordinates::sphrec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));
        break;
    case 30:
        result = smartastro::astrocore::spice_general_functions::subpnt_(qstringToCharP(args[0]), qstringToCharP(args[1]), qstringToDoublerealP(args[2]), qstringToCharP(args[3]), qstringToCharP(args[4]), qstringToCharP(args[5]), qstringToDoublerealP(args[6]), qstringToDoublerealP(args[7]), qstringToDoublerealP(args[8]), stringToFtnlen(args[9]), stringToFtnlen(args[10]), stringToFtnlen(args[11]), stringToFtnlen(args[12]), stringToFtnlen(args[13]));
        break;
    case 31:
        result = smartastro::astrocore::spice_general_functions::subslr_(qstringToCharP(args[0]), qstringToCharP(args[1]), qstringToDoublerealP(args[2]), qstringToCharP(args[3]), qstringToCharP(args[4]), qstringToCharP(args[5]), qstringToDoublerealP(args[6]), qstringToDoublerealP(args[7]), qstringToDoublerealP(args[8]), stringToFtnlen(args[9]), stringToFtnlen(args[10]), stringToFtnlen(args[11]), stringToFtnlen(args[12]), stringToFtnlen(args[13]));
        break;

        if(index == 18) {
            s = QString::number(result);
        }
        else {
            s = QString::number(doublerealResult);
        }
        return s;
    }
}


