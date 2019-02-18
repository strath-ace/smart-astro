#include "parameterinputcontroller.h"
#include "../../Astro-Core/spice_general_functions.h"
#include "../../Astro-Core/conversion_coordinates.h"


ParameterInputController::ParameterInputController()
{

}

ParameterInputController::~ParameterInputController()
{

}

void ParameterInputController::calculate(std::vector<QString> &args, int index)
{
    switch(index) {
    case 0:
    {
        double r = args[0].toDouble();
        double lonc = args[1].toDouble();
        double z = args[2].toDouble();
        double radius = args[3].toDouble();
        double lon = args[4].toDouble();
        double lat = args[5].toDouble();

        smartastro::astrocore::conversion_coordinates::cyllat(r, lonc, z, radius, lon, lat);

        args[0] = QString::number(r);
        args[1] = QString::number(lonc);
        args[2] = QString::number(z);
        args[3] = QString::number(radius);
        args[4] = QString::number(lon);
        args[5] = QString::number(lat);
    }
        break;
    case 1:
        /*result = smartastro::astrocore::conversion_coordinates::cylrec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 2:
    {
        double r = args[0].toDouble();
        double lonc = args[1].toDouble();
        double z = args[2].toDouble();
        double radius = args[3].toDouble();
        double colat = args[4].toDouble();
        double lon = args[5].toDouble();

        smartastro::astrocore::conversion_coordinates::cylsph(r, lonc, z, radius, colat, lon);

        args[0] = QString::number(r);
        args[1] = QString::number(lonc);
        args[2] = QString::number(z);
        args[3] = QString::number(radius);
        args[4] = QString::number(colat);
        args[5] = QString::number(lon);
    }
        break;
    case 3:
        /*result = smartastro::astrocore::conversion_coordinates::dcyldr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 4:
        /*result = smartastro::astrocore::conversion_coordinates::dgeodr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));*/
        break;
    case 5:
        /*result = smartastro::astrocore::conversion_coordinates::dlatdr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]),stringToDoublerealP(args[3]));*/
        break;
    case 7:
        /*result = smartastro::astrocore::conversion_coordinates::dpgrdr_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), stringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));*/
        break;
    case 8:
        /*result = smartastro::astrocore::conversion_coordinates::drdcyl_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 9:
        /*result = smartastro::astrocore::conversion_coordinates::drdgeo_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]),qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));*/
        break;
    case 10:
        /*result = smartastro::astrocore::conversion_coordinates::drdlat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 11:
        /*result = smartastro::astrocore::conversion_coordinates::drdpgr_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), stringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));*/
        break;
    case 12:
        /*result = smartastro::astrocore::conversion_coordinates::drdsph_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 13:
        /*result = smartastro::astrocore::conversion_coordinates::dsphdr_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), stringToDoublerealP(args[3]));*/
        break;
    case 14:
        /*result = smartastro::astrocore::conversion_coordinates::georec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));*/
        break;
    case 15:
    {
        double radius = args[0].toDouble();
        double lon = args[1].toDouble();
        double lat = args[2].toDouble();
        double r = args[3].toDouble();
        double lonc = args[4].toDouble();
        double z = args[5].toDouble();
        
        smartastro::astrocore::conversion_coordinates::latcyl(radius, lon, lat, r, lonc, z);

        args[0] = QString::number(radius);
        args[1] = QString::number(lon);
        args[2] = QString::number(lat);
        args[3] = QString::number(r);
        args[4] = QString::number(lonc);
        args[5] = QString::number(z);
    }
        break;
    case 16:
        /*result = smartastro::astrocore::conversion_coordinates::latrec_(stringToDoublerealP(args[0]),qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 17:
    {
        double radius = args[0].toDouble();
        double lon = args[1].toDouble();
        double lat = args[2].toDouble();
        double rho = args[3].toDouble();
        double colat = args[4].toDouble();
        double lons = args[5].toDouble();

        smartastro::astrocore::conversion_coordinates::latsph(radius, lon, lat, rho, colat, lons);

        args[0] = QString::number(radius);
        args[1] = QString::number(lon);
        args[2] = QString::number(lat);
        args[3] = QString::number(rho);
        args[4] = QString::number(colat);
        args[5] = QString::number(lons);
    }
        break;
    case 18:
    {
        const std::string body = args[0].toStdString();
        double et = args[1].toDouble();
        const std::string abcorr = args[2].toStdString();
        double result;
        
        result = smartastro::astrocore::spice_general_functions::lspcn(body, et, abcorr);

        args[0] = QString::fromStdString(body);
        args[1] = QString::number(et);
        args[2] = QString::fromStdString(abcorr);
        args.push_back(QString::number(result));
    }
        break;
    case 19:
        /*result = smartastro::astrocore::conversion_coordinates::pgrrec_(qstringToCharP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));*/
        break;
    case 20:
        /*result = smartastro::astrocore::conversion_coordinates::radrec_(stringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 21:
        /*result = smartastro::astrocore::conversion_coordinates::reccyl_(stringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 22:
        /*result = smartastro::astrocore::conversion_coordinates::recgeo_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]));*/
        break;
    case 23:
        /*result = smartastro::astrocore::conversion_coordinates::reclat_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 24:
        /*result = smartastro::astrocore::conversion_coordinates::recpgr_(qstringToCharP(args[0]),  qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]), qstringToDoublerealP(args[4]), qstringToDoublerealP(args[5]), qstringToDoublerealP(args[6]), qstringToFtnlen(args[7]));*/
        break;
    case 25:
        /*result = smartastro::astrocore::conversion_coordinates::recrad_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 26:
        /*result = smartastro::astrocore::conversion_coordinates::recsph_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 27:
    {
        double radius = args[0].toDouble();
        double colat = args[1].toDouble();
        double slon = args[2].toDouble();
        double r = args[3].toDouble();
        double lon = args[4].toDouble();
        double z = args[5].toDouble();

        smartastro::astrocore::conversion_coordinates::sphcyl(radius, colat, slon, r, lon, z);

        args[0] = QString::number(radius);
        args[1] = QString::number(colat);
        args[2] = QString::number(slon);
        args[3] = QString::number(r);
        args[4] = QString::number(lon);
        args[5] = QString::number(z);
    }
        break;
    case 28:
    {
        double r = args[0].toDouble();
        double colat = args[1].toDouble();
        double lons = args[2].toDouble();
        double radius = args[3].toDouble();
        double lon = args[4].toDouble();
        double lat = args[5].toDouble();

        smartastro::astrocore::conversion_coordinates::sphlat(r, colat, lons, radius, lon, lat);

        args[0] = QString::number(r);
        args[1] = QString::number(colat);
        args[2] = QString::number(lons);
        args[3] = QString::number(radius);
        args[4] = QString::number(lon);
        args[5] = QString::number(lat);
    }
        break;
    case 29:
        /*result = smartastro::astrocore::conversion_coordinates::sphrec_(qstringToDoublerealP(args[0]), qstringToDoublerealP(args[1]), qstringToDoublerealP(args[2]), qstringToDoublerealP(args[3]));*/
        break;
    case 30:
        /*result = smartastro::astrocore::spice_general_functions::subpnt_(qstringToCharP(args[0]), qstringToCharP(args[1]), qstringToDoublerealP(args[2]), qstringToCharP(args[3]), qstringToCharP(args[4]), qstringToCharP(args[5]), qstringToDoublerealP(args[6]), qstringToDoublerealP(args[7]), qstringToDoublerealP(args[8]), stringToFtnlen(args[9]), stringToFtnlen(args[10]), stringToFtnlen(args[11]), stringToFtnlen(args[12]), stringToFtnlen(args[13]));*/
        break;
    case 31:
        /*result = smartastro::astrocore::spice_general_functions::subslr_(qstringToCharP(args[0]), qstringToCharP(args[1]), qstringToDoublerealP(args[2]), qstringToCharP(args[3]), qstringToCharP(args[4]), qstringToCharP(args[5]), qstringToDoublerealP(args[6]), qstringToDoublerealP(args[7]), qstringToDoublerealP(args[8]), stringToFtnlen(args[9]), stringToFtnlen(args[10]), stringToFtnlen(args[11]), stringToFtnlen(args[12]), stringToFtnlen(args[13]));*/
        break;
   }
}


