/******************************************************************************
 *                               CONVERSION_TIME_H                            *
 *            List of conversion functions of the SMART-SIM toolbox           *
 ******************************************************************************/

#ifndef SMART_CONVERSION_TIME_H
#define SMART_CONVERSION_TIME_H

#include <vector>
#include <cmath>
#include "LinearAlgebra/Eigen/Dense"
#include "../exception.h"
#include "../constants.h"

namespace smartastro{
    namespace astrocore{

        /**
         * @brief The conversion_space class implements a set of static methods for calendar conversion.
         */
        class conversion_time{
        public:

            /******************************************************************************
             *                       CALENDAR CONVERSION FUNCTIONS                        *
             ******************************************************************************/

            /**
             * @brief hms2fracday - Hours, minutes, and seconds to a fraction of day.
             *
             * Function to convert hours, minutes, and seconds into a fraction of day.
             * @param[in] hms vector of 3 doubles.
             *            - hms[0] = number of hours (24 hours format)
             *            - hms[1] = number of minutes
             *            - hms[2] = number of seconds
             * @param[out] fracday double representing the time expressed in day
             * @return Error code
             *
             * @author Luca Masi 2008, Massimiliano Vasile 2008, Federico Zuiani 2008
             */
            static bool hms2fracday(const std::vector<double> &hms, double &fracday);

            /**
             * @brief fracday2hms - Fraction of day to hours, minutes, and seconds.
             *
             * Function to convert a fraction of day into hours, minutes, and seconds.
             * @param[in] fracday double representing the time expressed in day
             * @param[out] hms vector of 3 doubles.
             *            - hms[0] = number of hours (24 hours format)
             *            - hms[1] = number of minutes
             *            - hms[2] = number of seconds
             * @return Error code
             *
             * @author Massimiliano Vasile 2008, Federico Zuiani 2008
             */
            static bool fracday2hms(const double &fracday, std::vector<double> &hms);

            /**
             * @brief jd2mjd - Modified Julian date from Julian date.
             *
             * Function to convert Julian date to Modified Julian date
             * @param[in] jd double containing the date in Julian Day.
             *               The JD (Julian day) count is from 0 at 12:00
             *               noon, 1 January -4712 (4713 BC), Julian
             *               proleptic calendar. The corresponding date in
             *               Gregorian calendar is 12:00 noon, 24 November
             *               -4713.
             *
             * @param[out] mjd double containing the date in Modified Julian Day.
             *               The MJD count is from 00:00 midnight at the
             *               beginning of Wednesday November 17, 1858.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008.
             *
             */
            static bool jd2mjd(const double &jd, double &mjd);

            /**
             * @brief mjd2jd - Julian date from Modified Julian date.
             *
             * Function to convert Modified Julian date to Julian date.
             * @param[in] mjd double containing the date in Modified Julian Day.
             *            The MJD count is from 00:00 midnight at the
             *            beginning of Wednesday November 17, 1858.
             * @param[out] jd double containing the date in Julian Day.
             *         The JD (Julian day) count is from 0 at 12:00
             *         noon, 1 January -4712 (4713 BC), Julian
             *         proleptic calendar. The corresponding date in
             *         Gregorian calendar is 12:00 noon, 24 November
             *         -4713.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Massimiliano Vasile 2008
             */
            static bool mjd2jd(const double &mjd, double &jd);

            /**
             * @brief jd2mjd2000 - Modified Julian date 2000 from Julian date.
             *
             * Function to convert from Julian date to Modified Julian date 2000.
             * @param[in] jd double containing the date in Julian Day.
             *               The JD (Julian day) count is from 0 at 12:00
             *               noon, 1 January -4712 (4713 BC), Julian
             *               proleptic calendar. The corresponding date in
             *               Gregorian calendar is 12:00 noon, 24 November
             *               -4713.
             * @param[out] mjd2000 double containing the date in Modified Julian
             *         Day 2000.
             *         MJD2000 is defined as the number of days since
             *         01-01-2000, 12:00 noon.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008
             */
            static bool jd2mjd2000(const double &jd, double &mjd2000);

            /**
             * @brief mjd20002jd - Julian date from Modified Julian date 2000.
             *
             * Function to convert from Modified Julian date 2000 to Julian date.
             * @param[in] mjd2000  double containing the date in Modified Julian
             *                     Day 2000.
             *                     MJD2000 is defined as the number of days since
             *                     01-01-2000, 12:00 noon.
             * @param[out] jd double containing the date in Julian Day.
             *         The JD (Julian day) count is from 0 at 12:00
             *         noon, 1 January -4712 (4713 BC), Julian
             *         proleptic calendar. The corresponding date in
             *         Gregorian calendar is 12:00 noon, 24 November
             *         -4713.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Massimiliano Vasile 2008.
             */
            static bool mjd20002jd(const double &mjd2000, double &jd);

            /**
             * @brief mjd20002mjd - Julian date from Modified Julian date 2000.
             *
             * Function to convert from Modified Julian date 2000 to Julian date.
             * @param[in] mjd2000 double containing the date in Modified
             *                    Julian Day 2000.
             *                    MJD2000 is defined as the number of days
             *                    since 01-01-2000, 12:00 noon.
             * @param[out] mjd double containing the date in Modified
             *         Julian Day. The MJD count is from 00:00 midnight at the
             *         beginning of Wednesday November 17, 1858.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Massimiliano Vasile 2008.
             */
            static bool mjd20002mjd(const double &mjd2000, double &mjd);

            /**
             * @brief mjd2mjd2000 - Modfied Julian date 2000 from Modified Julian date.
             *
             * Function to convert from Modified Julian date to Modfied Julian date 2000.
             * @param[in] mjd double containing the date in Modified
             *                Julian Day.
             *                The MJD count is from 00:00 midnight at the
             *                beginning of Wednesday November 17, 1858.
             * @param[out] mjd2000 double containing the date in Modified
             *         Julian Day 2000.
             *         MJD2000 is defined as the number of days
             *         since 01-01-2000, 12:00 noon.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008.
             */
            static bool mjd2mjd2000(const double &mjd, double &mjd2000);

            /**
             * @brief date2jd - Julian date from Gregorian date.
             *
             * Function to convert from Gregorian date to Julian date.
             * This function computes the Julian date of the given date
             * (Gregorian calendar) plus a fractional part depending on the time of day
             *
             *  WARNINGS:
             *      - The function is valid for the whole range of dates since 12:00
             *         noon 24 November -4713, Gregorian calendar. (This bound is set
             *         in order to have symmetry with the inverse function JD2DATE)
             *      - The inputs must be feasible (i.e. the date must exist!). If an
             *         unfeasible date is inputed, wrong results are given because no
             *         check is done on that.
             *      - For dates before 1582, the resulting date components are valid
             *         only in the Gregorian proleptic calendar. This is based on the
             *         Gregorian calendar but extended to cover dates before its
             *         introduction.
             * @param[in] date vector of 6 doubles containing the
             *                 date in the Gregorian calendar.
             *                 date = [year, month, day, hour, minute, second]
             * @param[out] jd double containing the date in Julian Day.
             *                The JD (Julian day) count is from 0 at 12:00
             *                noon, 1 January -4712 (4713 BC), Julian
             *                proleptic calendar. The corresponding date in
             *                Gregorian calendar is 12:00 noon, 24 November
             *                -4713.
             * @return Error code
             *
             * @see Formula from http://scienceworld.wolfram.com/astronomy/JulianDate.html
             *    Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
             *    were found
             *
             * @author Federico Zuiani 2008, Luca Masi 2008
             */
            static bool date2jd(const std::vector<double> &date, double &jd);

            /**
             * @brief jd2date - Gregorian date from Julian date.
             *
             * This function computes the given date (Gregorian calendar) from the
             *  Julian date.
             *
             *  WARNINGS:
             *      - jd must be a non-negative real. This means that the function is
             *         valid for the whole range of dates since 12:00 noon 24 November
             *         -4713, Gregorian calendar.
             *      - For dates before 1582, the resulting date components are valid
             *         only in the Gregorian proleptic calendar. This is based on the
             *         Gregorian calendar but extended to cover dates before its
             *         introduction.
             * @param[in] jd double containing the date in Julian Day.
             *                The JD (Julian day) count is from 0 at 12:00
             *                noon, 1 January -4712 (4713 BC), Julian
             *                proleptic calendar. The corresponding date in
             *                Gregorian calendar is 12:00 noon, 24 November
             *                -4713.
             * @param[out] date vector of 6 doubles containing the
             *                  date in the Gregorian calendar.
             *                  date = [year, month, day, hour, minute, second]
             *                  date must be pre-allocated.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008.
             */
            static bool jd2date(const double &jd, std::vector<double> &date);


            /**
             * @brief date2mjd - Modified Julian date from Gregorian date.
             *
             * This function computes the modified Julian date of the given date
             *  (Gregorian calendar) plus a fractional part depending on the time of day
             *
             *  WARNINGS:
             *      - The function is valid for the whole range of dates since 12:00
             *         noon 24 November -4713, Gregorian calendar.
             *      - The inputs must be feasible (i.e. the date must exist!). If an
             *         unfeasible date is inputed, wrong results are given because no
             *         check is done on that.
             *      - For dates before 1582, the resulting date components are valid
             *         only in the Gregorian proleptic calendar. This is based on the
             *         Gregorian calendar but extended to cover dates before its
             *         introduction.
             * @param[in] date vector of 6 doubles containing the
             *                 date in the Gregorian calendar.
             *                 date = [year, month, day, hour, minute, second]
             * @param[out] mjd double containing the date in Modified Julian
             *                 Day.
             *                 The MJD count is from 00:00 midnight at the
             *                 beginning of Wednesday November 17, 1858.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008.
             */
            static bool date2mjd(const std::vector<double> &date, double &mjd);

            /**
             * @brief mjd2date - Gregorian date from modified Julian date.
             *
             *  This function computes the given date (Gregorian calendar) from the
             *  modified Julian date.
             * @param[in] mjd double containing the date in Modified Julian
             *                Day.
             *                The MJD count is from 00:00 midnight at the
             *                beginning of Wednesday November 17, 1858.
             * @param[out] date vector of 6 doubles containing the
             *                 date in the Gregorian calendar.
             *                 date = [year, month, day, hour, minute, second]
             * @return Error code
             *
             * @author Federico Zuiani 2008, Massimiliano Vasile 2008.
             */
            static bool mjd2date(const double &mjd, std::vector<double> &date);

            /**
             * @brief date2mjd2000 - Modified Julian date 2000 from Gregorian date.
             *
             * This function computes the modified Julian date 2000 of the given
             *  date (Gregorian calendar) plus a fractional part depending on the time
             *  of day.
             *  WARNINGS:
             *      - The function is valid for the whole range of dates since 12:00
             *         noon 24 November -4713, Gregorian calendar.
             *      - The inputs must be feasible (i.e. the date must exist!). If an
             *         unfeasible date is inputed, wrong results are given because no
             *         check is done on that.
             *      - For dates before 1582, the resulting date components are valid
             *         only in the Gregorian proleptic calendar. This is based on the
             *         Gregorian calendar but extended to cover dates before its
             *         introduction.
             *
             * @param[in] date vector of 6 doubles containing the
             *                 date in the Gregorian calendar.
             *                 date = [year, month, day, hour, minute, second]
             * @param[out] mjd2000 double containing the date in Modified Julian
             *                Day 2000.
             *                MJD2000 is defined as the number of days since
             *                01-01-2000, 12:00 noon.
             * @return Error code
             *
             * @author Federico Zuiani 2008, Luca Masi 2008.
             */
            static bool date2mjd2000(const std::vector<double> &date, double &mjd2000);

            /**
             * @brief mjd20002date - Gregorian date from modified Julian date 2000.
             *
             * This function computes the given date (Gregorian calendar) from the
             *  modified Julian date 2000.
             *
             *  WARNINGS:
             *      - mjd2000 must be a valid (ie the corresponding jd>=0). This means
             *         that the function is valid for the whole range of dates since
             *         12:00 noon 24 November -4713, Gregorian calendar.
             *      - For dates before 1582, the resulting date components are valid
             *         only in the Gregorian proleptic calendar. This is based on the
             *         Gregorian calendar but extended to cover dates before its
             *         introduction.
             * @param[in] mjd2000 double containing the date in Modified Julian
             *                    Day 2000.
             *                    MJD2000 is defined as the number of days since
             *                    01-01-2000, 12:00 noon.
             * @param[out] date vector of 6 doubles containing the
             *                 date in the Gregorian calendar.
             *                 date = [year, month, day, hour, minute, second]
             * @return Error code
             *
             * @author Federico Zuiani 2008, Massimiliano Vasile 2008.
             */
            static bool mjd20002date(const double &mjd2000, std::vector<double> &date);

            /**
             * @brief mean2eccentric_anomaly - solves Kepler equation
             *
             * This function computes the eccentric anomaly corresponding to a given mean anomaly by solving Kepler equation
             *
             * @param[in] M mean anomaly in radians
             * @param[in] eccentricity 
             * @param[out] E eccentric anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool mean2eccentric_anomaly(const double &M, const double &eccentricity, double &E);

            /**
             * @brief eccentric2mean_anomaly - converts eccentric to mean anomaly
             *
             * This function computes the mean anomaly corresponding to a given eccentric anomaly
             *
             * @param[in] E eccentric anomaly in radians
             * @param[in] eccentricity 
             * @param[out] M mean anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool eccentric2mean_anomaly(const double &E, const double &eccentricity, double &M);


            /**
             * @brief true2eccentric_anomaly - converts true to eccentric anomaly
             *
             * This function computes the eccentric anomaly corresponding to a given true anomaly
             *
             * @param[in] theta true anomaly in radians
             * @param[in] eccentricity 
             * @param[out] E eccentric anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool true2eccentric_anomaly(const double &theta, const double &eccentricity, double &E);

            /**
             * @brief eccentric2true_anomaly - converts eccentric to true anomaly
             *
             * This function computes the true anomaly corresponding to a given eccentric anomaly
             *
             * @param[in] E eccentric anomaly in radians
             * @param[in] eccentricity 
             * @param[out] theta true anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool eccentric2true_anomaly(const double &E, const double &eccentricity, double &theta);

            /**
             * @brief true2mean_anomaly - converts true to mean anomaly
             *
             * This function computes the true anomaly corresponding to a given mean anomaly
             *
             * @param[in] theta true anomaly in radians
             * @param[in] eccentricity 
             * @param[out] M mean anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool true2mean_anomaly(const double &theta, const double &eccentricity, double &M);

            /**
             * @brief mean2true_anomaly - converts mean to true anomaly
             *
             * This function computes the eccentric anomaly corresponding to a given true anomaly
             *
             * @param[in] M mean anomaly in radians
             * @param[in] eccentricity 
             * @param[out] theta true anomaly in radians
             * @return Error code
             *
             * @author Romain Serra 2017
             */
            static bool mean2true_anomaly(const double &M, const double &eccentricity, double &theta);


        };

    }
}


#endif /* SMART_CONVERSION_TIME_H */

