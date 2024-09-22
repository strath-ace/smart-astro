"""
    hms2fracday(hrs, mn, sec)
    converts hours, minutes, and seconds to the fraction of day

    INPUT
    hrs        Number of hours as integer greater or equal to 0 and lower
                or equal to 23.
    mn         Number of minutes as integer greater or equal to 0 and
                lower or equal to 59.
    sec        Number of seconds as a real greater or equal to 0 and
                strictly lower than 60.

    OUTPUT
    fracDay    A single real greater or equal to 0 and strictly lower than
                1.

    See also FRACDAY2HMS.

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    - Nicolas Croisard - 16/02/2008
    - Revised by Camilla Colombo - 29/02/2008

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language

"""
function hms2fracday(hrs, mn, sec)

    # Check the inputs
   
    # Check the hours
    if typeof(hrs)!=Int64
        Base.error("HMS2FRACDAY:incorrectInput : the hours must be given as a single integer");
    elseif floor(hrs)!=hrs || hrs<0 || hrs>23
        Base.error("HMS2FRACDAY:incorrectInput: the hours must be given as a float greater or equal to 0 and lower or equal to 23");
    end
    
    # Check the minutes
    if typeof(mn)!=Int64
        Base.error("HMS2FRACDAY:incorrectInput : The minutes must be given as a single integer");
    elseif floor(mn)!=mn || mn<0 || mn>59
        Base.error("HMS2FRACDAY:incorrectInput: The minutes must be given as an integer greater or equal to 0 and lower or equal to 59");
    end
    
    # Check the seconds
    if typeof(sec)!=Float64
        Base.error("HMS2FRACDAY:incorrectInput: The seconds must be given as a single float number");
    elseif sec<0 || sec>60
        Base.error("HMS2FRACDAY:incorrectInput: The seconds must be given as an real greater or equal to 0 and strictly lower then 60");
    end
    
    # Compute the fraction of day
    fracDay = (Float64(hrs) + (Float64(mn) + sec/60.0)/60.0) / 24.0;
    
    return fracDay
end

"""
    fracday2hms(fracDay)
    Converts a fraction of day into hours, minutes, and
	seconds.

    INPUT:
    fracDay[1] A single real greater or equal to 0 and strictly lower than
                1.

    OUTPUT:
    hrs[1]     Number of hours as integer greater or equal to 0 and lower
                or equal to 23.
    mn[1]      Number of minutes as integer greater or equal to 0 and
                lower or equal to 59.
    sec[1]     Number of seconds as a real greater or equal to 0 and
                strictly lower than 60.

    See also hms2fracday.

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    Nicolas Croisard, 16/02/2008, MATLAB, fracday2hms.m

    CHANGELOG:
    21/02/2008, REVISION, Matteo Ceriotti
    22/04/2010, Camilla Colombo: Header and function name in accordance
        with guidlines.

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language

"""
function fracday2hms(fracDay)   
    # Check the input
    if typeof(fracDay)!=Float64
        Base.error("FRACDAY2HMS:incorrectInput: the input should a single element");
    end
    if fracDay<0 || fracDay>=1
        error("FRACDAY2HMS:incorrectInput: the input should be real greater or equal to 0 and strictly lower than 1");
    end
    
    temp = fracDay*24;
    hrs   = temp
    hrs = Int(sign(hrs)*floor(abs(hrs)));
    mn = (temp-hrs)*60; 
    mn = Int(sign(mn)*floor(abs(mn))); 
    sec = (temp-hrs-mn/60)*3600;
    
    return  hrs, mn, sec;
end

"""
    date2jd(date)
    Returns the Julian day number of the given date (Gregorian calendar)
    plus a fractional part depending on the time of day.
    Note: The function is valid for the whole range of dates since 12:00
        noon 24 November -4713, Gregorian calendar. (This bound is set in
        order to have symmetry with the inverse function jd2date.m)
    Note: The inputs must be feasible (i.e. the date must exist!). If an
        unfeasible date is inputed, wrong results are given because no
        check is done on that.

    INPUT:
    date[6]     Date in the Gregorian calendar, as a 6-elements vector
                [year, month, day, hour, minute, second]. For dates before
                1582, the resulting date components are valid only in the
                Gregorian proleptic calendar. This is based on the
                Gregorian calendar but extended to cover dates before its
                introduction. Date must be after 12:00 noon, 24 November
                -4713.

    OUTPUT
    jd          Date in Julian Day. The JD (Julian day) count is from 0 at
                12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
                calendar. The corresponding date in Gregorian calendar is
                12:00 noon, 24 November -4713.

    REFERENCES:
    Formula from http://scienceworld.wolfram.com/astronomy/JulianDate.html
    (last visited 15/02/2008)
    Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
    were found

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    Nicolas Croisard, 16/02/2008, MATLAB, date2jd.m
    CHANGELOG:
    03/03/2008, REVISION: Camilla Colombo
    22/04/2010, Camilla Colombo: Header and function name in accordance
        with guidelines.

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
"""
function date2jd(date)
    # Check the input
    if size(date)[1]!=6
        Base.error("wrong input");
    end
    
    #  Manage the input
    Y   = Float64(date[1]);
    M   = Float64(date[2]);
    D   = Float64(date[3]);
    hrs = Int(date[4]);
    mn  = Int(date[5]);
    sec = Float64(date[6]);
    
    # Check the inputs
    if Y<-4713 || (Y==-4713 && (M<11 || (M==11 && D<24 || (D==24 && hrs<12))))
        Base.error("DATE2JD:incorrectInput: The function is valid for dates after since 12:00 noon 24 November -4713, Gregorian calendar");
    end
    
    # Formula converting Gregorian date into JD
    jd = 367.0*Y - floor(7.0*(Y+floor((M+9.0)/12.0))/4.0) - floor(3.0*floor((Y+(M-9.0)/7.0)/100.0+1.0)/4.0) + floor(275.0*M/9.0) + D + 1721028.5 + hms2fracday(hrs,mn,sec);

    return jd;

end

"""
    jd2date(jd)
    returns the Gregorian calendar date (year, month, day,
    hour, minute, and second) corresponding to the Julian day number JD.

    Note: jd must be a non-negative real. This means that the function is
        valid for the whole range of dates since 12:00 noon 24 November
        -4713, Gregorian calendar. (This bound is due to the function
        FRACDAY2HMS that does not work for negative inputs)

    INPUT
    jd          Date in Julian Day. The JD (Julian day) count is from 0 at
                12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
                calendar. The corresponding date in Gregorian calendar is
                12:00 noon, 24 November -4713. It must be a non-negative
                real.

    OUTPUT
    date        Date in the Gregorian calendar, as a 6-element vector
                [year, month, day, hour, minute, second]. For dates before
                1582, the resulting date components are valid only in the
                Gregorian proleptic calendar. This is based on the
                Gregorian calendar but extended to cover dates before its
                introduction.

    REFERENCES
    Formula from http://en.wikipedia.org/wiki/Julian_day
    (last visited 16/02/2008)
    Compared to http://pdc.ro.nu/mjd.cgi for a few dates, the same results
    were found


    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    - Nicolas Croisard - 16/02/2008
    - Revised by Matteo Ceriotti - 21/02/2008 - Validated with:
    - Fliegel, Van Flandern "A machine Algorithm for Processing Calendar
    Dates", Communications of the ACM, 1968. Also on Wertz, "Space
        Mission Analysis and Design".
    - A revised version of the algorithm on Vallado, "Fundamentals of
        Astrodynamics and Applications", third edition, for dates from year
        1900 to year 2100.

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
"""
function jd2date(jd)

    if jd < 0
        BASE.error("JD2DATE:jdLessThanZero: The input jd value cannot be negative");
    end
    
    # Adding 0.5 to JD and taking FLOOR ensures that the date is correct.
    j = floor(jd+0.5) + 32044.0;
    g = floor(j/146097.0);
    dg = mod(j,146097.0);
    c = floor((floor(dg/36524.0)+1.0) * 3.0/4.0);
    dc = dg - c*36524.0;
    b = floor(dc/1461.0);
    db = mod(dc,1461.0);
    a = floor((floor(db/365.0)+1.0) * 3.0/4.0);
    da = db - a*365.0;
    y = g*400.0 + c*100.0 + b*4.0 + a;
    m = floor((da*5.0 + 308.0)/153.0) - 2.0;
    d = da - floor((m+4.0)*153.0/5.0) + 122.0;
    
    # Year, Month and Day
    Y = y-4800.0 + floor((m+2.0)/12.0);
    M = mod((m+2.0),12.0) + 1.0;
    D = floor(d+1.0);
    
    # Hour, Minute and Second
    hrs, mn, sec = fracday2hms(mod(jd+0.5,floor(jd+0.5)));
    
    # Prepare output
    date = [Y, M, D, Int(hrs), Int(mn), sec];
    
    
    return date
end

"""
    date2mjd2000(date)
    Returns the modified Julian day 2000 number of the given date
    (Gregorian calendar) plus a fractional part depending on the time of
    day.
    Note: The function is valid for the whole range of dates since 12:00
    noon 24 November -4713, Gregorian calendar. (This bound is set in
    order to have symmetry with the inverse function jd2date)
    Note: The inputs must be feasible (i.e. the date must exist!). If an
    unfeasible date is inputed, wrong results are given because no
    check is done on that.

    INPUT:
    date[6] Date in the Gregorian calendar, as a 6-element vector
            [year, month, day, hour, minute, second]. For dates before
            1582, the resulting date components are valid only in the
            Gregorian proleptic calendar. This is based on the
            Gregorian calendar but extended to cover dates before its
            introduction. date must be after 12:00 noon, 24 November
            -4713.

    OUTPUT:
    mjd2000[1]  Date in MJD 2000. MJD2000 is defined as the number of days
                since 01-01-2000, 12:00 noon.


    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    AUTHOR:
    Nicolas Croisard, 16/02/2008, MATLAB, date2mjd2000.m

    CHANGELOG:
    03/03/2008, REVISION, Camilla Colombo
    22/04/2010, Camilla Colombo: Header and function name in accordance
    with guidlines.

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
"""
function date2mjd2000(date)

    jd = date2jd(date);
    mjd2000 = jd - 2400000.5 - 51544.5;
    return mjd2000;

end

"""
    mjd20002date(mjd2000)
    returns the Gregorian calendar date
    (year, month, day, hour, minute, and second) corresponding to the
    modified Julian day 2000 number.

    Note: Since this function calls jd2date, MJD2000 cannot be less than
        -2451545, that is 24 November -4713, 12:00 noon, Gregorian
        proleptic calendar.

    INPUT
    mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
                since 01-01-2000, 12:00 noon. It must be a real greater or
                equal than -2451545.

    OUTPUT
    date        Date in the Gregorian calendar, as a 6-element vector
                [year, month, day, hour, minute, second]. For dates before
                1582, the resulting date components are valid only in the
                Gregorian proleptic calendar. This is based on the
                Gregorian calendar but extended to cover dates before its
                introduction.

    MATLAB VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%
    - Nicolas Croisard - 16/02/2008
    - Revised by Matteo Ceriotti - 21/02/2008

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    08/2023 : Thomas Mcilwraith : translation in Julia language
"""
function mjd20002date(mjd2000)

    jd = mjd2000 + 2400000.5 + 51544.5 ;
    date = jd2date(jd);
    return date;

end

"""
    jd2greenwichsiderealtime(jd)
    returns Greenwich sidereal time 

    INPUT
    mjd2000     modified julian day 2000
    type        = 1 : output in radians
                = 2 : output in seconds
                = 3 : output in degree

    OUTPUT
    GST         Greenwich sidereal time

    JULIA VERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    01/2024   : Irene Cavallari

    REF
    Vallado,1997

"""
function mjd20002greenwichsiderealtime(mjd2000,type)

    # rad
    mjd2000_0 = floor(mjd2000);
    T = (mjd2000)/36525.0;
    nEarth = astroconstants(3)[6];
    ut = mjd2000-mjd2000_0;
    
    # rad
    if type == 1 || type == 3
        GST = 1.753368560  + 628.3319706889 * T + 6.7707*1e-6 * T^2 - 4.51e-10 * T^3;
        GST = mod(GST + nEarth*ut,2.0*pi);
        if type == 3
            GST = GST*pi/180.0;
        end
    else
     # seconds
        GST = 24110.54841 + 8640184.812866 * T + 0.093104 * T^2 - 0.0000062 * T^3;
        GST = GST + (nEarth*ut)*4.0*60.0*180.0/pi;
    end

    return GST;
end

"""
        MJD20002JD Julian day number from Modified Julian day 2000 number.

        JD = MJD20002JD(MJD2000) returns the Julian day number corresponding to
        the given modified Julian day 2000 number.

        INPUT
            mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
                    since 01-01-2000, 12:00 noon. 

        OUTPUT
        jd          Date in Julian Day. The JD (Julian day) count is from 0 at
                    12:00 noon, 1 January -4712 (4713 BC), Julian proleptic
                    calendar. The corresponding date in Gregorian calendar is
                    12:00 noon, 24 November -4713.

        See also JD2MJD2000.

        FUNCTIONS CALLED
        none

        MATLAB VERSION
        - Nicolas Croisard - 16/02/2008
        - Revised by Matteo Ceriotti - 20/02/2008
        
        ---- JULIA TRANSLATION
        - Irene Cavallari - 01/2024
"""
function mjd20002jd(mjd2000)
    mjd = mjd2000 + 51544.5;
    jd  = mjd + 2400000.5;
    return jd;
end


"""
        jd2mjd2000
        returns the modified Julian day 2000 number corresponding to the given Julian day number .

        INPUT
                jd = Date in Julian Day. The JD (Julian day) count is from 0 at 12:00 noon,
                1 January -4712 (4713 BC), Julian proleptic calendar. 
                The corresponding date in Gregorian calendar is  12:00 noon, 24 November -4713. 

        OUTPUT
                mjd2000     Date in MJD 2000. MJD2000 is defined as the number of days
                since 01-01-2000, 12:00 noon.

        See also MJD20002JD.

        FUNCTIONS CALLED
        none

        MATLAB VERSION
        - Nicolas Croisard - 16/02/2008
        - Revised by Matteo Ceriotti - 20/02/2008
        
        ---- JULIA TRANSLATION
        - Irene Cavallari - 01/2024
"""
function jd2mjd2000(jd)
    mjd = jd - 2400000.5;
    mjd2000 = mjd - 51544.5;
    return mjd2000;
end