"""
   exponential_atm_model(h,A)
   Function to read the exponential atmospheric model 
   given in input and determines the density, 
   
   INPUT:
         h = current altitude
         A = matrix containing atmospheric model
               in each row: [min altitude [km], density [kg/m^3], scale height [km]]

   OUTPUT:
         rho = density

"""
function exponential_atm_model(h,A)
    
   ier = 1;
   run = true;
   n = size(A,1);
   i = 1;

   #loop to find relevant position in atmosphere matrix from altitude
   while run == true
   
        if h < A[n,1] 
           n=n-1;
        else
           i = n;
           run = false;   
        end
   
        if n == 0
           ier = -1;
           run = false;
        end

   end
   
   
   rho0 = A[i,2];
   href = A[i,1];
   Href = A[i,3];

   rho = rho0*exp(-(h-href)/Href);
   return rho;
end

"""
   getExponentialAtmosphericDensityModel()
   exponential density model for the Earth
   ref Vallado, Table 7-4 (8-4 in 4th ed)

   INPUT:
   OUTPUT:  A = matrix containing atmospheric model
            in each row: [min altitude [km], density [kg/m^3], scale height [km]]
"""
function getExponentialAtmosphericDensityModelForEarth()

   #         [h_ellp, rho0, H], Vallado, Table 7-4 (8-4 in 4th ed)
   A  = [ 0      25    0  1.225       7.249;
         25     30   25  3.899e-2    6.349;
         30     40   30  1.774e-2    6.682;
         40     50   40  3.972e-3    7.554;
         50     60   50  1.057e-3    8.382;
         60     70   60  3.206e-4    7.714;
         70     80   70  8.770e-5    6.549;
         80     90   80  1.905e-5    5.799;
         90    100   90  3.396e-6    5.382;
         100    110  100  5.297e-7    5.877;
         110    120  110  9.661e-8    7.263;
         120    130  120  2.438e-8    9.473;
         130    140  130  8.484e-9   12.636;
         140    150  140  3.845e-9   16.149;
         150    180  150  2.070e-9   22.523;
         180    200  180  5.464e-10  29.740;
         200    250  200  2.789e-10  37.105;
         250    300  250  7.248e-11  45.546;
         300    350  300  2.418e-11  53.628;
         350    400  350  9.518e-12  53.298;
         400    450  400  3.725e-12  58.515;
         450    500  450  1.585e-12  60.828;
         500    600  500  6.967e-13  63.822;
         600    700  600  1.454e-13  71.835;
         700    800  700  3.614e-14  88.667;
         800    900  800  1.170e-14  124.64;
         900   1000  900  5.245e-15  181.05;
         1000   1e5  1000  3.019e-15  268.00];
      A = A[:,3:end]
   return A;
end

### it contains some approximations 
function getExponentialAtmosphericDensityModelForEarthOLD()

   #         [h_ellp, rho0, H], Vallado, Table 7-4 (8-4 in 4th ed)
   A  = [0	1.225	   7.249;
   25	   3.90e-2	   6.349;
   30	   1.77e-2	   6.682;
   40	   3.97e-3	   7.554;
   50	   1.06e-3	   8.382;
   60	   3.21e-4	   7.714;
   70	   8.77e-5	   6.549;
   80	   1.91e-5	   5.799;
   90	   3.40e-6	   5.382;
   100	5.2970e-7	5.877;
   110	9.6610e-8	7.263;
   120	2.4380e-8	9.473;
   130	8.4840e-9	12.636;
   140	3.8450e-9	16.149;
   150	2.0700e-9	22.523;
   180	5.4640e-10	29.74;
   200	2.7890e-10	37.105;
   250	7.2480e-11	45.546;
   300	2.4180e-11	53.628;
   350	9.1580e-11	53.298;
   400	3.7250e-12	58.515;
   450	1.5850e-12	60.828;
   500	6.9670e-13	63.822;
   600	1.45e-13	   71.835;
   700	3.61e-14	   88.667;
   800	1.17e-14	   124.64;
   900	5.25e-15	   181.05;
   1000	3.02e-15	   268;];
   return A;
end





"""
      GETATMOSPHEREPROPERTIES_NRLMSISE00 
      returns atmospheric characteristic
      according to the nrlmsise00 model

      INPUT:
            mjd2000        :    epoch in modified julian day 2000
            rV             :    position vector with respect to
                              an equatorial earth centered reference frame
                              with the x-axis towards the vernal equinox
                              [km]
            solarfluxdata  :    table with solar flux data with format:
                              [jd-AP1-AP2-AP3-AP4-AP5-AP6-AP7-AP8-AP9-AP_AVG-f107_OBS-f107_OBS_CENTER81]
                              AP1	Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
                              AP2	Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
                              AP3	Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
                              AP4	Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
                              AP5	Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
                              AP6	Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
                              AP7	Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
                              AP8	Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
                              AP_AVG	Arithmetic average of the 8 Ap indices for the day.
                              F10.7_OBS	Observed 10.7-cm Solar Radio Flux (F10.7). Measured at Ottawa at 1700 UT daily from 1947 Feb 14 until 1991 May 31 and measured at Penticton at 2000 UT from 1991 Jun 01 on. Expressed in units of 10-22 W/m2/Hz.
                              F10.7_OBS_CENTER81	Centered 81-day arithmetic average of F10.7 (observed).

      OUTPUT:
               rhoinf  :   density [kg/m^3]
               Tinf    :   temperature [K]
               ainf    :   sound velocity [m/s]
               ndens   :   number densities of He, O, N2, O2, Ar, H, N

      Irene Cavallari, 01/2024

"""
function getatmosphereproperties_nrlmsise00(mjd2000,rV,solarfluxdata)
   
   # # altitude, longitude, latitude;
   # hh,lat,long = EECIposition2hlonglat(rV,mjd2000,1);
   # hh = hh*1000.0;
   # # println(hh," ",lat, " ",long)

   # epoch
   jd = mjd20002jd(mjd2000);

   # eci 2 ecef
   RM = r_eci_to_ecef(J2000(), PEF(),jd)
   geod = ecef_to_geodetic(RM*rV*1000)
   lat   = geod[1]
   long  = geod[2]
   hh    = geod[3]
   
   # solar flux and geomagnetic index info
   f107,f107A90,Ap = getSolarFluxAndGeomagneticAp(jd,solarfluxdata);

   # molar mass
   MN  = 14.0067;
   MN2 = 28.0134;
   MO  = 15.9994;
   MO2 = 31.9988;
   MH  = 1.00849;
   MHe = 4.002602;
   MAr = 39.948;
   R   = 8.314472;
   Mj = [MHe,MO,MN2,MO2,MAr,MH,MN]/1000; # kg/mol
   
   # nrlmsise00 atmospheare
   atm = AtmosphericModels.nrlmsise00(jd,hh,lat,long,f107A90,f107,Ap,include_anomalous_oxygen=false);
   
   # temperature
   Tinf  = atm.temperature;

   # exospheric temperature
   Te = atm.exospheric_temperature
   
   # density
   rhoinf = atm.total_density;
   
   # number densities 
   ndens = [atm.He_number_density,atm.O_number_density,atm.N2_number_density,atm.O2_number_density, atm.Ar_number_density,atm.H_number_density,atm.N_number_density];
   
   # sound velocity
   # M = (atm.He_number_density*MHe + atm.O_number_density*MO + atm.N2_number_density*MN2 + atm.O2_number_density * MO2 + atm.Ar_number_density * MAr + atm.H_number_density * MH + atm.N_number_density*MN)/1000.0;
   # M = M/sum(ndens)
   M = sum(ndens.*Mj)/sum(ndens);
   # println(Tinf)
   ainf = sqrt(2*R*Tinf/M);

   return rhoinf,Tinf,ainf,ndens,Te;   
   
end
   
"""
      getSolarFluxAndGeomagneticAp extract solar flux and geomagnetic field dataSince

         INPUT
            jd            : epoch in julian day
            solarfluxdata : table with solar flux data with format:
                              [jd-AP1-AP2-AP3-AP4-AP5-AP6-AP7-AP8-AP9-AP_AVG-f107_OBS-f107_OBS_CENTER81]
                              AP1	Planetary Equivalent Amplitude (Ap) for 0000-0300 UT.
                              AP2	Planetary Equivalent Amplitude (Ap) for 0300-0600 UT.
                              AP3	Planetary Equivalent Amplitude (Ap) for 0600-0900 UT.
                              AP4	Planetary Equivalent Amplitude (Ap) for 0900-1200 UT.
                              AP5	Planetary Equivalent Amplitude (Ap) for 1200-1500 UT.
                              AP6	Planetary Equivalent Amplitude (Ap) for 1500-1800 UT.
                              AP7	Planetary Equivalent Amplitude (Ap) for 1800-2100 UT.
                              AP8	Planetary Equivalent Amplitude (Ap) for 2100-0000 UT.
                              AP_AVG	Arithmetic average of the 8 Ap indices for the day.
                              F10.7_OBS	Observed 10.7-cm Solar Radio Flux (F10.7). Measured at Ottawa at 1700 UT daily from 1947 Feb 14 until 1991 May 31 and measured at Penticton at 2000 UT from 1991 Jun 01 on. Expressed in units of 10-22 W/m2/Hz.
                              F10.7_OBS_CENTER81	Centered 81-day arithmetic average of F10.7 (observed).
         OUTPUT
               f107        : 10.7-cm daily solar flux
               f107A81     : 10.7-cm averaged solar flux, 81-day centered on input time.
               f107A90     : 10.7-cm averaged solar flux, 90-day centered on input time.
               Ap          : daily magnetic index.

         Irene Cavallari, 01/2024

"""
function getSolarFluxAndGeomagneticAp(jd,solarfluxdata)
   jd0 = floor(jd);
    
   if jd-jd0>=0.5
      jd0 = jd0+1;
   end
   
   idx = solarfluxdata[:,1] .== jd0;
   info = solarfluxdata[idx,:];
   f107 = info[11];
   # f107A81 = info[12];
   f107A90 = info[13];
   Ap      = info[10]; 
   #sum(info[2:9])/8;

   # jd0Min = jd0-44.5;
   # jd0Max = jd0+44.5;
   
   # if jd0Min > minimum(solarfluxdata[:,1]) && jd0Max < maximum(solarfluxdata[:,1])
   #    idx = solarfluxdata[:,1].>=jd0Min .&& solarfluxdata[:,1].<=jd0Max
   #    f107int = solarfluxdata[idx,11];
   #    f107A90 = sum(f107int)/length(f107int);        
   # else
   #    Base.error("out of range: the data do not include the selected epoch");
   # end

   return f107,f107A90,Ap;
end

