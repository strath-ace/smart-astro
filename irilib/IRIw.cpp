/* A wrapper program to call the IRI-2012 fortran code */
#include "./IRI_wrapper.h"
#include <iostream>
#include <map>

namespace IRI
{
    int JF[50];
    float OAR[100];
    int iut;
    int jmag;

    bool hasKey(std::map<std::string, bool> mp, std::string key){
        return (mp.find(key) != mp.end());
    }

    void initializeOAR(){
        for(int i = 0; i < 100; i++){
            OAR[i] = -1.0;
        }
    }

    void initializeConst2(){
        // initialize common block /consts/
        // is this really necessary?
        const2_.icall = 0;
        const2_.nmono = -1;
        const2_.iyearo = -1;
        const2_.idaynro = -1;
        const2_.rzino = false;
        const2_.igino = false;
        const2_.ut0 = -1.0;
    }

    int initializeIRI(std::map<std::string, bool> options) {
        initializeConst2();

        read_ig_rz_();
        readapf107_();

        initializeOAR();

        for(int i = 0; i < 50; i++){
            JF[i] = 1;
        }
        // default JF parameter
        JF[3]  = 0;
        JF[4]  = 0;
        JF[5]  = 0;
        JF[20] = 0;
        JF[21] = 0; // give the actual ion density
        JF[22] = 0;
        JF[27] = 0;
        JF[28] = 0;
        JF[29] = 0;
        JF[32] = 0; // auroral boundary model
        JF[33] = 0; // message off
        JF[34] = 0;

        // check options
        std::string key = "universal_time";
        if(hasKey(options, key) && options[key]) {
            iut = 1;
        } else {
            iut = 0;
        }
        key = "geomagnetic";
        if(hasKey(options, key) && options[key]) {
            jmag = 1;
        } else {
            jmag = 0;
        }
        key = "auroral_boundary";
        if(hasKey(options, key) && options[key]) {
            JF[32] = 1;
        }
        key = "foE_storm";
        if(hasKey(options, key) && options[key]) {
            JF[34] = 1;
        }
        key = "hmF2_foF2_storm";
        if(hasKey(options, key) && options[key]) {
            JF[35] = 0; // true/false => without / with
        }
        key = "topside_foF2_storm";
        if(hasKey(options, key) && options[key]) {
            JF[36] = 0; // true/false => without / with
        }

        return 0;
    }

    void checkInputs(float height, float latitude, float longitude, int year, int mmdd, float dhour){
        if(height < 80.0) {
            std::cout << "[WARNING] height must be higher than 80km." << std::endl;
        }
        if(latitude < -90.0 || latitude > 90.0) {
            std::cout << "[WARNING] latitude must be in the range from -90 deg. to 90 deg." << std::endl;
        }
        if(longitude < -180.0 || longitude > 360) {
            std::cout << "[WARNING] longitude must be in the range from 0 deg. to 360 deg. (If -180 to 0, it will be converted to 0 to 360)" << std::endl;
        }
        if((year > 2019) && (height < 2000.0)){
            std::cout << "[WARNING] in LEO (alt. < 2000km) year must be lower than 2020." << std::endl;
        }
        if(mmdd > 1231 || mmdd < -365.0) {
            std::cout << "[WARNING] mmdd must be in the range from -365 to 1231." << std::endl;
        }
        if(dhour > 24.0) {
            std::cout << "[WARNING] dhour must be lower than 24.0." << std::endl;
        }
    }

    // IRIresult getIRIProfileAt(float height, float latitude, float longitude, int year, int mmdd, float dhour, int iut, int jmag) {
    IRIresult getIRIProfileAt(float height, float latitude, float longitude, int year, int mmdd, float dhour, std::map<std::string, bool> options) {
        static bool initialized = initializeIRI(options);
        float xhour = dhour + 25.0 * iut;
        float OUTF[1000][20];
        float height_step  = 1.0;

        // checking the values
        checkInputs(height, latitude, longitude, year, mmdd, dhour);
        IRIresult res;

        if(height < 2000.0f) {
            iri_sub_(JF, &jmag, &latitude, &longitude, &year, &mmdd, &xhour, &height, &height, &height_step, OUTF, OAR);
            res.Ne = OUTF[0][0]; // electron density
            res.Te = TemperatureToElectronVolts(OUTF[0][3]); // electron temperature
            res.Ti = TemperatureToElectronVolts(OUTF[0][2]); // ion temperature
            res.NO = OUTF[0][4]; // O+ density
            res.NH = OUTF[0][5]; // H+ density
            res.NHe = OUTF[0][6]; // He+ density
            res.NO2 = OUTF[0][7]; // O2+ density
        } else {
            // if in MEO region, iri_sub_ must not be called.
            res.Ne = 0.0;
            res.Te = 0.0;
            res.Ti = 0.0;
            res.NO = 0.0;
            res.NH = 0.0;
            res.NHe = 0.0;
            res.NO2 = 0.0;
        }

        // save the input infomations
        res.dhour = dhour;
        res.iut = iut;
        res.height = height;
        res.longi = longitude;
        res.lati = latitude;

        return res;
    }

    float TemperatureToElectronVolts(float temp){
        const float kB = 1.3806488e-23;
        const float e = 1.60217657e-19;

        return temp * kB/e;
    }

    void IRIresult::print(){
        std::cout << "Te = " << Te << std::endl;
        std::cout << "Ti = " << Ti << std::endl;
        std::cout << "Ne = " << Ne << std::endl;
        std::cout << "NH+ = " << NH << std::endl;
        std::cout << "NO = " << NO << std::endl;
        std::cout << "NHe = " << NHe << std::endl;
        std::cout << "NO2 = " << NO2 << std::endl;
    }
}
