#ifndef __IRI_WRAPPER_INCLUDED__
#define __IRI_WRAPPER_INCLUDED__
#include <map>
#include <string>

extern "C" {
    void iri_sub_(
        int JF[],
        int *JMAG,
        float *ALATI,
        float *ALONG,
        int *IYYYY,
        int *MMDD,
        float *DHOUR,
        float *HEIBEG,
        float *HEIEND,
        float *HEISTP,
        float OUTF[][20],
        float OAR[]
    );

    void read_ig_rz_(void);
    void readapf107_(void);

    // IRI common blocks
    extern struct {
        int icall;
        int nmono;
        int iyearo;
        int idaynro;
        int rzino;
        int igino;
        float ut0;
    } const2_;
}

namespace IRI
{
    struct IRIresult {
        float Ne;
        float Te;
        float Ti;
        float NO;
        float NH;
        float NHe;
        float NO2;

        // to compute the Jph
        float dhour;
        int iut;
        float height;
        float longi;
        float lati;

        void print();
    };

    //int initializeIRI();
    // IRIresult getIRIProfileAt(float height, float latitude, float longitude, int year, int mmdd, float dhour, int iut, int jmag);
    IRIresult getIRIProfileAt(float height, float latitude, float longitude, int year, int mmdd, float dhour, std::map<std::string, bool> options = {});
    float TemperatureToElectronVolts(float temp);
}
#endif
