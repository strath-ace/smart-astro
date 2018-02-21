#ifndef __SCCAL_INCLUDED__
#define __SCCAL_INCLUDED__
#include "IRI_wrapper.h"
#include <iostream>
#include <fstream>
#include <math.h>
#define _USE_MATH_DEFINES

/* template library for computing spacecraft charging
 * author: Kento Hoshi
 * date  : 2016/10/07
 */

namespace SCCal {
    const double kB = 1.3806488e-23;
    const double e = 1.60217657e-19;
    const double me = 9.10938356e-31;
    const double mp = 1836.152 * me;
    const double mO = 16.0 * mp;
    const double mHe = 4.0 * mp;
    const double mO2 = 2.0 * mO;

    enum IonType {
        PROTON,
        OXYGEN,
        HELIUM,
        OXYGEN2,
    };

    template<typename T>
    struct AuroralCurrent {
        T latitude_from = 58.0;
        T latitude_to = 80.0;
        T density = 1.0e8;
        T drift_velocity = 2000.0; // eV
        T temperature = 20.0; // eV
        bool inAuroralRegion = false;
    };

    template<typename T>
    struct Spacecraft {
        T surface_area;
        T shape_factor;
        T potential;
        T cos_angle_Sun;
        T cos_angle_vel;
        T velocity; // to compute ion currents
    };

    //! utility functions
    //! @{
    template<typename T>
    T convertUTtoLT(T ut, T longitude){
        T lt = ut + longitude/15.0;
        if(lt > 24.0) lt -= 24.0;
        return lt;
    }

    // convert from "eV (electron volts)" to "velocity (m/s)"
    template<typename T>
    T convertEVtoVelocity(T ev, T mass){
        return sqrt(2.0 * ev * e/mass);
    }

    // compute true projected area for one-directional drift current
    template<typename T>
    T getSurfaceAreaForDrift(Spacecraft<T> obj){
        T area = obj.surface_area;

        if (obj.shape_factor == 1.0)
            area /= 4.0;
        if (obj.shape_factor == 0.5)
            area *= sqrt(1.0 - obj.cos_angle_vel * obj.cos_angle_vel) / M_PI;
        if (obj.shape_factor == 0.0)
            area *= obj.cos_angle_vel / 2.0;

        return area;
    }

    // compute true projected area for photoelectron current
    template<typename T>
    T getSurfaceAreaForPhoto(Spacecraft<T> obj){
        T area = obj.surface_area;

        if (obj.shape_factor == 1.0)
            area /= 4.0;
        if (obj.shape_factor == 0.5)
            area *= sqrt(1.0 - obj.cos_angle_Sun * obj.cos_angle_Sun) / M_PI;
        if (obj.shape_factor == 0.0)
            area *= obj.cos_angle_Sun / 2.0;

        return area;
    }

    //! functions for auroral current
    //! @{
    template<typename T>
    T getAuroralCurrentZero(Spacecraft<T> obj, AuroralCurrent<T> auro, T surface_area_for_aurora){
        T velocity = convertEVtoVelocity<T>(auro.drift_velocity, me);

        T Ia1 = sqrt(auro.temperature * e/(2.0 * M_PI * me)) * exp(- me * pow(velocity, 2) / (2.0 * e * auro.temperature));
        T Ia2 = (velocity / 2.0) * (1.0 + erf(sqrt(me/(e * auro.temperature)) * velocity));

        return e * auro.density * surface_area_for_aurora * (Ia1 + Ia2);
    }

    template<typename T>
    T getAuroralCurrent(Spacecraft<T> obj, AuroralCurrent<T> auro, T surface_area_for_aurora){
        T fact = 0.0;

        if(obj.potential < 0.0) {
            // res.Te is in the unit of "eV"
            // so (1 + eV/kB*Te) => (1 + V/Te) is enough
            fact = exp(obj.potential / auro.temperature);
        } else {
            fact = pow(1.0 + obj.potential / auro.temperature, obj.shape_factor);
        }

        return getAuroralCurrentZero<T>(obj, auro, surface_area_for_aurora) * fact;
    }

    template<typename T>
    T getAuroralCurrentPrime(Spacecraft<T> obj, AuroralCurrent<T> auro, T surface_area_for_aurora){
        T Ia0 = getAuroralCurrentZero<T>(obj, auro, surface_area_for_aurora);
        T e_kBTa = 1.0/auro.temperature; // because of res.Te is in eV

        if(obj.potential > 0.0) {
            if(obj.shape_factor == 1.0){
                return e_kBTa * Ia0;
            } else if (obj.shape_factor == 0.0) {
                return 0.0;
            } else {
                return obj.shape_factor * e_kBTa * Ia0 * pow(1.0 + obj.potential/auro.temperature, obj.shape_factor - 1.0);
            }
        } else {
            return e_kBTa * getAuroralCurrent<T>(obj, auro, surface_area_for_aurora);
        }
    }
    // ! @}

    //! functions for photoelectron current
    //! @{
    template<typename T>
    T getPhotoelectronCurrent(IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction){
        const T Jph = 3.0e-5; // current density, A/m^2 from Geotail paper

        return Jph * getSurfaceAreaForPhoto(obj) * ShadowFunction;
    }

    template<typename T>
    T getPhotoelectronCurrentPrime(IRI::IRIresult res, Spacecraft<T> obj){
        return 0.0;
    }
    //! @}

    template<typename T>
    T getElectronCurrent(IRI::IRIresult res, Spacecraft<T> obj){
        T v   = convertEVtoVelocity<T>(res.Te, me);

        // the coeff. 1/(2 * sqrt(M_PI)) is for a sphere, from the paper of J. E. Allen 1992
        // it is 2/sqrt(M_PI) for a cylinder from the same paper,
        // but there is no description for a plate
        // For consistency, we should use no coefficient for the original OML theory by Mott-Smith

        T Ie0 = obj.surface_area * res.Ne * e * v / (sqrt(M_PI) * 2.0);
        // T Ie0 = obj.surface_area * res.Ne * e * v;

        T fact = 0.0;

        if(obj.potential < 0.0) {
            // res.Te is in the unit of "eV"
            // so (1 + eV/kB*Te) => (1 + V/Te) is enough
            fact = exp(obj.potential / res.Te);
        } else {
            fact = pow(1.0 + obj.potential / res.Te, obj.shape_factor);
        }

        return Ie0 * fact;
    }

    template<typename T>
    T getElectronCurrentPrime(IRI::IRIresult res, Spacecraft<T> obj){
        T v   = convertEVtoVelocity<T>(res.Te, me);

        // T Ie0 = obj.surface_area * res.Ne * e * v;
        T Ie0 = obj.surface_area * res.Ne * e * v / (sqrt(M_PI) * 2.0);
        T e_kBTe = 1.0/res.Te; // because of res.Te is in eV

        if(obj.potential > 0.0) {
            if(obj.shape_factor == 1.0){
                return e_kBTe * Ie0;
            } else if (obj.shape_factor == 0.0) {
                return 0.0;
            } else {
                return obj.shape_factor * e_kBTe * Ie0 * pow(1.0 + obj.potential/res.Te, obj.shape_factor - 1.0);
            }
        } else {
            return e_kBTe * getElectronCurrent<T>(res, obj);
        }
    }

    template<typename T>
    T getIonDensity(IRI::IRIresult res, IonType ionType){
        switch(ionType){
            case PROTON:
                return res.NH;
                break;
            case OXYGEN:
                return res.NO;
                break;
            case HELIUM:
                return res.NHe;
                break;
            case OXYGEN2:
                return res.NO2;
                break;
            default:
                return res.NH;
                break;
        }
    }

    template<typename T>
    T getIonMass(IonType ionType){
        switch(ionType){
            case PROTON:
                return mp;
                break;
            case OXYGEN:
                return mO;
                break;
            case HELIUM:
                return mHe;
                break;
            case OXYGEN2:
                return mO2;
                break;
            default:
                return mp;
                break;
        }
    }

    template<typename T>
    T getIonThermalCurrent(IRI::IRIresult res, Spacecraft<T> obj, IonType ionType){
        T mass = getIonMass<T>(ionType);
        T v   = convertEVtoVelocity<T>(res.Ti, mass);

        // the coeff. 1/(2 * sqrt(M_PI)) is for a sphere, from the paper of J. E. Allen 1992
        // it is 2/sqrt(M_PI) for a cylinder from the same paper,
        // but there is no description for a plate
        // For consistency, we should use no coefficient for the original OML theory by Mott-Smith

        T Ii0 = obj.surface_area * getIonDensity<T>(res, ionType) * e * v / (sqrt(M_PI) * 2.0);
        // T Ii0 = obj.surface_area * getIonDensity<T>(res, ionType) * e * v;

        T fact = 0.0;

        if(obj.potential > 0.0) {
            fact = exp(-obj.potential / res.Ti);
        } else {
            fact = pow(1.0 - obj.potential / res.Ti, obj.shape_factor);
        }

        return Ii0 * fact;
    }

    template<typename T>
    T getIonThermalCurrentPrime(IRI::IRIresult res, Spacecraft<T> obj, IonType ionType){
        T mass = getIonMass<T>(ionType);
        T v   = convertEVtoVelocity<T>(res.Ti, mass);
        T Ii0 = obj.surface_area * getIonDensity<T>(res, ionType) * e * v / (sqrt(M_PI) * 2.0);
        T e_kBTi = 1.0/res.Ti; // because of res.Te is in eV

        if(obj.potential < 0.0) {
            if (obj.shape_factor == 1.0) {
                return -Ii0 * e_kBTi;
            } else if (obj.shape_factor == 0.0) {
                return 0.0;
            } else {
                return -obj.shape_factor * e_kBTi * Ii0 * pow(1.0 - obj.potential/res.Ti, obj.shape_factor - 1.0);
            }
        } else {
            return -e_kBTi * getIonThermalCurrent<T>(res, obj, ionType);
        }
    }

    template<typename T>
    T getIonDriftCurrent(IRI::IRIresult res, Spacecraft<T> obj, IonType ionType) {
        return getSurfaceAreaForDrift<T>(obj) * getIonDensity<T>(res, ionType) * e * obj.velocity;
    }

    template<typename T>
    T getIonCurrentLEO(IRI::IRIresult res, Spacecraft<T> obj, IonType ionType){
        T drift_current = getIonDriftCurrent(res, obj, ionType);
        T thermal_current = getIonThermalCurrent(res, obj, ionType);
        return drift_current + thermal_current;
    }

    template<typename T>
    T getIonCurrentPrimeLEO(IRI::IRIresult res, Spacecraft<T> obj, IonType ionType){
        const double drift_current_prime = 0.0;
        return drift_current_prime + getIonThermalCurrentPrime(res, obj, ionType);
    }

    template<typename T>
    void plotCurrent(std::string filename, IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction, AuroralCurrent<T> auro) {
        static bool isFirstCall = true;
        std::ofstream ofs;

        if(isFirstCall) {
            ofs.open(filename, std::ios::trunc);
            ofs << "## electron, proton, oxygen, helium, o2";

            if(auro.inAuroralRegion) ofs << ", auroral";

            ofs << ", photoelectron" << std::endl;
            isFirstCall = false;
        } else {
            ofs.open(filename, std::ios::app);
        }

        // positive: inflow currents (ions, photoelectron)
        // negative: outflow currents (electron)

        // electron current
        ofs << -1.0 * getElectronCurrent<T>(res, obj) << " ";

        // proton current
        ofs << getIonCurrentLEO<T>(res, obj, PROTON) << " ";

        // oxygen current
        ofs << getIonCurrentLEO<T>(res, obj, OXYGEN) << " ";

        // helium current
        ofs << getIonCurrentLEO<T>(res, obj, HELIUM) << " ";

        // o2 current
        ofs << getIonCurrentLEO<T>(res, obj, OXYGEN2) << " ";

        if (auro.inAuroralRegion)
        {
            ofs << -1.0 * getAuroralCurrent<T>(obj, auro, obj.surface_area / 4.0) << " ";
        }

        // photoelectron current
        ofs << getPhotoelectronCurrent<T>(res, obj, ShadowFunction) << " ";

        ofs << std::endl;
    }

    template<typename T>
    T solveCurrentBalance(IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction, bool considerAuroral = false, T eps = 1e-7, int LIMIT = 10000){
        T prev = obj.potential + 100.0; // initialize
        int itr = 0;
        T fx = 0.0; T fx_prime = 0.0;
        AuroralCurrent<T> auro;

        auro.inAuroralRegion = (considerAuroral) && (fabs(res.lati) >= auro.latitude_from) && (fabs(res.lati) <= auro.latitude_to);

        while( (fabs(prev - obj.potential) > eps && itr < LIMIT) ){
            // electron current
            fx = getElectronCurrent<T>(res, obj);
            fx_prime = getElectronCurrentPrime<T>(res, obj);

            // proton current
            if(res.NH > 0) {
                fx -= getIonThermalCurrent<T>(res, obj, PROTON);
                fx_prime -= getIonThermalCurrentPrime<T>(res, obj, PROTON);
            }

            // oxygen current
            if(res.NO > 0) {
                fx -= getIonThermalCurrent<T>(res, obj, OXYGEN);
                fx_prime -= getIonThermalCurrentPrime<T>(res, obj, OXYGEN);
            }

            // helium current
            if(res.NHe > 0) {
                fx -= getIonThermalCurrent<T>(res, obj, HELIUM);
                fx_prime -= getIonThermalCurrentPrime<T>(res, obj, HELIUM);
            }

            // o2 current
            if(res.NO > 0) {
                fx -= getIonThermalCurrent<T>(res, obj, OXYGEN2);
                fx_prime -= getIonThermalCurrentPrime<T>(res, obj, OXYGEN2);
            }

            if(auro.inAuroralRegion){
                // for sphere
                fx += getAuroralCurrent<T>(obj, auro, obj.surface_area/4.0);
                fx_prime += getAuroralCurrentPrime<T>(obj, auro, obj.surface_area/4.0);
            }

            // photoelectron current
            fx -= getPhotoelectronCurrent<T>(res, obj, ShadowFunction);
            fx_prime -= getPhotoelectronCurrentPrime<T>(res, obj);

            prev = obj.potential;
            obj.potential = obj.potential - (fx/fx_prime);
            itr++;
        }

        return obj.potential;
    }

    template<typename T>
    T solveCurrentBalanceMEO(IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction){
        // ref. L. K. Sarno-Smith et at., 2016
        // http://dx.doi.org/10.1002/2015SW001345
        double light = 5.0;
        double no_light = -5.0;

        return no_light + ShadowFunction*(light-no_light);
    }

    template<typename T>
    T solveCurrentBalanceLEO(IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction, bool printCurrentData = false, bool considerAuroral = false, T eps = 1e-7, int LIMIT = 10000){
        T prev = obj.potential + 100.0; // initialize
        int itr = 0;
        T fx = 0.0; T fx_prime = 0.0;
        AuroralCurrent<T> auro;

        auro.inAuroralRegion = (considerAuroral) && (fabs(res.lati) >= auro.latitude_from) && (fabs(res.lati) <= auro.latitude_to);

        while( (fabs(prev - obj.potential) > eps && itr < LIMIT) ){
            // electron current
            fx = getElectronCurrent<T>(res, obj);
            fx_prime = getElectronCurrentPrime<T>(res, obj);

            // proton current
            if(res.NH > 0) {
                fx -= getIonCurrentLEO<T>(res, obj, PROTON);
                fx_prime -= getIonCurrentPrimeLEO<T>(res, obj, PROTON);
            }

            // oxygen current
            if(res.NO > 0) {
                fx -= getIonCurrentLEO<T>(res, obj, OXYGEN);
                fx_prime -= getIonCurrentPrimeLEO<T>(res, obj, OXYGEN);
            }

            // helium current
            if(res.NHe > 0) {
                fx -= getIonCurrentLEO<T>(res, obj, HELIUM);
                fx_prime -= getIonCurrentPrimeLEO<T>(res, obj, HELIUM);
            }

            // o2 current
            if(res.NO > 0) {
                fx -= getIonCurrentLEO<T>(res, obj, OXYGEN2);
                fx_prime -= getIonCurrentPrimeLEO<T>(res, obj, OXYGEN2);
            }

            if(auro.inAuroralRegion){
                // for sphere
                fx += getAuroralCurrent<T>(obj, auro, obj.surface_area/4.0);
                fx_prime += getAuroralCurrentPrime<T>(obj, auro, obj.surface_area/4.0);
            }

            // photoelectron current
            fx -= getPhotoelectronCurrent<T>(res, obj, ShadowFunction);
            fx_prime -= getPhotoelectronCurrentPrime<T>(res, obj);

            prev = obj.potential;
            obj.potential = obj.potential - (fx/fx_prime);
            itr++;
        }

        if(printCurrentData) plotCurrent<T>("current_history.txt", res, obj, ShadowFunction, auro);

        return obj.potential;
    }

    template<typename T>
    T solve(IRI::IRIresult res, Spacecraft<T> obj, double ShadowFunction, bool plotCurrentData = false, bool considerAuroral = false, T eps = 1e-7, int LIMIT = 10000){
        if(res.height > 2000.0f) {
            return solveCurrentBalanceMEO(res, obj, ShadowFunction);
        } else {
            return solveCurrentBalanceLEO(res, obj, ShadowFunction, plotCurrentData, considerAuroral, eps, LIMIT);
        }
    }

    template<typename T>
    T solveCurrentBalanceAnalytical(IRI::IRIresult res){
        return res.Te * log(sqrt(me/mp) * sqrt(res.Ti/res.Te));
    }

}
#endif
