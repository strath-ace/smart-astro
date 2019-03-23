/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
-------Copyright (C) 2019 University of Strathclyde and Authors-------
--------- Author: Scott Hurley (GitHub: Scott_James_Hurley)-----------
-------- e-mail: scott.james.hurley.97@gmail.com ---------------------
*/

#ifndef SMARTASTO_CELESTIAL_OBJECT_H
#define SMARTASTO_CELESTIAL_OBJECT_H

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies
	{
		class Celestial_Object : public Astro_Body
		{

			protected:

				std::vector<double> positn;

				double mu;

			public:

				/**
         			* Constructor for providing a name and a spiceEphemeris
         			*
         			* @param givenName Name of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
         			*/

				Celestial_Object(std::string givenName, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn,  double givenMu);

				/**
         			* Constructor for providing a name and an integratedEphemeris
         			*
         			* @param givenName Name of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
         			*/

				Celestial_Object(std::string givenName, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn,  double givenMu);

				/**
         			* Constructor for providing an ID and a spiceEphemeris
         			*
         			* @param givenId ID of astro body
				* @param spiceEphemeridesParams Parameters to create spiceEphemeris. Format: referenceFrame, kernelsToLoad, target, abberrationCorrection
         			*/

				Celestial_Object(int givenId, std::vector<std::string> spiceEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				/**
         			* Constructor for providing an ID and an integratedEphemeris
         			*
         			* @param givenId ID of astro body
				* @param referenceFrame Reference frame for integratedEphemeris
				* @param integratedEphemeridesFunction Integration function for integratedEphemeris
				* @param integratedEphemeridesParams  Parameters to create integratedEphemeris. Format: ti, xi, hstep
         			*/

				Celestial_Object(int givenId, std::string referenceFrame, std::function<int(double,double,int,std::vector<double>,std::vector<double>)> integratedEphemeridesFunction, std::vector<double> integratedEphemeridesParams, std::vector<double> &Givenpositn, double givenMu);

				virtual ~Celestial_Object() = default;

				/**
				* @brief Map array of planetocentric longitude/latitude coordinate pairs 
   				* to surface points on a specified target body. The surface of the target 
				* body may be represented by a triaxial ellipsoid or by topographic data 
				* provided by DSK files. 
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] npts Number of coordinate pairs in input array.
				* @param[in] lonlat Array of longitude/latitude coordinate pairs.
				* @param[out] srfpts Array of surface points.
				* @return void
				*
				*/

				void latsrf(const std::string &method, const std::string &target, double &et, const std::string &fixref, const int &npts, const std::vector<double> &lonlat, std::vector<double> &srfpts);

				/**
				* @brief Map array of surface points on a specified target body to 
   				* the corresponding unit length outward surface normal vectors. 
 				* The surface of the target body may be represented by a triaxial 
   				* ellipsoid or by topographic data provided by DSK files. 
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] npts Number of coordinate pairs in input array.
				* @param[in] srfpts Array of surface points.
				* @param[out] normls Array of outward, unit length normal vectors.
				* @param[param] SPICE_DSKTOL_PTMEMM Default point-surface membership margin.
				* @return void
				*
				*/

				void srfnrm(const std::string &method, const std::string &target, double &et, const std::string &fixref, const int &npts, const std::vector<double> &srfpts, std::vector<double> &normls);

				/**
				* @brief This routine locates the point on the surface of an ellipsoid
   				* that is nearest to a specified position. It also returns the
   				* altitude of the position above the ellipsoid. 
				*
				* @param[in] positn Position of a point in bodyfixed frame.
				* @param[in] a Length of semi-axis parallel to x-axis.
				* @param[in] b Length of semi-axis parallel to y-axis.
				* @param[in] c Length on semi-axis parallel to z-axis.
				* @param[out] npoint Point on the ellipsoid closest to positn.
				* @param[out] alt Altitude of positn above the ellipsoid.
				* @return void
				*
				*/

				void nearpt(const double &a, const double &b, const double &c, std::vector<double> &npoint, double &alt);

				/**
				* @brief Determine the intersection of a line-of-sight vector with the 
   				* surface of an ellipsoid.
				*
				* @param[in] positn Position of the observer in body-fixed frame. 
				* @param[in] u Vector from the observer in some direction.
				* @param[in] a Length of the ellipsoid semi-axis along the x-axis.
				* @param[in] b Length of the ellipsoid semi-axis along the y-axis.
				* @param[in] c Length of the ellipsoid semi-axis along the z-axis.
				* @param[out] point Point on the ellipsoid pointed to by u.
				* @param[out] found Flag indicating if u points at the ellipsoid.
				* @return void
				*
				*/

				void surfpt(const std::vector<double> &u, const double &a, const double &b, const double &c, std::vector<double> &point, int &found);

				/**
				* @brief This routine computes the outward-pointing, unit normal vector 
   				* from a point on the surface of an ellipsoid.
				*
				* @param[in] a Length of the ellipsoid semi-axis along the x-axis.
				* @param[in] b Length of the ellipsoid semi-axis along the y-axis.
				* @param[in] c Length of the ellipsoid semi-axis along the z-axis.
				* @param[in] point Body-fixed coordinates of a point on the ellipsoid 
				* @param[out] normal Outward pointing unit normal to ellipsoid at point
				* @return void
				*
				*/

				void surfnm(const double &a, const double &b, const double &c, const std::vector<double> &point, std::vector<double> &normal);

				/**
				* @brief Find the limb of a triaxial ellipsoid, viewed from a specified 
   				* point
				*
				* @param[in] a Length of ellipsoid semi-axis lying on the x-axis. 
				* @param[in] b Length of ellipsoid semi-axis lying on the y-axis.
				* @param[in] c Length of ellipsoid semi-axis lying on the z-axis.
				* @param[in] viewpt Location of viewing point. 
				* @param[out] limb Limb of ellipsoid as seen from viewing point.
				* @return void
				*
				*/

				void edlimb(const double &a, const double &b, const double &c, const std::vector<double> &viewpt,  std::vector<double> &limb);

				/**
				* @brief Find the intersection of an ellipse and a plane.
				*
				* @param[in] ellips A CSPICE ellipse.
				* @param[in] plane A CSPICE plane. 
				* @param[in] spicePlaneType Type of SpicePlane. 1 for Normal vector and constant, 
				* 2 for Normal vector and point, 3 for Point and spanning vectors to plane.
				* @param[out] nxpts Number of intersection points of plane and ellipse.
				* @param[out] lxpt1 Intersection points. 
				* @param[out] lxpt2 Intersection points. 
				* @return void
				*
				*/

				void inelpl(const std::vector<double> &ellips, const std::vector<double> &plane, const int spicePlaneType, int &nxpts, std::vector<double> &xpt1, std::vector<double> &xpt2);

				/**
				* @brief Find nearest point on a triaxial ellipsoid to a specified line, 
   				* and the distance from the ellipsoid to the line.
				*
				* @param[in] a Length of ellipsoid's semi-axis in the x direction
				* @param[in] b Length of ellipsoid's semi-axis in the y direction
				* @param[in] c Length of ellipsoid's semi-axis in the z direction
				* @param[in] linept Point on line 
				* @param[in] linedr Direction vector of line 
				* @param[out] pnear Nearest point on ellipsoid to line 
				* @param[out] dist Distance of ellipsoid from line 
				* @return void
				*
				*/

				void npedln(const double &a, const double &b, const double &c, const std::vector<double> &linept, const std::vector<double> &linedr, std::vector<double> &pnear, double &dist);

				/**
				* @brief Find the intersection of a ray and a plane. 
				*
				* @param[in] vertex Vertex vector of ray. 
				* @param[in] dir Direction vector of ray. 
				* @param[in] plane A CSPICE plane.
				* @param[in] spicePlaneType Type of SpicePlane. 1 for Normal vector and constant, 
				* 2 for Normal vector and point, 3 for Point and spanning vectors to plane.
				* @param[out] nxpts Number of intersection points of ray and plane.
				* @param[out] xpt Intersection point, if nxpts = 1. 
				* @return void
				*
				*/

				void inrypl(const std::vector<double> &vertex, const std::vector<double> &dir, const std::vector<double> &plane, const int spicePlaneType, int &nxpts, std::vector<double> &xpt);

				/**
				* @brief Project an ellipse onto a plane, orthogonally. 
				*
				* @param[in] elin A CSPICE ellipse to be projected.
				* @param[in] plane A plane onto which elin is to be projected.
				* @param[out] elout A CSPICE ellipse resulting from the projection. 
				* @return void
				*
				*/

				void pjelpl(const std::vector<double> &elin, const std::vector<double> &plane, const int spicePlaneType, std::vector<double> &elout);

				/**
				* @brief Find semi-axis vectors of an ellipse generated by two arbitrary 
   				* three-dimensional vectors.
				*
				* @param[in] vec1 Vector used to generate an ellipse.
				* @param[in] vec2 Vector used to generate an ellipse. 
				* @param[out] smajor Semi-major axis of ellipse. 
				* @param[out] sminor Semi-minor axis of ellipse.
				* @return void
				*
				*/

				void saelgv(const std::vector<double> &vec1, const std::vector<double> &vec2,  std::vector<double> &smajor, std::vector<double> &sminor);

				/**
				* @brief Find the nearest point on a line to a specified point, and find 
   				* the distance between the two points. 
				*
				* @param[in] linpt Point on a line.
				* @param[in] lindir Line's direction vector.
				* @param[in] point A second point. 
				* @param[out] pnear Nearest point on the line to point. 
				* @param[out] dist Distance between point and pnear. 
				* @return void
				*
				*/

				void nplnpt(const std::vector<double> &linpt, const std::vector<double> &lindir, const std::vector<double> &point, std::vector<double> &pnear, double &dist);

				/**
				* @brief Given an observer and a direction vector defining a ray, compute 
   				* the surface intercept of the ray on a target body at a specified 
   				* epoch, optionally corrected for light time and stellar aberration. 
 				* The surface of the target body may be represented by a triaxial 
  				* ellipsoid or by topographic data provided by DSK files. 
 				* This routine supersedes srfxpt_c. 
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag. 
				* @param[in] dref Reference frame of ray's direction vector.
				* @param[in] dvec Ray's direction vector. 
				* @param[out] spoint Surface intercept point on the target body.
				* @param[out] trgepc Intercept epoch. 
				* @param[out] srfvec Vector from observer to intercept point. 
				* @param[out] found Flag indicating whether intercept was found.
				* @return void
				*
				*/

				void sincpt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &dref, const std::vector<double> &dvec,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, int &found);

				/**
				* @brief Compute ray-surface intercepts for a set of rays, using data 
   				* provided by multiple loaded DSK segments.
				*
				* @param[in] pri Data prioritization flag. 
				* @param[in] target Name of target body.
				* @param[in] nsurf Number of surface IDs in list.
				* @param[in] srflst Surface ID list. 
				* @param[in] et Epoch, expressed as seconds past J2000 TDB. 
				* @param[in] fixref Name of target body-fixed reference frame.
				* @param[in] nrays Number of rays. 
				* @param[in] vtxarr Array of vertices of rays.
				* @param[in] dirarr Array of direction vectors of rays.
				* @param[out] xptarr Intercept point array.
				* @param[out] fndarr Found flag array. 
				* @return void
				*
				*/

				void dskxv(const int &pri, const std::string &target, const int &nsurf, const std::vector<int> &srflst, const double &et, const std::string &fixref, const int &nrays, std::vector<double> &vtxarr, const std::vector<double> &dirarr, std::vector<double> &xptarr, std::vector<int> &fndarr);

				/**
				* @brief Compute a ray-surface intercept using data provided by  
   				* multiple loaded DSK segments. Return information about  
   				* the source of the data defining the surface on which the 
   				* intercept was found: DSK handle, DLA and DSK descriptors, 
  				* and DSK data type-dependent parameters.
				*
				* @param[in] pri Data prioritization flag. 
				* @param[in] target Name of target body.
				* @param[in] nsurf Number of surface IDs in list.
				* @param[in] srflst Surface ID list. 
				* @param[in] et Epoch, expressed as seconds past J2000 TDB. 
				* @param[in] fixref Name of target body-fixed reference frame.
				* @param[in] vertex Vertex of ray. 
				* @param[in] raydir Direction vector of ray.
				* @param[in] maxd Size of DC array.
				* @param[in] maxi Size of IC array.
				* @param[out] xpt Intercept point. 
				* @param[out] handle Handle of segment contributing surface data. 
				* @param[out] dladsc DLA descriptor of segment.
				* @param[out] dskdsc DSK descriptor of segment.
				* @param[out] dc Double precision component of source info.  
				* @param[out] ic Integer component of source info. 
				* @param[out] found Found flag.
				* @param[param] SPICE_DSKXSI_DCSIZE Required size of DC array.
				* @param[param] SPICE_DSKXSI_ICSIZE Required size of IC array.
				* @return void
				*
				*/

				void dskxsi(const int &pri, const std::string &target, const int &nsurf, const std::vector<int> &srflst, const double &et, const std::string &fixref, const std::vector<double> &vertex, const std::vector<double> &raydir, const int &maxd, const int &maxi, std::vector<double> &xpt, int &handle, std::vector<int> &dladsc, std::vector<double> &dskdsc, std::vector<double> &dc, std::vector<int> &ic, int &found);

				/**
				* @brief Find the illumination angles (phase, solar incidence, and
   				* emission) at a specified surface point of a target body.
 				* This routine supersedes illum_c.
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag. 
				* @param[in] spoint Body-fixed coordinates of a target surface point.
				* @param[out] trgepc Target surface point epoch.
				* @param[out] srfvec Vector from observer to target surface point.
				* @param[out] phase Phase angle at the surface point.  
				* @param[out] incdnc Solar incidence angle at the surface point. 
				* @param[out] emissn Emission angle at the surface point.
				* @return void
				*
				*/

				void ilumin(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn);

				/**
				* @brief Find the illumination angles (phase, incidence, and 
  				* emission) at a specified surface point of a target body. 
 				* The surface of the target body may be represented by a triaxial 
				* ellipsoid or by topographic data provided by DSK files. 
 				* The illumination source is a specified ephemeris object.
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] ilusrc Name of illumination source.
				* @param[in] et Epoch in TDB seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag. 
				* @param[in] spoint Body-fixed coordinates of a target surface point.
				* @param[out] trgepc Target surface point epoch.
				* @param[out] srfvec Vector from observer to target surface point.
				* @param[out] phase Phase angle at the surface point.  
				* @param[out] incdnc Solar incidence angle at the surface point. 
				* @param[out] emissn Emission angle at the surface point.
				* @return void
				*
				*/

				void illumg(const std::string &method, const std::string &target, const std::string &ilusrc, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn);

				/**
				* @brief Compute the illumination angles---phase, incidence, and
   				* emission---at a specified point on a target body. Return logical
   				* flags indicating whether the surface point is visible from
   				* the observer's position and whether the surface point is
   				* illuminated. The target body's surface is represented using 
				* topographic data provided by DSK files, or by a reference ellipsoid.
				* The illumination source is a specified ephemeris object.
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] ilusrc Name of illumination source.
				* @param[in] et Epoch in TDB seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction flag. 
				* @param[in] spoint Body-fixed coordinates of a target surface point.
				* @param[out] trgepc Target surface point epoch.
				* @param[out] srfvec Vector from observer to target surface point.
				* @param[out] phase Phase angle at the surface point.  
				* @param[out] incdnc Solar incidence angle at the surface point. 
				* @param[out] emissn Emission angle at the surface point.
				* @param[out] visibl Visibility flag (SPICETRUE == visible).
				* @param[out] lit Illumination flag (SPICETRUE == illuminated).
				* @return void
				*
				*/

				void illumf(const std::string &method, const std::string &target, const std::string &ilusrc, const double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn, int &visibl, int &lit);

				/**
				* @brief Find limb points on a target body. The limb is the set of points 
   				* of tangency on the target of rays emanating from the observer. 
   				* The caller specifies half-planes bounded by the observer-target 
   				* center vector in which to search for limb points. 
   				* The surface of the target body may be represented either by a 
   				* triaxial ellipsoid or by topographic data.
				*
				* @param[in] method Computation method.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction.
				* @param[in] corloc Aberration correction locus. 
				* @param[in] refvec Reference vector for cutting half-planes.
				* @param[in] rolstp Roll angular step for cutting half-planes.
				* @param[in] ncuts Number of cutting half-planes. 
				* @param[in] schstp Angular step size for searching.
				* @param[in] soltol Solution convergence tolerance.
				* @param[in] maxn Maximum number of entries in output arrays.
				* @param[out] npts Counts of limb points corresponding to cuts.
				* @param[out] points Limb points.
				* @param[out] epochs Times associated with limb points. 
				* @param[out] tangts Tangent vectors emanating from the observer. 
				* @return void
				*
				*/

				void limbpt(const std::string &method, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &corloc, const std::vector<double> &refvec, const double &rolstp, const int &ncuts, const double &schstp, const double &soltol, const int &maxn, std::vector<int> &npts, std::vector<double> &points, std::vector<double> &epochs, std::vector<double> &tangts);

				/**
				* @brief Find terminator points on a target body. The caller specifies
				* half-planes, bounded by the illumination source center-target center
				* vector, in which to search for terminator points. 
				* The terminator can be either umbral or penumbral. The umbral 
				* terminator is the boundary of the region on the target surface 
				* where no light from the source is visible. The penumbral 
				* terminator is the boundary of the region on the target surface 
				* where none of the light from the source is blocked by the target 
				* itself. 
				* The surface of the target body may be represented either by a 
				* triaxial ellipsoid or by topographic data.
				*
				* @param[in] method Computation method.
				* @param[in] ilusrc Illumination source.
				* @param[in] target Name of target body.
				* @param[in] et Epoch in ephemeris seconds past J2000 TDB.
				* @param[in] fixref Body-fixed, body-centered target body frame.
				* @param[in] abcorr Aberration correction.
				* @param[in] corloc Aberration correction locus. 
				* @param[in] refvec Reference vector for cutting half-planes.
				* @param[in] rolstp Roll angular step for cutting half-planes.
				* @param[in] ncuts Number of cutting planes. 
				* @param[in] schstp Angular step size for searching.
				* @param[in] soltol Solution convergence tolerance.
				* @param[in] maxn Maximum number of entries in output arrays.
				* @param[out] npts Counts of limb points corresponding to cuts.
				* @param[out] points Terminator points.
				* @param[out] epochs Times associated with terminator points. 
				* @param[out] trmvcs Terminator vectors emanating from the observer. 
				* @return void
				*
				*/

				void termpt(const std::string &method, const std::string &ilusrc, const std::string &target, const double &et, const std::string &fixref, const std::string &abcorr, const std::string &corloc, const std::vector<double> &refvec, const double &rolstp, const int &ncuts, const double &schstp, const double &soltol, const int &maxn,  std::vector<int> &npts, std::vector<double> &points, std::vector<double> &epochs, std::vector<double> &trmvcs);

				/**
				* @brief Determine the set of osculating conic orbital elements that
   				* corresponds to the state (position, velocity) of a body at
   				* some epoch.
				*
				* @param[in] elts Conic elements.
				* @param[in] et Input time.
				* @param[out] state State of orbiting body at et.
				* @return void
				*
				*/

				void conics(const std::vector<double> &elts, const double &et, std::vector<double> &state);

				/**
				* @brief Determine the state (position, velocity) of an orbiting body
   				* from a set of elliptic, hyperbolic, or parabolic orbital elements.
				*
				* @param[in] state tate of body at epoch of elements.
				* @param[in] et Epoch of elements.
				* @param[out] elts Equivalent conic elements
				* @return void
				*
				*/

				void oscelt(const std::vector<double> &state, const double &et, std::vector<double> &elts);

				/**
				* @brief Determines the occultation condition (not occulted, partially,
   				* etc.) of one target relative to another target as seen by
   				* an observer at a given time.
   				* The surfaces of the target bodies may be represented by triaxial
   				* ellipsoids or by topographic data provided by DSK files.
				*
				* @param[in] targ1 Name or ID of first target.
				* @param[in] shape1 Type of shape model used for first target.
				* @param[in] frame1 Body-fixed, body-centered frame for first body.
				* @param[in] targ2 Name or ID of second target.
				* @param[in] shape2 Type of shape model used for second target.
				* @param[in] frame2 Body-fixed, body-centered frame for second body.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] et Time of the observation (seconds past J2000).
				* @param[out] ocltid  Occultation identification code.
				* @return void
				*
				*/

				void occult(const std::string &targ1, const std::string &shape1, const std::string &frame1, const std::string &targ2, const std::string &shape2, const std::string &frame2, const std::string &abcorr, double &et, int &ocltid);

				/**
				* @brief Determine time intervals when an observer sees one target occulted
   				* by, or in transit across, another. The surfaces of the target bodies may 
				* be represented by triaxial ellipsoids or by topographic data provided by DSK files.
				*
				* @param[in] occtyp Type of occultation. 
				* @param[in] front Name of body occulting the other. 
				* @param[in] fshape Type of shape model used for front body. 
				* @param[in] fframe Body-fixed, body-centered frame for front body. 
				* @param[in] back Name of body occulted by the other.
				* @param[in] bshape Type of shape model used for back body. 
				* @param[in] bframe Body-fixed, body-centered frame for back body.
				* @param[in] abcorr Aberration correction flag.
				* @param[in] step Step size in seconds for finding occultation.
				* @param[in/out] cnfine SPICE window to which the search is restricted.
				* @param[out] result SPICE window containing results.
				* @param[param] SPICE_GF_CNVTOL Convergence tolerance. 
				* @return void
				*
				*/

				void gfoclt(const std::string &occtyp, const std::string &front, const std::string &fshape, const std::string &fframe, const std::string &back, const std::string &bshape, const std::string &bframe, const std::string &abcorr, const double &step, std::vector<double> &cnfine, std::vector<double> &result);

				/**
				* @brief Translate a surface ID code, together with a body string, to the 
   				* corresponding surface name. If no such surface name exists, 
   				* return a string representation of the surface ID code.
				*
				* @param[in] code Integer surface ID code to translate to a string. 
				* @param[in] srflen Length of output string `srfstr'.
				* @param[out] srfstr tring corresponding to surface ID code.
				* @param[out] isname Flag indicating whether output is a surface name.
				* @param[param] SPICE_SRF_SFNMLN Maximum length of surface name.
				*
				*/

				void srfcss(const int &code, const int srflen, std::string &srfstr, int &isname);

				/**
				* @brief Translate a surface string, together with a body string, to the 
   				* corresponding surface ID code. The input strings may contain 
   				* names or integer ID codes. 
				*
				* @param[in] srfstr Surface name or ID string. 
				* @param[out] code Integer surface ID code.
				* @param[out] found Flag indicating whether surface ID was found. 
				*
				*/

				void srfs2c(const std::string &srfstr, int &code, int &found);

				/**
				* @brief Translate a surface ID code, together with a body ID code, to the 
   				* corresponding surface name. If no such name exists, return a 
   				* string representation of the surface ID code. 
				*
				* @param[in] code Integer surface ID code to translate to a string.
				* @param[in] srflen Length of output string `srfstr'.
				* @param[out] srfstr String corresponding to surface ID code. 
				* @param[out] isname Logical flag indicating output is a surface name.
				* @param[param] SPICE_SRF_SFNMLN Maximum length of surface name.
				*
				*/

				void srfc2s(const int &code, const int srflen, std::string &srfstr, int &isname);

				/**
				* @brief Translate a surface string, together with a body ID code, to the 
   				* corresponding surface ID code. The input surface string may 
   				* contain a name or an integer ID code. 
				*
				* @param[in] srfstr Surface name or ID string.
				* @param[out] code Integer surface ID code.
				* @param[out] found Flag indicating whether surface ID was found.
				*
				*/

				void srfscc(const std::string &srfstr, int &code, int &found);

				/**
				* @brief Convert planetocentric latitude and longitude of a surface 
   				* point on a specified body to rectangular coordinates.
				*
				* @param[in] longitude Longitude of point in radians.
				* @param[in] latitude Latitude of point in radians.
				* @param[out] rectan Rectangular coordinates of the point. 
				*
				*/

				void srfrec(const double &longitude, const double &latitude,  std::vector<double> &rectan);
		};
	}
}

#endif // SMARTASTRO_CELESTIAL_OBJECT_H
