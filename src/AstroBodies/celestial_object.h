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
			public:
				Celestial_Object( std::string givenName, int givenId,  std::vector<double> &Givenpositn,  double givenMu);

				~Celestial_Object();

				/*void latsrf_(const std::string &method, const std::string &target, double &et, const std::string &fixref, int &npts, const std::vector<std::vector::<double>> &lonlat, std::vector<std::vector::<double>> &srfpts);

				int srfnrm_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &srfpts,  std::vector<double> &normls, int method_len, int target_len, int fixref_len);*/

				void nearpt(double &a, double &b, double &c, std::vector<double> &npoint, double &alt);

				void surfpt(const std::vector<double> &u, double &a, double &b, double &c,  std::vector<double> &point, int &found);

				void surfnm(double &a, double &b, double &c, const std::vector<double> &point,  std::vector<double> &normal);

				/*int edlimb_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &viewpt,  std::vector<double> &limb);

				int inelpl_( std::vector<double> &ellips,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt1,  std::vector<double> &xpt2);*/

				void npedln(double &a, double &b, double &c, const std::vector<double> &linept, const std::vector<double> &linedr, std::vector<double> &pnear, double &dist);

				/*int inrypl_( std::vector<double> &vertex,  std::vector<double> &dir,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt);

				int pjelpl_( std::vector<double> &elin,  std::vector<double> &plane,  std::vector<double> &elout);*/

				void saelgv(const std::vector<double> &vec1, const std::vector<double> &vec2,  std::vector<double> &smajor,  std::vector<double> &sminor);

				void nplnpt(const std::vector<double> &linpt, const std::vector<double> &lindir, const std::vector<double> &point,  std::vector<double> &pnear, double &dist);

				void sincpt(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::string &dref, const std::vector<double> &dvec,  std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, int &found);

				/*int dskxv_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &nrays,  std::vector<double> &vtxarr,  std::vector<double> &dirarr,  std::vector<double> &xptarr,  std::vector<int> &fndarr, int target_len, int fixref_len);

				int dskxsi_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<double> &vertex,  std::vector<double> &raydir,  std::vector<int> &maxd,  std::vector<int> &maxi,  std::vector<double> &xpt,  std::vector<int> &handle,  std::vector<int> &dladsc,  std::vector<double> &dskdsc,  std::vector<double> &dc,  std::vector<int> &ic,  std::vector<int> &found, int target_len, int fixref_len);*/

				void ilumin(const std::string &method, const std::string &target, double &et, const std::string &fixref, const std::string &abcorr, const std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn);

				void illumg(const std::string &method, const std::string &target, const std::string &ilusrc, double &et, const std::string &fixref, const std::string &abcorr, std::vector<double> &spoint, double &trgepc,  std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn);

				void illumf(const std::string &method, const std::string &target, const std::string &ilusrc, double &et, const std::string &fixref, const std::string &abcorr, std::vector<double> &spoint, double &trgepc, std::vector<double> &srfvec, double &phase, double &incdnc, double &emissn, int &visibl, int &lit);

				/*int limbpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,    std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &tangts, int method_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len);

				int termpt_( std::string &method,  std::string &ilusrc,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,    std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &trmvcs, int method_len, int ilusrc_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len);*/

				void conics(const std::vector<double> &elts, double &et,  std::vector<double> &state);

				void oscelt(const std::vector<double> &state, double &et, std::vector<double> &elts);

				void occult(const std::string &targ1, const std::string &shape1, const std::string &frame1, const std::string &targ2, const std::string &shape2, const std::string &frame2, const std::string &abcorr, double &et, int &ocltid);

				/*int gfoclt_( std::string &occtyp,  std::string &front,  std::string &fshape,  std::string &fframe,  std::string &back,  std::string &bshape,  std::string &bframe,  std::string &abcorr,    std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int occtyp_len, int front_len, int fshape_len, int fframe_len, int back_len, int bshape_len, int bframe_len, int abcorr_len, int obsrvr_len);*/

				void srfcss(int &code, int srflen, std::string &srfstr, int &isname);

				void srfs2c(const std::string &srfstr, int &code, int &found);

				void srfc2s(int &code, int srflen, std::string &srfstr, int &isname);

				void srfscc(const std::string &srfstr, int &code, int &found);

				void bodc2n(int &lenout, int &found);

				void srfrec(double &longitude, double &latitude,  std::vector<double> &rectan);

			protected:
				std::vector<double> positn;

				double mu;
		};
	}
}

#endif // SMARTASTRO_CELESTIAL_OBJECT_H
