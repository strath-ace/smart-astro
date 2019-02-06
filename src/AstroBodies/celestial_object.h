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

extern "C" {
#include <cspice/SpiceUsr.h>
#include <cspice/SpiceZfc.h>
}

#include "astro_body.h"

namespace smartastro
{
	namespace astrobodies 
	{
		class Celestial_Object : public Astro_Body
		{
			public:
				Celestial_Object( std::string givenName, int givenId,  std::vector<double> &Givenpositn,  std::vector<double> &givenMu);

				~Celestial_Object();

				int latsrf_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &lonlat,  std::vector<double> &srfpts, int method_len, int target_len, int fixref_len);

				int srfnrm_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &npts,  std::vector<double> &srfpts,  std::vector<double> &normls, int method_len, int target_len, int fixref_len);

				int nearpt_( std::vector<double> &positn,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &npoint,  std::vector<double> &alt);

				int dnearp_( std::vector<double> &state,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &dnear,  std::vector<double> &dalt,  std::vector<int> &found);

				int surfpt_( std::vector<double> &positn,  std::vector<double> &u,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &point,  std::vector<int> &found);

				int surfnm_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &point,  std::vector<double> &normal);

				int ednmpt_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &normal,  std::vector<double> &point);

				int edpnt_( std::vector<double> &p,  std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &ep);

				int edlimb_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &viewpt,  std::vector<double> &limb);

				int inelpl_( std::vector<double> &ellips,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt1,  std::vector<double> &xpt2);

				int npedln_( std::vector<double> &a,  std::vector<double> &b,  std::vector<double> &c__,  std::vector<double> &linept,  std::vector<double> &linedr,  std::vector<double> &pnear,  std::vector<double> &dist);

				int inrypl_( std::vector<double> &vertex,  std::vector<double> &dir,  std::vector<double> &plane,  std::vector<int> &nxpts,  std::vector<double> &xpt);

				int pjelpl_( std::vector<double> &elin,  std::vector<double> &plane,  std::vector<double> &elout);

				int saelgv_( std::vector<double> &vec1,  std::vector<double> &vec2,  std::vector<double> &smajor,  std::vector<double> &sminor);

				int nplnpt_( std::vector<double> &linpt,  std::vector<double> &lindir,  std::vector<double> &point,  std::vector<double> &pnear,  std::vector<double> &dist);

				int npsgpt_( std::vector<double> &ep1,  std::vector<double> &ep2,  std::vector<double> &point,  std::vector<double> &pnear,  std::vector<double> &dist);

				int sincpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::string &dref,  std::vector<double> &dvec,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<int> &found, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len, int dref_len);

				int dskxv_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<int> &nrays,  std::vector<double> &vtxarr,  std::vector<double> &dirarr,  std::vector<double> &xptarr,  std::vector<int> &fndarr, int target_len, int fixref_len);

				int dskxsi_( std::vector<int> &pri,  std::string &target,  std::vector<int> &nsurf,  std::vector<int> &srflst,  std::vector<double> &et,  std::string &fixref,  std::vector<double> &vertex,  std::vector<double> &raydir,  std::vector<int> &maxd,  std::vector<int> &maxi,  std::vector<double> &xpt,  std::vector<int> &handle,  std::vector<int> &dladsc,  std::vector<double> &dskdsc,  std::vector<double> &dc,  std::vector<int> &ic,  std::vector<int> &found, int target_len, int fixref_len);

				int ilumin_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &solar,  std::vector<double> &emissn, int method_len, int target_len, int fixref_len, int abcorr_len, int obsrvr_len);

				int illumg_( std::string &method,  std::string &target,  std::string &ilusrc,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &incdnc,  std::vector<double> &emissn, int method_len, int target_len, int ilusrc_len, int fixref_len, int abcorr_len, int obsrvr_len);

				int illumf_( std::string &method,  std::string &target,  std::string &ilusrc,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &spoint,  std::vector<double> &trgepc,  std::vector<double> &srfvec,  std::vector<double> &phase,  std::vector<double> &incdnc,  std::vector<double> &emissn,  std::vector<int> &visibl,  std::vector<int> &lit, int method_len, int target_len, int ilusrc_len, int fixref_len, int abcorr_len, int obsrvr_len);

				int limbpt_( std::string &method,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::string &obsrvr,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &tangts, int method_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len);

				int termpt_( std::string &method,  std::string &ilusrc,  std::string &target,  std::vector<double> &et,  std::string &fixref,  std::string &abcorr,  std::string &corloc,  std::string &obsrvr,  std::vector<double> &refvec,  std::vector<double> &rolstp,  std::vector<int> &ncuts,  std::vector<double> &schstp,  std::vector<double> &soltol,  std::vector<int> &maxn,  std::vector<int> &npts,  std::vector<double> &points,  std::vector<double> &epochs,  std::vector<double> &trmvcs, int method_len, int ilusrc_len, int target_len, int fixref_len, int abcorr_len, int corloc_len, int obsrvr_len);

				int conics_( std::vector<double> &elts,  std::vector<double> &et,  std::vector<double> &state);

				int oscelt_( std::vector<double> &state,  std::vector<double> &et,  std::vector<double> &mu,  std::vector<double> &elts);

				int occult_( std::string &targ1,  std::string &shape1,  std::string &frame1,  std::string &targ2,  std::string &shape2,  std::string &frame2,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &et,  std::vector<int> &ocltid, int targ1_len, int shape1_len, int frame1_len, int targ2_len, int shape2_len, int frame2_len, int abcorr_len, int obsrvr_len);

				int gfoclt_( std::string &occtyp,  std::string &front,  std::string &fshape,  std::string &fframe,  std::string &back,  std::string &bshape,  std::string &bframe,  std::string &abcorr,  std::string &obsrvr,  std::vector<double> &step,  std::vector<double> &cnfine,  std::vector<double> &result, int occtyp_len, int front_len, int fshape_len, int fframe_len, int back_len, int bshape_len, int bframe_len, int abcorr_len, int obsrvr_len);

				int srfcss_( std::vector<int> &code,  std::string &bodstr,  std::string &srfstr,  std::vector<int> &isname, int bodstr_len, int srfstr_len);

				int srfs2c_( std::string &srfstr,  std::string &bodstr,  std::vector<int> &code,  std::vector<int> &found, int srfstr_len, int bodstr_len);

				int srfc2s_( std::vector<int> &code,  std::vector<int> &bodyid,  std::string &srfstr,  std::vector<int> &isname, int srfstr_len);

				int srfscc_( std::string &srfstr,  std::vector<int> &bodyid,  std::vector<int> &code,  std::vector<int> &found, int srfstr_len);

				int bodc2n_( std::vector<int> &code,  std::string &name__,  std::vector<int> &found, int name_len);

				int srfrec_( std::vector<int> &body,  std::vector<double> &long__,  std::vector<double> &lat,  std::vector<double> &rectan);

			protected:
				std:vector<double> positn;

				std:vector<double> mu;
		};
	}
}

#endif // SMARTASTRO_CELESTIAL_OBJECT_H
