/*
  This file is part of the TOPPE development environment for platform-independent MR pulse programming.

  TOPPE is free software: you can redistribute it and/or modify
  it under the terms of the GNU Library General Public License as published by
  the Free Software Foundation version 2.0 of the License.

  TOPPE is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Library General Public License for more details.

  You should have received a copy of the GNU Library General Public License
  along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
 
  (c) 2016 The Regents of the University of Michigan
  Jon-Fredrik Nielsen, jfnielse@umich.edu

  $Id: jfn_rf_files2.h,v 1.1 2018/10/04 19:01:20 jfnielse Exp $
 */

#ifndef JFN_RF_FILES2_H
#define JFN_RF_FILES2_H

typedef struct {
	char    fname[80];          /* rf file name */
	short   ncoils;             /* number of coils/channels */
	short   npre;               /* number of points before start of rf waveform */
	short   rfres;              /* number of points in rf waveform */
	short   res;                /* number of points in gradient waveform (>= rfres) */
	short   npulses;            /* number of different waveforms in .wav file */
	float   b1max;              /* Gauss */
	float   gmax;               /* Gauss/cm (4.0 on um3t) */
	short   dataoffset;         /* total header size (including ASCII and binary parts) */
	short   nparamsint16;       /* # of int16 parameters */
	short   nparamsfloat;       /* # of float parameters */
	short   paramsint16[32];    /* int16 parameters */
	float   paramsfloat[32];    /* float parameters */
	short*** rho;
	short*** theta;
	short**  gx;
	short**  gy;
	short**  gz;
} rfstruct;

int jfn_rf_getfilename(char *fname, int index, const char* directory);
int jfn_rf_readheader(char *fname, rfstruct *rfinfo);
int jfn_rf_allocatemem(rfstruct *rfinfo);
int jfn_rf_readwaveforms(rfstruct *rfinfo, int ishard);
int jfn_rf_freemem(rfstruct *rfinfo);
int readshort(short* i, int n, FILE* fid);

#endif
