/*@Start***********************************************************
 This file is the GE interpreter for the TOPPE development environment 
 for platform-independent MR pulse programming.

 (c) 2016 The Regents of the University of Michigan
 Jon-Fredrik Nielsen, jfnielse@umich.edu
 
 $Id: toppev2b.e,v 1.2 2018/10/26 16:02:06 jfnielse Exp $
 $SOurce: $
*@End*************************************************************/

@inline epic.h
@inline intwave1.h

@global
#include <math.h>
#ifndef HW_IO
#include <stdio.h>
#include <stdlib.h>
#else /* HW_IO */
#include <stdioLib.h>
#endif /* !HW_IO */
#include <string.h>

#include "stddef_ep.h"
#include "epicconf.h"
#include "pulsegen.h"
#include "em_psd_ermes.in"
#include "epic_error.h"
#include "support_func.h"
#include "filter.h"
#include "epicfuns.h"
#include "epic_loadcvs.h"
#include "InitAdvisories.h"
#include "psdiopt.h"
#include "epic_iopt_util.h"

#include "grad_rf_sprlio.globals.h"

#include "jfn_rf_files2.h"   /* functions for loading external waveforms from file JFN Jan 13 */
#include "cores.h"   /* functions for reading .coredef file */
#include "jfn_globaldefs.h"   

#define TSP             2.0us   /* slow rcvr sampling 500kHz */
#define IRINT(f) ((((f)-(int)(f)) < 0.5) ? ((int) (f)) : (1+(int)(f)))
#define MAX(a,b) ( (a > b) ? a : b )
#define MIN(a,b) ( (a < b) ? a : b )
#define ABS(a) ( (a < 0) ? -a : a )

#define MAXCORES 20 /* max number of different cores in modules.txt. This number was determined empirically; may change in the future. */

#define GAMMA_H1 26754          /* in (rad/s)/Gauss */

#define GAMMA 4258

/* begin change for DV26.0 -- hack */
int debug = 0;
/* begin change for DV26.0 -- hack */

@inline Prescan.e PSglobal
int debugstate = 1;

@cv
@inline loadrheader.e rheadercv
@inline Prescan.e PScvs

int nextra = 50 with {0,,,,"number of disdaqs",};
int nframes = 1 with {1,,,,"number of time frames",};
int total_views;	/*  total shots to gather per slice */
int gating = TRIG_INTERN with {,,,,"1=line,3=ecg,5=aux,7=intern,9=nline",};
int psdseqtime;     /* sequence repetition time */
int timessi=100us with {0,,400us,INVIS,"time from eos to ssi in intern trig",};
float gslew = 150.0 with {0.0,3000.0,1,VIS, "readout gradient max slew, mT/m/ms",};
float rtimescale = 2.1 with {0.1,10.0,1,VIS, "grad risetime scale relative to cfsrmode",};
int daqdel = 128us; /*  gradient delay vs readout */
int minte;
int endtime = 500ms with {0,,,,"time at end of seq in rt mode",};

int ncores;
int nstartseq; 
int maxslice;
int maxecho;
int maxview;

int slord;

int filter_echo1;
float tsp = 4.0us with {1.0,32,,, "A/D sampling interval",};
float bandwidth = 125.0 with {2.0,250.0,1,VIS, "CERD rec low pass freq, kHz",};
float decimation = 1.0;

int filter_aps2 = 1;
int psfrsize = 512 with {0,,,INVIS, "num data acq during APS2",};
int queue_size = 1024 with {0,,,INVIS, "Rsp queue size",};

int obl_debug = 0; /* for obloptimize, I guess */
int obl_method = 1 with {0,1,0,INVIS,
"On(=1) to optimize the targets based on actual rotation matrices",};

int 	xres=64;            /* used in place of opxres to fix a simulation error */

float target;
int   rtime, ftime;

/* some variables to adjust timing - leftovers from signa */
float 	slwid180 = 1.2 with {0.0,4.0,,,"180 slice width as fraction of 90",}; /*this is so that the 180 stuff will work */
int	tpre=0;
int	tdel=0;
int	mytpre = 0;
int 	tcs;  /* duration of the chem sat pulses */
int	cyc_rf1 = 4; /* this seems like it was used by the stanford code, but was missing from declarations */
int 	echoshift = 0us with {,,,,"180 advance for spin echo T2* wting",};

/* tip-down pulse: loaded from tipdown.mod; used fill in rfpulse struct */
float flip_rfd;  
int pw_rfd, res_rfd, ia_rfd;
int wg_rfd = TYPRHO1 with {0, WF_MAX_PROCESSORS*2-1,TYPRHO1, VIS, , };
float a_rfd;
float thk_rfd = 5;


int npre = 0 ;   /* start of acquisition in readout.wav */
int ndaq = 0 ;   /* rhfrsize (if using external readout .wav file) */

int rftype    = 0 with {0,2,,, "Select RF pulse used in imaging pulse train. (0) external RF pulse (in rffiles.txt), (1) min-phase pulse (rf3d.rho), (2) tagging pulse.", };
float gsliceselect = 4.0;   /* [G/cm] Needed for slice offset calculation. */
int ishard    = 0 with {0,,,, "Play out hard pulse?", };
/* int yres    = 64 with {1,512,,, "# pixels along y",}; */

/* CVs: Velocity Encoding  */
float vres;

/* pulse sequence timing */
int myrfdel    = 94us ;   /* measured by JFN Nov-21-2012, from looking at rotation of rectangular tip-up excitation */
int start_core = 224us;    /* earliest start time for pulses in a core. If gradients only, can set as low as 24us. */
int timetrwait = 64us;    /* time in SEQLEN must be this much longer than end of RF waveforms, determined empirically JFN */

/* B1 scaling */
float xmtaddScan;




@ipgexport
@inline Prescan.e PSipgexport
long savrot[TRIG_ROT_MAX][9];   /* copy of rotation matrices */
RF_PULSE_INFO rfpulseInfo[RF_FREE];

int looparr[9000000];    /* scan loop definition, loaded from scanloop.txt */
int hasRF[MAXCORES];     /* I couldn't get dynamically allocated array (on host side) to work in scan() */
int hasDAQ[MAXCORES];
int coredur[MAXCORES];

float nom_fa;          /* nominal (design) flip angle in external RF (.wav) file */
char filenamerfd[100]; 


@host
#include "stdio.h"

/* changes for DV26 */
/* #include "sar_pm.h" */
#include "psdopt.h"
#include "sar_burst_api.h"
#include "sar_display_api.h"
#include "sar_limit_api.h"
#include "sar_pm.h"
/* end changes for DV26 */

#include "grad_rf_sprlio.h"
#include "fudgetargets.c"

static char supfailfmt[] = "Support routine %s exploded!";
FILTER_INFO echo1_filt;
FILTER_INFO aps2_filt;

@inline Prescan.e PShostVars

/* begin change for DV26 */
@inline loadrheader.e rheaderhost
/* end change for DV26 */

abstract("Pulseq interpreter module for GE");
psdname("toppe");  /* The End Of Pulse Programming */

int cvinit()
{
	int temp;

	EpicConf();
	inittargets(&loggrd, &phygrd);
	fudgetargets(&loggrd, &phygrd, rtimescale);
	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
					exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE)
		return FAILURE;
	
	/* rcvn_flag = 2;  */  /* commented out for DV26 */  /* Necessary change around DV23 or so. rcvn_flag is declared in Prescan.e. JFN Dec 2015 */
@inline Prescan.e PScvinit
#include "cvinit.in"	/* Runs the code generated by macros in preproc.*/
	
	pitrnub = 0;
   temp = _optr.fixedflag;
	optr = 6000ms;
	_optr.fixedflag = temp;

	pixresnub = 0;
	
	/* restrict sequence type to gradient echo */
	cvdef(oppseq,2);
	cvmin(oppseq,2);
	cvmax(oppseq,2);
	
	/* don't make the user even worry about selecting the number of echoes
	 or NEX  */
	piechnub = 0;
	pinexnub = 0;
	
	/* turn off variable bandwidth button */
	pircbnub = 0;
	pircb2nub = 0;
	pircbval2 = 15.525;
	pircbval3 = 31.25;
	pircbval3 = 62.5;
	pircbval3 = 125.0;
	cvmin(oprbw, 125);
	cvmax(oprbw, 125);
	cvdef(oprbw, 125);
	oprbw = 125;
	
	return SUCCESS;
	
}

int cveval()
{
	int tmptr, entry, temp;
	rfstruct rfdinfo;   /* contains rf header info, and the waveforms themselves */
	/* int dephasetime; */
	float dt = 4e-3;       /* msec */
	float gamma = 4.2575;  /* kHz/Gauss */
	char tmpfname[100];
	
	fprintf(stderr,"\nCveval stuff:");
	if (_psd_rf_wait.fixedflag == 0)  { /* sets psd_grd_wait and psd_rf_wait */
		if (setsysparms() == FAILURE)  {
			epic_error(use_ermes,"Support routine setsysparams failed",
					   EM_PSD_SUPPORT_FAILURE,1, STRING_ARG,"setsysparms");
			return FAILURE;
		}
	}
	
	if (obloptimize(&loggrd, &phygrd, scan_info, exist(opslquant), exist(opplane),
					exist(opcoax), obl_method, obl_debug, &opnewgeo, cfsrmode)==FAILURE)
		return FAILURE;
	
	if (existcv(opcgate) && (opcgate == PSD_ON)) {
		tmptr = RUP_GRD((int)((float)(exist(ophrep))
							  *(60.0/exist(ophrate))* 1e6));
		pitrnub = 0;
	}
	else {
		tmptr = optr;
		pitrnub = 0;
		if (oppseq == 2) { /* gradient echo */
			pitrval2 = 80ms;
			pitrval3 = 250ms;
			pitrval4 = 500ms;
			pitrval5 = 1s;
			pitrval6 = 2s;
		}
		else {
			pitrval2 = 300ms;
			pitrval3 = 500ms;
			pitrval4 = 1s;
			pitrval5 = 1.5s;
			pitrval6 = 2s;
		}
	}
	pisctim1 = nframes*tmptr;
	pisctim2 = pinexval2*pisctim1;
	pisctim3 = pinexval3*pisctim1;
	pisctim4 = pinexval4*pisctim1;
	pisctim5 = pinexval5*pisctim1;
	pisctim6 = pinexval6*pisctim1;
	
	opflip = 90.0;
	cvdef(opflip, 90.0);
	if (oppseq == 2)
	{
		pifanub = 6;
		pifaval2 = 10;
		pifaval3 = 30;
		pifaval4 = 40;
		pifaval5 = 60;
		pifaval6 = 90;
	}
	else /* spin echo */
	{
		pifanub = 0; /* turn off flip angle buttons */
	}
	
	/* label TE buttons */
	pite1nub = 63; /* apparently a bit mask */
	pite1val2 = PSD_MINIMUMTE;
	pite1val3 = 20ms;
	pite1val4 = 30ms;
	pite1val5 = 40ms;
	pite1val6 = 50ms;
	
	/* label FOV buttons */
	cvdef(opfov, 200);
	opfov = 200;
	pifovval2 = 200;
	pifovval3 = 220;
	pifovval4 = 240;
	pifovval5 = 360;
	pifovval6 = 480;
	
	slord = TYPNORMORDER;

	piamnub = 7;
	pixresnub = 0;
	cvmin(opxres, 16);
	cvdef(opxres, 64);
    temp = _opxres.fixedflag;
	opxres = 64;
	_opxres.fixedflag = temp;
	pixresval3 = 96;
	pixresval4 = 128;
	piyresnub = 0;

  /* TOPPE: turn off non-relevant pinubs */
  pite1nub = 0;
  temp = _opte.fixedflag;
  opte = 5ms;
  _opte.fixedflag = temp;

  pifanub = 0;
/*
  temp = _opflip.fixedflag;
  opflip = 10;
  _opflip.fixedflag = temp;
*/

  pifovnub = 0;
  temp = _opfov.fixedflag;
  opfov = 240;
  _opfov.fixedflag = temp;

  pislquant = 0;
  temp = _opslquant.fixedflag;
  opslquant = 1;
  _opslquant.fixedflag = temp;


  /* CVEVAL: fill in rfpulse struct values based on tipdown.mod */
  fprintf(stderr,"cveval(): reading rf file headers\n");
  jfn_rf_readheader("tipdown.mod", &rfdinfo);
  pw_rfd  = 4*rfdinfo.rfres;
  res_rfd = rfdinfo.rfres;
  fprintf(stderr,"cveval(): rfdinfo.paramsint[0],[1] = %d, %d\n", rfdinfo.paramsint16[0], rfdinfo.paramsint16[1]);
  rfpulse[RFD_SLOT].abswidth      = rfdinfo.paramsfloat[1];
  rfpulse[RFD_SLOT].effwidth      = rfdinfo.paramsfloat[2];
  rfpulse[RFD_SLOT].area          = rfdinfo.paramsfloat[3];
  rfpulse[RFD_SLOT].dtycyc        = rfdinfo.paramsfloat[4];
  rfpulse[RFD_SLOT].maxpw         = rfdinfo.paramsfloat[5];
  rfpulse[RFD_SLOT].num           = rfdinfo.paramsfloat[6];
  rfpulse[RFD_SLOT].max_b1        = rfdinfo.paramsfloat[7];
  rfpulse[RFD_SLOT].max_int_b1_sq = rfdinfo.paramsfloat[8];
  rfpulse[RFD_SLOT].max_rms_b1    = rfdinfo.paramsfloat[9];
  rfpulse[RFD_SLOT].nom_fa        = rfdinfo.paramsfloat[10];
  rfpulse[RFD_SLOT].nom_pw        = rfdinfo.paramsfloat[11];
  rfpulse[RFD_SLOT].nom_bw        = rfdinfo.paramsfloat[12];

  nom_fa = rfpulse[RFD_SLOT].nom_fa; 

  /* Set opflip to nom_fa of tipdown.mod */
  temp = _opflip.fixedflag;
  _opflip.fixedflag = 0;
  opflip = rfpulse[RFD_SLOT].nom_fa;
  _opflip.fixedflag = temp;
 
  /* set tipdown flip angle to opflip */
  a_rfd = 1.0; /* opflip/nom_fa; */
  flip_rfd = opflip;

	/* RF Scaling: Scale Pulses to the peak B1 in whole seq */
/*
	maxB1Seq = 0.0;
	for (entry=0; entry < MAX_ENTRY_POINTS; entry++) {
		if (peakB1(&maxB1[entry], entry, RF_FREE, rfpulse) == FAILURE) {
			epic_error(use_ermes,"peakB1 failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"peakB1");
			return FAILURE;
		}
		if (maxB1[entry] > maxB1Seq)
			maxB1Seq = maxB1[entry];
	}
*/

	if( findMaxB1Seq(&maxB1Seq, maxB1, MAX_ENTRY_POINTS, rfpulse, RF_FREE) == FAILURE ) {
		epic_error(use_ermes,supfailfmt,EM_PSD_SUPPORT_FAILURE,EE_ARGS(1),STRING_ARG,
		"findMaxB1Seq");
		return FAILURE;
	}

	xmtaddScan = 0;
	extraScale = 1.0;

	if (setScale(L_SCAN, RF_FREE, rfpulse, maxB1[L_SCAN], extraScale) == FAILURE) { 
		epic_error(0, "setScale() failed in predownload.",0,0);
		return FAILURE; 
	}

	ia_rfd    = (int)(max_pg_iamp*(*rfpulse[RFD_SLOT].amp));

	entry_point_table[L_SCAN].epxmtadd = (short) rint((double)xmtaddScan);
	
	// LHG 6/29/12:   this may not be necessary
	rhtype1 = 0;
	
	/* get the specs of the z gradient */
        gettarget(&target, ZGRAD, &loggrd);
        getramptime(&rtime, &ftime, ZGRAD, &loggrd);


	/* predownload(): new z and y phase-encode calculations JFN */
	gettarget(&target, YGRAD, &loggrd);
	getramptime(&rtime, &ftime, YGRAD, &loggrd);

@inline Prescan.e PScveval

	fprintf(stderr,"\nEnded cveval");
	
	return SUCCESS;
}

/* begin change for DV26 */
void getAPxParam(optval *min,optval *max,optdelta *delta,optfix *fix,float coverage,int algorithm)
{
/* Need to be filled when APx is supported in this PSD */
}

int getAPxAlgorithm(optparam *optflag, int *algorithm)
{
return APX_CORE_NONE;
}
/* end change for DV26 */

int cvcheck()
{
	return SUCCESS;
}

int predownload()
{
	int pdi, pdj;
	int i,j,jj,isl;
	float max_rbw;
	FILE *fpin;
	rfstruct roinfo;            /* contains contents of readout.mod file */
	corestruct coredefinfo;  /* contain list of .mod files for array of cores and related info. See cores.h */
	int loophdr[5];
	int scandur, tmpmin, tmpsec;

  /* velocity-encoding local variables */
/*
  float maxTa, risecoeff;
  float veTa, veT, tempfloat;
  float b, c, d;
  float specified_1st_moment, actual_1st_moment;
  float vfov;
  float maxvfov, minvfov, deltakvmin, deltakvmax;
  float tempampl;
  int k;
*/

	fprintf(stderr,"\nPredownload stuff:");

	opnecho = 1;
	cvdef(opnecho, 1);
	cvmin(opnecho, 1);
	cvmax(opnecho, 1);
	
	/* image header variables set for correct annotation */
	ihflip = opflip;
	ihnex = opnex;
	ihtr = optr;
	
	/*  make the end dead time = TR */
	endtime = optr;
	
	tmin = RUP_GRD(4ms);

	psdseqtime =  RUP_GRD(tmin) ;
	
	/* adjust grads, tsp, to reflect bandwidth  */
	calcvalidrbw((double)oprbw,&bandwidth,&max_rbw,&decimation,OVERWRITE_OPRBW,0);
	tsp = TSP*decimation;
	
	/* readout gradients - calculate 'em here  */
	
	/* get ndaq */
	if (jfn_rf_readheader("readout.mod", &roinfo) == JFN_FAILURE) {
		epic_error(use_ermes, "error in jfn_rf_readheader when reading readout.mod.",0,0);
		return FAILURE;
	}
	npre = roinfo.npre;
	ndaq = roinfo.rfres;

	/* core array stuff (TOPPE) */

	/* get scan loop header info */
	if ( cores_getloophdr("scanloop.txt", loophdr) == JFN_FAILURE ) {
		epic_error(use_ermes, "cores_getloophdr call failed",0,0);
		return FAILURE;
	}
	nstartseq = loophdr[0];
	maxslice  = loophdr[1];
	maxecho   = loophdr[2];
	maxview   = loophdr[3];
	scandur   = loophdr[4];     /* microseconds (int) */
	fprintf(stderr, "predownload(): nstartseq = %d, maxslice/maxecho/maxview = %d/%d/%d \n", nstartseq, maxslice, maxecho, maxview);
	tmpmin = scandur/(60*1000000);               /* C integer division acts as floor */
	tmpsec = scandur/1000000 - tmpmin*60;  
	fprintf(stderr, "predownload(): scan duration: %d min %d sec (scandur = %d)\n", tmpmin, tmpsec, scandur);

	/* get total number of unique cores, and fill hasDAQ and hasRF arrays */
	if (cores_getinfo("modules.txt", &coredefinfo)== JFN_FAILURE) {
		epic_error(use_ermes, "Error in cores_getinfo when reading modules.txt",0,0);
		return FAILURE;
	};
	ncores = coredefinfo.ncores;
	if (ncores  > MAXCORES ) {
		epic_error(use_ermes, "number of cores in modules.txt exceeds MAXCORES",0,0);
		return FAILURE;
	}
	for (j=0; j<ncores; j++) {
		hasRF[j] = coredefinfo.hasRF[j];
		hasDAQ[j] = coredefinfo.hasDAQ[j];
		fprintf(stderr, "hasRF[%d] = %d, hasDAQ[%d] = %d\n", j, hasRF[j], j, hasDAQ[j]);
		if (hasRF[j] && hasDAQ[j]) {
			epic_error(use_ermes, "Core can only use either ADC or RF excitation.",0,0);
			return FAILURE;
		}
	}

	/* fill looparr (needed in scan()) */
	if ( cores_readloop("scanloop.txt", looparr) == JFN_FAILURE ) {
		epic_error(use_ermes, "cores readloop call failed",0,0);
		return FAILURE;
	}

	/* if(gram_duty() == FAILURE) return FAILURE; */

#include "predownload.in"
	
	/* set up clock */
	if ((exist(opcgate) == PSD_ON) && existcv(opcgate))
	{
		pidmode = PSD_CLOCK_CARDIAC;
		piviews = nextra+nframes;
		piclckcnt = ophrep;
		pitscan = (float)(optr)*(nframes);
		pitslice = psdseqtime;
	}
	else
	{
		pidmode = PSD_CLOCK_NORM;
		pitslice = psdseqtime;
		pitscan = (float)(optr)*(nframes);
	}

	pidmode = PSD_CLOCK_NORM;
	/* pitslice = (float) (scandur/1000000.0/(float)opslquant); */   /* don't know if this is needed */        
	pitscan = scandur;    /* us */         /* (float)(scandur)/1000000.0; */      /* sec */
	fprintf(stderr, "pitscan = %d\n", pitscan);
	
	minte = 4us;
	cvmin(opte, minte);
	cvdef(opte, minte);
	if ((exist(opautote) == PSD_MINTE)||(exist(opautote) == PSD_MINTEFULL))
		opte = minte;
	ihte1 = opte;
	
@inline loadrheader.e rheaderinit
	
	rhimsize = 64;
	while(rhimsize < opxres)  rhimsize *= 2;
	rhxoff = rhimsize*scan_info[0].oprloc/opfov;
	rhyoff = rhimsize*scan_info[0].opphasoff/opfov;

	rhbline = 0;
	rhnecho = maxecho+1;
	rhnslices = maxslice+1;
	rhrcctrl = 1;    		/* lx wants to create images.  */
	rhexecctrl = 2;  		/* just save the raw data */
	/* autolock = 1; */
	/* rhdacqctrl += 8192; */    /* CHANGE for DV26 */
	
	if(opuser19==1.0)           /* turn on/off RDS for realtime */
		rhtype1 |= (0x00080000 | 0x00100000);
	else
		rhtype1 &=  0x0007FFFF;
	
	/* acqs = 1; */
	slquant1 = opslquant;
	
	/* for straight sequential order of slices. */
	if (!orderslice(slord, opslquant, slquant1, gating))
		epic_error(use_ermes,"orderslice call failed",0,0);

	/* needed for conversion to DV23 */
	for (isl=0;isl<opslquant;isl++) {
        for (jj=0;jj<9;jj++) {
            rsprot[isl][jj] = hostToRspRotMat(scan_info[isl].oprot[jj]);
        }
    }
	
	/* initialize copy of original rotation matrices */
	for (pdi = 0; pdi < opslquant; pdi++)
		for (pdj = 0; pdj < 9; pdj++)
			savrot[pdi][pdj] = rsprot[pdi][pdj];
	scalerotmats(rsprot, &loggrd, &phygrd, opslquant, 0);
	
	/* save stuff for maxwell correction */
	rhmaxcoef1a = rsprot[0][0]/(float)cfxfull;	/* save x rotator */
	rhmaxcoef1b = rsprot[0][1]/(float)cfxfull;
	rhmaxcoef2a = rsprot[0][3]/(float)cfyfull; /* y  */
	rhmaxcoef2b = rsprot[0][4]/(float)cfyfull;
	rhmaxcoef3a = rsprot[0][6]/(float)cfzfull; /* z  */
	rhmaxcoef3b = rsprot[0][7]/(float)cfzfull;
	
	rhdab0s = cfrecvst;
	rhdab0e = cfrecvend;
	
	if(entrytabinit(entry_point_table, (int)ENTRY_POINT_MAX)
	   == FAILURE) {
		epic_error(use_ermes,"Can't initialize entry point table.",0,0);
		return FAILURE;
	}
	
	/* set up receiver */
	
	initfilter();
	
	cvmax(rhfrsize, 32768);		/* for now  */
	rhfrsize = ndaq*4us/tsp;      /* num points sampled */
	rhfrsize = 4*(rhfrsize/4);          /* wants to be divisible by 4 */
	rhnframes =2*((maxview+1)/2);  /* has to be an even number */
	rhuser4 = rhnframes;   

	/* the following is for backward compatibility with 'host' .e file */
	/* yres = 64; */
	
	if (calcfilter( &echo1_filt,bandwidth,rhfrsize,OVERWRITE_OPRBW ) == FAILURE) {
		epic_error(use_ermes,"%s failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"calcfilter");
		return FAILURE;
	}
	setfilter( &echo1_filt, SCAN );
	filter_echo1 = echo1_filt.fslot;
	
	if (calcfilter( &aps2_filt,bandwidth,psfrsize,OVERWRITE_OPRBW ) == FAILURE) {
		epic_error(use_ermes,"%s failed",EM_PSD_SUPPORT_FAILURE,1,STRING_ARG,"calcfilter");
		return FAILURE;
	}
	setfilter( &aps2_filt, PRESCAN );
	filter_aps2 = aps2_filt.fslot;
	
	rhrawsize = 2*rhptsize*rhfrsize*(rhnframes+1)*rhnslices*rhnecho;
	numrecv = rhdab0e-rhdab0s+1;
	if((float)rhrawsize*numrecv > cftpssize-2e7)  {      /* reserve 20 MB */
		epic_error(use_ermes,"Oops! Tooo much memory requested.\n",0,0);
		return FAILURE;
	}
	
	daqdel = psd_grd_wait + 0.5;
	
	strcpy(entry_point_table[L_SCAN].epname, "scan");
	entry_point_table[L_SCAN].epfastrec = 0;
	entry_point_table[L_SCAN].epstartrec = rhdab0s;
	entry_point_table[L_SCAN].ependrec = rhdab0e;
	entry_point_table[L_SCAN].epfilter = (unsigned char)echo1_filt.fslot;
	entry_point_table[L_SCAN].epprexres = rhfrsize;
	entry_point_table[L_SCAN].epxmtadd = txCoilInfo[getTxIndex(coilInfo[0])].coilAtten;
	entry_point_table[L_APS2] =
	entry_point_table[L_MPS2] =
	entry_point_table[L_SCAN];      /* copy scan into APS2 & MPS2 */
	strcpy(entry_point_table[L_APS2].epname,"aps2");
	strcpy(entry_point_table[L_MPS2].epname,"mps2");
	entry_point_table[L_APS2].epfilter = (unsigned char)aps2_filt.fslot;
	entry_point_table[L_MPS2].epfilter = (unsigned char)aps2_filt.fslot;
	entry_point_table[L_APS2].epprexres = psfrsize;
	entry_point_table[L_MPS2].epprexres = psfrsize;
	
	/*  Some prescan stuff  */
	
	pislquant = opslquant;
	
@inline Prescan.e PSfilter
@inline Prescan.e PSpredownload
	
	phys_record_flag = opuser18;        /* put here so it isn't overwritten in predownload */
	phys_record_channelsel = 14;        /* PG wave,trig & RESP */

	eepf = 0;
	oepf = 0;
	eeff = 1;
	oeff = 1; 
	set_echo_flip(&rhdacqctrl, &chksum_rhdacqctrl, eepf, oepf, eeff, oeff);

	return SUCCESS;
} /* End-Of-Predownload */

@inline Prescan.e PShost



@rsp

int pre = 0; /* prescan flag */
short thamp;

CHAR *entry_name_list[ENTRY_POINT_MAX] = { "scan", "mps2", "aps2",
@inline Prescan.e PSeplist
};

int *fxmit, *frec;

/* LHG 9/26/12/ Tables for kz encoding (grad amplitudes and phases) */ 

@rspvar

int iv, ifr, isl, it, ikv;
int tmp;
int i, k;
int trig, dtype;

short trigonpkt[3] = {0, SSPOC+DREG, SSPD+DSSPD4};
short trigoffpkt[3] = {0, SSPOC+DREG, SSPD};
short trigonwd, trigoffwd;

/* variables needed for prescan */
short chopamp;
int view, slice, dabop, excitation, seqCount, ec, rf1Phase, seqCount;
int rspent, rspdda, rspbas, rspvus, rspgy1, rspasl;
int rspesl, rspchp, rspnex, rspslq, rspsct;
int dabmask;

@inline Prescan.e PSrspvar	/* For Prescan */
extern PSD_EXIT_ARG psdexitarg;



@pg

#include <epic_loadcvs.h>
int dur_tipdowncore, dur_seqcore;
int dur_rfcore, dur_spoilcore, dur_astscore;

long 	deadtime_astcore;
long	deadtime_controlcore;
long	deadtime_tdelaycore;
long	deadtime_tadjustcore;
long	deadtime_preBScore;


WF_PULSE rfd, thetad, gxrfd, gyrfd, gzrfd; /* tip-down pulses */

/* stuff for array of cores (TOPPE) */
WF_PULSE initpulse = INITPULSE;
WF_PULSE* rhocores;
WF_PULSE* thetacores;
WF_PULSE* gxcores;
WF_PULSE* gycores;
WF_PULSE* gzcores;
WF_PULSE* waitcores;
WF_PULSE* cores;
SEQUENCE_ENTRIES* off_cores;
#ifndef IPG
  int* idx_cores;
#endif 
WF_PULSE* echocores;

/* pulses and pointers needed for real-time waveform switching */
WF_PULSE** gxwav;
WF_PULSE** gywav;
WF_PULSE** gzwav;
WF_PULSE** rhowav;
WF_PULSE** thetawav;
WF_HW_WAVEFORM_PTR** gxwavp;
WF_HW_WAVEFORM_PTR** gywavp;
WF_HW_WAVEFORM_PTR** gzwavp;
WF_HW_WAVEFORM_PTR** rhowavp;
WF_HW_WAVEFORM_PTR** thetawavp;

WF_PULSE echo1 = INITPULSE; 

/* Create x/y/z/rho/theta pulse from 'rfstruct' structure (defined in jfn_rf_files2.h). */
WF_PULSE jfn_rf_makepulse(WF_PROCESSOR wfp, char *pname, rfstruct *rfinfo, int tbeg)
{
   WF_PULSE *ep;
   WF_PULSE proto = INITPULSE;
   short* wave;
   int res;

   switch (wfp) {
      case TYPXGRAD:
      case TYPYGRAD:
      case TYPZGRAD:
         res = rfinfo->res;
         break;
      case TYPRHO1:
      case TYPTHETA:
         res = rfinfo->rfres;
         break;
   }

   ep = (WF_PULSE *) AllocNode(sizeof(WF_PULSE));
   memcpy((char*) ep, (char*) &proto, sizeof(WF_PULSE));
   pulsename(ep, pname);
   createreserve(ep, wfp, res);

   switch (wfp) {
      case TYPXGRAD:
         wave = rfinfo->gx[0];
         break;
      case TYPYGRAD:
         wave = rfinfo->gy[0];
         break;
      case TYPZGRAD:
         wave = rfinfo->gz[0];
         break;
      case TYPRHO1:
         wave = &(rfinfo->rho[0][0][rfinfo->npre]);
         break;
      case TYPTHETA:
         wave = &(rfinfo->theta[0][0][rfinfo->npre]);
         break;
      default:
         break;
   }

   /* make sure EOS bit of last point is set */
   if (!wave[res-1] % 2)
      wave[res-1]++;

   movewaveimm(wave, ep, 0, res, TOHARDWARE);
   createinstr(ep, tbeg, 4*res, MAX_PG_IAMP);
   if (wfp == TYPRHO1) {
       /* addrfbits(ep, 0, tbeg, 4*res); */
       fastAddrfbits(ep,0,tbeg,4*res, 50us);
   }

   return(*ep);
}

/* create wait pulse (for adjusting TR) */
WF_PULSE jfn_wait(WF_PROCESSOR wfp, char *pname, int tbeg)
{
   WF_PULSE *ep;
   WF_PULSE proto = INITPULSE;
   int res = 1;

   ep = (WF_PULSE *) AllocNode(sizeof(WF_PULSE));
   memcpy((char*) ep, (char*) &proto, sizeof(WF_PULSE));
   pulsename(ep, pname);
   createreserve(ep, wfp, res);

   createconst(ep, wfp, 4*res, 0); 
   createinstr(ep, tbeg, 4*res, 0);
   
   return(*ep);
}


STATUS pulsegen(void)
{
	int j,i; 
	char tstr[40];
	rfstruct rfdinfo;   /* contains rf header info, and the waveforms themselves  JFN Jan 2013 */
	rfstruct* coresinfo;  /* array of rfstruct (.mod file info) for array of cores */
	corestruct coredefinfo;  /* contain list of .mod files for array of cores and related info. See cores.h */
	char name[16];
	char rhoname[16];
	char thetaname[16];
	char gxname[16];
	char gyname[16];
	char gzname[16];
	char waitname[16];
	int  mindur;

	sspinit(psd_board_type);

	/****************************************************************/
	/* create an array of cores                                     */
	/****************************************************************/
	cores_getinfo("modules.txt", &coredefinfo);
	ncores = coredefinfo.ncores;
	fprintf(stderr, "\npulsegen(): ncores = %d \n", ncores);
	fprintf(stderr, "pulsegen(): nstartseq = %d, maxslice/maxecho/maxview = %d/%d/%d \n", nstartseq, maxslice, maxecho, maxview);

	cores     = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	off_cores = (SEQUENCE_ENTRIES*) AllocNode(ncores*sizeof(SEQUENCE_ENTRIES));
#ifndef IPG
	idx_cores = (int*) AllocNode(ncores*sizeof(int));
#endif 
	gxcores   = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	gycores   = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	gzcores   = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	rhocores  = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	thetacores  = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));
	waitcores  = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));

	coresinfo = (rfstruct *) AllocNode(ncores*sizeof(rfstruct));
	echocores = (WF_PULSE *) AllocNode(ncores*sizeof(WF_PULSE));

	gxwav  = (WF_PULSE **) AllocNode(ncores*sizeof(WF_PULSE *));
	gywav  = (WF_PULSE **) AllocNode(ncores*sizeof(WF_PULSE *));
	gzwav  = (WF_PULSE **) AllocNode(ncores*sizeof(WF_PULSE *));
	rhowav    = (WF_PULSE **) AllocNode(ncores*sizeof(WF_PULSE *));
	thetawav  = (WF_PULSE **) AllocNode(ncores*sizeof(WF_PULSE *));
	gxwavp  = (WF_HW_WAVEFORM_PTR **) AllocNode(ncores*sizeof(WF_HW_WAVEFORM_PTR *));
	gywavp  = (WF_HW_WAVEFORM_PTR **) AllocNode(ncores*sizeof(WF_HW_WAVEFORM_PTR *));
	gzwavp  = (WF_HW_WAVEFORM_PTR **) AllocNode(ncores*sizeof(WF_HW_WAVEFORM_PTR *));
	rhowavp = (WF_HW_WAVEFORM_PTR **) AllocNode(ncores*sizeof(WF_HW_WAVEFORM_PTR *));
	thetawavp  = (WF_HW_WAVEFORM_PTR **) AllocNode(ncores*sizeof(WF_HW_WAVEFORM_PTR *));

	for ( j=0; j<ncores; j++ ) {

		/* read and load .mod file for core j */
		jfn_rf_readheader(coredefinfo.wavfiles[j], &coresinfo[j]);
		jfn_rf_allocatemem(&coresinfo[j]);     
		jfn_rf_readwaveforms(&coresinfo[j], ishard);

		/* Gradient waveforms */
    	sprintf(gxname, "gxcore%d", j);
		gxcores[j] = jfn_rf_makepulse(TYPXGRAD, gxname, &coresinfo[j], RUP_GRD(start_core)); /* TODO */
    	sprintf(gyname, "gycore%d", j);
		gycores[j] = jfn_rf_makepulse(TYPYGRAD, gyname, &coresinfo[j], RUP_GRD(start_core));
    	sprintf(gzname, "gzcore%d", j);
		gzcores[j] = jfn_rf_makepulse(TYPZGRAD, gzname, &coresinfo[j], RUP_GRD(start_core));

		/* RF pulse (if present) */
		myrfdel = psd_rf_wait;  /* psd_rf_wait has been a bit off in the past, so check this. TODO */
		if (hasRF[j] == 1) {
    		sprintf(rhoname, "rhocore%d", j);
			rhocores[j] = jfn_rf_makepulse(TYPRHO1, rhoname, &coresinfo[j], RUP_GRD(start_core+myrfdel));
    		sprintf(thetaname, "thetacore%d", j);
			thetacores[j] = jfn_rf_makepulse(TYPTHETA, thetaname, &coresinfo[j], pbeg(&rhocores[j], rhoname, 0)); 
		}

		/* data acquisition pulse (if present) */
		if (hasDAQ[j] == 1) {
			memcpy(&echocores[j], &initpulse, sizeof(WF_PULSE)); 
			npre = coresinfo[j].npre;
			ndaq = coresinfo[j].rfres;
			sprintf(name, "echocores%d", j);
			pulsename(&echocores[j], name);
			acqq(&echocores[j], RUP_GRD(pbeg(&gxcores[j], gxname, 0)+npre*4us+daqdel), 
				(long)(DEFAULTPOS), (long)(DEFAULTPOS), (long)filter_echo1, (TYPDAB_PACKETS)DABNORM);
		}

		/* create pointer to all waveforms in this module (core) */
		gxwav[j] = (WF_PULSE *) AllocNode(coresinfo[j].npulses*sizeof(WF_PULSE));
		gywav[j] = (WF_PULSE *) AllocNode(coresinfo[j].npulses*sizeof(WF_PULSE));
		gzwav[j] = (WF_PULSE *) AllocNode(coresinfo[j].npulses*sizeof(WF_PULSE));
		gxwavp[j] = (WF_HW_WAVEFORM_PTR *) AllocNode(coresinfo[j].npulses*sizeof(WF_HW_WAVEFORM_PTR));
		gywavp[j] = (WF_HW_WAVEFORM_PTR *) AllocNode(coresinfo[j].npulses*sizeof(WF_HW_WAVEFORM_PTR));
		gzwavp[j] = (WF_HW_WAVEFORM_PTR *) AllocNode(coresinfo[j].npulses*sizeof(WF_HW_WAVEFORM_PTR));
		if (hasRF[j] == 1) {
			rhowav[j] = (WF_PULSE *) AllocNode(coresinfo[j].npulses*sizeof(WF_PULSE));
			thetawav[j]  = (WF_PULSE *) AllocNode(coresinfo[j].npulses*sizeof(WF_PULSE));
			rhowavp[j] = (WF_HW_WAVEFORM_PTR *) AllocNode(coresinfo[j].npulses*sizeof(WF_HW_WAVEFORM_PTR));
			thetawavp[j]  = (WF_HW_WAVEFORM_PTR *) AllocNode(coresinfo[j].npulses*sizeof(WF_HW_WAVEFORM_PTR));
		}
		for (i = 0; i < coresinfo[j].npulses; i++) {
			sprintf(tstr, "gxwav_%d_%d", j, i);
			pulsename(&gxwav[j][i], tstr);
			createreserve(&(gxwav[j][i]), TYPXGRAD, coresinfo[j].res);
			movewaveimm(coresinfo[j].gx[i], &(gxwav[j][i]), (int) 0, coresinfo[j].res, TOHARDWARE);
			gxwavp[j][i] = gxwav[j][i].wave_addr;

			sprintf(tstr, "gywav_%d_%d", j, i);
			pulsename(&gywav[j][i], tstr);
			createreserve(&(gywav[j][i]), TYPYGRAD, coresinfo[j].res);
			movewaveimm(coresinfo[j].gy[i], &(gywav[j][i]), (int) 0, coresinfo[j].res, TOHARDWARE);
			gywavp[j][i] = gywav[j][i].wave_addr;

			sprintf(tstr, "gzwav_%d_%d", j, i);
			pulsename(&gzwav[j][i], tstr);
			createreserve(&(gzwav[j][i]), TYPZGRAD, coresinfo[j].res);
			movewaveimm(coresinfo[j].gz[i], &(gzwav[j][i]), (int) 0, coresinfo[j].res, TOHARDWARE);
			gzwavp[j][i] = gzwav[j][i].wave_addr;

			if (hasRF[j] == 1) {
				sprintf(tstr, "rhowav_%d_%d", j, i);
				pulsename(&rhowav[j][i], tstr);
				createreserve(&(rhowav[j][i]), TYPRHO1, coresinfo[j].res);
				movewaveimm(coresinfo[j].rho[i][0], &(rhowav[j][i]), (int) 0, coresinfo[j].res, TOHARDWARE);
				rhowavp[j][i] = rhowav[j][i].wave_addr;

				sprintf(tstr, "thetawav_%d_%d", j, i);
				pulsename(&thetawav[j][i], tstr);
				createreserve(&(thetawav[j][i]), TYPTHETA, coresinfo[j].res);
				movewaveimm(coresinfo[j].theta[i][0], &(thetawav[j][i]), (int) 0, coresinfo[j].res, TOHARDWARE);
				thetawavp[j][i] = thetawav[j][i].wave_addr;
			}
		}

		/* initial load TODO: remove, no need*/
		setwave(gxwavp[j][0], &(gxcores[j]), 0);    
		setwave(gywavp[j][0], &(gycores[j]), 0);
		setwave(gzwavp[j][0], &(gzcores[j]), 0);
		if (hasRF[j] == 1) {
			setwave(rhowavp[j][0], &(rhocores[j]), 0);
			setwave(thetawavp[j][0], &(thetacores[j]), 0);
		}

		/* set core duration; add wait pulse for real-time TR adjustment */
		mindur = RUP_GRD(pend(&gxcores[j],gxname,0) + timetrwait);  /* temporary calculation */
		if (hasRF[j] == 1) {
			mindur += myrfdel;
		}
		if (hasDAQ[j] == 1) {
			mindur += daqdel;
		}

    	sprintf(waitname, "waitcore%d", j);
		waitcores[j] = jfn_wait(SSP, waitname, RUP_GRD(mindur+4));
		mindur += 12 + timessi;             /* final value */

		coredur[j] = (coredefinfo.dur[j] < mindur) ? mindur : coredefinfo.dur[j] ;

		/* build core */
		memcpy(&cores[j], &initpulse, sizeof(WF_PULSE));  /* like WF_PULSE echo1 = INITPULSE */
    	sprintf(name, "core%d", j);
    	pulsename(&cores[j], name);
		createseq(&cores[j], RUP_GRD(coredur[j]-timessi), off_cores[j]);

#ifndef IPG
		updateIndex(&idx_cores[j]);
#endif 
	}

	FreeNode(coresinfo);

	/*  wait_for_scan_to_start sequence */
	WAIT(SSP, waitStart, 24us, 1ms);
	SEQLENGTH(waitpass, 10ms, waitpass);
	
	/*  wait for scanend for real time */
	WAIT(SSP, waitEnd, 24us, 2us);
	SEQLENGTH(waitend, endtime, waitend);

	/*  pass packet sequence (pass).  */
	PASSPACK(endpass, 49ms);
	SEQLENGTH(pass, 50ms, pass);
	
	fprintf(stderr, "\n inline Prescan.e PSpulsegen." );
@inline Prescan.e PSpulsegen	/*  prescan sequences  */
	fprintf(stderr, "\t ...  done .\n" );
	
	fprintf(stderr, "buildinstr() ...\n" );
	buildinstr();              	/* load the sequencer memory */
	fprintf(stderr, "\t ...  done .\n" );
	return SUCCESS;
	
} /* end of pulsegen */

@inline Prescan.e PSipg
/* end of @pg */



@rsp

STATUS scancore(void);

@inline Prescan.e PScore

/*manual prescan */
int mps2() {
	pre = 2;
	scancore();
	rspexit();
	return SUCCESS;
}
/*auto prescan */
int aps2() {
	pre = 1;
	scancore();
	rspexit();
	return SUCCESS;
}
/*Actual scan*/
int scan()
{
	pre = 0;
	scancore();
	rspexit();
	return SUCCESS;
}


STATUS scancore()
{
	int counter;
#define RESERVED_IPG_MEMORY (0xbfff00)
	int i, j, k, n;
	int myxmitfreq, myrecfreq;
/*	int 	*pcasl_iphase; */
	int icore;
	int rotmatx[1][9];
	float phi, cphi, sphi;
	
	printf("\nEntering scancore   pre = %d\n", pre);

	setrfconfig(ENBL_RHO1 + ENBL_THETA);
	setssitime((LONG)timessi/GRAD_UPDATE_TIME);
	rspqueueinit(queue_size);
/*
	scopeon(&seqcore);
	syncon(&seqcore);
*/
	syncoff(&pass);
	setrotatearray((SHORT)opslquant, rsprot[0]);
	settriggerarray((SHORT)opslquant, rsptrigger);

	dabmask = PSD_LOAD_DAB_ALL;
	if(pre)  {
		fprintf(stderr,"\nSetting filters for prescan...");
		for ( j=0; j<ncores; j++ ) {
			if (hasDAQ[j]) {
				setrfltrs(filter_aps2, &echocores[j]);
			}
		}
	}
	else  {
		fprintf(stderr,"\nSetting filters...");
		for ( j=0; j<ncores; j++ ) {
			if (hasDAQ[j]) {
				setrfltrs(filter_echo1, &echocores[j]);
			}
		}
	}
	fprintf(stderr,"\tdone");
	
	/* fix the clock for aux trig */
	
	if (!pre && gating==TRIG_AUX)  {
		setscantimemanual();
		setscantimestop();
		setscantimeimm(pidmode, pitscan, piviews,
					   pitslice, opslicecnt);
	}

	/* Allocate memory for RF pulse param tables  (legacy code) */ 
	fxmit = 	(int *) AllocNode(opslquant*sizeof(int));
	frec  = 	(int *) AllocNode(opslquant*sizeof(int));

	fprintf(stderr,"\nSetting up PCASL pulse Phase table  ...");

	setupslices(fxmit, rsp_info, opslquant, gsliceselect, 1.0, opfov, TYPTRANSMIT);

	for (isl = 0; isl < opslquant; isl++)
		rsp_info[isl].rsprloc = 0;
	setupslices(frec, rsp_info, opslquant, 0.0, 1.0, 2.0, TYPREC);

	/* set transmit and receive frequencies */
	myxmitfreq = (int)((fxmit[opslquant/2] + fxmit[opslquant/2-1])/2);
	myrecfreq = (int)((frec[opslquant/2] + frec[opslquant/2-1])/2);
	fprintf(stderr,"\nSetting tx/rx freq for core array...");
	for ( j=0; j<ncores; j++ ) {
		if (hasRF[j]) {
			setfrequency(myxmitfreq, &rhocores[j], 0); 
		}
		if (hasDAQ[j]) {
			setfrequency(myrecfreq , &echocores[j], 0);  
		}
	}
	
	dabmask = PSD_LOAD_DAB_ALL;
	dabop = DABSTORE;

	settrigger(TRIG_INTERN,0);

	for (j = 0; j < nstartseq; j++) {	
		/* get core index */
		icore = looparr[j*NL+0] - 1;

		/* set gradient waveforms and waveform amplitudes*/
		setwave(gxwavp[icore][looparr[j*NL+15]-1], &(gxcores[icore]), 0);
		setwave(gywavp[icore][looparr[j*NL+15]-1], &(gycores[icore]), 0);
		setwave(gzwavp[icore][looparr[j*NL+15]-1], &(gzcores[icore]), 0);
		setiamp(looparr[j*NL+3], &gxcores[icore], 0);
		setiamp(looparr[j*NL+4], &gycores[icore], 0);
		setiamp(looparr[j*NL+5], &gzcores[icore], 0);
	
		/* set RF waveform, amplitude, phase, and frequency */
		if (hasRF[icore]) {
			setwave(rhowavp[icore][looparr[j*NL+15]-1], &(rhocores[icore]), 0);
			setwave(thetawavp[icore][looparr[j*NL+15]-1], &(thetacores[icore]), 0);
			setiamp(looparr[j*NL+1], &rhocores[icore], 0);
			setiamp(looparr[j*NL+2], &thetacores[icore], 0);
			setiphase(looparr[j*NL+11], &rhocores[icore], 0);
			setfrequency((int)(looparr[j*NL+14]/TARDIS_FREQ_RES), &rhocores[icore], 0); 
		}

		/* set up data acquisition */
		if (hasDAQ[icore]) {
			dtype = looparr[j*NL+9] ? DABON : DABOFF;
			loaddab(&echocores[icore], looparr[j*NL+6], looparr[j*NL+7], dabop, looparr[j*NL+8], dtype, dabmask);
			setiphase(looparr[j*NL+12], &echocores[icore], 0);
		}

		/* set (in-plane) rotation */
		phi = (float) (1.0*looparr[j*NL+10]/max_pg_iamp*M_PI);
		cphi = cos(phi); sphi = sin(phi);
		for (k = 0; k < 9; k += 3)
		{
			rotmatx[0][k] = IRINT(cphi*savrot[0][k]+sphi*savrot[0][k+1]);
			rotmatx[0][k+1] = IRINT(-sphi*savrot[0][k]+cphi*savrot[0][k+1]);
			rotmatx[0][k+2] = savrot[0][k+2];
		}
		scalerotmats(rotmatx, &loggrd, &phygrd, 1, 0);
		setrotate(rotmatx[0],0);

		/* set wait pulse duration (for extending TR) */
		setperiod(MAX(looparr[j*NL+13],12), &waitcores[icore], 0);

		/* play module */
		boffset(off_cores[icore]); 
		startseq(0, MAY_PAUSE);  
	}
	
	if(opuser19 == 1.0) {      /* wait around for grecon to finish */
		boffset(off_waitend);
		settrigger(TRIG_INTERN, 0);
		startseq(0, MAY_PAUSE);   /* fat lady sings */
	}
	
	/* tell 'em it's over */
	boffset(off_pass);
	setwamp(SSPD+DABPASS+DABSCAN, &endpass, 2);
	settrigger(TRIG_INTERN, 0);
	startseq(0, MAY_PAUSE);	/* fat lady sings */

	FreeNode(cores);
	FreeNode(off_cores);
	FreeNode(rhocores);
	FreeNode(gxcores);
	FreeNode(gycores);
	FreeNode(gzcores);
	FreeNode(waitcores);
	FreeNode(echocores);
#ifndef IPG
	FreeNode(idx_cores);
#endif 
	
	return SUCCESS;
}


@pg
/********************************************
 * dummylinks
 *
 * This routine just pulls in routines from
 * the archive files by making a dummy call.
 ********************************************/
void dummylinks()
{
	epic_loadcvs("thefile"); /* for downloading CVs */
}


