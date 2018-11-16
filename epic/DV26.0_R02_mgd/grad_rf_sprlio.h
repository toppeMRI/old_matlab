/* *************************************
 * grad_rf_sprlio.h		1/30/2008 	ghg
 * This structure is used to track the 
 * rf heating, SAR heating, Grad coil heating,
 * grad amplifier heating.
 * ********************************** */

/* only do this once in any given compilation.*/
#ifndef  grad_rf_sprlio__INCL
#define  grad_rf_sprlio_INCL

RF_PULSE rfpulse[MAX_RFPULSE] = {
  /* external tip-down pulse JFN Jan 2013*/
  /* Dummy values, will be overwritten in .e file */
  {  (int *)&pw_rfd,   /* pw */
     (float *)&a_rfd,  /* amp */
     0.5767,           /* abswidth */
     0.4291,           /* effwidth */
     0.5767,           /* area */
     1.0,              /* dtycyc */
     1.0,              /* maxpw */
     1,                /* num */
     0.117851,         /* maxb1 */
     0.00286052,       /* max_int_b1_sq */
     0.0771972,        /* max_rms_b1 */
     50.0,             /* nom_fa */
     &flip_rfd,        /* act_fa */
     480.0,            /* nom_pw */
     3125.0,           /* nom_bw */
     PSD_APS2_ON + PSD_MPS2_ON + PSD_SCAN_ON,      /* activity */
     0,                /* reference (not used) */
     0,                /* isodelay */
     1.0,              /* scale */
     (int*)&res_rfd,   /* res (not used) */
     0,                /* extgradflag (not used) */
     (int *)&wg_rfd    /* waveform generator (e.g., TYPRHO1 or TYPRHO2) */
  },
#include "rf_Prescan.h"
};

#define MAX_ENTRY_POINTS 15
float maxB1[MAX_ENTRY_POINTS], maxB1Seq;

GRAD_PULSE gradx[MAX_GRADX] = {
  {
  }
};

GRAD_PULSE grady[MAX_GRADY] = {
  {
  }
};

GRAD_PULSE gradz[MAX_GRADZ] = {
  {
  }
};

#endif  
