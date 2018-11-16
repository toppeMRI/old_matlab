/*@Start***********************************************************/
/* GEMSBG Include File
 * Copyright (C) 1995 The General Electric Company
 *
 *      Include File Name:  grad_rf_grass.globals
 *      Developer:              T. Hlaban        Original for 5.5
 *
 * $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/DV26.0_R02_mgd/grad_rf_sprlio.globals.h,v $
 * $Revision: 1.1 $  $Date: 2018/10/04 19:01:19 $
 *      prescan.globals.h        10/1/95 ghg
 */

/*@Synopsis
  This has global #defines for ipg & host
*/

/*@Description

*/

/*@End*********************************************************/

/* only do this once in any given compilation.*/
#ifndef  grad_rf_globals_sprlio_INCL
#define  grad_rf_globals_sprlio_INCL

#define MAX_RFPULSE 25
#define MAX_GRADX 26
#define MAX_GRADY 26
#define MAX_GRADZ 26

#define RFD_SLOT 0  /* bSSFP/STFR tip-down pulse JFN Jan 2013 */
#define RF_FREE1 1  /* changed from 0 to 2, JFN Jan 2013 */

#define GX_FREE 0

#define GY_FREE 0

#define GZ_FREE 0

#include "rf_Prescan.globals.h"

#endif 
