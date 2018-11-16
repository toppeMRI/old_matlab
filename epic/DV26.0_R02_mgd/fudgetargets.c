/* This routine scales rise/fall times by scale
   so that the PNS doesn't get you


*@End*********************************************************/


/* System header files */
#include <values.h>
#include <stdio.h>




/*
 *  fudgetargets
 *  
 *  Type: Public Function
 *  
 *
 *  Arguments:
 *    lgrad: logical gradient characteristics
 *    pgrad: physical gradient characteristics
 *  
 */
STATUS 
fudgetargets( LOG_GRAD *lgrad,
             PHYS_GRAD *pgrad, float tscale )
{


    pgrad->xrt *= tscale;
    pgrad->yrt *= tscale;
    pgrad->zrt *= tscale;
    pgrad->xft *= tscale;
    pgrad->yft *= tscale;
    pgrad->zft *= tscale;

    lgrad->xrt = lgrad->yrt = lgrad->zrt = IMax(3,cfrmp2xfs,cfrmp2yfs,cfrmp2zfs)*tscale;
    lgrad->xft = lgrad->yft = lgrad->zft = IMax(3,cffall2x0,cffall2y0,cffall2z0)*tscale;

    return SUCCESS;
}   
