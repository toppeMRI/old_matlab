#ifndef JFN_GLOBALDEFS_H_
#define JFN_GLOBALDEFS_H_

#define JFN_FAILURE -1
#define JFN_SUCCESS 0

/* NL = size of 2nd dimension of looparr (implemented as 1D array) */
/* NL = 10 for toppe.e (see README) */
/* NL = 11 for toppe2.e */
/* NL = 13 for toppe3.e, toppe4.e, toppe7.e, toppe7c.e */
/* NL = 14 for toppe5.e, toppe7d.e */
/* #define NL 15 for toppev1 */
#define NL 16    /*  toppev2 */

/* #define MAX_PG_IAMP 32767 */
#define JFN_GRADUPDATETIME 4

#define M_PI            3.14159265358979323846
#define TM_PI           6.28318530717958647692
#define JFN_GAMMA       4.25759   /* KHz/G */

#define JFN_MIN(a,b) ((a) < (b) ? (a) : (b))
#define JFN_MAX(a,b) ((a) > (b) ? (a) : (b))
#define JFN_ABS(x)  ((x) > 0 ? (x) : -(x))

#endif
