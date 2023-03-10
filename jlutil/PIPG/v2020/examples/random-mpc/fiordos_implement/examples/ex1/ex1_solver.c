/**
 * This file is generated by FiOrdOs, a program licensed under GPL
 * by copyright holder Automatic Control Laboratory, ETH Zurich.
 * 
 * If you are interested in using this file commercially,
 * please contact the copyright holder.
 */

#include "ex1_solver.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

/* <<< print >>> */
static void printArray(realtype *arr, int dim) {
  int i;
  printf("[\n"); ;
  for (i=0; i<dim; i++) {
    printf("%8.16e\n",*arr++);
  }
  printf("]\n\n");
}

static void printArrayNamed(realtype *arr, int dim, char *txt) {
  printf("%s = ",txt);
  printArray(arr,dim);
}

/* <<< copy vectors >>> */
static void copyvec10(realtype *dest, const realtype *src) {
    {
        int i;
        for (i=0; i<10; ++i) {
            dest[i] = src[i];
        }
    }
}

/* <<< projection >>> */
static void algo_project(realtype *thevec, ex1_Work *work, ex1_Params *params) {
    realtype *theveci;

    /* ====== EssBox X1 ====== */
    theveci = thevec;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X1.l[icoord]) {
            theveci[icoord] = work->Prob.X1.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X1.u[icoord]) {
            theveci[icoord] = work->Prob.X1.u[icoord];
          }
        }
    }
    /* ====== EssBox X2 ====== */
    theveci = thevec+1;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X2.l[icoord]) {
            theveci[icoord] = work->Prob.X2.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X2.u[icoord]) {
            theveci[icoord] = work->Prob.X2.u[icoord];
          }
        }
    }
    /* ====== EssBox X3 ====== */
    theveci = thevec+2;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X3.l[icoord]) {
            theveci[icoord] = work->Prob.X3.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X3.u[icoord]) {
            theveci[icoord] = work->Prob.X3.u[icoord];
          }
        }
    }
    /* ====== EssBox X4 ====== */
    theveci = thevec+3;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X4.l[icoord]) {
            theveci[icoord] = work->Prob.X4.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X4.u[icoord]) {
            theveci[icoord] = work->Prob.X4.u[icoord];
          }
        }
    }
    /* ====== EssBox X5 ====== */
    theveci = thevec+4;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X5.l[icoord]) {
            theveci[icoord] = work->Prob.X5.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X5.u[icoord]) {
            theveci[icoord] = work->Prob.X5.u[icoord];
          }
        }
    }
    /* ====== EssBox X6 ====== */
    theveci = thevec+5;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X6.l[icoord]) {
            theveci[icoord] = work->Prob.X6.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X6.u[icoord]) {
            theveci[icoord] = work->Prob.X6.u[icoord];
          }
        }
    }
    /* ====== EssBox X7 ====== */
    theveci = thevec+6;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X7.l[icoord]) {
            theveci[icoord] = work->Prob.X7.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X7.u[icoord]) {
            theveci[icoord] = work->Prob.X7.u[icoord];
          }
        }
    }
    /* ====== EssBox X8 ====== */
    theveci = thevec+7;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X8.l[icoord]) {
            theveci[icoord] = work->Prob.X8.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X8.u[icoord]) {
            theveci[icoord] = work->Prob.X8.u[icoord];
          }
        }
    }
    /* ====== EssBox X9 ====== */
    theveci = thevec+8;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X9.l[icoord]) {
            theveci[icoord] = work->Prob.X9.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X9.u[icoord]) {
            theveci[icoord] = work->Prob.X9.u[icoord];
          }
        }
    }
    /* ====== EssBox X10 ====== */
    theveci = thevec+9;
    /* ------ projection ------ */
    {
        int icoord;
        for (icoord=0; icoord<1; ++icoord) {
          if (theveci[icoord] < work->Prob.X10.l[icoord]) {
            theveci[icoord] = work->Prob.X10.l[icoord];
          }
          else if (theveci[icoord] > work->Prob.X10.u[icoord]) {
            theveci[icoord] = work->Prob.X10.u[icoord];
          }
        }
    }
}

/* <<< oracle >>> */
static void algo_oracle(ex1_Work *work, ex1_Params *params, realtype *z, realtype *fval, realtype *grad) {
    realtype hz[10];
    {
        /* (hz):= (work->Prob.H)^T*(z);   // (10x10)*(10x1) */
        int i,k;
        realtype tmp;
        for (i=0; i<10; ++i) {
            tmp = 0.0;
            for (k=0; k<10; ++k) {
              tmp += work->Prob.H[k+i*10]*z[k];
            }
            hz[i] = tmp;
        }
    }

    if (NULL!=fval) {
        int i;
        realtype tmp;
        tmp = 0.0;
        for (i=0; i<10; ++i) {
          tmp += (0.5*hz[i] + params->g[i]) * z[i];
        }
        fval[0] = tmp + params->c[0];
    }

    if (NULL!=grad) {
        int i;
        for (i=0; i<10; ++i) {
          grad[i] = hz[i] + params->g[i];
        }
    }
}

/* <<< SOLVE >>> */
void ex1_solve(ex1_Params *params, ex1_Settings *settings, ex1_Result *result, ex1_Work *work) {

    /* << prepare data >> */
    {
        /* === initial points === */
        copyvec10(work->algo.init,settings->algo.init);  /* (work->algo.init):=(settings->algo.init) */
        algo_project(work->algo.init,work,params);
        
    }

    /* << main algorithm >> */
    /* AlgoFgm */
    {
        int stopgNext;
        stopgNext = 1;
        work->algo.res_stopcode = 0;
        {
            realtype z[10];
            realtype zp[10];
            realtype y[10];
            realtype grad[10];
            int iter;
            int icoord;
            realtype alpha,alpha_p;
            realtype aamul,muLinv;
            realtype beta;
            
            iter =0;
            
            copyvec10(z,work->algo.init);  /* (z):=(work->algo.init) */
            copyvec10(y,z);  /* (y):=(z) */
            alpha = 1;
            muLinv = work->algo.glob_mu*work->algo.glob_Linv;
            copyvec10(zp,z);  /* (zp):=(z) */
            
            while (iter<settings->algo.maxit) {
                ++iter;
            
                algo_oracle(work,params,y,NULL,grad);
                for (icoord=0; icoord<10; ++icoord) {
                    zp[icoord] = y[icoord] - work->algo.glob_Linv*grad[icoord];
                }
                algo_project(zp,work,params);
            
                if (iter==stopgNext) {
                    int icoord;
                    realtype gradmapi;
                    realtype normsq;
                    normsq = 0;
                    for (icoord=0; icoord<10; ++icoord) {
                        gradmapi=(y[icoord]-zp[icoord]);
                        normsq += gradmapi * gradmapi;
                    }
                    normsq *= 1.0/(work->algo.glob_Linv*work->algo.glob_Linv);
                    if ( settings->algo.stopgEps>=0 && normsq <= settings->algo.stopgEps * settings->algo.stopgEps ) {
                        work->algo.res_stopcode = 2;
                        break;
                    }
                    stopgNext += settings->algo.stopgStride;
                }
            
                aamul = alpha*alpha - muLinv;
                alpha_p= -0.5 * aamul + sqrt( aamul*aamul * 0.25 + alpha*alpha);
                beta = alpha*(1-alpha)/(alpha*alpha+alpha_p);
            
                for (icoord=0; icoord<10; ++icoord) {
                    y[icoord] = zp[icoord] + beta*(zp[icoord]-z[icoord]);
                }
            
                copyvec10(z,zp);  /* (z):=(zp) */
                alpha = alpha_p;
            }
            copyvec10(work->algo.res_z,zp);  /* (work->algo.res_z):=(zp) */
            work->algo.res_iter = iter;
        }
    }
    
    
    /* << finalize result >> */
    copyvec10(result->x,work->algo.res_z);  /* (result->x):=(work->algo.res_z) */
    algo_oracle(work,params,work->algo.res_z,&(result->f),NULL);
    result->iter = work->algo.res_iter;
    result->exitflag = work->algo.res_stopcode;

}



/* <<< INIT >>> */
void ex1_init(ex1_Params *params, ex1_Settings *settings, ex1_Result *result, ex1_Work *work) {
work->Prob.H[0] = (realtype) 1.9199932803366923e+02;
work->Prob.H[1] = (realtype) 1.7139304576073437e+02;
work->Prob.H[2] = (realtype) 1.4717895799941547e+02;
work->Prob.H[3] = (realtype) 1.2030997616426666e+02;
work->Prob.H[4] = (realtype) 9.2495002055330303e+01;
work->Prob.H[5] = (realtype) 6.5318287958048415e+01;
work->Prob.H[6] = (realtype) 4.0351704824100011e+01;
work->Prob.H[7] = (realtype) 1.9265899226280002e+01;
work->Prob.H[8] = (realtype) 3.9462587508000011e+00;
work->Prob.H[9] = (realtype) -3.3801555360000002e+00;
work->Prob.H[10] = (realtype) 1.7139304576073437e+02;
work->Prob.H[11] = (realtype) 1.6709565179495237e+02;
work->Prob.H[12] = (realtype) 1.4565216984817323e+02;
work->Prob.H[13] = (realtype) 1.2099032909297509e+02;
work->Prob.H[14] = (realtype) 9.4226930013728420e+01;
work->Prob.H[15] = (realtype) 6.7283127681060009e+01;
work->Prob.H[16] = (realtype) 4.2016704227880005e+01;
work->Prob.H[17] = (realtype) 2.0350568758800005e+01;
work->Prob.H[18] = (realtype) 4.4043053640000007e+00;
work->Prob.H[19] = (realtype) -3.3642599999999998e+00;
work->Prob.H[20] = (realtype) 1.4717895799941547e+02;
work->Prob.H[21] = (realtype) 1.4565216984817323e+02;
work->Prob.H[22] = (realtype) 1.4048573714233189e+02;
work->Prob.H[23] = (realtype) 1.1857442124137643e+02;
work->Prob.H[24] = (realtype) 9.4014970788900001e+01;
work->Prob.H[25] = (realtype) 6.8143095971880015e+01;
work->Prob.H[26] = (realtype) 4.3162072726800005e+01;
work->Prob.H[27] = (realtype) 2.1291760644000004e+01;
work->Prob.H[28] = (realtype) 4.9181394000000012e+00;
work->Prob.H[29] = (realtype) -3.2487311999999999e+00;
work->Prob.H[30] = (realtype) 1.2030997616426666e+02;
work->Prob.H[31] = (realtype) 1.2099032909297509e+02;
work->Prob.H[32] = (realtype) 1.1857442124137643e+02;
work->Prob.H[33] = (realtype) 1.1292525974634002e+02;
work->Prob.H[34] = (realtype) 9.1109056010280014e+01;
work->Prob.H[35] = (realtype) 6.7444632598799998e+01;
work->Prob.H[36] = (realtype) 4.3553915844000002e+01;
work->Prob.H[37] = (realtype) 2.2007533800000001e+01;
work->Prob.H[38] = (realtype) 5.4951828000000011e+00;
work->Prob.H[39] = (realtype) -2.9980319999999998e+00;
work->Prob.H[40] = (realtype) 9.2495002055330303e+01;
work->Prob.H[41] = (realtype) 9.4226930013728420e+01;
work->Prob.H[42] = (realtype) 9.4014970788900001e+01;
work->Prob.H[43] = (realtype) 9.1109056010280014e+01;
work->Prob.H[44] = (realtype) 8.5542256150800014e+01;
work->Prob.H[45] = (realtype) 6.4601082564000009e+01;
work->Prob.H[46] = (realtype) 4.2888311399999999e+01;
work->Prob.H[47] = (realtype) 2.2390342800000003e+01;
work->Prob.H[48] = (realtype) 6.1439880000000020e+00;
work->Prob.H[49] = (realtype) -2.5665600000000000e+00;
work->Prob.H[50] = (realtype) 6.5318287958048415e+01;
work->Prob.H[51] = (realtype) 6.7283127681060009e+01;
work->Prob.H[52] = (realtype) 6.8143095971880015e+01;
work->Prob.H[53] = (realtype) 6.7444632598799998e+01;
work->Prob.H[54] = (realtype) 6.4601082564000009e+01;
work->Prob.H[55] = (realtype) 5.9857233000000008e+01;
work->Prob.H[56] = (realtype) 4.0772518800000000e+01;
work->Prob.H[57] = (realtype) 2.2300068000000000e+01;
work->Prob.H[58] = (realtype) 6.8744400000000017e+00;
work->Prob.H[59] = (realtype) -1.8959999999999999e+00;
work->Prob.H[60] = (realtype) 4.0351704824100011e+01;
work->Prob.H[61] = (realtype) 4.2016704227880005e+01;
work->Prob.H[62] = (realtype) 4.3162072726800005e+01;
work->Prob.H[63] = (realtype) 4.3553915844000002e+01;
work->Prob.H[64] = (realtype) 4.2888311399999999e+01;
work->Prob.H[65] = (realtype) 4.0772518800000000e+01;
work->Prob.H[66] = (realtype) 3.7701348000000003e+01;
work->Prob.H[67] = (realtype) 2.1555239999999998e+01;
work->Prob.H[68] = (realtype) 7.6980000000000004e+00;
work->Prob.H[69] = (realtype) -9.1199999999999948e-01;
work->Prob.H[70] = (realtype) 1.9265899226280002e+01;
work->Prob.H[71] = (realtype) 2.0350568758800005e+01;
work->Prob.H[72] = (realtype) 2.1291760644000004e+01;
work->Prob.H[73] = (realtype) 2.2007533800000001e+01;
work->Prob.H[74] = (realtype) 2.2390342800000003e+01;
work->Prob.H[75] = (realtype) 2.2300068000000000e+01;
work->Prob.H[76] = (realtype) 2.1555239999999998e+01;
work->Prob.H[77] = (realtype) 2.0922000000000001e+01;
work->Prob.H[78] = (realtype) 8.6280000000000001e+00;
work->Prob.H[79] = (realtype) 4.8000000000000043e-01;
work->Prob.H[80] = (realtype) 3.9462587508000015e+00;
work->Prob.H[81] = (realtype) 4.4043053640000016e+00;
work->Prob.H[82] = (realtype) 4.9181394000000012e+00;
work->Prob.H[83] = (realtype) 5.4951828000000011e+00;
work->Prob.H[84] = (realtype) 6.1439880000000020e+00;
work->Prob.H[85] = (realtype) 6.8744400000000017e+00;
work->Prob.H[86] = (realtype) 7.6980000000000004e+00;
work->Prob.H[87] = (realtype) 8.6280000000000001e+00;
work->Prob.H[88] = (realtype) 1.0680000000000000e+01;
work->Prob.H[89] = (realtype) 2.4000000000000004e+00;
work->Prob.H[90] = (realtype) -3.3801555360000002e+00;
work->Prob.H[91] = (realtype) -3.3642599999999998e+00;
work->Prob.H[92] = (realtype) -3.2487311999999999e+00;
work->Prob.H[93] = (realtype) -2.9980319999999998e+00;
work->Prob.H[94] = (realtype) -2.5665600000000000e+00;
work->Prob.H[95] = (realtype) -1.8959999999999999e+00;
work->Prob.H[96] = (realtype) -9.1199999999999948e-01;
work->Prob.H[97] = (realtype) 4.8000000000000043e-01;
work->Prob.H[98] = (realtype) 2.4000000000000004e+00;
work->Prob.H[99] = (realtype) 6.0000000000000000e+00;
work->Prob.X1.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X1.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X2.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X2.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X3.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X3.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X4.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X4.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X5.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X5.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X6.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X6.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X7.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X7.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X8.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X8.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X9.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X9.u[0] = (realtype) 7.0000000000000000e+00;
work->Prob.X10.l[0] = (realtype) -7.0000000000000000e+00;
work->Prob.X10.u[0] = (realtype) 7.0000000000000000e+00;
work->algo.glob_Linv = (realtype) 1.4494698026030653e-03;
work->algo.glob_mu = (realtype) 3.4691582875472093e+00;
{int i; for (i=0;i<10;++i) {settings->algo.init[i] = (realtype) 0.0;}}
settings->algo.maxit = 5000;
settings->algo.stopgEps = (realtype) 1.0000000000000000e-03;
settings->algo.stopgStride = 1;
}

