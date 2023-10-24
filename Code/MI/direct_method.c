#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Needed for building mex file:                                           */
#include <mex.h>
#include <matrix.h>

void panzeri_treves_96(double *C_ptr, mwSize Xtot, double N, double *bias);

/*
 *   Copyright (C) 2010 Cesare Magri
 *   Version: 6a
 */

/*
 * -------
 * LICENSE
 * -------
 * This software is distributed free under the condition that:
 *
 * 1. it shall not be incorporated in software that is subsequently sold;
 *
 * 2. the authorship of the software shall be acknowledged and the following
 *    article shall be properly cited in any publication that uses results
 *    generated by the software:
 *
 *      Magri C, Whittingstall K, Singh V, Logothetis NK, Panzeri S: A
 *      toolbox for the fast information analysis of multiple-site LFP, EEG
 *      and spike train recordings. BMC Neuroscience 2009 10(1):81;
 *
 * 3.  this notice shall remain in place in each source file.
 *
 */

/* ----------
 * DISCLAIMER
 * ----------
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


/* Global variables:                                                       */
static bool initialized = false;
static mxArray *ClogC_mxArray = NULL;


/* ========================================================================*/
/* VECTORIAL KRONECKER'S PRODUCT                                           */
/* ========================================================================*/
void vecKronProd(
                 double  *xPtr,      /* input: pointer to array x          */
                 mwSize   xLen,      /* input: length of x                 */
                 double  *yPtr,      /* input: pointer to array y          */
                 mwSize   yLen,      /* input: length of y                 */
                 double **K_ptr_ptr, /* output                             */
                 mwSize*  KLen_ptr   /* output                             */
                )
{   
    double *K_ptr;
    mwIndex i, j, indx;
    
    /* Where to store the final array:                                     */
    K_ptr = *K_ptr_ptr;
    
    /* Length of final array:                                              */
    *KLen_ptr = xLen*yLen;
    
    K_ptr = mxCalloc(*KLen_ptr, sizeof(double));
    
    /* In the main funciton we map bins into responses as follows:         */
    /*     r = X[0,t,s] + X[1,t,s]*Nb[0] + X[2,t,s]*Nb[0]*Nb[1] + ...      */
    /* this explains why we first loop over j and then over i in the       */
    /* following nested loop:                                              */
    indx = 0;
    for(j=0; j<yLen; j++)
        for(i=0; i<xLen; i++)
            K_ptr[indx++] = xPtr[i] * yPtr[j];
    
    /* Storing result in the right place:                                  */
    *K_ptr_ptr = K_ptr;
}



/* ========================================================================*/
/* DIVIDE AND CONQUER                                                      */
/* ========================================================================*/
void splitProds(
                double *X_ptr, /* pointer to data matrix                   */
                mwSize M,      /* number of lines of X                     */
                mwSize N,      /* number of columns of X                   */
                mwSize* n_ptr, /* pointer to array which specifies how many*/
                               /* elements in each row of X are really     */
                               /* occupied by data                         */
                
                double** K_ptr_ptr,
                mwSize*  KLen_ptr
               )
{      
    double *K_ptr;
    

    /* If the number of columns of X is greater than two the fucntion calls*/
    /* itself recursively. The following variable are thus used to break   */
    /* the problem in two parts.                                           */
    double *K1_ptr, *K2_ptr;
    mwSize  K1Len, K2Len;
    mwSize  N1, N2;
    
    K_ptr = *K_ptr_ptr;
    
    /* If N=2 compute the product directly                                 */
    if(N==2) {
        vecKronProd(&X_ptr[0], n_ptr[0], &X_ptr[M], n_ptr[1], &K_ptr, KLen_ptr);
        
    /* Otherwise split X in two parts and work on each half independently  */
    } else {
        /* First haf ------------------------------------------------------*/
        N1 = N/2; /* note: floored int ratio */
        if(N1==1) {
            K1_ptr = X_ptr;
            K1Len = n_ptr[0];
        } else {
            splitProds(X_ptr, M, N1, n_ptr, &K1_ptr, &K1Len);
        }
        
        /* Second half ----------------------------------------------------*/
        N2 = N - N1;
        splitProds(&X_ptr[N1*M], M, N2, &n_ptr[N1], &K2_ptr, &K2Len);
        
        /* Final product --------------------------------------------------*/
        vecKronProd(K1_ptr, K1Len, K2_ptr, K2Len, &K_ptr, KLen_ptr);
        
        if(N1>1) /* Freeing K1 when N1=1 would erase part of X */
            mxFree(K1_ptr);

        mxFree(K2_ptr);
        
    }
    
    *K_ptr_ptr = K_ptr;
}



/* ========================================================================*/
/* CREATE ClogC                                                            */
/* ========================================================================*/
double *create_ClogC(mwSize Ns, mwSize Nt)
{
    double *ClogC;
    
    ClogC_mxArray = mxCreateDoubleMatrix(Ns*Nt+100, 1, mxREAL);
    mexMakeArrayPersistent(ClogC_mxArray);
    ClogC = mxGetPr(ClogC_mxArray);
    
    return ClogC;
}



/* ========================================================================*/
/* CLEANUP FUNCTION                                                        */
/* ========================================================================*/
void cleanup(void) {
    mexPrintf("MEX-file is terminating, destroying MEX arrays\n");
    mxDestroyArray(ClogC_mxArray);
}



/* ========================================================================*/
/* INCREMENT                                                               */
/* ========================================================================*/
double increment(double C, double *ClogC)
{
    double out;
    mwSize indx;
    
    indx = (mwIndex) C;

    if(ClogC[indx+1] == 0)
        ClogC[indx] = (C) * log(C);

    out = ClogC[indx] - ClogC[indx-1];
        
    return out;
}



/* ========================================================================*/
/* ENTROPY COMPUTATION (GATEWAY FUNCTION)                                  */
/* ========================================================================*/
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

    /* Input --------------------------------------------------------------*/
    
    double  *X;         /* response matrix                                 */
    double  *Nt;        /* number of trials per stimulus                   */
    int     biasCorrNum;

    bool     doHX;      /* compute H(X)       (flag)                       */
    bool     doHXY;     /* compute H(X|Y)     (flag)                       */
    bool     doHlX;     /* compute H_lin(X)   (flag)                       */
    bool     doHlXY;    /* compute H_lin(X|Y) (flag)                       */
    bool     doHiX;     /* compute H_ind(X)   (flag)                       */
    bool     doHiXY;    /* compute H_ind(X|Y) (flag)                       */
    bool     doChiX;    /* compute Chi(X)     (flag)                       */
    bool     doHshX;    /* compute H_sh(X)    (flag)                       */
    bool     doHshXY;   /* compute H_sh(X|Y)  (flag)                       */
    
    bool     testMode;   /* specifies if working in test condition  (flag)  */
    
    /* Internal variables -------------------------------------------------*/
    
    mxArray *PrToInputMxArray;  /* generic pointer to mxArray, used to read*/
                                /* several input                           */
    
    /* Other useful flags:                                                 */
    bool     doMap;     /* map multi-dim response into scalar (flag)       */
    bool     doCr;      /* compute C(r)       (flag)                       */
    bool     doPir;     /* compute P_ind(r)   (flag)                       */
    bool     doCirs;    /* compute C_ind(r|s) (flag)                       */
    bool     doCrcs;    /* compute C(r_c|s)   (flag)                       */
    bool     doSh;      /* do the shuffling   (flag)                       */
    
    /* Dimensions ---------------------------------------------------------*/
    mwSize   L;         /* number of cells                                 */
    mwSize   Ns;        /* number of stimuli                               */
    mwSize  *Nb;        /* number of bins per cell                         */
    mwSize   Nr;        /* number of possible responses                    */
    
    mwSize   maxNb;     /* max number of bins: maxNb=max(Nb)               */
    mwSize  *cumprodNb; /* cumulative product of Nb:                       */
                        /*     cumprodNb(n) = Nb(n)*Nb(n-1)*...*Nb(1)      */
    mwSize  *base;      /* base used for mapping each response from an     */
                        /* L-dimensional array to single number:           */
                        /*     base = cumprodNb ./ Nb(1)                   */
    
    /* Storing useful dimension-related quantities:                        */
    mwSize   maxNt;             /* max number of trials: maxNt=max(Nt)     */
    double   inv_Nt;            /* inverse of Nt: inv_Nt=1/Nt              */
    double   NtpowL;            /* Nt power L: NtpowL=Nt^L                 */
    double   inv_NtpowL;        /* inverse of NtpowL: inv_NtpowL=1/(Nt^L)  */
    double   inv_Ntpow_Lminus1; /* inverse of Nt power (L - 1):            */
                                /*   inv_Ntpow_Lminus1 = 1/(Nt^(L-1))      */
    
    double   totNt;     /* total numner of trials: totNt=sum(Nt)           */
    double   inv_totNt; /* inverse of totNt: inv_totNt=1/totNt             */
    
    /* Other useful quantities                                             */
    double   log_2 = log(2);

    /* Indexes                                                             */
    mwIndex  c;         /* cell index                                      */
    mwIndex  t;         /* trial index                                     */
    mwIndex  s;         /* stimulus index                                  */
    mwIndex  b;         /* bin index                                       */
    mwIndex  r;         /* response index                                  */

    mwIndex  i;         /* counter                                         */
    mwIndex  j;         /* counter                                         */
    mwIndex  indx;      /* used used to store index to elements of a matrix*/

    /* Variables used for computing cell-shuffled quantities:              */
    double  *Xsh;       /* shuffled response matrix                        */
    mwIndex  rsh;       /* shuffled response index                         */
    mwIndex  bsh;       /* shuffled bin index;                             */
    mwIndex  tsh;       /* shuffled trial index                            */
    mwIndex  indxSh;    /* used used to store index to an element of a     */
                        /* shuffled matrix                                 */

    /* Count and probability matrices:                                     */
    double  *Cr;        /* pointer to C(r) values                          */
    double  *Crs;       /* pointer to C(r|s) values                        */
    double  *Crc;       /* pointer to C(rc) values                         */
    double  *Crcs;      /* pointer to C(rc|s) values                       */
    double  *Pir;       /* pointer to C_ind(r) values                      */
    double   Cirs;      /* pointer to C_ind(rc|s) values                   */
    double  *Cshr;      /* pointer to C_sh(r) values                       */
    double  *Cshrs;     /* pointer to C_sh(r|s) values                     */
    
    /* C*log(C) function:                                                  */
    double  *ClogC;
    mwSize   ClogCLen;
    
    /* Variables needed for calls to splitProd                             */
    double  *X_ptr;
    mwSize   M;
    mwSize   N1;
    mwSize   N2;
    mwSize  *n_ptr;
    
    double  *K1_ptr;
    double  *K2_ptr;
    
    mwSize   K1Len;
    mwSize   K2Len;
    
    /* Bias Values                                                         */
    double  *biasPtr;
    double   biasHlXY;
    double   biasHlX;
    
    /* Output -------------------------------------------------------------*/
    
    double  *HX;
    double  *HXY;
    double  *HlX;
    double  *HlXY;
    double  *HiX;
    double  *HiXY;
    double  *ChiX;
    double  *HshX;
    double  *HshXY;
    
    mwIndex k;
    
    
    /* Reading input ------------------------------------------------------*/
    X  = mxGetPr(prhs[0]);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "Nt");
    Nt = mxGetPr(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "biasCorrNum");
    biasCorrNum = (int) *mxGetPr(PrToInputMxArray);
        
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHX");
    doHX = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHXY");
    doHXY = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHlX");
    doHlX = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHlXY");
    doHlXY = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHiX");
    doHiX = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHiXY");
    doHiXY = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doChiX");
    doChiX = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHshX");
    doHshX = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "doHshXY");
    doHshXY = *mxGetLogicals(PrToInputMxArray);
    
    PrToInputMxArray = mxGetField(prhs[1], 0, "testMode");
    testMode = *mxGetLogicals(PrToInputMxArray);
    
    /* Other useful flags:                                                 */
    doMap   = doHX || doHXY || doChiX;
    doCr    = doHX || doChiX;
    doPir   = doHiX || doChiX;
    doCirs  = doPir || doHiXY;
    doCrcs  = doHlXY || doCirs;
    doSh    = doHshX || doHshXY;

    /* Computing dimensions -----------------------------------------------*/
    
    L     = mxGetDimensions(prhs[0])[0];
    maxNt = mxGetDimensions(prhs[0])[1];
    /* If the users inputs a response matrix with a single stimulus then   */
    /* the number of dimensions for the response matrix will be 2 and      */
    /* asking for the number of elements in the third dimension will return*/
    /* a meaningless value. Thus first we need to check that there is      */
    /* indeed a third dimension. We do this using MXGETNUMBEROFDIMENSIONS  */
    if(mxGetNumberOfDimensions(prhs[0])>2)
        Ns = mxGetDimensions(prhs[0])[2];
    else
        Ns = 1;
    
    Nb         = mxCalloc(L, sizeof(mwSize));
    base       = mxCalloc(L, sizeof(mwSize));
    cumprodNb  = mxCalloc(L, sizeof(mwSize));
    
    if(doSh)
        Xsh = mxCalloc(L * maxNt * Ns, sizeof(double));
    
    maxNb = 0;
    Nr    = 1;
    for(c=0; c<L; c++) {
            
        for(s=0; s<Ns; s++) {
            
            for(t=0; t<Nt[s]; t++) {

                indx = c + t*L + s*L*maxNt;

                /* Find max(X(c,:,:))                                      */
                if(X[indx] > Nb[c]) {
                    Nb[c] = X[indx];
                }

                /* Shuffling                                               */
                /* ---------                                               */
                /* The shuffling routine is based on the following simple  */
                /* algorithm for permuting an array (i.e., any ordering    */
                /* of the array is equally probable or any element has     */
                /* equal chance of being in any position):                 */
                /*                                                         */
                /*     For i=1 to n-1:                                     */
                /*     - Let j=r(i+1)                                      */
                /*     - Swap a[i] and a[j]                                */
                /*                                                         */
                /* where r[n] is a random number generated between 0 and   */
                /* n-1.                                                    */

                if(doSh) {
                    /* In test-mode we do not want to shuffle              */
                    if(testMode) {
                        Xsh[indx] = X[indx];
                    }
                    else {
                    /* Generating random number between 0 and t            */
                    tsh = rand() % (t+1);
                    indxSh = c + tsh*L + s*L*maxNt;

                    /* Swapping X[c, t, s] and X[c, tsh, s]                */
                    Xsh[indx]   = Xsh[indxSh];
                    Xsh[indxSh] = X[indx];
                    }
                }
            }
        }
        
        /* The response is assumed to start from zero, the number of bins  */
        /* must thus be increased by one                                   */
        Nb[c]++;
    
        /* Finding max number of bins for cell c                           */
        if(Nb[c] > maxNb) {
            maxNb = Nb[c];
        }
        
        /* Building cumulative product                                     */
        if(c==0) {
            cumprodNb[c] = Nb[c];
            base[c]      = 1;
        } else {
            cumprodNb[c] = cumprodNb[c-1] * Nb[c];
            base[c]      = cumprodNb[c-1];
        }
        
        /* Number of possible responses                                    */
        Nr *= Nb[c];
    }

	/* Checking Nr size ---------------------------------------------------*/
	if (Nr == 0) {
		mexErrMsgIdAndTxt("MI_DirectMethod:ResponseDimErr",
						  "Number of possible responses for neuron/population is zero. \n" \
						  "This is likely due to a too high number of bins, leading to "\
						  "a number of possible responses that cannot be represented by "\
						  "an integer number.\nAborting now.");
	}

    /* Allocating memory --------------------------------------------------*/
    if(doCr)    Cr    = mxCalloc(Nr             , sizeof(double));
    if(doHXY)   Crs   = mxCalloc(Nr * Ns        , sizeof(double));
    if(doHlX)   Crc   = mxCalloc(maxNb * L      , sizeof(double));
    if(doCrcs)  Crcs  = mxCalloc(maxNb * L * Ns , sizeof(double));
    if(doPir)   Pir   = mxCalloc(cumprodNb[L-1] , sizeof(double));
    if(doHshX)  Cshr  = mxCalloc(Nr             , sizeof(double));
    if(doHshXY) Cshrs = mxCalloc(Nr * Ns        , sizeof(double));
   

    /* HX                                                                  */
    plhs[0] = mxCreateDoubleScalar(0);
    HX = mxGetPr(plhs[0]);
    
    /* HXY                                                                 */
    if(doHXY) plhs[1] = mxCreateDoubleMatrix(Ns,1,mxREAL);
    else      plhs[1] = mxCreateDoubleScalar(0);
    HXY = mxGetPr(plhs[1]);

    /* HlX                                                                 */
    plhs[2] = mxCreateDoubleScalar(0);
    HlX = mxGetPr(plhs[2]);
    
    /* HlXY                                                                */
    if(doHlXY) plhs[3] = mxCreateDoubleMatrix(Ns,1,mxREAL);    
    else       plhs[3] = mxCreateDoubleScalar(0);
    HlXY = mxGetPr(plhs[3]);
    
    /* HiX                                                                 */
    plhs[4] = mxCreateDoubleScalar(0);
    HiX = mxGetPr(plhs[4]);
    
    /* HiXY                                                                */
    if(doHiXY) plhs[5] = mxCreateDoubleMatrix(Ns,1,mxREAL);
    else       plhs[5] = mxCreateDoubleScalar(0);
    HiXY = mxGetPr(plhs[5]);    

    /* ChiX                                                                */
    plhs[6] = mxCreateDoubleScalar(0);
    ChiX = mxGetPr(plhs[6]);
    
    /* HshX                                                                */
    plhs[7] = mxCreateDoubleScalar(0);
    HshX = mxGetPr(plhs[7]);
    
    /* HshXY                                                               */
    if(doHshXY) plhs[8] = mxCreateDoubleMatrix(Ns,1,mxREAL);
    else        plhs[8] = mxCreateDoubleScalar(0);
    HshXY = mxGetPr(plhs[8]);
    
    biasPtr = mxCalloc(1, sizeof(double));
    
    /* Initializing ClogC -------------------------------------------------*/
    
    /* If ClogC has not been already initialized:                          */
    if (!initialized) {
        /* Create persistent array and register its cleanup:               */
        ClogC = create_ClogC(Ns, maxNt);
        mexAtExit(cleanup);
        initialized = true;
    }
    else {
        /* If the size of ClogC is too small for the new input matrix we   */
        /* re-initialize it.                                               */
        ClogCLen = mxGetNumberOfElements(ClogC_mxArray);
        
        if(ClogCLen<Ns*maxNt) {
            ClogC = mxGetPr(ClogC_mxArray);
            mxFree(ClogC);
            ClogC = create_ClogC(Ns, maxNt);
        }
        else
            ClogC = mxGetPr(ClogC_mxArray);
    }
   
    /* Computing entropy --------------------------------------------------*/
    totNt = 0;
    biasHlX  = 0;
    for(s=0; s<Ns; s++) {
        
        totNt += Nt[s];
        inv_Nt = 1 / Nt[s];
        
            NtpowL        = pow(Nt[s],L);
        inv_NtpowL        = 1 / NtpowL;
        inv_Ntpow_Lminus1 = 1 / pow(Nt[s],L-1);
        
        biasHlXY = 0;

        for(t=0;t<Nt[s];t++) {
            
            r = 0;
            rsh = 0;
            
            for(c=0;c<L;c++) {
                
                b = X[c + t*L + s*L*maxNt];
                if(doMap)
                    r += b * base[c];
                
                if(doSh) {
                    bsh = Xsh[c + t*L + s*L*maxNt];            
                    rsh += bsh * base[c];
                }
            
                if(doCrcs) {
                    indx = b + c*maxNb + s*maxNb*L;
                    ++Crcs[indx];
                    
                    if(doHlXY) {
                        HlXY[s] += increment(Crcs[indx], ClogC);
                    }
                }

                if(doHlX) {
                    indx = b + c*maxNb;
                    ++Crc[indx];
                    
                    *HlX += increment(Crc[indx], ClogC);
                }
            } /* End of loop on c ---------------------------------------- */
            
            if(doCr) {
                ++Cr[r];
                                
                if(doHX)
                    HX[0] += increment(Cr[r], ClogC);

            }
            
            if(doHXY) {
                indx = r + s*Nr;
                ++Crs[indx];
                
                HXY[s] += increment(Crs[indx], ClogC);
            }
            
            if(doHshX) {
                ++Cshr[rsh];
                
                *HshX += increment(Cshr[rsh], ClogC);
            }
            
            if(doHshXY) {
                indx = rsh + s*Nr;
                ++Cshrs[indx];
                
                HshXY[s] += increment(Cshrs[indx], ClogC);
            }
        } /* End of loop on t -------------------------------------------- */
        
        if(doCirs) {

            /* This part is essentially identical to part of splitProds.   */
            /* However, since we want to perform the last product in the   */
            /* main function we need to repeat this piece of code (this has*/
            /* the advantage of greatly reducing the number of variables   */
            /* that need to be passed to the subfunctions).                */
            
            /* If L=2 we just do the product directly                      */
            if(L==2) {
                K1_ptr = &Crcs[s*maxNb*L];
                K2_ptr = &Crcs[s*maxNb*L + maxNb];
                
                K1Len = Nb[0];
                K2Len = Nb[1];
            } else {
                
                /* First half                                              */
                N1 = L/2; /* floored ratio                                 */
                
                if(N1 == 1) {
                    /* If N1=1 we it must be L=3, thus we split the        */
                    /* products as                                         */
                    /*      P1 * (P2 * P3)                                 */

                    /* In this case for K1 we point directly to the        */
                    /* beginning of Crcs and its length is Nb[0]           */
                    K1_ptr = &Crcs[s*maxNb*L];
                    K1Len = Nb[0];
                } else {
                    X_ptr = &Crcs[s*maxNb*L];
                    M = maxNb;

                    n_ptr = &Nb[0];
                    
                    splitProds(X_ptr, M, N1, n_ptr, &K1_ptr, &K1Len);
                }
                                
                X_ptr = &Crcs[s*maxNb*L + N1*maxNb];
                M = maxNb;
                N2 = L - N1;
                n_ptr = &Nb[N1];
                
                splitProds(X_ptr, M, N2, n_ptr, &K2_ptr, &K2Len);
            }

            /* We perform the last Kronecker product here                  */
            indx = 0;
            for(j=0; j<K2Len; j++) {
                for(i=0; i<K1Len; i++) {
                    
                    Cirs = K1_ptr[i] * K2_ptr[j];
                    
                    if(Cirs>0) {
                        if(doHiXY)
                            HiXY[s] += Cirs * log(Cirs);

                        if(doPir)
                            Pir[indx] += inv_Ntpow_Lminus1 * Cirs;
                    }

                    if(s==Ns-1) {
                        if(doHiX && Pir[indx]>0)
                            *HiX += Pir[indx] * log(Pir[indx]);

                        if(doChiX && Pir[indx]>0)
                            *ChiX +=  Cr[indx] * log(Pir[indx]);
                        
                    }
                    
                    ++indx;
                }
            }
            
            /* Freeing K1 and K2 when L=2 would erase Crcs                 */
            if(L> 2) {
                /* Freeing K1 when N1=1 would erase part of Crcs           */
                if(N1>1) {
                    mxFree(K1_ptr);
                }
                mxFree(K2_ptr);
            }
            
        }
        
        if(doHXY) {
            HXY[s] = (log(Nt[s]) - inv_Nt * HXY[s]) / log_2;
        }
        
        if(doHshXY) {
            HshXY[s] = (log(Nt[s]) - inv_Nt * HshXY[s]) / log_2;
        }
        
        if(doHlXY) {
            HlXY[s] = (log(Nt[s]) * L - inv_Nt * HlXY[s]) / log_2;
        }
        
        if(doHiXY)
            HiXY[s] = (log(NtpowL) - inv_NtpowL * HiXY[s]) / log_2;

    } /* End of loop on s ------------------------------------------------ */
    
    inv_totNt  = 1 / totNt;
    
    if(doHX) {
        HX[0] = (log(totNt) - inv_totNt * HX[0]) / log_2;
    }

    if(doHshX) {
        HshX[0] = (log(totNt) - inv_totNt * HshX[0]) / log_2;
    }

    if(doHlX) {
        HlX[0]  = (log(totNt) * L - inv_totNt * HlX[0]) / log_2;
    }

    if(doHiX)
        HiX[0]  = (log(totNt) - inv_totNt * HiX[0]) / log_2;

    if(doChiX)
        ChiX[0] = (log(totNt) - inv_totNt * ChiX[0]) / log_2;

    /* Freeing memory                                                      */
    if(doHX)    mxFree(Cr);
    if(doHXY)   mxFree(Crs);
    if(doHlX)   mxFree(Crc);
    if(doCrcs)  mxFree(Crcs);
    if(doPir)   mxFree(Pir);
    if(doHshX)  mxFree(Cshr);
    if(doHshXY) mxFree(Cshrs);
}
