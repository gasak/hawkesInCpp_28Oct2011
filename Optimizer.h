//
//  Optimizer.h
//  hawkesInCpp
//
//  Created by Siva Kumar Gorantla on 9/26/11.
//  Copyright 2011 University of Illinois. All rights reserved.
//

#ifndef hawkesInCpp_Optimizer_h
#define hawkesInCpp_Optimizer_h
#include "Matrix.h"
#include "MyVector.h"
#include "math.h"
#define big 1.0e+35
#include "OptimFn.h"

class Optimizer{
public:
        
    void performNMmin(Matrix<double> initParamsMat, int dimension, Hawkes &h, const SubSpikeStruct& ss, int numParams, double* X, double abstol, double intol, double alpha, double bet, double gamm, int trace, int *fncount, int maxit,double *Fmin,int *fail);
    
    void nmmin(int n, double *Bvec, double *X, double *Fmin, OptimFn* fminfn, int *fail, double abstol, double intol,
               double alpha, double bet, double gamm, int trace,
               int *fncount, int maxit);
    
    void performBFGS(Matrix<double> initParamsMat, int dimension, Hawkes &h, const SubSpikeStruct& ss, int numParams, double* X, double abstol, double intol, double alpha, double bet, double gamm, int trace, int *fncount, int maxit,double *Fmin,int *fail);
    
    void vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr, int maxit, int trace, int *mask, double abstol, double reltol, int nREPORT, void *ex, int *fncount, int *grcount, int *fail);

    

private:
    
    bool isFinite(double f){
        bool ret;
        if ((f > 2.22507e-308) && (f<1.79769e+308)) {
            ret = true; //cout << "f is finite " << f << endl;
        }
        else{
            ret = false; cout << "f is not finite " << f << endl; }
        return ret;
    }

};



/* Nelder-Mead */
void Optimizer::nmmin(int n, double *Bvec, double *X, double *Fmin, OptimFn* fminfn,
           int *fail, double abstol, double intol,
           double alpha, double bet, double gamm, int trace,
           int *fncount, int maxit)
{
    char action[50];
    int C;
    bool calcvert;
    double convtol, f;
    int funcount=0, H, i, j, L=0;
    int n1=0;
    double oldsize;
    double **P;
    double size, step, temp, trystep;
    char tstr[6];
    double VH, VL, VR;
    
    if (maxit <= 0) {
        /* for(int iter1 = 0 ; iter1 < n; iter1++){
            initParamsMat.Set_Element(dimension,iter1,Bvec[i]);
        }
        *Fmin = h.computeNegativeLogLikelihood(ss, initParamsMat, dimension);*/
        //*Fmin = fminfn->getFunVal(n, Bvec);
        *Fmin = fminfn->getFunVal(n, Bvec); 
        *fncount = 0;
        *fail = 0;
        return;
    }
    if (trace)
        cout << "  Nelder-Mead direct search function minimizer" << endl;
    //P = matrix(n, n+1);
    P = new double*[n+1]; 
    for(int iter1 = 0; iter1 < n+1; iter1++){
        P[iter1] = new double[n+2];
    }
    *fail = false;
    f = fminfn->getFunVal(n, Bvec); 
    //f = fminfn(n, Bvec, ex);
    if (!isFinite(f)) {
        cout << "function cannot be evaluated at initial parameters" << endl;
        *fail = true;
    } else {
        if (trace) cout << "function value for initial parameters = " << f << endl;
        funcount = 1;
        convtol = intol * (fabs(f) + intol);
        if (trace) cout << "  Scaled convergence tolerance is " << convtol << endl;
        n1 = n + 1;
        C = n + 2;
        P[n1 - 1][0] = f;
        for (i = 0; i < n; i++)
            P[i][0] = Bvec[i];
        
        L = 1;
        size = 0.0;
        
        step = 0.0;
        for (i = 0; i < n; i++) {
            if (0.1 * fabs(Bvec[i]) > step)
                step = 0.1 * fabs(Bvec[i]);
        }
        if (step == 0.0) step = 0.1;
        if (trace) cout << "Stepsize computed as "<< step << "\n";
        for (j = 2; j <= n1; j++) {
            strcpy(action, "BUILD          ");
            for (i = 0; i < n; i++)
                P[i][j - 1] = Bvec[i];
            
            trystep = step;
            while (P[j - 2][j - 1] == Bvec[j - 2]) {
                P[j - 2][j - 1] = Bvec[j - 2] + trystep;
                trystep *= 10;
            }
            size += trystep;
        }
        oldsize = size;
        calcvert = true;
        do {
            if(trace) {
                cout << "Iter " << funcount << " ;Params =";
                for(int junknumber = 0 ; junknumber < n; junknumber++){
                    cout << Bvec[junknumber] << ",";
                }
                cout << "\n";
            }
            if (calcvert) {
                for (j = 0; j < n1; j++) {
                    if (j + 1 != L) {
                        for (i = 0; i < n; i++)
                            Bvec[i] = P[i][j];
                        f = fminfn->getFunVal(n, Bvec);   
                        //f = fminfn(n, Bvec, ex);
                        if (!isFinite(f)) f = big;
                        funcount++;
                        P[n1 - 1][j] = f;
                    }
                }
                calcvert = false;
            }
            
            VL = P[n1 - 1][L - 1];
            VH = VL;
            H = L;
            
            for (j = 1; j <= n1; j++) {
                if (j != L) {
                    f = P[n1 - 1][j - 1];
                    if (f < VL) {
                        L = j;
                        VL = f;
                    }
                    if (f > VH) {
                        H = j;
                        VH = f;
                    }
                }
            }
            
            if (VH <= VL + convtol || VL <= abstol) break;
            //if (VH <= VL + convtol) break;

            sprintf(tstr, "%5d", funcount);
            if (trace) cout << action << tstr << " " << VH << " " << VL << endl;
            //Rprintf("%s%s %f %f\n", action, tstr, VH, VL);
            
            for (i = 0; i < n; i++) {
                temp = -P[i][H - 1];
                for (j = 0; j < n1; j++)
                    temp += P[i][j];
                P[i][C - 1] = temp / n;
            }
            for (i = 0; i < n; i++)
                Bvec[i] = (1.0 + alpha) * P[i][C - 1] - alpha * P[i][H - 1];
            f = fminfn->getFunVal(n, Bvec);  
            //f = fminfn(n, Bvec, ex);
            if (!isFinite(f)) f = big;
            funcount++;
            strcpy(action, "REFLECTION     ");
            VR = f;
            if (VR < VL) {
                P[n1 - 1][C - 1] = f;
                for (i = 0; i < n; i++) {
                    f = gamm * Bvec[i] + (1 - gamm) * P[i][C - 1];
                    P[i][C - 1] = Bvec[i];
                    Bvec[i] = f;
                }
                f = fminfn->getFunVal(n, Bvec);   
                //f = fminfn(n, Bvec, ex);
                if (!isFinite(f)) f = big;
                funcount++;
                if (f < VR) {
                    for (i = 0; i < n; i++)
                        P[i][H - 1] = Bvec[i];
                    P[n1 - 1][H - 1] = f;
                    strcpy(action, "EXTENSION      ");
                } else {
                    for (i = 0; i < n; i++)
                        P[i][H - 1] = P[i][C - 1];
                    P[n1 - 1][H - 1] = VR;
                }
            } else {
                strcpy(action, "HI-REDUCTION   ");
                if (VR < VH) {
                    for (i = 0; i < n; i++)
                        P[i][H - 1] = Bvec[i];
                    P[n1 - 1][H - 1] = VR;
                    strcpy(action, "LO-REDUCTION   ");
                }
                
                for (i = 0; i < n; i++)
                    Bvec[i] = (1 - bet) * P[i][H - 1] + bet * P[i][C - 1];
                f = fminfn->getFunVal(n, Bvec);   
                //f = fminfn(n, Bvec, ex);
                if (!isFinite(f)) f = big;
                funcount++;
                
                if (f < P[n1 - 1][H - 1]) {
                    for (i = 0; i < n; i++)
                        P[i][H - 1] = Bvec[i];
                    P[n1 - 1][H - 1] = f;
                } else {
                    if (VR >= VH) {
                        strcpy(action, "SHRINK         ");
                        calcvert = true;
                        size = 0.0;
                        for (j = 0; j < n1; j++) {
                            if (j + 1 != L) {
                                for (i = 0; i < n; i++) {
                                    P[i][j] = bet * (P[i][j] - P[i][L - 1])
                                    + P[i][L - 1];
                                    size += fabs(P[i][j] - P[i][L - 1]);
                                }
                            }
                        }
                        if (size < oldsize) {
                            oldsize = size;
                        } else {
                            if (trace)
                                cout << "Polytope size measure not decreased in shrink\n";
                            *fail = 10;
                            break;
                        }
                    }
                }
            }
            
        } while (funcount <= maxit);
        
    }
    
    if (trace) {
        cout << "Exiting from Nelder Mead minimizer\n";
        cout << funcount << " function evaluations used\n";
    }
    *Fmin = P[n1 - 1][L - 1];
    for (i = 0; i < n; i++) X[i] = P[i][L - 1];
    if (funcount > maxit) *fail = 1;
    *fncount = funcount;
}

void Optimizer::performNMmin(Matrix<double> initParamsMat, int dimension, Hawkes &h, const SubSpikeStruct& sOpt, int numParams, double* X, double abstol, double intol, double alpha, double bet, double gamm, int trace, int *fncount, int maxit, double *Fmin, int *fail){
    
    h.setData(sOpt);
    h.setParams(initParamsMat);
    h.setCurrentDimension(dimension);
    h.setTrace(trace); //Don't trace
    
    OptimFn * fptr = &h;
    
    double* Bvec; Bvec = new double[numParams]; // initial parameters
    for(int i = 0; i < numParams; i++){
        Bvec[i] = initParamsMat.Get_Element(dimension, i);
    }
    
    //cout << "n is " << numParams << " Bvec[0] = " << Bvec[0] << " Bvec[1] = " << Bvec[1] << endl;
    nmmin(numParams, Bvec,X, Fmin, fptr, fail, abstol, intol, 
            alpha,  bet,  gamm,  trace,
            fncount,  maxit);
    
}



/*  BFGS variable-metric method, based on Pascal code
 in J.C. Nash, `Compact Numerical Methods for Computers', 2nd edition,
 converted by p2c then re-crafted by B.D. Ripley */

void
vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
      int maxit, int trace, int *mask,
      double abstol, double reltol, int nREPORT, void *ex,
      int *fncount, int *grcount, int *fail)
{
    Rboolean accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l;
    
    if (maxit <= 0) {
        *fail = 0;
        *Fmin = fminfn(n0, b, ex);
        *fncount = *grcount = 0;
        return;
    }
    
    if (nREPORT <= 0)
        error(_("REPORT must be > 0 (method = \"BFGS\")"));
    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;
    g = vect(n0);
    t = vect(n);
    X = vect(n);
    c = vect(n);
    B = Lmatrix(n);
    f = fminfn(n0, b, ex);
    if (!R_FINITE(f))
        error(_("initial value in 'vmmin' is not finite"));
    if (trace) Rprintf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    fmingr(n0, b, g, ex);
    iter++;
    ilast = gradcount;
    
    do {
        if (ilast == gradcount) {
            for (i = 0; i < n; i++) {
                for (j = 0; j < i; j++) B[i][j] = 0.0;
                B[i][i] = 1.0;
            }
        }
        for (i = 0; i < n; i++) {
            X[i] = b[l[i]];
            c[i] = g[l[i]];
        }
        gradproj = 0.0;
        for (i = 0; i < n; i++) {
            s = 0.0;
            for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
            for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
            t[i] = s;
            gradproj += s * g[l[i]];
        }
        
        if (gradproj < 0.0) {	/* search direction is downhill */
            steplength = 1.0;
            accpoint = FALSE;
            do {
                count = 0;
                for (i = 0; i < n; i++) {
                    b[l[i]] = X[i] + steplength * t[i];
                    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
                        count++;
                }
                if (count < n) {
                    f = fminfn(n0, b, ex);
                    funcount++;
                    accpoint = R_FINITE(f) &&
                    (f <= *Fmin + gradproj * steplength * acctol);
                    if (!accpoint) {
                        steplength *= stepredn;
                    }
                }
            } while (!(count == n || accpoint));
            enough = (f > abstol) &&
            fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);
            /* stop if value if small or if relative change is low */
            if (!enough) {
                count = n;
                *Fmin = f;
            }
            if (count < n) {/* making progress */
                *Fmin = f;
                fmingr(n0, b, g, ex);
                gradcount++;
                iter++;
                D1 = 0.0;
                for (i = 0; i < n; i++) {
                    t[i] = steplength * t[i];
                    c[i] = g[l[i]] - c[i];
                    D1 += t[i] * c[i];
                }
                if (D1 > 0) {
                    D2 = 0.0;
                    for (i = 0; i < n; i++) {
                        s = 0.0;
                        for (j = 0; j <= i; j++)
                            s += B[i][j] * c[j];
                        for (j = i + 1; j < n; j++)
                            s += B[j][i] * c[j];
                        X[i] = s;
                        D2 += s * c[i];
                    }
                    D2 = 1.0 + D2 / D1;
                    for (i = 0; i < n; i++) {
                        for (j = 0; j <= i; j++)
                            B[i][j] += (D2 * t[i] * t[j]
                                        - X[i] * t[j] - t[i] * X[j]) / D1;
                    }
                } else {	/* D1 < 0 */
                    ilast = gradcount;
                }
            } else {	/* no progress */
                if (ilast < gradcount) {
                    count = 0;
                    ilast = gradcount;
                }
            }
        } else {		/* uphill search */
            count = 0;
            if (ilast == gradcount) count = n;
            else ilast = gradcount;
            /* Resets unless has just been reset */
        }
        if (trace && (iter % nREPORT == 0))
            Rprintf("iter%4d value %f\n", iter, f);
        if (iter >= maxit) break;
        if (gradcount - ilast > 2 * n)
            ilast = gradcount;	/* periodic restart */
    } while (count != n || ilast != gradcount);
    if (trace) {
        Rprintf("final  value %f \n", *Fmin);
        if (iter < maxit) Rprintf("converged\n");
        else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
}


#endif
