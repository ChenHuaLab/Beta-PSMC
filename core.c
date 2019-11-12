#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h> // add by liujf 2019-07-17
#include <float.h> // add by liujf 2019-07-17
#include "psmc.h"

int n_dis; // add by liujf 2019-07-17

// add by liujf 2019-07-17
void error2 (char * message)
{ printf("\nError: %s.\n", message); exit(-1); }

// add by liujf 2019-07-17
long factorial (int n)
{
   long f=1, i;
   if (n>11) error2("n>10 in factorial");
   for (i=2; i<=(long)n; i++) f *= i;
   return (f);
}

// add by liujf 2019-07-17
double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x;

   if((double)nx==x && nx>=0 && nx<=11)
      lng = log((double)factorial(nx-1));
   else {
      if(x<=0) {
         printf("LnGamma(%.6f) not implemented", x);
         if((int)x-x==0) { puts("lnGamma undefined"); return(-1); }
         for (fneg=1; x<0; x++) fneg /= x;
         if(fneg<0) 
            error2("strange!! check lngamma");
         fneg = log(fneg);
      }
      if (x<7) {
         f = 1;
         z = x-1;
         while (++z<7)  
            f *= z;
         x = z;   
         f = -log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}

// add by liujf 2019-07-17
double CDFBeta (double x, double pin, double qin, double lnbeta)
{
/* Returns distribution function of the standard form of the beta distribution, 
   that is, the incomplete beta ratio I_x(p,q).

   This is also known as the incomplete beta function ratio I_x(p, q)

   lnbeta is log of the complete beta function; provide it if known,
   and otherwise use 0.

   This is called from QuantileBeta() in a root-finding loop.

    This routine is a translation into C of a Fortran subroutine
    by W. Fullerton of Los Alamos Scientific Laboratory.
    Bosten and Battiste (1974).
    Remark on Algorithm 179, CACM 17, p153, (1974).
*/
   double ans, c, finsum, p, ps, p1, q, term, xb, xi, y, small=1e-15;
   int n, i, ib;
   static double eps = 0, alneps = 0, sml = 0, alnsml = 0;

   if(x<small)        return 0;
   else if(x>1-small) return 1;
   if(pin<=0 || qin<=0)  { 
      printf("p=%.4f q=%.4f: parameter outside range in CDFBeta",pin,qin); 
      return (-1); 
   }

   if (eps == 0) {/* initialize machine constants ONCE */
      eps = pow((double)FLT_RADIX, -(double)DBL_MANT_DIG);
      alneps = log(eps);
      sml = DBL_MIN;
      alnsml = log(sml);
   }
   y = x;  p = pin;  q = qin;

    /* swap tails if x is greater than the mean */
   if (p / (p + q) < x) {
      y = 1 - y;
      p = qin;
      q = pin;
   }

   if(lnbeta==0) lnbeta = LnGamma(p) + LnGamma(q) - LnGamma(p+q);

   if ((p + q) * y / (p + 1) < eps) {  /* tail approximation */
      ans = 0;
      xb = p * log(max2(y, sml)) - log(p) - lnbeta;
      if (xb > alnsml && y != 0)
         ans = exp(xb);
      if (y != x || p != pin)
      ans = 1 - ans;
   }
   else {
      /* evaluate the infinite sum first.  term will equal */
      /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
      ps = q - floor(q);
      if (ps == 0)
         ps = 1;

      xb=LnGamma(ps)+LnGamma(p)-LnGamma(ps+p);
      xb = p * log(y) - xb - log(p);

      ans = 0;
      if (xb >= alnsml) {
         ans = exp(xb);
         term = ans * p;
         if (ps != 1) {
            n = (int)max2(alneps/log(y), 4.0);
         for(i=1 ; i<= n ; i++) {
            xi = i;
            term = term * (xi - ps) * y / xi;
            ans = ans + term / (p + xi);
         }
      }
   }

   /* evaluate the finite sum. */
   if (q > 1) {
      xb = p * log(y) + q * log(1 - y) - lnbeta - log(q);
      ib = (int) (xb/alnsml);  if(ib<0) ib=0;
      term = exp(xb - ib * alnsml);
      c = 1 / (1 - y);
      p1 = q * c / (p + q - 1);

      finsum = 0;
      n = (int) q;
      if (q == (double)n)
         n = n - 1;
         for(i=1 ; i<=n ; i++) {
            if (p1 <= 1 && term / eps <= finsum)
               break;
            xi = i;
            term = (q - xi + 1) * c * term / (p + q - xi);
            if (term > 1) {
               ib = ib - 1;
               term = term * sml;
            }
            if (ib == 0)
               finsum = finsum + term;
         }
         ans = ans + finsum;
      }
      if (y != x || p != pin)
         ans = 1 - ans;
      if(ans>1) ans=1;
      if(ans<0) ans=0;
   }
   return ans;
}

// add by liujf 2019-07-17
FLOAT psmc_beta(FLOAT p, FLOAT q, int j, int k)
{
	FLOAT a0, a1, y;
	a0 = (FLOAT) j/k;
	a1 = (FLOAT) (j+1)/k;
	y = CDFBeta(a1, p, q, 0) - CDFBeta(a0, p, q, 0);
	return y;
}

void psmc_update_intv(int n, FLOAT t[], FLOAT max_t, FLOAT alpha, FLOAT *inp_ti)
{
	int k;
	FLOAT *t_tmp, s = 0.0; // add by liujf 2019-07-12
	int j; // add by liujf 2019-07-12
	t_tmp = (FLOAT*)malloc(sizeof(FLOAT) * (n + 2)); // \t_tmp_k add by liujf 2019-07-12
	if (inp_ti == 0) {
		FLOAT beta;
		beta = log(1.0 + max_t / alpha) / (n / n_dis); // beta controls the sizes of intervals change n to n/n_dis by liujf 2019-07-18
		//for (k = 0; k < n; ++k) // comment by liujf 2019-07-12
			//t[k] = alpha * (exp(beta * k) - 1); // comment by liujf 2019-07-12
		// change the way to compute the value of pd->t by liujf 2019-07-12, 2019-07-17
		if (n_dis > 1) {
			for (k = 0; k < (n / n_dis); ++k) {
				t_tmp[k] = alpha * (exp(beta * k) - 1);
				if (k == (n / n_dis) - 1) {
					s = max_t - t_tmp[k];
				} else {
					s = alpha * (exp(beta * (k + 1)) - 1) - t_tmp[k];
				}
				for (j = 0; j < n_dis; ++j) {
					t[k*n_dis + j] = t_tmp[k] + s * j / n_dis;
				}
			}
		} else {
			for (k = 0; k < n; ++k) {
				t[k] = alpha * (exp(beta * k) - 1);
			}
		}
		
		t[n] = max_t; t[n+1] = PSMC_T_INF; // the infinity: exp(PSMC_T_INF) > 1e310 = inf
	} else {
		memcpy(t, inp_ti, (n+1) * sizeof(FLOAT));
		t[n+1] = PSMC_T_INF;
	}
	free(t_tmp); // add by liujf 2019-07-12
}

psmc_data_t *psmc_new_data(psmc_par_t *pp)
{
	psmc_data_t *pd;
	int k, n = pp->n;
	pd = (psmc_data_t*)calloc(1, sizeof(psmc_data_t));
	//pd->n_params = pp->n_free + PSMC_N_PARAMS + ((pp->flag & PSMC_F_DIVERG)? 1 : 0); // one addition parameter for the divergence model // comment by liujf 2019-07-10
	// change the value of pd->n_params by liujf 2019-07-10
	if (pp->n_discrete > 1) {
		pd->n_params = (pp->n_free) * 3 + PSMC_N_PARAMS + ((pp->flag & PSMC_F_DIVERG)? 1 : 0); // one addition parameter for the divergence model
	} else {
		pd->n_params = pp->n_free + PSMC_N_PARAMS + ((pp->flag & PSMC_F_DIVERG)? 1 : 0); // one addition parameter for the divergence model
	}
	pd->hp = hmm_new_par(2, n + 1);
	// initialize
	pd->sigma = (FLOAT*)calloc(n+1, sizeof(FLOAT));
	pd->post_sigma = (FLOAT*)calloc(n+1, sizeof(FLOAT));
	pd->t = (FLOAT*)malloc(sizeof(FLOAT) * (n + 2)); // $t_0,\ldots,t_{n+1}$
	pd->params = (FLOAT*)calloc(pd->n_params, sizeof(FLOAT)); // free lambdas + theta + rho
	// initialize t[] and params[]
	if (pp->inp_pa) { // pameters are loaded from a file
		memcpy(pd->params, pp->inp_pa, sizeof(FLOAT) * pd->n_params); // FIXME: not working for the divergence model
	} else {
		FLOAT theta;
		// initialize psmc_data_t::params[]
		theta = -log(1.0 - (FLOAT)pp->sum_n / pp->sum_L);
		pd->params[0] = theta; // \theta_0
		pd->params[1] = theta / pp->tr_ratio; // \rho_0
		pd->params[2] = pp->max_t;
		for (k = PSMC_N_PARAMS; k != pp->n_free + PSMC_N_PARAMS; ++k) {
			// \lambda_k
			pd->params[k] = 1.0 + (drand48() * 2.0 - 1.0) * pp->ran_init;
			if (pd->params[k] < 0.1) pd->params[k] = 0.1;
			// two shape parameters for Beta distribution by liujf 2019-07-10
			if (pp->n_discrete > 1) {
				pd->params[k + pp->n_free] = 1.0 + (drand48() * 2.0 - 1.0) * pp->ran_init;
				pd->params[k + (pp->n_free) * 2] = 1.0 + (drand48() * 2.0 - 1.0) * pp->ran_init;
			}
		}
		if (pp->flag & PSMC_F_DIVERG) pd->params[pd->n_params - 1] = pp->dt0;
	}
	psmc_update_hmm(pp, pd);
	return pd;
}
void psmc_delete_data(psmc_data_t *pd)
{
	if (pd == 0) return;
	free(pd->sigma); free(pd->post_sigma);
	free(pd->t); free(pd->params);
	hmm_delete_par(pd->hp);
	free(pd);
}
void psmc_update_hmm(const psmc_par_t *pp, psmc_data_t *pd) // calculate the a_{kl} and e_k(b)
{
	FLOAT *q, tmp, sum_t, max_t, *alpha, *beta, *q_aux, *lambda, theta, rho, *t, *tau, dt = 0;
	hmm_par_t *hp = pd->hp;
	int k, l, n = pp->n;
	FLOAT *lambda_tmp; // add by liujf 2019-07-10
	int j; //add by liujf 2019-07-10
	t = pd->t;
	lambda_tmp = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \lambda_k add by liujf 2019-07-12
	lambda = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \lambda_k
	alpha = (FLOAT*)malloc(sizeof(FLOAT) * (n + 2)); // \alpha_k
	beta = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \beta_k
	q_aux = (FLOAT*)malloc(sizeof(FLOAT) * n); // for acceleration
	q = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // q_{kl}
	tau = (FLOAT*)malloc(sizeof(FLOAT) * (n + 1)); // \tau_k
	// calculate population parameters: \theta_0, \rho_0, \lambda_k and max_t
	theta = pd->params[0]; rho = pd->params[1]; max_t = pd->params[2];
	//for (k = 0; k <= n; ++k) // comment by liujf 2019-07-10
		//lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS]; // comment by liujf 2019-07-10
	//change the way to compute the value of lambda_k by liujf 2019-07-10
	if (pp->n_discrete > 1) {
		for (k = 0; k <= (n / pp->n_discrete); ++k) {
			lambda_tmp[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
			if (k == (n / pp->n_discrete)) {
				lambda[k*(pp->n_discrete)] = lambda_tmp[k];
			} else {
				for (j = 0; j < pp->n_discrete; ++j) {
					lambda[k*(pp->n_discrete) + j] = lambda_tmp[k] * (pp->n_discrete) * psmc_beta(pd->params[pp->par_map[k] + PSMC_N_PARAMS + pp->n_free], pd->params[pp->par_map[k] + PSMC_N_PARAMS + (pp->n_free)*2], j, pp->n_discrete);
				}
			}
		}		
	} else {
		for (k = 0; k <= n; ++k) {
			lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
		}
	}
	n_dis = pp->n_discrete; // add by liujf 2019-07-17
	psmc_update_intv(pp->n, pd->t, max_t, pp->alpha, pp->inp_ti);
	// set the divergence time parameter if necessary
	if (pp->flag & PSMC_F_DIVERG) {
		dt = pd->params[pd->n_params - 1];
		if (dt < 0) dt = 0;
	}
	// calculate \tau_k
	for (k = 0; k <= n; ++k) tau[k] = t[k+1] - t[k];
	// calculate \alpha
	for (k = 1, alpha[0] = 1.0; k <= n; ++k)
		alpha[k] = alpha[k-1] * exp(-tau[k-1] / lambda[k-1]);
	alpha[k] = 0.0;
	// calculate \beta
	for (k = 1, beta[0] = 0.0; k <= n; ++k)
		beta[k] = beta[k-1] + lambda[k-1] * (1.0 / alpha[k] - 1.0 / alpha[k-1]);
	// calculate q_aux
	for (l = 0; l < n; ++l)
		q_aux[l] = (alpha[l] - alpha[l+1]) * (beta[l] - lambda[l] / alpha[l]) + tau[l];
	// calculate C_pi and C_sigma
	for (l = 0, pd->C_pi = 0.0; l <= n; ++l)
		pd->C_pi += lambda[l] * (alpha[l] - alpha[l+1]);
	pd->C_sigma = 1.0 / (pd->C_pi * rho) + 0.5;
	// calculate all the rest
	for (k = 0, sum_t = 0.0; k <= n; ++k) {
		FLOAT *aa, avg_t, ak1, lak, pik, cpik;
		ak1 = alpha[k] - alpha[k+1]; lak = lambda[k]; // just for convenient
		// calculate $\pi_k$, $\sigma_k$ and Lak
		cpik = ak1 * (sum_t + lak) - alpha[k+1] * tau[k];
		pik = cpik / pd->C_pi;
		pd->sigma[k] = (ak1 / (pd->C_pi * rho) + pik / 2.0) / pd->C_sigma;
		// calculate avg_t, the average time point where mutation happens
		avg_t = - log(1.0 - pik / (pd->C_sigma*pd->sigma[k])) / rho;
		if (isnan(avg_t) || avg_t < sum_t || avg_t > sum_t + tau[k]) // in case something bad happens
			avg_t = sum_t + (lak - tau[k] * alpha[k+1] / (alpha[k] - alpha[k+1]));
		// calculate q_{kl}
		tmp = ak1 / cpik;
		for (l = 0; l < k; ++l) q[l] = tmp * q_aux[l]; // q_{kl}, l<k
		q[l++] = (ak1 * ak1 * (beta[k] - lak/alpha[k]) + 2*lak*ak1 - 2*alpha[k+1]*tau[k]) / cpik; // q_{kk}
		if (k < n) {
			tmp = q_aux[k] / cpik;
			for (; l <= n; ++l) q[l] = (alpha[l] - alpha[l+1]) * tmp; // q_{kl}, l>k
		}
		// calculate p_{kl} and e_k(b)
		tmp = pik / (pd->C_sigma * pd->sigma[k]);
		for (aa = hp->a[k], l = 0; l <= n; ++l) aa[l] = tmp * q[l];
		aa[k] = tmp * q[k] + (1.0 - tmp);
		hp->a0[k] = pd->sigma[k];
		hp->e[0][k] = exp(-theta * (avg_t + dt));
		hp->e[1][k] = 1.0 - hp->e[0][k];
		// update sum_lt
		sum_t += tau[k];
		// for (l = 0, tmp = 0.0; l <= n; ++l) tmp += q[l]; fprintf(stderr, "%d\t%lf\n", k, tmp); // for testing only
	}
	// for (l = 0, tmp = 0.0; l <= n; ++l) tmp += hp->a0[l]; fprintf(stderr, "%lf\n", tmp); // for testing only
	// free
	free(q); free(alpha); free(beta); free(q_aux); free(lambda); free(tau);
	free(lambda_tmp); //add by liujf 2019-07-10
}
// compute the average time for each interval; a strip-down version of psmc_update_hmm()
void psmc_avg_t(const psmc_par_t *pp, const psmc_data_t *pd, double *avg_t)
{
	FLOAT sum_t, *alpha, *lambda, rho = pd->params[1], *tau, dt = 0;
	int k, n = pp->n;
	FLOAT *lambda_tmp; // add by liujf 2019-07-12
	int j; //add by liujf 2019-07-12
	lambda_tmp = (FLOAT*)alloca(sizeof(FLOAT) * (n + 1)); // \lambda_tmp_k	add by liujf 2019-07-12
	lambda = (FLOAT*)alloca(sizeof(FLOAT) * (n + 1)); // \lambda_k
	alpha = (FLOAT*)alloca(sizeof(FLOAT) * (n + 2)); // \alpha_k
	tau = (FLOAT*)alloca(sizeof(FLOAT) * (n + 1)); // \tau_k
	//for (k = 0; k <= n; ++k) // comment by liujf 2019-07-12
		//lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS]; // comment by liujf 2019-07-12
	//change the way to compute the value of lambda_k by liujf 2019-07-12
	if (pp->n_discrete > 1) {
		for (k = 0; k <= (n / pp->n_discrete); ++k) {
			lambda_tmp[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
			if (k == (n / pp->n_discrete)) {
				lambda[k*(pp->n_discrete)] = lambda_tmp[k];
			} else {
				for (j = 0; j < pp->n_discrete; ++j) {
					lambda[k*(pp->n_discrete) + j] = lambda_tmp[k] * (pp->n_discrete) * psmc_beta(pd->params[pp->par_map[k] + PSMC_N_PARAMS + pp->n_free], pd->params[pp->par_map[k] + PSMC_N_PARAMS + (pp->n_free)*2], j, pp->n_discrete);
				}
			}
		}		
	} else {
		for (k = 0; k <= n; ++k) {
			lambda[k] = pd->params[pp->par_map[k] + PSMC_N_PARAMS];
		}
	}
	if (pp->flag & PSMC_F_DIVERG) {
		dt = pd->params[pd->n_params - 1];
		if (dt < 0) dt = 0;
	}
	for (k = 0; k <= n; ++k) tau[k] = pd->t[k+1] - pd->t[k];
	for (k = 1, alpha[0] = 1.0; k <= n; ++k)
		alpha[k] = alpha[k-1] * exp(-tau[k-1] / lambda[k-1]);
	alpha[k] = 0.0;
	for (k = 0, sum_t = 0.0; k <= n; ++k) {
		FLOAT ak1, lak, pik;
		ak1 = alpha[k] - alpha[k+1]; lak = lambda[k]; // just for convenient
		pik = (ak1 * (sum_t + lak) - alpha[k+1] * tau[k]) / pd->C_pi;
		avg_t[k] = - log(1.0 - pik / (pd->C_sigma*pd->sigma[k])) / rho;
		if (isnan(avg_t[k]) || avg_t[k] < sum_t || avg_t[k] > sum_t + tau[k]) // in case something bad happens
			avg_t[k] = sum_t + (lak - tau[k] * alpha[k+1] / (alpha[k] - alpha[k+1]));
		avg_t[k] += dt;
		sum_t += tau[k];
	}
	free(lambda_tmp); //add by liujf 2019-07-12
}
