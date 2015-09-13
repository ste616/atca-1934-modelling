#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include "cpgplot.h"

#define BUFSIZE 1024
#define NEXTRA 2
#define NFIT_SAULT 4
#define NFIT_REYNOLDS 4
#define NFIT_STEVENS 3
#define NFIT_OLD 4

void minmax(int narr, float *arr, float *minarr, float *maxarr);
void linfit_order(int order, int n, float *x, float *y, float *w, float **params);

void fitstring(float *c, int n, char *s) {
  int i;
  char tl[BUFSIZE];
  for (i=0; i<n; i++) {
    if (i == 0) {
      sprintf(tl, " %.4f", c[i]);
      strcat(s, tl);
    } else {
      sprintf(tl, " + %.4f x\\u%d\\d", c[i], i);
      strcat(s, tl);
    }
  }

}

int main(int argc, char *argv[]) {
  float *x=NULL,*y=NULL,minx,maxx,miny,maxy,cx, *oparams=NULL,*nparams=NULL;
  float *rx=NULL,*ry=NULL,*w=NULL,*wparams=NULL,*wx=NULL,*wy=NULL,*ww=NULL;
  float *y_sault_fit=NULL, *x_fit=NULL, *y_new_fit=NULL, *y_reynolds_fit=NULL;
  float *y_stevens_fit=NULL,*y_whole_fit=NULL,*y_old_fit=NULL;
  int i,j,n=0,n_fit=100,new_fit_order=2,whole_fit_order=5,nr=0,nw=0;
  /* float extra_x[NEXTRA]={ 93, 95 }, extra_y[NEXTRA] = { 0.1223, 0.1168 }; */
  float extra_x[NEXTRA]={ 93, 95 }, extra_y[NEXTRA] = { 0.1116, 0.1056 };
  float extra_u[NEXTRA]={ 0.01356, 0.01399 };
  float fitp_sault[NFIT_SAULT]={ -202.6259, 149.7321, -36.4943, 2.9372 };
  float fitp_reynolds[NFIT_REYNOLDS]={ -30.7667, 26.4908, -7.0977, 0.605334 };
  float fitp_stevens[NFIT_STEVENS]={ -1.237160, 2.005317, -0.400622 };
  float fitp_old[NFIT_OLD]={ -23.839, 19.569, -4.8168, 0.35836 };
  float *ratio_reynolds_fit=NULL,*ratio_stevens_fit=NULL,*ratio_sault_fit=NULL;
  float *ratio_new_fit=NULL;
  float vpx1, vpx2, vpy1, vpy2, vpy3, lx, ly, dly;
  char fitlabel[BUFSIZE];

  /* Generate the cm fit points. */
  for (cx=1.0; cx<10.0; cx+=0.1) {
    nr++;
    rx = realloc(rx, nr * sizeof(float));
    ry = realloc(ry, nr * sizeof(float));
    rx[nr-1] = log10f(cx * 1000);
    ry[nr-1] = 0.0;
    for (i=0; i<NFIT_REYNOLDS; i++) {
      ry[nr-1] += fitp_reynolds[i] * powf(rx[n-1], (float)i);
    }
  }

  /* Generate the 15mm fit points. */
  for (cx=10.0; cx<=24.0; cx+=0.128) {
    n++;
    x = realloc(x, n * sizeof(float));
    y = realloc(y, n * sizeof(float));
    w = realloc(w, n * sizeof(float));
    x[n-1] = log10f(cx * 1000);
    y[n-1] = 0.0;
    for (i=0; i<NFIT_REYNOLDS; i++) {
      y[n-1] += fitp_sault[i] * powf(x[n-1], (float)i);
    }
    w[n-1] = 1.0/0.1;
  }

  /* Do the fit. */
  linfit_order(NFIT_SAULT, n, x, y, w, &oparams);
  for (i=0; i<NFIT_SAULT; i++) {
    printf("i = %d c[i] = %.4f\n", i, oparams[i]);
  }

  /* Add the 3mm flux points. */
  for (i=0; i<NEXTRA; i++) {
    n++;
    x = realloc(x, n * sizeof(float));
    y = realloc(y, n * sizeof(float));
    w = realloc(w, n * sizeof(float));
    x[n-1] = log10f(extra_x[i] * 1000);
    y[n-1] = log10f(extra_y[i]);
    w[n-1] = 1/extra_u[i];
  }

  /* Do another fit. */
  linfit_order(new_fit_order, n, x, y, w, &nparams);
  for (i=0; i<new_fit_order; i++) {
    printf("i = %d nc[i] = %.4f\n", i, nparams[i]);
  }

  /* Generate the whole range fit points. */
  minx=log10f(900);
  maxx=log10f(100000);
  miny=-2;
  maxy=log10f(20);
  x_fit = malloc(n_fit * sizeof(float));
  for (i=0; i<n_fit; i++) {
    x_fit[i] = minx + i * ((maxx - minx)/(float)n_fit);
    nw++;
    wx = realloc(wx, nw * sizeof(float));
    wy = realloc(wy, nw * sizeof(float));
    ww = realloc(ww, nw * sizeof(float));
    wx[nw-1] = x_fit[i];
    wy[nw-1] = 0.0;
    ww[nw-1] = 1;
    if (x_fit[i] < log10f(11143)) {
      /* Use the Reynolds fit. */
      for (j=0; j<NFIT_REYNOLDS; j++) {
	wy[nw-1] += fitp_reynolds[j] * powf(x_fit[i], (float)j);
      }
    } else {
      /* Use the new fit. */
      for (j=0; j<new_fit_order; j++) {
	wy[nw-1] += nparams[j] * powf(x_fit[i], (float)j);
      }
    }
  }

  /* Do a whole-range fit. */
  linfit_order(whole_fit_order, nw, wx, wy, ww, &wparams);
  for (i=0; i<whole_fit_order; i++) {
    printf("i = %d wc[i] = %.4f\n", i, wparams[i]);
  }

  // minmax(n, x, &minx, &maxx);
  // minmax(n, y, &miny, &maxy);

  y_sault_fit = malloc(n_fit * sizeof(float));
  y_reynolds_fit = malloc(n_fit * sizeof(float));
  y_stevens_fit = malloc(n_fit * sizeof(float));
  y_new_fit = malloc(n_fit * sizeof(float));
  y_whole_fit = malloc(n_fit * sizeof(float));
  y_old_fit = malloc(n_fit * sizeof(float));
  ratio_reynolds_fit = malloc(n_fit * sizeof(float));
  ratio_stevens_fit = malloc(n_fit * sizeof(float));
  ratio_sault_fit = malloc(n_fit * sizeof(float));
  ratio_new_fit = malloc(n_fit * sizeof(float));
  /* minx=log10f(50); */
  minx=log10f(1000);
  /* maxx=log10f(500000); */
  maxx=log10f(110000);
  for (i=0; i<n_fit; i++) {
    x_fit[i] = minx + i * ((maxx - minx)/(float)n_fit);
    y_sault_fit[i] = 0.0;
    y_reynolds_fit[i] = 0.0;
    y_new_fit[i] = 0.0;
    y_stevens_fit[i] = 0.0;
    y_whole_fit[i] = 0.0;
    y_old_fit[i] = 0.0;
    for (j=0; j<NFIT_SAULT; j++) {
      y_sault_fit[i] += fitp_sault[j] * powf(x_fit[i], (float)j);
    }
    for (j=0; j<NFIT_REYNOLDS; j++) {
      y_reynolds_fit[i] += fitp_reynolds[j] * powf(x_fit[i], (float)j);
    }
    for (j=0; j<NFIT_STEVENS; j++) {
      y_stevens_fit[i] += fitp_stevens[j] * powf(x_fit[i], (float)j);
    }
    for (j=0; j<new_fit_order; j++) {    
      y_new_fit[i] += nparams[j] * powf(x_fit[i], (float)j);
    }
    for (j=0; j<whole_fit_order; j++) {
      y_whole_fit[i] += wparams[j] * powf(x_fit[i], (float)j);
    }
    for (j=0; j<NFIT_OLD; j++) {
      y_old_fit[i] += fitp_old[j] * powf(x_fit[i], (float)j);
    }
    ratio_reynolds_fit[i] = powf(10, (y_reynolds_fit[i] - y_whole_fit[i]));
    ratio_stevens_fit[i] = powf(10, (y_stevens_fit[i] - y_whole_fit[i]));
    ratio_sault_fit[i] = powf(10, (y_sault_fit[i] - y_whole_fit[i]));
    ratio_new_fit[i] = powf(10, (y_new_fit[i] - y_whole_fit[i]));
  }

  /* cpgopen("11/xs"); */
  cpgopen("1934-638_models.ps/cps");
  /* cpgopen("1934-638_models.png/png"); */
  cpgqvp(0, &vpx1, &vpx2, &vpy1, &vpy2);
  vpy3 = vpy1 + (vpy2 - vpy1) / 5.0;
  /* cpgsvp(vpx1, vpx2, vpy3, vpy2); */
  cpgswin(minx, maxx, miny, maxy);
  lx = minx + (maxx - minx) / 9.0;
  ly = miny + (maxy - miny) / 3.0;
  dly = (maxy - miny) / 20.0;
  cpgsch(1.0);
  cpgbox("BCLNTS",0,0,"BCLNTS",0,0);
  cpglab("Frequency (MHz)", "Flux Density (Jy)", "1934-638 Model Comparison");
  cpgsch(0.8);
  cpgpt(n, x, y, 4);
  /* cpgpt(nw, wx, wy, 4); */
  cpgsci(2);
  /* cpgpt(nr, rx, ry, 4); */
  cpgline(n_fit, x_fit, y_sault_fit);
  strcpy(fitlabel, "Sault: ");
  fitstring(fitp_sault, NFIT_SAULT, fitlabel);
  cpgtext(lx, ly, fitlabel);
  cpgsci(3);
  cpgline(n_fit, x_fit, y_new_fit);
  strcpy(fitlabel, "Stevens (linear): ");
  fitstring(nparams, new_fit_order, fitlabel);
  ly -= dly;
  cpgtext(lx, ly, fitlabel);
  cpgsci(4);
  cpgline(n_fit, x_fit, y_reynolds_fit);
  strcpy(fitlabel, "Reynolds: ");
  fitstring(fitp_reynolds, NFIT_REYNOLDS, fitlabel);
  ly -= dly;
  cpgtext(lx, ly, fitlabel);
  cpgsci(5);
  cpgline(n_fit, x_fit, y_stevens_fit);
  strcpy(fitlabel, "Stevens (Miriad): ");
  fitstring(fitp_stevens, NFIT_STEVENS, fitlabel);
  ly -= dly;
  cpgtext(lx, ly, fitlabel);
  cpgsci(6);
  cpgline(n_fit, x_fit, y_old_fit);
  strcpy(fitlabel, "Pre-1994: ");
  fitstring(fitp_old, NFIT_OLD, fitlabel);
  ly -= dly;
  cpgtext(lx, ly, fitlabel);
  /* cpgsci(6); */
  /* cpgline(n_fit, x_fit, y_whole_fit); */
  /* strcpy(fitlabel, "Stevens (New): "); */
  /* fitstring(wparams, whole_fit_order, fitlabel); */
  /* ly -= dly; */
  /* cpgtext(lx, ly, fitlabel); */
  /* cpgsvp(vpx1, vpx2, vpy1, vpy3); */
  /* cpgsci(1); */
  /* cpgswin(minx, maxx, 0.9, 1.1); */
  /* cpgsch(1.0); */
  /* cpgbox("BCLNTS",0,0,"BCMTS",0,0); */
  /* cpglab("Frequency (MHz)", "Model Ratio", ""); */
  /* cpgsci(2); */
  /* cpgline(n_fit, x_fit, ratio_sault_fit); */
  /* cpgsci(3); */
  /* cpgline(n_fit, x_fit, ratio_new_fit); */
  /* cpgsci(4); */
  /* cpgline(n_fit, x_fit, ratio_reynolds_fit); */
  /* cpgsci(5); */
  /* cpgline(n_fit, x_fit, ratio_stevens_fit); */
  cpgclos();

  exit(0);
}

void minmax(int narr, float *arr, float *minarr, float *maxarr) {
  int i;
  
  *minarr = arr[0];
  *maxarr = arr[1];
  for (i=1; i<narr; i++) {
    *minarr = (arr[i] < *minarr) ? arr[i] : *minarr;
    *maxarr = (arr[i] > *maxarr) ? arr[i] : *maxarr;
  }
}

void linfit_order(int order, int n, float *x, float *y, float *w, float **params) {
  int i,j;
  double chisq;
  gsl_matrix *X, *cov;
  gsl_vector *Y, *W, *c;
  gsl_multifit_linear_workspace *work;

  X = gsl_matrix_alloc(n, order);
  Y = gsl_vector_alloc(n);
  W = gsl_vector_alloc(n);
  
  c = gsl_vector_alloc(order);
  cov = gsl_matrix_alloc(order, order);

  for (i=0; i<n; i++) {
    for (j=0; j<order; j++) {
      gsl_matrix_set(X, i, j, pow((double)x[i], (double)j));
    }
    gsl_vector_set(Y, i, (double)y[i]);
    gsl_vector_set(W, i, (double)w[i]);
  }

  work = gsl_multifit_linear_alloc(n, order);
  gsl_multifit_wlinear(X, W, Y, c, cov, &chisq, work);

  *params = malloc(order * sizeof(float));
  for (j=0; j<order; j++) {
    (*params)[j] = (float)(gsl_vector_get(c, j));
  }

  gsl_matrix_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(W);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
}
