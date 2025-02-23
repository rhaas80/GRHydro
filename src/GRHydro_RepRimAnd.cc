#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <cmath>
#include <iostream>

#include "reprimand/con2prim_imhd.h"
#include "reprimand/eos_idealgas.h"
#include "reprimand/eos_idealgas.h"
#include "reprimand/eos_barotr_poly.h"

extern "C"
void CCTK_FCALL CCTK_FNAME(GRHydro_RPR_Con2Prim_pt)(
   CCTK_REAL *dens_in,  CCTK_REAL *scon1_in,
   CCTK_REAL *scon2_in,  CCTK_REAL *scon3_in,
   CCTK_REAL *tau_in, CCTK_REAL *g11_in,
   CCTK_REAL *g12_in, CCTK_REAL *g13_in,
   CCTK_REAL *g22_in, CCTK_REAL *g23_in,
   CCTK_REAL *g33_in, CCTK_REAL *rho,
  CCTK_REAL *eps, CCTK_REAL *press, CCTK_REAL *vel1,
  CCTK_REAL *vel2, CCTK_REAL *vel3, CCTK_REAL *error, CCTK_REAL *w_lorentz
) {
  DECLARE_CCTK_PARAMETERS

  using namespace EOS_Toolkit;

  // avoid possible hickups if real_t != CCTK_REAL exactly
  
  real_t g11 = *g11_in;
  real_t g12 = *g12_in;
  real_t g13 = *g13_in;
  real_t g22 = *g22_in;
  real_t g23 = *g23_in;
  real_t g33 = *g33_in;

  

  //Order is xx, yx, yy, zx, zy, zz
  sm_symt3l metric(sm_symt3l{g11, g12, g22, g13, g23, g33});
  sm_metric3 g(metric);

  // Reprimand's "evolved" type includes the sqrt(detg), so we do not have to
  // divide it out here
  // https://wokast.github.io/RePrimAnd/notation.html
  real_t dens = *dens_in ;
  real_t scon1 = *scon1_in ;
  real_t scon2 = *scon2_in ;
  real_t scon3 = *scon3_in ;
  real_t tau = *tau_in ;

  // Faber 2007 fix for max s2
  real_t s2 = (    g.lo(0,0) * scon1*scon1 +
               2 * g.lo(1,0) * scon2*scon1 +
               2 * g.lo(2,0) * scon3*scon1 +
                   g.lo(1,1) * scon2*scon2 +
               2 * g.lo(1,2) * scon2*scon3 +
                   g.lo(2,2) * scon3*scon3) * g.det;
  if(s2 >= tau*(tau + 2*dens)) {
    real_t corrfac = sqrt(0.98 * tau*(tau + 2*dens)/s2);
    *scon1_in = (scon1 *= corrfac);
    *scon2_in = (scon2 *= corrfac);
    *scon3_in = (scon3 *= corrfac);
  }

  real_t max_eps = 11.;
  real_t max_rho = 1e6;
  real_t adiab_ind = 1./(gl_gamma - 1.); // n = 1/(Gamma - 1)
  auto eos = make_eos_idealgas(adiab_ind, max_eps, max_rho);

  //Set up atmosphere
  // from docs (https://wokast.github.io/RePrimAnd/eos_barotr_available.html)
  // K = rho_p^(1-Gamma) => rho_p = K^(1/(1-Gamma)) = K^-n
  real_t adiab_ind_atmo = 1./(poly_gamma - 1.); // n = 1/(Gamma - 1)
  real_t K_atmo = poly_k; // poly_k from parfile
  real_t rmd_p_atmo = pow(K_atmo,-adiab_ind);
  auto eos_atmo = make_eos_barotr_poly(adiab_ind_atmo, rmd_p_atmo, max_rho);
  real_t atmo_rho = rho_abs_min; // the parameter, GRHydro_rho_min is a grid scalar nad harder to access
  real_t atmo_eps = eos_atmo.at_rho(atmo_rho).eps();
  real_t atmo_ye = 0.5;
  real_t atmo_cut = atmo_rho * (1.+GRHydro_atmo_tolerance);
  real_t atmo_p = eos_atmo.at_rho(atmo_rho).press();

  atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

  //Primitive recovery parameters det
  bool  ye_lenient = false;
  int max_iter = 30;
  real_t c2p_acc = 1e-8;
  real_t max_b = 10.;
  real_t max_z = 1e3;

  //Primitive recovery parameters 
  real_t rho_strict = 100 * atmo_rho;
  //Get a recovery function
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, 
                     atmo, c2p_acc, max_iter);
  real_t trye = 0.5*dens;
  
  //collect
  cons_vars_mhd evolved{dens, tau, trye, 
                        {scon1,scon2,scon3}, {0.,0.,0.}};    
  prim_vars_mhd primitives;
  con2prim_mhd::report rep;
  
  //recover
  cv2pv(primitives, evolved, g, rep);

  if(rep.failed()) {
    *error = 1;
    *rho = primitives.rho;
    *eps = primitives.eps;
    *press = primitives.press;
    *vel1 = primitives.vel(0);
    *vel2 = primitives.vel(1);
    *vel3 = primitives.vel(2);
    *w_lorentz = primitives.w_lor;
    *scon1_in = evolved.scon(0);
    *scon2_in = evolved.scon(1);
    *scon3_in = evolved.scon(2);
    *dens_in = evolved.dens;
    *tau_in = evolved.tau;
  } else {
    *error = 0;
    static FILE * fh;
    if (fh == NULL) {
      fh = fopen("failures.txt", "wb");
      setvbuf(fh, NULL, _IOLBF, BUFSIZ);
    }
    if(primitives.eps == 0.) {
      fprintf(fh, "%.18e %.18e %.18e %.18e %.18e ", *dens_in, *scon1_in, *scon2_in, *scon3_in, *tau_in);
      fprintf(fh, "%.18e %.18e %.18e %.18e %.18e %.18e ", g11, g12, g13, g22, g23, g33);
      fprintf(fh, "%d\n", int(rep.status));
      exit(1);
    }
    //write back primitive vars
    *rho = primitives.rho;
    *eps = primitives.eps;
    *press = primitives.press;
    *vel1 = primitives.vel(0);
    *vel2 = primitives.vel(1);
    *vel3 = primitives.vel(2);
    *w_lorentz = primitives.w_lor;
    if (rep.adjust_cons) {
      //write back corrected evolved vars to grid here
      *scon1_in = evolved.scon(0);
      *scon2_in = evolved.scon(1);
      *scon3_in = evolved.scon(2);
      *dens_in = evolved.dens;
      *tau_in = evolved.tau;
    }
  }

  return;
}
