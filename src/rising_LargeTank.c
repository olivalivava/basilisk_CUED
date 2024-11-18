#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"

#include "maxruntime.h"
#define RHOR 1000.
#define MUR 100.


// Bubble 19 of Cano-Lozano et al, P.R.Fluids, 2016
const double Ga = 100.8;
const double Bo = 4.;
const double MAXTIME = 82;

const double WIDTH = 120. [1];
const double Zi = 3.5;
int LEVEL = 12;

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  size (WIDTH);
  origin (-L0/2, 0, -L0/2);
  init_grid (128);

  rho1 = 1. [0];
  rho2 = rho1/RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;
  TOLERANCE = 1e-4 [*];
  run();
}

event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(x) + sq(y - Zi) + sq(z) - sq(0.75) < 0 && level < LEVEL);
    fraction (f, sq(x) + sq(y - Zi) + sq(z) - sq(.5));
  }
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}

event adapt (i++) {
  double uemax = 1e-2;
  adapt_wavelet ({f,u}, (double[]){0.01,uemax,uemax,uemax}, LEVEL, 5);
}

event movie (t = 21; t <= MAXTIME; t += 0.25)
{
  
  view (fov = 5.0278, quat = {-0.132839,0.513023,0.0748175,0.844727},
	tx = 0.00149469, ty = -0.355489, width = 300, height = 800);
  travelling (50, 60, fov = 2.07254, tx = 0.00524944, ty = -0.513744);
  travelling (60, 82, tx = 0.00898601, ty = -0.703841);
  clear();
  draw_vof ("f", fc = {0.13,0.47,0.77});
  #if dimension == 3
  scalar l2[];
  lambda2 (u, l2);
  isosurface ("l2", -0.0002);
  #endif
  save ("bubble.mp4");
}