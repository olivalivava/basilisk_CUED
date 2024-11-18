#include "grid/multigrid.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "vof.h"
#include "tension.h"
#include "view.h"

#define radius (0.1)
#define Re 2000.
#define vrho1 (1.)
#define vrho2 (1./1e3)
#define vmu1 (2.*radius/Re)
#define vmu2 (2.*radius/Re*(1.81e-5/1e-3))


#define xcc1 1.
#define xcc2 -1.

double runtime = 6.; //set runtime length

int main (int argc, char * argv[])
{
  init_grid (256);
  size (4.);
  origin (-L0/2., -L0/2.);
  rho1 = vrho1, rho2 = vrho2;
  mu1 = vmu1, mu2 = vmu2;

  run();
}

event init (t = 0) {
  fraction (f, max(sq(radius) - sq(x - xcc1) - sq(y),
                    sq(radius) - sq(x - xcc2) - sq(y)));
  foreach() {
    u.x[] = f[]/2.;
  }
  boundary ({f,u});
}

event movie (t += 0.04; t <= runtime)
{
  clear();
  squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");
  box();
  save ("movie.mp4");
}
