#include "grid/multigrid.h"

#define dimension 3

#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
#include "axi.h"

#ifndef LEVEL
# define LEVEL 8
#endif

double We = 61.4;
double Re = 296.5;
double B = 0.0;
double R = 0.000168;

#define rho1c 762.0 //tetradecane
#define rho2c 1.251 //nitrogen
#define sigmac 0.0276 //tetradecane
#define MUR 97.7


//double R1 = R; //set the radius for the left droplet (R1<1.)
//double R2 = R; //set the radius for the right (R2<1.)

const double Zi = 3.5;
//double uvel = 0.5; //set colliding speed (uniform velocity ??)
double runtime = 2.0; //set runtime length
//#include "two-phase-generic.h"

// Open boundaries at the left and right sides
u.n[left] = neumann(0.);    // Free flow: no gradient in normal velocity
u.t[left] = neumann(0.);    // Free slip: no gradient in tangential velocity
p[left] = dirichlet(0.);    // Pressure is fixed at zero

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);
p[right] = dirichlet(0.);

int main()
{
  size (15*R);
  origin (-L0/2., 0, -L0/2.);
  init_grid(64); // Higher grid resolution
  double uvelc = sqrt((We*sigmac)/(rho1c*2*R));
  rho1 = rho1c;
  rho2 = rho2c;
  mu1 = (rho1*uvelc*2*R)/Re;
  mu2 = mu1/MUR;
  f.sigma = sigmac;
  TOLERANCE = 1e-5 [*];
  run();
}

event init (t = 0)
{
  fraction (f, max (- (sq(x + 2*R) + sq(y) + sq(z)- sq(R)),
		    - (sq(x - 2*R) + sq(y) + sq(z) - sq(R))));
  foreach() {
      double uvelc = sqrt((We*sigmac)/(rho1c*2*R));
      u.x[] = - sign(x)*f[] * uvelc; //how to assign velocity to each droplet?
  }
}

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1;
}


event movie (t += 0.0004; t <= runtime)
{
  clear();
  view (width = 20*R, height = 10*R);
  //squares ("u.x", spread = -1, linear = true);
  draw_vof ("f");

  box();
  save ("movie3D.mp4");
}
