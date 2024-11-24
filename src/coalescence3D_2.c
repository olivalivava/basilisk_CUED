#include "grid/multigrid.h"
#define dimension 3

#include "navier-stokes/centered.h"
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
//#include "axi.h" // Not symmetrical about z-axis

#define LEVEL 8

double We = 70.8;
double Re = 327.7;
double B = 0.25;
double R = 0.178; // R1 = R/Rr


#define rho1c 762.0 //tetradecane
#define rho2c 1.251 //nitrogen
#define sigmac 0.0276 //tetradecane 0.0276
#define MUR 97.7
#define Rr 1000

double runtime = 10.0; //set runtime length

//double R1 = R; //set the radius for the left droplet (R1<1.)
//double R2 = R; //set the radius for the right (R2<1.)
//double uvel = 0.5; //set colliding speed (uniform velocity ??)
//#include "two-phase-generic.h"

// Open boundaries at the left and right sides
u.n[left] = neumann(0.);    // Free flow: no gradient in normal velocity
u.t[left] = neumann(0.);    // Free slip: no gradient in tangential velocity

u.n[right] = neumann(0.);
u.t[right] = neumann(0.);


int main()
{
  size (15.*R);
  init_grid(128); // Base resolution
  origin (-L0/2., -L0/2. , -L0/2.);     //changed the origin
  double uvel = sqrt((We*sigmac)/(rho1c*2*(R)));
  rho1 = rho1c;               //kg/m^3
  rho2 = rho2c;               //kg/m^3
  mu1 = (rho1*uvel*2*(R))/Re;  // (Pa * s)
  mu2 = mu1/MUR;              //Pa s
  f.sigma = sigmac;           //N/m
  TOLERANCE = 1e-4 [*];       //defult 1e04
  run();
}

event init (t = 0)
{
  double X = B*(2*(R));
  fraction (f, max (- (sq(x + R) + sq(y + X/2.) + sq(z)- sq(R)),
		                - (sq(x - R) + sq(y - X/2.) + sq(z) - sq(R))));
  foreach() {
      double uvel = sqrt((We*sigmac)/(rho1c*2*(R)));
      u.x[] = - sign(x)*f[] * uvel; //how to assign velocity to each droplet?
  }
}

// event adapt(i++) {
//   foreach() {
//     refine(f[] > 0.01 && f[] < 0.99 && level < LEVEL);
//   }
// }

// event acceleration (i++) {
//   double g = 9.81 * sq(1/Rr) //using Bo = (rho1 * g * (2R)**2)/sigma
//   face vector av = a;
//   foreach_face(y)
//     av.y[] -= g;   //gravity = 9.81
// }


// event adapt (i++) {
//   double uemax = 1e-2;
//   adapt_wavelet ({f,u}, (double[]){0.001,uemax,uemax,uemax}, LEVEL, 5);
// }


event movie (t += 0.20; t <= runtime)
{
  clear();
  //view (width = 20*R, height = 10*R);
  //squares ("u.x", spread = -1, linear = true); //removed this
  draw_vof ("f");

  box();
  save ("movie3D_2.mp4");
}
