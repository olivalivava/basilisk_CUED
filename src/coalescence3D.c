#include "grid/multigrid.h"
#define dimension 3

#include "navier-stokes/centered.h"
//#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"
#include "vof.h"
#include "tension.h"
#include "view.h"
//#include "axi.h" //disabled for the droplet_h_sigma=10 run

#define LEVEL 8

double We = 61.4;
double Re = 296.5;
double B = 0.06;
double R = 0.000168; // R1 = R/Rr, 0.000168

#define rho1c 1.251 //nitrogen (outside)
#define rho2c 762.0 //tetradecane (inside)
#define mu1c 1.787e-5//nitrogen 1.787e-5
#define mu2c 2.128e-3 //tetradecane
#define sigmac 0.0276 //tetradecane
#define RHOR 0.7 // 1.6417e-3, 0.06
#define MUR 0.0082*1.5 //8.186e-3, 1.25, 0.0123

//#define RHOR 0.06
//#define MUR 1.25
//#define Rr 1000
double runtime = 0.0003; //set runtime length
double uvel = 50.;

//double R1 = R; //set the radius for the left droplet (R1<1.)
//double R2 = R; //set the radius for the right (R2<1.)
//double uvel = 0.5; //set colliding speed (uniform velocity ??)
//#include "two-phase-generic.h"

// // Open boundaries at the left and right sides
// u.n[left] = neumann(0.);    // Free flow: no gradient in normal velocity
// u.t[left] = neumann(0.);    // Free slip: no gradient in tangential velocity

// u.n[right] = neumann(0.);
// u.t[right] = neumann(0.);


int main()
{
  size (10.*R);
  init_grid(64); // Base resolution
  origin (-L0/2., -L0/2. , -L0/2.);     //changed the origin
  rho2 = (We * sigmac)/(2*R*uvel*uvel);               //kg/m^3
//  rho1 = rho1c;
  // mu2 = mu2c;
  // mu1 = mu1c;
  rho1 = rho2*RHOR;      //used the parameters given on the website
//  double uvel = sqrt((We*sigmac)/(rho2*2*(R)));
  mu2 = rho2*uvel*(2*R)/Re;                      //Pa s
  mu1 = mu2*MUR;
  f.sigma = sigmac;           //N/m
  TOLERANCE = 1e-4 [*];       //defult 1e04
  run();
}

event init (t = 0)
{
  double X = B*(2*(R));
  fraction (f, max (- (sq(x + 1.3*R) + sq(y) + sq(z)- sq(R)),
		                - (sq(x - 1.3*R) + sq(y) + sq(z-X) - sq(R))));
  foreach() {
//      double uvel = sqrt((We*sigmac)/(rho2*2*(R)));
      u.x[] = - sign(x)*f[] * uvel; //how to assign velocity to each droplet?
  }
}

// event adapt(i++) {
//   foreach() {
//     refine(f[] > 0.01 && f[] < 0.99 && level < LEVEL);
//   }
// }

// event acceleration (i++) {
//   double g = 9.81 * sq(1/Rr); //using Bo = (rho1 * g * (2R)**2)/sigma
//   face vector av = a;
//   foreach_face(y)
//     av.y[] -= g;   //gravity = 9.81
// }


// event adapt (i++) {
//   double uemax = 1e-2;
//   adapt_wavelet ({f,u}, (double[]){0.001,uemax,uemax,uemax}, LEVEL, 5);
// }

event logfile (i++) {
//  double uvel = sqrt((We*sigmac)/(rho2*2*(R)));
  fprintf(stderr, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g", 
                    t, dt,rho1, rho2, mu1, mu2, f.sigma, uvel, perf.t, perf.speed);
  putchar ('\n');
  fflush (stdout);
}

event movie (t += 2.0e-6; t <= runtime)
{
  clear();
  //view (width = 20*R, height = 10*R);
  //squares ("u.x", spread = -1, linear = true); //removed this
  fprintf(stderr, "Generating frame at t = %g\n", t);
  draw_vof ("f");

  box();
  save ("movie3D.mp4");
}
