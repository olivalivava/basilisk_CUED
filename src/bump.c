#include "saint-venant.h"

event init (t = 0) {
  double a = 1., b = 200.;
  foreach()
    h[] = 0.1 + a*exp(- b*(x*x + y*y));
}

event graphs (i++) {
  stats s = statsf (h);
  fprintf (stderr, "%g %g %g\n", t, s.min, s.max);
}

event images (i++) {
  output_ppm(h);
}

event end (i = 10) {
  printf ("i = %d t = %g\n", i, t);
}

int main() {
  origin (-0.5, -0.5);
  init_grid (256);
  run();
}
