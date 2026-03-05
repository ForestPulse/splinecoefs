# splinecoefs

Spline coefficient fitting to FORCE ARD

(C) 2025-2026, David Klehr, David Frantz

```
Usage: splinecoefs -j cpus -i input-table -x mask-image -o output-image
          -k order -c control-points -l lambda -y target-year
          -r = RED-band -n = NIR-band

  -j = number of CPUs to use

  -i = input-table = csv file with input images
  -x = mask image
  -o = output file (.tif)

  -k = order of the spline
  -c = number of control points (equiv. to number of output coefficients)
  -l = lambda (smoothing parameter)
  -y = target year
  -w = maximum weight for past data points
  -r = RED band number
  -n = NIR band number
```
