# splinecoefs

Spline coefficient fitting to FORCE ARD

(C) 2025-2026, David Klehr, David Frantz


## Usage

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

## Parameters used for testing:

```
-k = 4
-c = 22
-l = 1000
-w = 0.2
```

## Input data

The process is invoked on an individual tile. 
You need to loop over tiles with your tool of choice.

FORCE ARD in the usual format are expected. 

The user needs to input the data via a csv-table:

- 3 unnamed columns, separated by ``,`` or ``;``
- column 1: BOA images (from or compatible with FORCE)
- column 2: QAI images (from or compatible with FORCE)
- column 3: THT images (see next chapter)
- full file paths or relative file paths (with respect to workspace when calling the program)

The BOA images are expected to be from a single sensor-family and need to include bands for NDVI calculation. You cannot mix sensors where bands differ, neither in ordering nor in the band number.

The products in each row need to be from the same date.
The date needs to be encoded in the filename, i.e., the first 8 digits of the basename in *YYYYMMDD*.

The rows should be time-ordered - earliest on top, latest on bottom.

There shouldn't be images later than the target year.

There can be multiple years before the target year to stabilize the spline fit in case there is not enough data in the target year. Each year is weighted in comparison to its similarity to the target year (``-y``), but with a maximum given weight (``-w``).

## THT images

The splines are fitted such that the x-axis represents the time, i.e., the day-of-year (0-365).
Instead of using the image acquisition date, cumulative thermal time (THT) is being used as a means to homogenize phenology across larger areas.
The THT images should be formatted similar to the FORCE data, and need to contain the cumulative THT for each pixel for each BOA observation.
The THT should be normalized to [0-365].


