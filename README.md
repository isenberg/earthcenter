# Geographic Center of Earth Calculator
Calculates the geographic center of all land surfaces on Earth.

Author and copyright 2003-2025:
Holger Isenberg, [@areoinfo](https://x.com/areoinfo), https://areo.info

Geometrically defined, the geographic center is the geometric median
of all land surfaces within the two dimensions of the spherical surface
of Earth. This location is also known as solution of the Facility Location Problem
in Operations Research. The Geoid would more exactly approximate Earth's outer shape,
but as the elevation differences on Earth in relation to its circumference,
are insignificant, the spherical surface is for this calculation purpose sufficient.

The geographic center as the geometric median is the location on the surface with the
minimum distance sum to all other locations on land. The distances are measured
on the sphere surface, which is two-dimensional without boundary,
and not three-dimensional as straight line through the Earth underground.

The distance definition, the shortest path between two locations
on the spherical surface, is the arc of a great circle (orthodrome).

This calculation uses a simple gradient descent search to find the center.

The center definition and its calculation result is due to its use of
the great circle distance independent from the type of map projection.
A transformation from the projection of the input map to geographic coordinates
and unit areas of sample points is done, but this implementation only contains
the transformation for an equidistant cylindrical projection.

An alternative center definition can be selected with the squaremode parameter.
That definition uses the minimum sum of the distance squares (power of 2)
and is equivalent to the center of gravity for the 2d spherical surface
(not the 3d shape), in geometry also known as centroid.

For better code-readability a single-threaded variant is available
as EarthCenter2025SingleThread.java. It will run about 10 times slower
on modern notebooks.

A current MacBook completes the calculation in less than 1 minute.

## References

### The Center of the Earth
Andrew J. Woods, San Diego, 1973

https://archive.org/details/centerofearth0000wood

https://www.icr.org/article/50  by Henry M. Morris, without tables and maps

 * Part I by Henry M. Morris provides a philosophical context.
 * Part II by Andrew J. Woods explains the calculation and presents the results.
 * Uses as center definition the minimum sum of great-circle distances.
 * In addition discusses an alternative calculation using the distance squares.
 * The calculations are done based on the current sea level
  and alternative levels -70m and +70m.

### Giza, Center of Earth
Holger Isenberg

1st edition, October 19, 2003: https://mars-news.de/pyramids/gizacenter.html

2nd edition, August 25, 2025: https://x.com/areoinfo/status/1960235682200531262 and https://areoinfo.substack.com/giza-center-of-earth

 * Article about the calculation history since 1864 and the results.
 * Uses as center definition the minimum sum of great-circle distances
   for today's sea level and in 10m steps up to +300m.

### A New Method for Finding Geographic Centers, with Application to U.S. States
Peter A. Rogerson, 2015

https://ui.adsabs.harvard.edu/abs/2015ProfG..67..686R/abstract

 * Uses as center definition the minimum sum of distance squares.
 * ```EarthCenter2025.java``` also offers that alternative calculation by adding ```-Dsquaremode=true```

## Usage Preparations
Download a Java 21 JDK for best performance, or if you already have Java
installed: at least version 11 JDK is needed.
For example available on: https://javaalmanac.io/jdk/21

Download a single GeoTIFF global map file into the same directory
where EarthCenter2025.java is located from
https://ncei.noaa.gov/products/etopo-global-relief-model
for example 60 Arc-Second Resolution Ice surface_elevation geotiff:
```ETOPO_2022_v1_60s_N90W180_surface.tif```

Decompress the map-TIFF to avoid Java ImageIO read failures
(Illegal value for Predictor in TIFF file).
using ImageMagick from https://imagemagick.org:
```
convert ETOPO_2022_v1_60s_N90W180_surface.tif -compress none ETOPO_2022_v1_60s_N90W180_surface_uncompressed.tif
```
or using GDAL from https://gdal.org:
```
gdal_translate ETOPO_2022_v1_60s_N90W180_surface.tif ETOPO_2022_v1_60s_N90W180_surface_uncompressed.tif -co COMPRESS=NONE
```

## Usage Examples
To calculate the center for today's sea level using the default
map ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif:
```
java -Dsealevel=0 EarthCenter2025.java
```

at a sea level of 178 meter above today's with a higher resolution map:
```
java -Dsealevel=178 -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
```

today's sea level, but find the minimum sum of the distance squares:
```
java -Dsquaremode=true -Dsealevel=0 EarthCenter2025.java
```

quick preview with reduced map resolution of 15' (28km) for today's sea level:
```
java -Dmapwidth=1440 -Dsealevel=0 EarthCenter2025.java
```

today's sea level, consider antarctic ice shelfs as land:
```
java -Dsealevel=0 -Dantarctic=0 EarthCenter2025.java
```

repeat Woods' 1973 calculation:
```
java -Dsealevel=0 -Dmapwidth=180 -Dcalcwidth=360 -Dantarctic=0 -Dstartlat=20 -Dstartlon=20 EarthCenter2025.java
```

repeat Wood's 1973 alternative sum of distance squares calculation, similar to Rogerson's:
```
java -Dsquaremode=true -Dsealevel=0 -Dmapwidth=180 -Dcalcwidth=360 -Dantarctic=0 -Dstartlat=20 -Dstartlon=20 EarthCenter2025.java
```

## Usage Notes
 * In case a java.lang.OutOfMemoryError occurs,
   for the ETOPO 60s map, adding -Xmx4G is sufficient,
   for the large ETOPO 30s map, -Xmx8G is needed.
   By default, Java uses at maximum 25% of the system RAM size.
 * To exit, close the window with the red icon or X icon or Command Q,
   or press Control C in the terminal.

## Known Issues
 * escapes map area for some starting points, or gets stuck at low resolution.
   As workaround add ```-Dstartlat=0 -Dstartlon=0``` or other starting location.

## Example Results
```
% java -Xmx8G -Dsealevel=0 -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distances
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=0m, antarctic threshold=70m, elevation max=8353m, map width: 43200px
land area: 171099298, sea area: 422943343, land ratio: 29%
elevation mean: 146.83m
elevation median above sealevel: 423.00m
geographic center: 40.357°N, 34.655°E, 0m sealevel, 7358.975km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 70m antarctic threshold
completed: 28s, 116 iterations
```

```
% java -Xmx8G -Dsealevel=178 -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distances
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=+178m, antarctic threshold=70m, elevation max=8175m, map width: 43200px
land area: 127589793, sea area: 466452848, land ratio: 21%
elevation mean: 118.49m
elevation median above sealevel: 460.00m
geographic center: 30.002°N, 31.535°E, +178m sealevel, 7339.747km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 70m antarctic threshold
completed: 43s, 122 iterations
```

```
% java -Xmx8G -Dantarctic=0 -Dsealevel=0 -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distances
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=0m, antarctic threshold=0m, elevation max=8353m, map width: 43200px
land area: 172467295, sea area: 421575346, land ratio: 29%
elevation mean: 146.89m
elevation median above sealevel: 418.00m
geographic center: 39.375°N, 34.395°E, 0m sealevel, 7414.960km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 0m antarctic threshold
completed: 29s, 114 iterations
```

```
% java -Xmx8G -Dsealevel=0 -Dsquaremode=true -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distance squares
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=0m, antarctic threshold=70m, elevation max=8353m, map width: 43200px
land area: 171099298, sea area: 422943343, land ratio: 29%
elevation mean: 146.83m
elevation median above sealevel: 423.00m
geographic center of distance squares: 33.867°N, 27.223°E, 0m sealevel, 10882.522km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 70m antarctic threshold
completed: 28s, 123 iterations
```

```
% java -Xmx8G -Dsealevel=178 -Dsquaremode=true -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distance squares
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=+178m, antarctic threshold=70m, elevation max=8175m, map width: 43200px
land area: 127589793, sea area: 466452848, land ratio: 21%
elevation mean: 118.49m
elevation median above sealevel: 460.00m
geographic center of distance squares: 28.173°N, 27.320°E, +178m sealevel, 10817.454km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 70m antarctic threshold
completed: 24s, 115 iterations
```

```
% java -Xmx8G -Dantarctic=0 -Dsealevel=0 -Dsquaremode=true -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
Geographic Center of Earth Calculator
(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info
version: sequential gradient descent 20250825
loading map: ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif
mode: minimum sum of greatcircle distance squares
display: 3840x1920
map raster: 43200x21600
calculation raster: 216000x108000
processors: 10
sea level=0m, antarctic threshold=0m, elevation max=8353m, map width: 43200px
land area: 172467295, sea area: 421575346, land ratio: 29%
elevation mean: 146.89m
elevation median above sealevel: 418.00m
geographic center of distance squares: 32.322°N, 27.043°E, 0m sealevel, 11029.671km average distance
parameters: 30arcsec (0.93km) map resolution, 6arcsec (0.19km) calculation resolution, 0m antarctic threshold
completed: 41s, 120 iterations
```
