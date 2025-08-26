/*
 * Geographic Center of Earth Calculator
 *
 * Calculates the geographic center of all land surfaces on Earth.
 *
 * Author and copyright: Holger Isenberg, @areoinfo, https://areo.info
 * this version 2025, first simple C version October 2003
 *
 * Geometrically defined, the geographic center is the geometric median
 * of all land surfaces within the two dimensions of the spherical surface
 * of Earth. The Geoid would more exactly approximate Earth's outer shape,
 * but as the elevation differences on Earth in relation to its circumference,
 * are insignificant, the spherical surface is for this calculation purpose sufficient.
 *
 * The geographic center as the geometric median is the location on the surface with the
 * minimum distance sum to all other locations on land. The distances are measured
 * on the sphere surface, which is two-dimensional without boundary,
 * and not three-dimensional as straight line through the Earth underground. 
 *
 * The distance definition, the shortest path between two locations
 * on the spherical surface, is the arc of a great circle (orthodrome).
 *
 * This calculation uses a simple gradient descent search to find the center.
 *
 * The center definition and its calculation result is due to its use of
 * the great circle distance independent from the type of map projection.
 * A transformation from the projection of the input map to geographic coordinates
 * and unit areas of sample points is done, but this implementation only contains
 * the transformation for an equidistant cylindrical projection.
 * 
 * References:
 *
 * The Center of the Earth. Andrew J. Woods, San Diego, 1973
 *   Part I by Henry M. Morris provides a philosophical context.
 *   Part II by Andrew J. Woods explains the calculation and presents the results.
 *   Uses as center definition the minimum sum of great-circle distances.
 *   In addition discusses an alternative calculation using the distance squares.
 *   The calculations are done based on the current sea level
 *   and alternative levels -70m and +70m.
 *   https://archive.org/details/centerofearth0000wood
 *   https://www.icr.org/article/50  (without tables and maps)
 *
 * Giza, Center of Earth. Holger Isenberg, October 19, 2003
 *   Article about the calculation history since 1864 and the results.
 *   Uses as center definition the minimum sum of great-circle distances
 *   for today's sea level and in 10m steps up to +300m.
 *   https://mars-news.de/pyramids/gizacenter.html
 *
 * A New Method for Finding Geographic Centers, with Application to U.S. States,
 *   Peter A. Rogerson, 2015
 *   Uses as center definition the minimum sum of distance squares.
 *   EarthCenter2025.java also offers that alternative calculation by adding -Dsquaremode=true
 *   https://ui.adsabs.harvard.edu/abs/2015ProfG..67..686R/abstract
 *
 * Usage Preparations:
 *  1. Download a Java 21 JDK for best performance, or if you already have Java
 *     installed: at least version 11 JDK is needed.
 *     For example available on: https://javaalmanac.io/jdk/21
 *  2. Download a single GeoTIFF global map file into the same directory
 *     where EarthCenter2025.java is located from
 *     https://ncei.noaa.gov/products/etopo-global-relief-model
 *     for example 60 Arc-Second Resolution Ice surface_elevation geotiff:
 *     ETOPO_2022_v1_60s_N90W180_surface.tif
 *  3. Decompress the map-TIFF to avoid Java ImageIO read failures
 *     (Illegal value for Predictor in TIFF file).
 *     using ImageMagick from https://imagemagick.org:
 *       convert ETOPO_2022_v1_60s_N90W180_surface.tif -compress none ETOPO_2022_v1_60s_N90W180_surface_uncompressed.tif
 *     or using GDAL from https://gdal.org:
 *       gdal_translate ETOPO_2022_v1_60s_N90W180_surface.tif ETOPO_2022_v1_60s_N90W180_surface_uncompressed.tif -co COMPRESS=NONE
 *
 * Usage Examples:
 *
 *  To calculate the center for today's sea level using the default
 *  map ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif:
 *  java -Dsealevel=0 EarthCenter2025.java
 *
 *  at a sea level of 178 meter above today's with a higher resolution map:
 *  java -Dsealevel=178 -Dmap=ETOPO_2022_v1_30s_N90W180_surface_uncompressed.tif EarthCenter2025.java
 *
 *  today's sea level, but find the minimum sum of the distance squares:
 *  java -Dsquaremode=true -Dsealevel=0 EarthCenter2025.java
 *
 *  quick preview with reduced map resolution of 15' (28km) for today's sea level:
 *  java -Dmapwidth=1440 -Dsealevel=0 EarthCenter2025.java
 *
 *  today's sea level, consider antarctic ice shelfs as land:
 *  java -Dsealevel=0 -Dantarctic=0 EarthCenter2025.java
 *
 *  repeat Woods' 1973:
 *  java -Dsealevel=0 -Dmapwidth=180 -Dcalcwidth=360 -Dantarctic=0 -Dstartlat=20 -Dstartlon=20 EarthCenter2025.java
 *
 *  repeat Wood's 1973 alternative square sum:
 *  java -Dsquaremode=true -Dsealevel=0 -Dmapwidth=180 -Dcalcwidth=360 -Dantarctic=0 -Dstartlat=20 -Dstartlon=20 EarthCenter2025.java
 *
 * Usage Notes:
 *  In case a java.lang.OutOfMemoryError occurs,
 *  for the ETOPO 60s map, adding -Xmx4G is sufficient,
 *  for the large ETOPO 30s map, -Xmx8G is needed.
 *  By default, Java uses at maximum 25% of the system RAM size.
 *  To exit, close the window with the red icon or X icon or Cmd Q,
 *  or press Control C in the terminal.
 *
 * Known Issues:
 *  * escapes map area for some starting points, or gets stuck at low resolution.
 *    As workaround add -Dstartlat=0 -Dstartlon=0 or other starting location.
 */

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.awt.image.Raster;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.Stack;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;
import java.util.stream.IntStream;

import javax.imageio.ImageIO;
import javax.imageio.ImageReader;
import javax.imageio.stream.ImageInputStream;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class EarthCenter2025 extends JPanel {
	static BufferedImage mapImg, tmpImg;
	static final long serialVersionUID = 20250825L;

	public static void main(String[] args) {
		// Program start.
		JFrame frame = new JFrame("Geographic Center of Earth Calculator    \u00A9 2003 - 2025 Holger Isenberg  @areoinfo  https://areo.info");
		try {
			// map in equidistant cylindrical projection (digital elevation model) in uncompressed geoTIFF format:
			String mapFilename = System.getProperty("map", "ETOPO_2022_v1_60s_N90W180_surface_uncompressed.tif");
			int level = Integer.getInteger("sealevel", 0); // sea level elevation in meter, example: 0 or 178
			String startLatStr = System.getProperty("startlat", "-40.0"); // calculation start latitude in positive north degree
			double startLat = Double.parseDouble(startLatStr);
			String startLonStr = System.getProperty("startlon", "-110.0"); // calculation start longitude in positive east degree
			double startLon = Double.parseDouble(startLonStr);
			String startStepStr = System.getProperty("startstep", "10.0"); // degrees, initial step size
			double startStep = Double.parseDouble(startStepStr);
			int mapWidthOpt = Integer.getInteger("mapwidth", 0); // pixels, scale down to this map width for calculation, default from provided map file
			int calcWidthOpt = Integer.getInteger("calcwidth", 0); // pixels, final near center calculation resolution, default mapWidth
			String name = System.getProperty("name", null); // optional output result filename
			int antarcticThreshold = Integer.getInteger("antarctic", 70); // elevation below that threshold south of 60°S is considered water
			boolean squareMode = Boolean.getBoolean("squaremode"); // if enabled, calculates sum of distance squares instead of plain sum

			// start message
			System.out.println("Geographic Center of Earth Calculator\n(c) 2003 - 2025 Holger Isenberg @areoinfo https://areo.info");
			System.out.println("version: sequential gradient descent " + serialVersionUID);
			System.out.println("loading map: " + mapFilename);

			// map loading
			File mapFile = new File(mapFilename);
			Raster tiffDEM = loadGeoTIFFDEM(mapFile);
			
			// display resizing to map size
			int mapWidth;
			int calcWidth;
			if (mapWidthOpt == 0) {
				mapWidth = tiffDEM.getWidth();
			} else {
				mapWidth = mapWidthOpt;
			}
			if (calcWidthOpt == 0 || calcWidthOpt < mapWidth) {
				calcWidth = mapWidth * 5;
			} else {
				calcWidth = calcWidthOpt;
			}
			Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
			//int width = Math.min(resPerc * land_w / 100, screenSize.width);
			int width = Math.min(mapWidth, screenSize.width);
			//frame.setContentPane(new EarthCenter2025());
			frame.setContentPane(new EarthCenter2025());
			frame.getContentPane().setLayout(null);
			frame.setSize(width, width/2 + 20);
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			//frame.pack();
			frame.setVisible(true);
			
			// display preparation
			mapImg = new BufferedImage(frame.getWidth(), frame.getHeight(), BufferedImage.TYPE_4BYTE_ABGR);
			tmpImg = new BufferedImage(frame.getWidth(), frame.getHeight(), BufferedImage.TYPE_4BYTE_ABGR);
			Graphics2D gmap = mapImg.createGraphics();
			Graphics2D gtmp = tmpImg.createGraphics();
			gmap.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			gtmp.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			
			// calculation start
			GradientCalc centerCalc = new GradientCalc(frame, gmap, gtmp);
			centerCalc.calculate(tiffDEM, level, startLat, startLon, startStep, mapWidth, calcWidth, name, antarcticThreshold, squareMode);
		} catch(Exception ex) {
			System.out.println(ex);
			frame.dispose();
		}
	}

	EarthCenter2025() {
		// make the content of the map window visible
		setOpaque(true);
	}

	protected void paintComponent(Graphics g) {
		// draw the map and calculation indicators
		super.paintComponent(g);
		g.drawImage(mapImg, 0, 0, null);
		g.drawImage(tmpImg, 0, 0, null);
	}

	static Raster loadGeoTIFFDEM(File mapFile) throws IOException {
		// Loads the GeoTIFF map into memory.
		// Doesn't support ZIP-compressed TIFFs, those will result in this error:
		// javax.imageio.IIOException: Illegal value for Predictor in TIFF file
		Iterator<ImageReader> readers = ImageIO.getImageReadersByFormatName("TIFF");
		if (!readers.hasNext()) {
			throw new IOException("No TIFF readers available in this JRE.");
		}
		ImageReader tiffReader = readers.next();
		ImageInputStream tiffStream = ImageIO.createImageInputStream(mapFile);
		// ignoreMetadata = true to avoid javax.imageio.IIOException: Unexpected count 15 for DateTime field
		tiffReader.setInput(tiffStream, false,  true);
		//int ni = tiffReader.getNumImages(true);
		BufferedImage img = tiffReader.read(0);
		Raster tiffDEM = img.getData();
		tiffReader.dispose();
		tiffStream.close();
		return tiffDEM;
	}

	static class GradientCalc {
		// The actual center calculation is contained in this class.
		JFrame frame;
		Graphics2D gmap, graphics;
		Raster tiffDEM;
		final double earthCircumference = 40075.0;
		final double antarcticMaxLat = -60.0;
		double screenScaleMap = 1.0;
		double calcScaleScreen = 1.0;
		int antarcticThreshold = 0;
		boolean squareMode = false;
		double piWidth, piRes, piCalcRes;
		int geoTiffWidth;
		int geoTiffHeight;
		int antarcticYmax;
		int seaLevel;

		GradientCalc(JFrame frame, Graphics2D gmap, Graphics2D graphics) {
			// Initializes the display window.
			// If the map raster resolution is larger than the screen width in pixels
			// the display window will scale down to the screen width.
			// That won't affect the result precision or the map resolution used for the calculation.
			this.frame = frame;
			this.gmap = gmap;
			this.graphics = graphics;
		}

		int getElevation(int x, int y) {
			// Returns the elevation in meter above the set sea level
			// at the given map raster coordinates.
			// Considers surfaces south of antarctic_ymax below antarcticThreshold meters
			// as water and returns 0 for those.
			int[] pixel = new int[3];
			tiffDEM.getPixel(x, y, pixel);
			if (y > antarcticYmax && pixel[0] < antarcticThreshold) {
				// consider this as floating ice shelf which needs to be seen as water
				return 0;
			} else {
				return pixel[0]; // all 3 tuple elements are identical for monochrome GeoTIFFs
			}
		}

		class MPoint {
			// point on map in raster coordinates
			// x==0: 180W longitude
			// x==mapWidth: 180E longitude
			// y==0: 90N latitude
			// y==mapHeight-1: 90S latitude
			public int x, y;
			public double sum;
			public MPoint(int x, int y, double sum) {
				this.x = x;
				this.y = y;
				this.sum = sum;
			}
			public MPoint() {
				this.x = Integer.MIN_VALUE;
				this.y = Integer.MIN_VALUE;
				this.sum = Double.NaN;
			}
		}

		class MVector {
			// vector on the map
			public double x, y, sum;
			MVector() {
				this.x = Double.NaN;
				this.y = Double.NaN;;
				this.sum = Double.NaN;
			}
			MVector(double x, double y, double sum) {
				this.x = x;
				this.y = y;
				this.sum = sum;
			}
		}

		class PathSegment {
			// path segment between two map points for the gradient search 
			public double sLat, sLon, eLat, eLon;
			public double gradient;
			PathSegment(double sLat, double sLon, double eLat, double eLon, double gradient) {
				this.sLat = sLat;
				this.sLon = sLon;
				this.eLat = eLat;
				this.eLon = eLon;
				this.gradient = gradient;
			}
		}

		public void calculate(Raster tiffDEM, int seaLevel, double startLat, double starLon, double startStep, int mapWidth, int calcWidth, String outputprefix, int antarcticLevel, boolean squareMode) throws IOException {
			// The geographic center calculation with map display.
			this.tiffDEM = tiffDEM;
			int displayWidth = frame.getWidth();
			int displayHeight = displayWidth / 2;
			int mapHeight = mapWidth / 2;
			int calcHeight = calcWidth / 2;
			Color oceanColor = new Color(0, 30, 50);
			Color landColor = Color.DARK_GRAY;
			String levelSign = (seaLevel>0?"+":"");
			gmap.setBackground(oceanColor);
			gmap.clearRect(0, 0, displayWidth, displayHeight);
			frame.repaint();

			// read map size
			geoTiffWidth = tiffDEM.getWidth();
			geoTiffHeight = tiffDEM.getHeight();
			antarcticYmax = (int)Math.round((90-antarcticMaxLat)/180*mapHeight);
			
			// unit conversions
			this.antarcticThreshold = antarcticLevel;
			this.squareMode = squareMode;
			this.piWidth=2*Math.PI/geoTiffWidth;
			this.piRes=2*Math.PI/mapWidth;
			this.piCalcRes=2*Math.PI/calcWidth;
			this.seaLevel = seaLevel;
			screenScaleMap = (double) mapWidth / displayWidth;
			calcScaleScreen =  (double)displayWidth / calcWidth;

			// optional output map graphics writer
			PrintWriter filePGM = null;
			if (outputprefix != null) {
				String namePGM = outputprefix + seaLevel + "m.pgm";
				FileWriter writerPGM = new FileWriter(namePGM);
				filePGM = new PrintWriter(writerPGM);
				System.out.println("PNM: output \"" + namePGM + "\" " + mapWidth + " x " + mapHeight);
			}

			// start message
			if (squareMode) {
				System.out.println("mode: minimum sum of greatcircle distance squares");
			} else {
				// default
				System.out.println("mode: minimum sum of greatcircle distances");
			}
			System.out.println("display: " + displayWidth + "x" + displayHeight);
			System.out.println("map raster: " + mapWidth + "x" + mapHeight);
			System.out.println("calculation raster: " + calcWidth + "x" + calcHeight);
			System.out.println("processors: " + Runtime.getRuntime().availableProcessors());

			// map display and elevation statistic calculation
			// those statistics are not used for the actual center calculation
			BufferedImage mapImage = new BufferedImage(displayWidth, displayHeight, BufferedImage.TYPE_INT_ARGB);
			int[] mapPixels = ((DataBufferInt) mapImage.getRaster().getDataBuffer()).getData();
			double piRes = 2 * Math.PI/mapWidth;
			ConcurrentHashMap<Integer, Double> elevationHistogram = new ConcurrentHashMap<Integer, Double>(10000);
			AtomicInteger elevationMaxAdder = new AtomicInteger(0);
			DoubleAdder landAreaAdder = new DoubleAdder();
			DoubleAdder seaAreaAdder = new DoubleAdder();
			DoubleAdder elevationSumFlatAdder = new DoubleAdder();
			DoubleAdder elevationSumAdder = new DoubleAdder();
			IntStream.range(0, mapHeight).parallel().forEachOrdered(yt -> {
				double area = Math.sin(piRes * yt);
				int[] rowDraw = new int[displayWidth];
				IntStream.range(0, mapWidth).parallel().forEach(xt -> {
					int elevation = getElevation(xt*geoTiffWidth/mapWidth, yt*geoTiffWidth/mapWidth) - seaLevel;
					if(elevation > 0) {
						int dispX = (int)(xt/screenScaleMap);
						if (xt % (int)screenScaleMap == 0) {
							rowDraw[dispX] = landColor.getRGB();
						}
						elevationSumFlatAdder.add(elevation);
						elevationSumAdder.add(area * elevation);
						landAreaAdder.add(area);
						elevationMaxAdder.accumulateAndGet(elevation, Math::max);
						elevationHistogram.merge(elevation, area, Double::sum);
					} else {
						seaAreaAdder.add(area);
					}
				});
				if (yt % (int)screenScaleMap == 0) {
					int dispY = (int)(yt / screenScaleMap);
					System.arraycopy(rowDraw, 0, mapPixels, dispY * displayWidth, rowDraw.length);
				}
				if (yt % (int)screenScaleMap == 0) {
					gmap.drawImage(mapImage, null, 0, 0);
					frame.repaint();
				}
			});
			gmap.drawImage(mapImage, null, 0, 0);
			frame.repaint();
			int elevationMax = elevationMaxAdder.get();
			double landArea = landAreaAdder.sum();
			double seaArea = seaAreaAdder.sum();
			double landRatio = landArea/(landArea+seaArea);
			double elevationSumFlat = elevationSumFlatAdder.sum();
			double elevationSum = elevationSumAdder.sum();
			double elevationMeanGlobeFlat = elevationSumFlat / (mapWidth * mapHeight);
			double elevationMeanGlobeHalfFlat = elevationSum / (mapWidth * mapHeight);
			double elevationMeanGlobe = elevationSum / (landArea + seaArea);
			double elevationMeanLand = elevationSum / landArea;
			double elevationMedianLand = median(elevationHistogram);
			elevationHistogram.put(0, seaArea);
			double elevationMedianGlobe = median(elevationHistogram);
			System.out.printf("sea level=%s%dm, antarctic threshold=%dm, elevation max=%dm, map width: %dpx\n", levelSign, seaLevel, antarcticLevel, elevationMax, mapWidth);
			System.out.printf("land area: %.0f, sea area: %.0f, land ratio: %.0f%%\n", landArea, seaArea, 100*landRatio); 
			System.out.printf("elevation mean flat: %.2fm\n", elevationMeanGlobeFlat);
			System.out.printf("elevation mean half-flat: %.2fm\n", elevationMeanGlobeHalfFlat);
			System.out.printf("elevation mean global: %.2fm\n", elevationMeanGlobe);
			System.out.printf("elevation mean above sealevel: %.2fm\n", elevationMeanLand);
			System.out.printf("elevation median global: %.2fm\n", elevationMedianGlobe);
			System.out.printf("elevation median above sealevel: %.2fm\n", elevationMedianLand);
			if (filePGM != null) {
				// write map with selected sea level for reference
				filePGM.println("P2");
				filePGM.println("# Earth with sea level " + levelSign + seaLevel + "m, antarctic threshold " + antarcticLevel + "m, "
						+ "land/water ratio " + (int)(landRatio*100.0) + "%");
				filePGM.println("# global elevation data ETOPO2 from NOAA, http://www.ngdc.noaa.gov/mgg/");
				filePGM.println(mapWidth + " " + mapHeight);
				filePGM.println("255");
				for (int yt = 0; yt < mapHeight; yt++) {
					for (int xt = 0; xt < mapWidth; xt++) {
						filePGM.printf("%4d", mapPixels[yt*mapWidth + xt] & 0xFF);
					}
					filePGM.println();
				}
				filePGM.close();
			}

			// center calculation
			Stack<MPoint> path = new Stack<MPoint>();
			List<PathSegment> drawPath = new ArrayList<PathSegment>();
			int iter = 0;
			int start_xt=(int)((180.0+starLon)/360.0*calcWidth+0.5);
			int start_yt=(int)((90.0-startLat)/360.0*calcWidth+0.5);
			int start_st=(int)(startStep/360.0*mapWidth+0.5);
			int mapStep = start_st;
			int prevMapStep = mapStep;
			int minMapStep = (int)Math.round(calcWidth/mapWidth);
			int calcStep = mapStep;
			int prevCalcStep = calcStep;
			boolean skipPrint = true;
			double printThreshold = 20;
			boolean completed = false;
			System.out.printf("pos:  %7.3fN, %8.3fE, [%5d, %5d]\n", 90.0-(double)(start_yt)*360.0/calcWidth, (double)(start_xt)*360.0/calcWidth-180.0, start_xt, start_yt);
			double sum = distsum(start_xt*piCalcRes, start_yt*piCalcRes, mapStep, Color.green);
			path.push(new MPoint(start_xt, start_yt, sum));
			Date timeStarted = new Date();
			while(!completed) {
				iter++;
				MPoint pp = path.peek();
				if (prevMapStep != mapStep) {
					// recalculation needed due to small numerical result differences at other mapstep size
					pp.sum = distsum(pp.x*piCalcRes, pp.y*piCalcRes, mapStep, Color.orange);
				}
				double mapRes = 360.0/mapWidth*calcStep;
				if (mapStep < printThreshold) {
					String sumStr = "";
					if (mapStep <= minMapStep) {
						sumStr = String.format(", %.4fkm", pp.sum);
					}
					System.out.printf("pos:  %7.3fN, %8.3fE, [%5d, %5d], %6.3f\u00B0 mapres, %.3f\u00B0 calcres%s\n",
							90.0-(double)(pp.y)*360.0/calcWidth, (double)(pp.x)*360.0/calcWidth-180.0, pp.x, pp.y, mapRes, 360.0/calcWidth*calcStep, sumStr);
				}
				MPoint np = nextPoint(pp, mapStep, calcStep);
				np.sum = distsum(np.x*piCalcRes, np.y*piCalcRes, mapStep, Color.orange);
				if(np.x != 0 || np.y != 0) {
					drawPath.add(new PathSegment(pp.y*calcScaleScreen, pp.x*calcScaleScreen, np.y*calcScaleScreen, np.x*calcScaleScreen, np.sum - pp.sum));
				}
				if(np.sum <= pp.sum) {
					// reached smaller distance sum
					path.push(new MPoint(np.x, np.y, np.sum));
					graphics.setColor(Color.green);
					graphics.fillOval((int)Math.round(np.x*calcScaleScreen) - 5,  (int)Math.round(np.y*calcScaleScreen) - 5,  10,  10);
					frame.repaint();
				} else {
					// landed at larger distance sum
					if (mapStep > 2*minMapStep) {
						prevMapStep = mapStep;
						mapStep = mapStep / 2;
					} else if (mapStep > minMapStep) {
						prevMapStep = mapStep;
						mapStep = mapStep - 1;
					} else if (calcStep > 4) {
						prevCalcStep = calcStep;
						calcStep = Math.max(3, calcStep / 2);
					} else if (calcStep > 1) {
						prevCalcStep = calcStep;
						calcStep--;
					} else {
						// previous iteration already reached minimum final calculation step size
						completed = true;
					}
				}
				graphics.setBackground(new Color(255, 255, 255, 0));
				graphics.clearRect(0, 0, displayWidth, displayHeight);
				for(PathSegment ps : drawPath) {
					if(ps.gradient > 0) {
						graphics.setColor(Color.orange);
					} else {
						graphics.setColor(Color.green);
					}
					graphics.drawLine((int)Math.round(ps.sLon), (int)Math.round(ps.sLat),  (int)Math.round(ps.eLon),  (int)Math.round(ps.eLat));
					graphics.setColor(Color.white);
					graphics.drawLine((int)Math.round(ps.sLon), (int)Math.round(ps.sLat),  (int)Math.round(ps.sLon),  (int)Math.round(ps.sLat));
				}
				frame.repaint();
				if (mapStep < printThreshold) {
					if (skipPrint) {
						System.out.println();
						skipPrint = false;
					}
					System.out.printf("%s: %7.3fN, %8.3fE, [%5d, %5d]\n",
							(np.sum < pp.sum) ? "next" : "skip",
									90.0-(double)(np.y)*360.0/calcWidth, (double)(np.x)*360.0/calcWidth-180.0, np.x, np.y);
				} else {
					System.out.print(".");
				}
			}
			Date timeFinished = new Date();
			MPoint c = path.peek();
			System.out.printf("\ngeographic center%s: %.3f\u00B0N, %.3f\u00B0E, %s%dm sealevel, %.3fkm average distance\n",
					squareMode?" of distance squares":"",
							90.0-(double)(c.y)*360.0/calcWidth, (double)(c.x)*360.0/calcWidth-180.0, levelSign, seaLevel, c.sum);
			System.out.printf("parameters: %.0farcsec (%.2fkm) map resolution, %.0farcsec (%.2fkm) calculation resolution, %dm antarctic threshold\n",
					360.0*60*60/mapWidth, earthCircumference/mapWidth, 360.0*60*60/calcWidth, earthCircumference/calcWidth, antarcticThreshold);
			System.out.printf("completed: %ds, %d iterations\n",
					Math.round((timeFinished.getTime() - timeStarted.getTime())/1000),
					iter);
		}

		MPoint nextPoint(MPoint currPoint, int mapStep, int step) {
			// Calculates the next map point for the gradient search
			// to check for the minimum sum and returns that next map point.
			if(step == 1) {
				return nextPoint(currPoint.x, currPoint.y);
			}
			// gradient calculation
			// select 3 probe locations around current location for higher precision than just 2
			MPoint gradtest[] = new MPoint[3];
			double dw = 2*Math.PI / gradtest.length;
			MPoint next;
			IntStream.range(0, gradtest.length).parallel().forEach(i -> {
				// parallel calculation of the sample triangle
				double w = 2*Math.PI/24 + i*dw;
				gradtest[i] = new MPoint();
				gradtest[i].x = (int)(currPoint.x + step/2*Math.sin(w));
				gradtest[i].y = (int)(currPoint.y - step/2*Math.cos(w));
				gradtest[i].sum = distsum(gradtest[i].x*piCalcRes, gradtest[i].y*piCalcRes, mapStep, Color.lightGray);
				//System.out.printf("gradtest[%d]: [%d, %d] %.3fkm\n", i, gradtest[i].x, gradtest[i].y, gradtest[i].sum);
			});
			// calculate gradient direction (steepest descent direction to minimum)
			MVector u = new MVector();
			u.x = gradtest[1].x - gradtest[0].x;
			u.y = gradtest[1].y - gradtest[0].y;
			u.sum = gradtest[1].sum - gradtest[0].sum;
			MVector v = new MVector();
			v.x = gradtest[2].x - gradtest[0].x;
			v.y = gradtest[2].y - gradtest[0].y;
			v.sum = gradtest[2].sum - gradtest[0].sum;
			MVector n = new MVector();
			n.x = u.y * v.sum - u.sum * v.y;
			n.y = u.sum * v.x - u.x * v.sum;
			n.sum = u.x * v.y - u.y * v.x;
			double len = Math.sqrt(n.x*n.x + n.y*n.y);
			next = new MPoint();
			next.x = (int)Math.round(currPoint.x + step * n.x/len);
			next.y = (int)Math.round(currPoint.y + step * n.y/len);
			//double direction = 90.0 - Math.atan2(-n.y, n.x) * 180.0 / Math.PI;
			//if (direction < 0.0) {
			//	direction = 360 + direction;
			//}
			//System.out.printf("grad: [%.2f, %.2f] %d°\n", n.x, n.y, Math.round(direction));
			return next;
		}

		MPoint nextPoint(int xt, int yt) {
			// Calculates the next map point for the final single-step neighbor search
			// to check for the minimum sum.
			// Returns the point with the minimum sum if that is smaller than the current sum,
			// else return infinite sum.
			MPoint test[] = new MPoint[8];
			test[0] = new MPoint(xt-1, yt-1, Integer.MAX_VALUE);
			test[1] = new MPoint(xt, yt-1, Integer.MAX_VALUE);
			test[2] = new MPoint(xt+1, yt-1, Integer.MAX_VALUE);
			test[3] = new MPoint(xt+1, yt, Integer.MAX_VALUE);
			test[4] = new MPoint(xt+1, yt+1, Integer.MAX_VALUE);
			test[5] = new MPoint(xt, yt+1, Integer.MAX_VALUE);
			test[6] = new MPoint(xt-1, yt+1, Integer.MAX_VALUE);		
			test[7] = new MPoint(xt-1, yt, Integer.MAX_VALUE);
			int min=-1;
			double sum=Integer.MAX_VALUE;
			IntStream.range(0, test.length).parallel().forEach(i -> {
				test[i].sum = distsum(test[i].x*piCalcRes, test[i].y*piCalcRes, 1, Color.lightGray);
			});
			for(int i=0; i<test.length; i++) {
				System.out.printf("  test[%d]: [%d, %d] %.4fkm\n", i, test[i].x, test[i].y, test[i].sum);
				if(test[i].sum < sum) {
					min = i;
					sum = test[i].sum;
				}
			}
			if (min >= 0) {
				return test[min];
			} else {
				return new MPoint(0,0,Integer.MAX_VALUE);
			}
		}

		double distsum(double xtrad, double ytrad, int st, Color color) {
			// Calculates the great-circle distance sum from this map point
			// to all other points on land and displays the calculation activity
			// Returns not the plain distance sum, but the average sum
			// to be able to approximately visually compare calculations
			// at different step sizes and scale it to km on Earth.
			// As precise comparison is not possible between different step sizes
			// the distance sum is recalculated in case the steps size changes. 
			graphics.setColor(color);
			graphics.drawOval((int)Math.round((xtrad/piCalcRes - st/2)*calcScaleScreen), (int)Math.round((ytrad/piCalcRes - st/2)*calcScaleScreen),  (int)Math.round(st*calcScaleScreen), (int)Math.round(st*calcScaleScreen));
			frame.repaint();
			//printf("\rxtrad=%3.0f<B0> ytrad=%3.0f<B0> xt=%5d yt=%5d", xtrad*180/PI, ytrad*180/PI, xt, yt);
			double[] total = IntStream.range(0, geoTiffHeight)
					.filter(y -> y % st == 0)
					.parallel()
					.mapToObj(y -> {
						double yrad = y * piWidth;
						double localSum = 0;
						double localCount = 0;
						for (int x = 0; x < geoTiffWidth; x += st) {
							if (getElevation(x, y) > seaLevel) {
								double xrad = x * piWidth;
								// great circle distance between locations (yrad, xrad) and (ytrad,xtrad)
								// with coordinates in radiant, origin at north pole, zero meridian at 180E
								double d = Math.acos(Math.cos(yrad) * Math.cos(ytrad) + Math.sin(yrad) * Math.sin(ytrad) * Math.cos(xrad - xtrad));
								// calculate sample point unit area on cylinder projection with equidistant linear latitude and longitude
								double area = Math.sin(yrad);
								if (squareMode) {
									// optional mode for sum of distance square calculation
									d = Math.pow(d, 2);
								}
								localSum += d * area;
								localCount += area;
							}
						}
						return new double[]{localSum, localCount};
					})
					.reduce(new double[]{0, 0}, (a, b) -> new double[]{a[0] + b[0], a[1] + b[1]});
			double sum = total[0];
			double count = total[1];
			graphics.setColor(color);
			graphics.fillOval((int)Math.round((xtrad/piCalcRes - st/2)*calcScaleScreen), (int)Math.round((ytrad/piCalcRes - st/2)*calcScaleScreen), (int)Math.round(st*calcScaleScreen), (int)Math.round(st*calcScaleScreen));
			frame.repaint();
			sum = sum * earthCircumference / count / (2*Math.PI);
			return sum;
		}

		double median(Map<Integer, Double> histogram) {
			// statistical median calculation for the statistical map data display
			// not used for the actual center calculation
			double n = (count(histogram) - 1d)/2 + 1;
			double d = n - Math.floor(n);
			SortedSet<Integer> bins = new ConcurrentSkipListSet<Integer>(histogram.keySet());
			int observationsBelowBinInclusive = 0;
			Integer lowBin = bins.first();
			Double valuePercentile = null;
			for (Integer highBin : bins)
			{
				observationsBelowBinInclusive += histogram.get(highBin);
				if (n <= observationsBelowBinInclusive)
				{
					if ((d == 0f) || (histogram.get(highBin) > 1L))
					{
						lowBin = highBin;
					}
					valuePercentile = lowBin.doubleValue() + ((highBin - lowBin) * d);
					break;
				}
				lowBin = highBin;
			}
			return valuePercentile;
		}

		long count(Map<Integer, Double> histogram) {
			// helper for the statistical map data display
			// not used for the actual center calculation
			long observations = 0L;
			for (Integer value : histogram.keySet())
			{
				observations += histogram.get(value);
			}
			return observations;
		}
	}
}
