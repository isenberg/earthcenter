/* Kopie der Dateiversion Oktober 2003 von https://mars-news.de/pyramids/geocenter.c */
/*
 * Berechnung des Schwerpunkts aller Landflächen der Erde
 * (geographisches Zentrum)
 *
 * Der Schwerpunkt sei hier die Punktemenge,
 * deren Summe aller Entfernungen zu allen Punkten
 * auf Landflächen minimal ist.
 * Die Entfernung wird dabei ueber den sphärischen
 * Grosskreis bestimmt.
 *
 * Algorithmus: Bergsteigen (hillclimbing)
 *
 * Autor:
 * Holger Isenberg 
 * <a href="mailto:web@mars-news.de">web@mars-news.de</a>
 * http://mars-news.de
 *
 * Data-Source:
 * NOAA Satellites & Data
 * ETOPO2 (Global 2-min gridded data)
 * http://www.ngdc.noaa.gov/mgg/global/global.html
 *
 * GIF-Reader from giftopnm.c by David Koblas (koblas@netcom.com)
 *
 */

#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <values.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI M_PI
#define MAXCOLORMAPSIZE  256

#define STR_NORESULT "NO result yet, calculation still running!\n"

#define CM_RED           0
#define CM_GREEN         1
#define CM_BLUE          2
#define TRUE 1
#define FALSE 0

#define MAX_LWZ_BITS     12

#define INTERLACE          0x40
#define LOCALCOLORMAP      0x80
#define BitSet(byte, bit)  (((byte) & (bit)) == (bit))

#define ReadOK(file,buffer,len) (fread(buffer, len, 1, file) != 0)
#define LM_to_uint(a,b)         (((b)<<8)|(a))

#define GRAYSCALE        1
#define COLOR            2

typedef unsigned char CMap[3][MAXCOLORMAPSIZE];

static struct
{
  unsigned int Width;
  unsigned int Height;
  CMap ColorMap;
  unsigned int BitPixel;
  unsigned int ColorResolution;
  unsigned int Background;
  unsigned int AspectRatio;
  /*
   **
   */
  int GrayScale;
} GifScreen;

static struct
{
  int transparent;
  int delayTime;
  int inputFlag;
  int disposal;
} Gif89 = { -1, -1, -1, 0 };

int  highest_used_index = 0;
static unsigned char   used_cmap[3][256];
static unsigned char   gimp_cmap[768];
int ZeroDataBlock = FALSE;

static int ReadColorMap (FILE *, int, CMap, int *);
static int DoExtension (FILE *, int);
static int GetDataBlock (FILE *, unsigned char *);
static int GetCode (FILE *, int, int);
static int LWZReadByte (FILE *, int, int);
unsigned char* ReadImage (FILE *fd, char *filename, int   len, int   height, CMap  cmap,
	   int   ncols, int   format, int   interlace, int   number, unsigned int   leftpos,
	   unsigned int   toppos,  unsigned int screenwidth, unsigned int screenheight);
short int* loadGIF(char* filename, int* width, int* height);
short int* loadDEM(char* filename, int width, int height);

void usage(char* arg0)
{
	printf("usage: %s -i inputfile -o outputfile\n", arg0);
	printf("\t\t-n northern latitude -e eastern longitude\n");
	printf("\t\t-r result-width -l water-level\n\n");
	printf("Inputfile is a 16bit DEM in host byte order or a grayscale GIF.\n");
	printf("Map projection type should be linear.\n");
	printf("Outputfile is the filename-prefix for the result-files.\n");
	printf("Result-width is the width of the output PGM-file.\n");
	printf("Lat and long, in decimal degrees, is the start position.\n");
	printf("Water-level is the gobal sea level in DEM-units [ETOPO2: m].\n");
}

int
main(int argc, char* argv[])
{
	double xtrad, ytrad,
		sum, max_sum, min_sum, min, old_sum, mass,
		xs, ys,
		piWidth, piRes,
		xrad, yrad,
		d;
	long x,y,xt,yt,xc,yc,xstep,ystep,step,i,w=0;
	double area_x=31.0, area_y=30.0;
	double area_xt, area_yt, ratio_l, ratio_w, max_x, max_y;
	int calc;
	long fx,fy;
	long lastx, lasty;
	char c, newcalc,newy;
	FILE* outfile;
	int fd;
	struct stat statbuf;
	char* inputfile=NULL;
	char* outputprefix=NULL;
	char outputfile[3][1024];
	char buffer[10000];
	double* result;
	short int* land;
	int level=0, res=1080, land_h=0, land_w=0,map,demsize;

	while (1)
	{
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		
		c = getopt(argc, argv, "i:w:e:n:s:o:r:l:h");
		if (c == -1)
			break; 
		
		switch (c)    
		{
		    case 0:
			    printf ("option");
			    if (optarg)
				    printf (" with arg %s", optarg);
			    printf ("0");
			    break;
			    break;
			    
		    case 'h':
			    usage(argv[0]);
			    exit(1);
			    break;

		    case 'e':
			    area_x=atof(optarg);
			    break;

		    case 'n':
			    area_y=atof(optarg);
			    break;

		    case 'i':
			    inputfile=strdup(optarg);
			    break;

		    case 'w':
			    land_w=atoi(optarg);
			    break;

		    case 'o':
			    outputprefix = strdup(optarg);
			    break;
			    
		    case 'r':
			    res=atoi(optarg);
			    break;

		    case 'l':
			    level=atoi(optarg);
			    break;
			    
		    default:
			    usage(argv[0]);
			    exit(1);
		}
	}
	
	if(res<=0 || inputfile==0 || !outputprefix) {
		usage(argv[0]);
		exit(1);
	}

	sprintf(outputfile[0], "%s_%d_%dm.pgm", outputprefix, res, level);
	sprintf(outputfile[1], "%s_%d_%dm.txt", outputprefix, res, level);
	sprintf(outputfile[2], "%s_%d_%dm_area.pgm", outputprefix, res, level);

	if(strstr(inputfile, ".gif") || strstr(inputfile, ".GIF")) {
		land = loadGIF(inputfile, &land_w, &land_h);
	} else {
		fd=open(inputfile, 0);
		if (fd<0) {
			fprintf(stderr, "can't open \"%s\"\n", inputfile);
			exit(1);
		}
		fstat(fd, &statbuf);
		demsize = statbuf.st_size;
		land_w=(int)sqrt(demsize);
		land_h=land_w/2;
		printf("DEM: %s %dx%d\n", inputfile, land_w, land_h);
		if( (land = (short int*)mmap(NULL, demsize, PROT_READ, MAP_FILE | MAP_SHARED, fd, 0)) <0) {
			fprintf(stderr, "can't mmap() %d bytes!\n", demsize);
			exit(1);
		}				
	}
	if(!land) {
		fprintf(stderr, "could not load GIF\n");
		exit(1);
	}
	if((result=(double*)malloc(res*res/2*sizeof(double)))==NULL) {
		printf("cannot allocate %ld bytes\n", res*res/2*sizeof(double));
		exit(1);
	}

	//step=land_w/res;
	step=1;
	piWidth=2*M_PI/land_w;
	piRes=2*M_PI/res;
	old_sum=MAXDOUBLE;
	min_sum=MAXDOUBLE-1.0;
	max_sum=0.0;
	ystep=step;
	xstep=step;
	area_xt=res*(180.0+area_x)/360.0;
	area_yt=res*(90.0-area_y)/360.0;

	printf("sea level: %dm\noutput map: %s %dx%d ", level, outputfile[0], res, res/2);
	fflush(0);
	for(y=0,ratio_w=.0,ratio_l=.0; y<land_h; y++) {
		yrad=(double)(y)*piWidth;
		for(x=0; x<land_w; x++)
			if(land[x+y*land_w]>level) 
				ratio_l+=sin(yrad);
			else
				ratio_w+=sin(yrad);
	}
	printf("%d%% land/sea ratio\n",(int)(ratio_l/(ratio_l+ratio_w)*100.0));
	printf("start: %.4fN %.4fE\n", area_y, area_x);

	// Weltkarte mit angegebenen Meeresspiegel speichern
	if(!(outfile=fopen(outputfile[0], "w"))) {
		fprintf(stderr, "cannot write to %s", outputfile[0]);
		exit(2);
	}
	fprintf(outfile, "P2\n");
	fprintf(outfile, "# Earth with sea level %s%dm, land/water ratio %d%%\n",
		level>=0?"+":"", level, (int)(ratio_l/(ratio_l+ratio_w)*100.0));
	fprintf(outfile, "# global elevation data ETOPO2 from NOAA, http://www.ngdc.noaa.gov/mgg/\n");
	fprintf(outfile, "%d %d\n255\n", res, res/2);
	for(yt=0; yt<res/2; yt++) {
		for(xt=0; xt<res; xt++) {
			result[xt+yt*res]=-1.0;
			if(land[xt*land_w/res+yt*land_w/res*land_w]>level) {
				fprintf(outfile, "%4d", 255);
			} else {
				fprintf(outfile, "%4d", 0);
			}
		}
		fprintf(outfile, "\n");
	}
	fclose(outfile);
		
	xc=area_xt;
	yc=area_yt;
	do {
		old_sum=min_sum;
		area_xt = xc;
		area_yt = yc;
		for(yt=area_yt-1; yt<=(area_yt+1); yt++) {
			for(xt=area_xt-1, ytrad=yt*piRes; xt<=(area_xt+1); xt++) {
				if (result[xt+yt*res]<.0) {
					printf("%8.4fN %8.4fE: ", 90.0-(double)(yt)*360.0/res, (double)(xt)*360.0/res-180.0);
					fflush(0);
					for(y=0, sum=.0, mass=.0, xtrad=xt*piRes; y<land_h; y++)
						for(x=0, yrad=y*piWidth; x<land_w; x++)
							if(land[x+y*land_w]>level) {
								xrad=x*piWidth;
								
								/* Entfernung zwischen den Punkten (yrad, xrad) und (ytrad,xtrad)
								 * über den sphärischen Großkreis. Die Koordinaten sind im Bogenmaß
								 * angegeben: Nullpunkt im Nordpol, Null-Meridian ist bei 180W.
								 */
								d=acos(cos(yrad)*cos(ytrad)+sin(yrad)*sin(ytrad)*cos(xrad-xtrad));
								
								/* Transformation der Zylinderprojektion
								 * (die verwendeten Daten stammen aus einer Zyl.projektion
								 * mit äquidistanten Längen- und Breitengeraden)
								 */
								sum+=d*sin(yrad);
								mass+=sin(yrad);
							}
					result[xt+yt*res]=sum;
					printf("%.3fkm %.3fradsum", (sum/mass)*40075.0/(2*PI), sum);
					if(sum>.0 && sum<min_sum) {
						min_sum=sum;
						xc=xt;
						yc=yt;
						printf(" -\n");
					} else if(sum==min_sum) {
						min_sum=sum;
						xc=xt;
						yc=yt;
						printf(" ==\n");
					} else if (sum>max_sum) {
						max_sum = sum;
						printf(" +\n");
					} else {
						printf(" +\n");
					}

					if(min_sum==sum) {						
						sprintf(buffer, "Description geographic center of all land areas on Earth calculated with geocenter\n\
Author Holger Isenberg, H.Isenberg@ping.de, http://mars-news.de\n\
\"input data\" global elevation data ETOPO2 from NOAA, http://www.ngdc.noaa.gov/mgg/\n\
\"land/water ratio\" %d%%\n\
\"mean distance\" %.3fkm\n\
\"global sea level\" %s%dm\n\
\"geographic pixel center\" x=%dpx, y=%dpx, map: %dx%d\n\
resolution %.3fkm, %.4fdeg\n\
\"geographic center\" %.4f N, %.4f E\n\
%s",
							(int)(ratio_l/(ratio_l+ratio_w)*100.0),
							(sum/mass)*40075.0/(2*PI),
							level>=0?"+":"", level,
							xc, yc, res, res/2,
							40075.0/res, 360.0/res,
							90.0-(double)(yc)*360.0/res, (double)(xc)*360.0/res-180.0,
							min_sum<old_sum?STR_NORESULT:"");
						if(!(outfile=fopen(outputfile[1], "w"))) {
							fprintf(stderr, "cannot write to %s", outputfile[1]);
							exit(2);
						}
						fprintf(outfile, buffer);
						fclose(outfile);
					}
				}
			}
		}
	} while (min_sum < old_sum);
	if(!(outfile=fopen(outputfile[1], "w"))) {
		fprintf(stderr, "cannot write to %s", outputfile[1]);
		exit(2);
	}
	buffer[strlen(buffer)-strlen(STR_NORESULT)]='\0';
	fprintf(outfile, buffer);
	fclose(outfile);
	printf("\ngeographic center:\n %dx, %dy [pixel]\n %.4fN, %.4fE\n",
	       xc, yc,
	       90.0-(double)(yc)*360.0/res, (double)(xc)*360.0/res-180.0);
	if(fd>0) {
		munmap(land, demsize);
		close(fd);
	}
	return 0;
}

short int*
loadGIF(char* filename, int* width, int* height)
{
  FILE *fd;
  char * name_buf;
  unsigned char buf[16];
  unsigned char c;
  CMap localColorMap;
  int grayScale;
  int useGlobalColormap;
  int bitPixel;
  int imageCount = 0;
  char version[4];
  unsigned char* pixmap;
  short int* dem;
  long i;

  fd = fopen (filename, "rb");
  if (!fd) {
	  fprintf(stderr, "can't open \"%s\"\n", filename);
	  return 0;
  }

  name_buf = (char*)malloc (strlen (filename) + 11);

  if (!ReadOK (fd, buf, 6)) {
	  fprintf(stderr, "error reading magic number\n");
	  return 0;
  }

  if (strncmp ((char *) buf, "GIF", 3) != 0) {
	  fprintf(stderr, "not a GIF file\n");
	  return 0;
    }

  strncpy (version, (char *) buf + 3, 3);
  version[3] = '\0';

  if ((strcmp (version, "87a") != 0) && (strcmp (version, "89a") != 0)) {
	  fprintf(stderr, "bad version number, not '87a' or '89a'\n");
	  return 0;
  }
  
  if (!ReadOK (fd, buf, 7)) {
	  fprintf(stderr, "failed to read screen descriptor\n");
	  return 0;
  }

  GifScreen.Width = LM_to_uint (buf[0], buf[1]);
  GifScreen.Height = LM_to_uint (buf[2], buf[3]);
  GifScreen.BitPixel = 2 << (buf[4] & 0x07);
  GifScreen.ColorResolution = (((buf[4] & 0x70) >> 3) + 1);
  GifScreen.Background = buf[5];
  GifScreen.AspectRatio = buf[6];

  *width =   GifScreen.Width;
  *height = GifScreen.Height;

  if (BitSet (buf[4], LOCALCOLORMAP)) {
	  /* Global Colormap */
	  if (ReadColorMap (fd, GifScreen.BitPixel, GifScreen.ColorMap, &GifScreen.GrayScale)) {
		  fprintf(stderr, "error reading global colormap\n");
		  return 0;
	  }
  }

  if (GifScreen.AspectRatio != 0 && GifScreen.AspectRatio != 49)
    {
      fprintf(stderr, "GIF: warning - non-square pixels\n");
    }


  highest_used_index = 0;
  pixmap=0;

  for(;;) {
  if (!ReadOK (fd, &c, 1))
  {
	  fprintf(stderr, "GIF: EOF / read error on image data\n");
	  return 0;
  }
  
  if (c == ';')
  {
	  /* GIF terminator */
	  fprintf(stderr, "GIF: got ;\n");
	  //return pixmap;
	  return 0;
  }

  if (c == '!')
  {
	  /* Extension */
	  if (!ReadOK (fd, &c, 1))
	  {
		  fprintf(stderr, "GIF: OF / read error on extention function code\n");
		  return 0;
	  }
	  DoExtension (fd, c);
	  continue;
  }
  
  if (c != ',')
  {
	  /* Not a valid start character */
	  fprintf(stderr, "GIF: bogus character 0x%02x, ignoring\n", (int) c);
	  continue;
  }
  
  ++imageCount;
  
  if (!ReadOK (fd, buf, 9))
  {
	  fprintf(stderr, "GIF: couldn't read left/top/width/height\n");
	  return 0;
  }
  
  useGlobalColormap = !BitSet (buf[8], LOCALCOLORMAP);
  
  bitPixel = 1 << ((buf[8] & 0x07) + 1);

  pixmap = ReadImage (fd, filename, LM_to_uint (buf[4], buf[5]),
		      LM_to_uint (buf[6], buf[7]),
		      GifScreen.ColorMap, GifScreen.BitPixel,
		      GifScreen.GrayScale,
		      BitSet (buf[8], INTERLACE), imageCount,
		      (unsigned int) LM_to_uint (buf[0], buf[1]),
		      (unsigned int) LM_to_uint (buf[2], buf[3]),
		      GifScreen.Width,
		      GifScreen.Height
		      );
  dem = malloc(GifScreen.Width * GifScreen.Height*2);
  if(dem==NULL)
	  return NULL;
  
  for(i=0; i<GifScreen.Width * GifScreen.Height; i++)
	  dem[i] = pixmap[i] - 128;
  free(pixmap);
  return dem;
  }
}

int
ReadColorMap (FILE *fd,
	      int   number,
	      CMap  buffer,
	      int  *format)
{
  int i;
  unsigned char rgb[3];
  int flag;

  flag = TRUE;

  for (i = 0; i < number; ++i)
    {
      if (!ReadOK (fd, rgb, sizeof (rgb)))
	{
	  fprintf(stderr, "GIF: bad colormap\n");
	  return TRUE;
	}

      buffer[CM_RED][i] = rgb[0];
      buffer[CM_GREEN][i] = rgb[1];
      buffer[CM_BLUE][i] = rgb[2];

      flag &= (rgb[0] == rgb[1] && rgb[1] == rgb[2]);
    }

  *format = (flag) ? GRAYSCALE : COLOR;

  return FALSE;
}

static int
GetDataBlock (FILE          *fd,
	      unsigned char *buf)
{
  unsigned char count;

  if (!ReadOK (fd, &count, 1))
    {
      fprintf(stderr, "GIF: error in getting DataBlock size\n");
      return -1;
    }

  ZeroDataBlock = count == 0;

  if ((count != 0) && (!ReadOK (fd, buf, count)))
    {
      fprintf(stderr, "GIF: error in reading DataBlock\n");
      return -1;
    }

  return count;
}

static int
GetCode (FILE *fd,
	 int   code_size,
	 int   flag)
{
  static unsigned char buf[280];
  static int curbit, lastbit, done, last_byte;
  int i, j, ret;
  unsigned char count;

  if (flag)
    {
      curbit = 0;
      lastbit = 0;
      done = FALSE;
      return 0;
    }

  if ((curbit + code_size) >= lastbit)
    {
      if (done)
	{
	  if (curbit >= lastbit)
	    {
	      fprintf(stderr, "GIF: ran off the end of by bits\n");
	      exit(2);
	    }
	  return -1;
	}
      buf[0] = buf[last_byte - 2];
      buf[1] = buf[last_byte - 1];

      if ((count = GetDataBlock (fd, &buf[2])) == 0)
	done = TRUE;

      last_byte = 2 + count;
      curbit = (curbit - lastbit) + 16;
      lastbit = (2 + count) * 8;
    }

  ret = 0;
  for (i = curbit, j = 0; j < code_size; ++i, ++j)
    ret |= ((buf[i / 8] & (1 << (i % 8))) != 0) << j;

  curbit += code_size;

  return ret;
}

static int
LWZReadByte (FILE *fd,
	     int   flag,
	     int   input_code_size)
{
  static int fresh = FALSE;
  int code, incode;
  static int code_size, set_code_size;
  static int max_code, max_code_size;
  static int firstcode, oldcode;
  static int clear_code, end_code;
  static int table[2][(1 << MAX_LWZ_BITS)];
  static int stack[(1 << (MAX_LWZ_BITS)) * 2], *sp;
  register int i;

  if (flag)
    {
      set_code_size = input_code_size;
      code_size = set_code_size + 1;
      clear_code = 1 << set_code_size;
      end_code = clear_code + 1;
      max_code_size = 2 * clear_code;
      max_code = clear_code + 2;

      GetCode (fd, 0, TRUE);

      fresh = TRUE;

      for (i = 0; i < clear_code; ++i)
	{
	  table[0][i] = 0;
	  table[1][i] = i;
	}
      for (; i < (1 << MAX_LWZ_BITS); ++i)
	table[0][i] = table[1][0] = 0;

      sp = stack;

      return 0;
    }
  else if (fresh)
    {
      fresh = FALSE;
      do
	{
	  firstcode = oldcode =
	    GetCode (fd, code_size, FALSE);
	}
      while (firstcode == clear_code);
      return firstcode;
    }

  if (sp > stack)
    return *--sp;

  while ((code = GetCode (fd, code_size, FALSE)) >= 0)
    {
      if (code == clear_code)
	{
	  for (i = 0; i < clear_code; ++i)
	    {
	      table[0][i] = 0;
	      table[1][i] = i;
	    }
	  for (; i < (1 << MAX_LWZ_BITS); ++i)
	    table[0][i] = table[1][i] = 0;
	  code_size = set_code_size + 1;
	  max_code_size = 2 * clear_code;
	  max_code = clear_code + 2;
	  sp = stack;
	  firstcode = oldcode =
	    GetCode (fd, code_size, FALSE);
	  return firstcode;
	}
      else if (code == end_code)
	{
	  int count;
	  unsigned char buf[260];

	  if (ZeroDataBlock)
	    return -2;

	  while ((count = GetDataBlock (fd, buf)) > 0)
	    ;

	  if (count != 0)
	    fprintf (stderr, "GIF: missing EOD in data stream (common occurence)");
	  return -2;
	}

      incode = code;

      if (code >= max_code)
	{
	  *sp++ = firstcode;
	  code = oldcode;
	}

      while (code >= clear_code)
	{
	  *sp++ = table[1][code];
	  if (code == table[0][code])
	    {
	      fprintf(stderr, "GIF: circular table entry BIG ERROR\n");
	      exit(2);
	    }
	  code = table[0][code];
	}

      *sp++ = firstcode = table[1][code];

      if ((code = max_code) < (1 << MAX_LWZ_BITS))
	{
	  table[0][code] = oldcode;
	  table[1][code] = firstcode;
	  ++max_code;
	  if ((max_code >= max_code_size) &&
	      (max_code_size < (1 << MAX_LWZ_BITS)))
	    {
	      max_code_size *= 2;
	      ++code_size;
	    }
	}

      oldcode = incode;

      if (sp > stack)
	return *--sp;
    }
  return code;
}

unsigned char *
ReadImage (FILE *fd,
	   char *filename,
	   int   len,
	   int   height,
	   CMap  cmap,
	   int   ncols,
	   int   format,
	   int   interlace,
	   int   number,
	   unsigned int   leftpos,
	   unsigned int   toppos,
	   unsigned int screenwidth,
	   unsigned int screenheight)
{
  static int frame_number = 1;

  unsigned char *dest, *temp;
  unsigned char c;
  int xpos = 0, ypos = 0, pass = 0;
  int cur_progress, max_progress;
  int v;
  int i, j;
  char framename[200]; /* FIXME */
  char alpha_frame = FALSE;
  int nreturn_vals;
  static int previous_disposal;

  /*
   **  Initialize the Compression routines
   */
  if (!ReadOK (fd, &c, 1))
    {
      fprintf(stderr, "GIF: EOF / read error on image data\n");
      return 0;
    }

  if (LWZReadByte (fd, TRUE, c) < 0)
    {
      fprintf(stderr, "GIF: error while reading\n");
      return 0;
    }

  for (i = 0, j = 0; i < ncols; i++)
	{
		used_cmap[0][i] = gimp_cmap[j++] = cmap[0][i];
		used_cmap[1][i] = gimp_cmap[j++] = cmap[1][i];
		used_cmap[2][i] = gimp_cmap[j++] = cmap[2][i];
	}
  
  if (Gif89.delayTime < 0)
	  strcpy(framename, "Background");
  else
	  sprintf(framename, "Background (%dms)", 10*Gif89.delayTime);
  
  previous_disposal = Gif89.disposal;
  
  if (Gif89.delayTime < 0)
	  sprintf(framename, "Frame %d", frame_number);
  else
	  sprintf(framename, "Frame %d (%dms)",
		  frame_number, 10*Gif89.delayTime);
  
  switch (previous_disposal)
  {
      case 0x00: break; /* 'don't care' */
      case 0x01: strcat(framename," (combine)"); break;
      case 0x02: strcat(framename," (replace)"); break;
      case 0x03: strcat(framename," (combine)"); break;
      case 0x04:
      case 0x05:
      case 0x06:
      case 0x07:
	      strcat(framename," (unknown disposal)");
	      fprintf(stderr, "GIF: Hmm... please forward this GIF to the "
		      "GIF plugin author!\n  (adam@foxbox.org)\n");
	      break;
      default: fprintf(stderr, "GIF: Something got corrupted.\n"); break;
  }
  
  previous_disposal = Gif89.disposal;
  
  dest = (unsigned char *) malloc (len * height);
  if(!dest) {
	  fprintf(stderr, "cannot allocate %d bytes\n", len*height);
	  return 0;
  }
  fprintf (stderr, "GIF: reading %d by %d%s GIF image, ncols=%d\n",
	   len, height, interlace ? " interlaced" : "", ncols);

  while ((v = LWZReadByte (fd, FALSE, c)) >= 0) {
	  if (((unsigned char)v > highest_used_index) && !(v == Gif89.transparent))
		  highest_used_index = (unsigned char)v;
	  
	  temp = dest + (ypos * len) + xpos;
	  *(temp  ) = (unsigned char) cmap[0][v];
	  ++xpos;
	  if (xpos == len)
	  {
		  xpos = 0;
		  if (interlace)
		  {
			  switch (pass)
			  {
			      case 0:
			      case 1:
				      ypos += 8;
				      break;
			      case 2:
				      ypos += 4;
				      break;
			      case 3:
				      ypos += 2;
				      break;
			  }

			  if (ypos >= height)
			  {
				  ++pass;
				  switch (pass)
				  {
				      case 1:
					      ypos = 4;
					      break;
				      case 2:
					      ypos = 2;
					      break;
				      case 3:
					      ypos = 1;
					      break;
				      default:
					      goto fini;
				  }
			  }
		  }
		  else
		  {
			  ++ypos;
		  }
	  }
	  if (ypos >= height)
		  break;
  }

  fini:
  if (LWZReadByte (fd, FALSE, c) >= 0)
	  fprintf(stderr, "GIF: too much input data, ignoring extra...\n");

  return dest;
}

static int
DoExtension (FILE *fd,
	     int   label)
{
  static unsigned char buf[256];
  char *str;

  switch (label)
    {
    case 0x01:			/* Plain Text Extension */
      str = "Plain Text Extension";
#ifdef notdef
      if (GetDataBlock (fd, (unsigned char *) buf) == 0)
	;

      lpos = LM_to_uint (buf[0], buf[1]);
      tpos = LM_to_uint (buf[2], buf[3]);
      width = LM_to_uint (buf[4], buf[5]);
      height = LM_to_uint (buf[6], buf[7]);
      cellw = buf[8];
      cellh = buf[9];
      foreground = buf[10];
      background = buf[11];

      while (GetDataBlock (fd, (unsigned char *) buf) != 0)
	{
	  PPM_ASSIGN (image[ypos][xpos],
		      cmap[CM_RED][v],
		      cmap[CM_GREEN][v],
		      cmap[CM_BLUE][v]);
	  ++index;
	}

      return FALSE;
#else
      break;
#endif
    case 0xff:			/* Application Extension */
      str = "Application Extension";
      break;
    case 0xfe:			/* Comment Extension */
      str = "Comment Extension";
      while (GetDataBlock (fd, (unsigned char *) buf) != 0) {}
      return FALSE;
    case 0xf9:			/* Graphic Control Extension */
      str = "Graphic Control Extension";
      (void) GetDataBlock (fd, (unsigned char *) buf);
      Gif89.disposal = (buf[0] >> 2) & 0x7;
      Gif89.inputFlag = (buf[0] >> 1) & 0x1;
      Gif89.delayTime = LM_to_uint (buf[1], buf[2]);
      if ((buf[0] & 0x1) != 0)
	Gif89.transparent = buf[3];
      else
	Gif89.transparent = -1;

      while (GetDataBlock (fd, (unsigned char *) buf) != 0)
	;
      return FALSE;
    default:
      str = (char *)buf;
      sprintf ((char *)buf, "UNKNOWN (0x%02x)", label);
      break;
    }

  fprintf(stderr, "GIF: got a '%s' extension\n", str);

  while (GetDataBlock (fd, (unsigned char *) buf) != 0)
    ;

  return FALSE;
}

