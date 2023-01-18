/*
	ECE 4310
	Lab 8
	Roderick "Rance" White
	
	This program takes in an image and segments an object in it.
	
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SQR(x) ((x)*(x))
#define MAX_QUEUE 10000	/* max perimeter size (pixels) of border wavefront */
#define WINDOW_SIZE 5
#define ANGLE_THRESH 1.04
#define IMAGE_THRESHOLD 150
#define PIXEL_DISTANCE 3


/* ----------------------------------------------------------------------------------- */
/* Basic Image Functions */

/* This function serves to read and open the image, given it's name. 
 * It will return the values within the file, the number of rows, the number of columns, and
 * the number of bytes.
 */
unsigned char *Image_Read(char *FileName, char *header, int *r, int *c, int *b)
{
	FILE						*fpt;
	unsigned char		*image;
	int							ROWS,COLS,BYTES;
	
	/* read image */
	if ((fpt=fopen(FileName,"rb")) == NULL)
	{
		printf("Unable to open %s for reading\n", FileName);
		exit(0);
	}

	/* read image header (simple 8-bit greyscale PPM only) */
	fscanf(fpt,"%s %d %d %d",header,&COLS,&ROWS,&BYTES);	
	if (strcmp(header,"P5") != 0 || BYTES != 255)
	{
		printf("Not a greyscale 8-bit PPM image\n");
		exit(0);
	}
	
	/* allocate dynamic memory for image */
	image=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
	header[0]=fgetc(fpt);	/* read white-space character that separates header */
	fread(image,1,COLS*ROWS,fpt);
	fclose(fpt);
	
	*r = ROWS;		//Return number of rows for the image
	*c = COLS;		//Return number of columns for the image
	*b = BYTES;
	
	return image;
}

/* This function serves to free the array of image pixels 
 * Despite being a very simple function to perform, I like the readability
 */
void Image_Free(unsigned char *image)
{
	free(image);
}

/* This function serves to write the image to the appropriate file for more concise code */
void Image_Write(unsigned char *image, char *FileName, int ROWS, int COLS)
{
	FILE						*fpt;
	fpt=fopen(FileName,"w");
	fprintf(fpt,"P5 %d %d 255\n",COLS,ROWS);
	fwrite(image,COLS*ROWS,1,fpt);
	fclose(fpt);
}


/* Function to threshold an image 
 * It will leave behind all parts of the image that exceed the threshold, setting all others to 0
 */
unsigned char *Threshold_Image(unsigned char *image, int ROWS, int COLS)
{
	unsigned char *Matches;
	int r, c;
	Matches=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
	
	//Check all pixels to see which ones exceed the threshold
	for(r=0; r<ROWS; r++)
	{
		for(c=0; c<COLS; c++)
		{
			if(image[r*COLS+c]<=IMAGE_THRESHOLD)
				Matches[r*ROWS+c] = image[r*COLS+c];
			else
				Matches[r*COLS+c] = 255;
		}
	}
	return Matches;
}

/* ----------------------------------------------------------------------------------- */
/* Basic Coordinate Functions */


/* This function allocates space for 3D coordinates based on the number of rows and columns
 */
double **Coordinate_Allocate(int ROWS, int COLS)
{
	double **Coordinates;

	/* Allocate space for the coordinates */
	//Allocate space for the X, Y, and Z
	Coordinates =(double **)calloc(3,sizeof(double *));

	//Allocate the space for each individual coordinate value
	Coordinates[0]=(double *)calloc(ROWS*COLS,sizeof(double));		//X coordinate
	Coordinates[1]=(double *)calloc(ROWS*COLS,sizeof(double));		//Y coordinate
	Coordinates[2]=(double *)calloc(ROWS*COLS,sizeof(double));		//Z coordinate

	return Coordinates;
}


/* This function frees the allocated space for the 3D coordinates.
 */
void Coordinate_Free(double **Coordinates)
{
	/* Free space for the coordinates */
	//Free the space for each individual coordinate value
	free(Coordinates[0]);		//X coordinate
	free(Coordinates[1]);		//Y coordinate
	free(Coordinates[2]);		//Z coordinate

	//Free space for the entire array
	free(Coordinates);

}


/* This function takes in the X, Y, and Z array and writes it to a file.
 */
void Coordinate_Write(char* Filename, double **Coordinates, int ROWS, int COLS)
{
	FILE	*fpt;
	int i;
	char	Outfile[160];

	/* Write to a file purely for checking results */
	sprintf(Outfile,"%s-coords.coords",Filename);
	fpt=fopen(Outfile,"w");
	fwrite(Coordinates[0],8,ROWS*COLS,fpt);
	fwrite(Coordinates[1],8,ROWS*COLS,fpt);
	fwrite(Coordinates[2],8,ROWS*COLS,fpt);
	fclose(fpt);
	
	/* CSV file for project */
	sprintf(Outfile,"%s-coords.csv",Filename);
	fpt=fopen(Outfile,"w");
	fprintf(fpt, "X,Y,Z\n");
	for(i=0;i<(ROWS*COLS);i++)
		fprintf(fpt,"%lf,%lf,%lf\n", Coordinates[0][i],Coordinates[1][i],Coordinates[2][i]);
	fclose(fpt);
	
}


/*
**	This routine converts the data in an Odetics range image into 3D
**	cartesian coordinate data. The range image is 8-bit, and comes
**	already separated from the intensity image.
*/
//void odetics_to_coords(unsigned char *RangeImage, int ROWS, int COLS, double ***Coordinates)
double **odetics_to_coords(unsigned char *RangeImage, int ROWS, int COLS)
{
	int	r,c;
	double	cp[7];
	double	xangle,yangle,dist;
	double	ScanDirectionFlag,SlantCorrection;
	double **P;
	int ImageTypeFlag;	
	
//	double 	*X,*Y,*Z;										// Coordinates for the image
	
	/* Allocate space for the coordinates */
	//Allocate space for the X, Y, and Z
	P=Coordinate_Allocate(ROWS,COLS);
/*	
	X=(double *)calloc(ROWS*COLS,sizeof(double));
	Y=(double *)calloc(ROWS*COLS,sizeof(double));
	Z=(double *)calloc(ROWS*COLS,sizeof(double));
*/

/*
	printf("Up(-1), Down(1) or Neither(0)? ");
	scanf("%d",&ImageTypeFlag);
*/
	ImageTypeFlag=0;				//Always assume scan is directionally downward <- image type is 0
	//This could have been entirely removed, but I wanted to change the code very little.

	cp[0]=1220.7;									/* horizontal mirror angular velocity in rpm */
	cp[1]=32.0;										/* scan time per single pixel in microseconds */
	cp[2]=(COLS/2)-0.5;						/* middle value of columns */
	cp[3]=1220.7/192.0;						/* vertical mirror angular velocity in rpm */
	cp[4]=6.14;										/* scan time (with retrace) per line in milliseconds */
	cp[5]=(ROWS/2)-0.5;						/* middle value of rows */
	cp[6]=10.0;										/* standoff distance in range units (3.66cm per r.u.) */

	cp[0]=cp[0]*3.1415927/30.0;		/* convert rpm to rad/sec */
	cp[3]=cp[3]*3.1415927/30.0;		/* convert rpm to rad/sec */
	cp[0]=2.0*cp[0];							/* beam ang. vel. is twice mirror ang. vel. */
	cp[3]=2.0*cp[3];							/* beam ang. vel. is twice mirror ang. vel. */
	cp[1]/=1000000.0;							/* units are microseconds : 10^-6 */
	cp[4]/=1000.0;								/* units are milliseconds : 10^-3 */

	switch(ImageTypeFlag)
		{
		case 1:		/* Odetics image -- scan direction upward */
			ScanDirectionFlag=-1;
			break;
		case 0:		/* Odetics image -- scan direction downward */
			ScanDirectionFlag=1;
			break;
		default:		/* in case we want to do this on synthetic model */
			ScanDirectionFlag=0;
			break;
		}

		/* start with semi-spherical coordinates from laser-range-finder: */
		/*			(r,c,RangeImage[r*COLS+c])			*/
		/* convert those to axis-independant spherical coordinates:		*/
		/*			(xangle,yangle,dist)				*/
		/* then convert the spherical coordinates to cartesian:					 */
		/*			(P => X[] Y[] Z[])				*/

	if (ImageTypeFlag != 3)
	{
		for (r=0; r<ROWS; r++)
		{
			for (c=0; c<COLS; c++)
			{
				SlantCorrection=cp[3]*cp[1]*((double)c-cp[2]);
				xangle=cp[0]*cp[1]*((double)c-cp[2]);
				yangle=(cp[3]*cp[4]*(cp[5]-(double)r))+		/* Standard Transform Part */
					SlantCorrection*ScanDirectionFlag;			/*	+ slant correction */
				dist=(double)RangeImage[r*COLS+c]+cp[6];
				P[2][r*COLS+c]=sqrt((dist*dist)/(1.0+(tan(xangle)*tan(xangle))
					+(tan(yangle)*tan(yangle))));
				P[0][r*COLS+c]=tan(xangle)*P[2][r*COLS+c];
				P[1][r*COLS+c]=tan(yangle)*P[2][r*COLS+c];
			}
		}
	}
	
	return P;
}


/* Function to find the surface normals of the coordinates */
double **Surface_Normals(unsigned char* image, double **P, int ROWS, int COLS)
{
	double **SurfaceNormal;
	int r, c;	
	double A[3], B[3], AP[3], BP[3];
	
	//Allocate space for the surface normal coordinates
	SurfaceNormal=Coordinate_Allocate(ROWS,COLS);
	
	/* Loop to calculate the surface normal for each coordinate */
	for (r=0; r<ROWS; r++)
	{
		for (c=0; c<COLS; c++)
		{
			if(!(((r+PIXEL_DISTANCE)>=ROWS)||((c+PIXEL_DISTANCE)>=COLS)||image[r*COLS+c]==255))
			{
				/* Collect the coordinates that are at the pixel distance from the current pixel */
				//Collect coordinates for pixel A (pixel that is shifted along columns for the pixel distance
				A[0]=P[0][r*COLS + (c-PIXEL_DISTANCE)];
				A[1]=P[1][r*COLS + (c-PIXEL_DISTANCE)];
				A[2]=P[2][r*COLS + (c-PIXEL_DISTANCE)];
				
				//Collect coordinates for pixel B (pixel that is shifted along rows for the pixel distance
				B[0]=P[0][(r-PIXEL_DISTANCE)*COLS + c];
				B[1]=P[1][(r-PIXEL_DISTANCE)*COLS + c];
				B[2]=P[2][(r-PIXEL_DISTANCE)*COLS + c];
				
				/* Calculate the coordinate difference for each pixel from the current one */
				//Distance of P from A
				AP[0]=A[0]-P[0][r*COLS + c];
				AP[1]=A[1]-P[1][r*COLS + c];
				AP[2]=A[2]-P[2][r*COLS + c];
				
				//Distance of P from B
				BP[0]=B[0]-P[0][r*COLS + c];
				BP[1]=B[1]-P[1][r*COLS + c];
				BP[2]=B[2]-P[2][r*COLS + c];
				
				/* Calculate the cross product of the A distance and B distance to find the surface normal */
				SurfaceNormal[0][r*COLS+c] = (AP[1]*BP[2])-(AP[2]*BP[1]);			//Surface normal X
				SurfaceNormal[1][r*COLS+c] = (AP[2]*BP[0])-(AP[0]*BP[2]);			//Surface normal Y
				SurfaceNormal[2][r*COLS+c] = (AP[0]*BP[1])-(AP[1]*BP[0]);			//Surface normal Z
			
			}
		}
	}
	return SurfaceNormal;
}


	/*
	** Given an image, a starting point, and a label, this routine
	** paint-fills (8-connected) the area with the given new label
	** according to the given criteria (pixels close to the average
	** intensity of the growing region are allowed to join).
	*/
void RegionGrow(unsigned char *image,		/* image data */
		unsigned char *labels,							/* segmentation labels */
		double **P, 												/* Coordinates of the object */
		int ROWS,int COLS,									/* size of image */
		int r,int c,												/* pixel to paint from */
		int paint_over_label,								/* image label to paint over */
		int new_label,											/* image label for painting */
		int *indices,												/* output: indices of pixels painted */
		int *count,													/* output: count of pixels painted */
		double *avg)
{
	int	r2,c2;
	int	queue[MAX_QUEUE],qh,qt;
	int Index;
	double DotProduct, a_dist, b_dist, angle;

	*count=0;

	if (indices != NULL)
		indices[0]=r*COLS+c;
	queue[0]=r*COLS+c;
	qh=1;	/* queue head */
	qt=0;	/* queue tail */
	(*count)=1;

	/* Initial average normals */
	avg[0] = P[0][(queue[qt]/COLS)*COLS+queue[qt]%COLS];
	avg[1] = P[1][(queue[qt]/COLS)*COLS+queue[qt]%COLS];
	avg[2] = P[2][(queue[qt]/COLS)*COLS+queue[qt]%COLS];

	while (qt != qh)
	{
		for (r2=-1; r2<=1; r2++)
			for (c2=-1; c2<=1; c2++)
			{
				/* This is a variable to hold the index for readability and my sanity */
				//Change the index to have the corresponding qt every loop
				Index = (queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2;
				
				if (labels[Index]!=paint_over_label)
					continue;
				
				if((P[0][Index]!=0) && (P[1][Index]!=0) && (P[2][Index]!=0))
				{
					/* test criteria to join region */
					/* Find the angle */
					avg[0] = (avg[0] + P[0][(r+r2)*COLS+(c+c2)])/2.0;
					avg[1] = (avg[1] + P[1][(r+r2)*COLS+(c+c2)])/2.0;
					avg[2] = (avg[2] + P[2][(r+r2)*COLS+(c+c2)])/2.0;
					
					DotProduct = (avg[0]*P[0][Index]) + (avg[1]*P[1][Index]) + (avg[2]*P[2][Index]);
					a_dist = sqrt(pow(avg[0], 2) + pow(avg[1], 2) + pow(avg[2], 2));
					b_dist = sqrt(pow(P[0][Index], 2) + pow(P[1][Index], 2) + pow(P[2][Index], 2));
					angle = acos(DotProduct/(a_dist*b_dist));

					if(angle > ANGLE_THRESH)
						continue;
					labels[Index]=new_label;
					if (indices != NULL)
						indices[*count]=Index;
					(*count)++;
					queue[qh]=Index;
					qh=(qh+1)%MAX_QUEUE;
					if (qh == qt)
					{
						printf("Max queue size exceeded\n");
						exit(0);
					}					
				}
			}
		qt=(qt+1)%MAX_QUEUE;
	}
}




	/*
	** Reads a greyscale image from the comand line, and creates
	** a segmentation of regions based upon similar greyscale areas.
	** This version demonstrates how to use a paint-fill technique
	** along with region criteria to do segmentation.
	**
	** Good demonstration on targets.ppm.
	*/

unsigned char *RG_Queue(unsigned char *image, double **P, int ROWS, int COLS)

{
	unsigned char	*labels;
	FILE		*fpt;
	int		r,c,r2,c2;
	int		*indices,i;
	int		RegionSize,TotalRegions;
	double avgSN[3];
	int Wr=WINDOW_SIZE/2;
	int IsLabeled;

	/* segmentation image = labels; calloc initializes all labels to 0 */
	labels=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
	
	/* used to quickly erase small grown regions */
	indices=(int *)calloc(ROWS*COLS*3,sizeof(int));
	
	if ((fpt=fopen("chair-seg-regions.csv","wb")) == NULL)
  {
  	printf("Unable to open file for writing\n");
  	exit(0);
  }
  
  /* Write region and average normals to a file */
  //Write header
  fprintf(fpt,"Number of Regions,Number of Pixels,Avg Normal X,Avg Normal Y,Avg Normal Z\n");
  
	TotalRegions=0;
	for (r=Wr; r<ROWS-Wr; r++)
	{
		for (c=Wr; c<COLS-Wr; c++)
		{
			/* Check if any pixel in the 5 x 5 window is already labeled */
			IsLabeled=0;	//Reset the label check;
			for (r2=-Wr; r2<=Wr; r2++)
				for (c2=-Wr; c2<=Wr; c2++)
					if (labels[(r+r2)*COLS+(c+c2)] != 0)
						IsLabeled=1;

      /* condition for seeding a new region is that the pixel is within threshold and not labeled */
			if (IsLabeled!=1 && image[r*COLS+c]!=255)
			{
					TotalRegions++;
					RegionGrow(image,labels,P,ROWS,COLS,r,c,0,TotalRegions,indices,&RegionSize,avgSN);
					if (RegionSize < 40)
					{	
						/* erase region (relabel pixels back to 0) */
						for (i=0; i<RegionSize; i++)
							labels[indices[i]]=0;
						TotalRegions--;
					}
					else
						fprintf(fpt,"%d,%d,%lf,%lf,%lf\n",TotalRegions,RegionSize,avgSN[0],avgSN[1],avgSN[2]);
			}
		}
	}
  fclose(fpt);
	printf("%d total regions were found\n",TotalRegions);
	free(indices);
	  
  return(labels);
}

/* Function to make the labels into a visible image */
unsigned char *ConvertLabelsToShades(unsigned char *labels, int ROWS, int COLS)
{
	int i;
	unsigned char *segmented_image;
	segmented_image=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char));
	for(i=0; i<ROWS*COLS;i++)
		segmented_image[i] = (labels[i]*50)%255;
	return segmented_image;
}



int main()
{
	unsigned char		*InputImage, *ThreshImage, *LabelImage, *SegmentImage;
	double					**P, **SN;

//	unsigned char 	*NormalizedImage;
	char						header[320];
	int							ROWS,COLS,BYTES;						// The information for the input image
//	int							MSFROWS,MSFCOLS,MSFBYTES;		// The information for the MSF image
//	int 						i=0, filesize=0;				// The boundaries of the range
//	int 						*ContRow, *ContCol;
//	char 						c;
	
	/* Read in the images */
	InputImage = Image_Read("chair-range.ppm", header, &ROWS, &COLS, &BYTES);

	/* Threshold the image */
	ThreshImage = Threshold_Image(InputImage, ROWS, COLS);
	Image_Write(ThreshImage, "chair-thresholded.ppm", ROWS, COLS);
	
	/* Convert image to coordinates */
	P = odetics_to_coords(InputImage, ROWS, COLS);
	Coordinate_Write("chair-range", P, ROWS, COLS);
	
	/* Find the surface normals of the coordinates */
	SN = Surface_Normals(ThreshImage, P, ROWS, COLS);
	Coordinate_Write("chair-surfacenormal", SN, ROWS, COLS);
	
	/* Find the segmented image */
	LabelImage = RG_Queue(ThreshImage, SN, ROWS, COLS);
	SegmentImage = ConvertLabelsToShades(LabelImage, ROWS, COLS);	//Give the labels visible shades
	Image_Write(SegmentImage, "chair-segmented.ppm", ROWS, COLS);

	
	/* Free space at the end */
	Image_Free(InputImage);
	Image_Free(ThreshImage);
	Image_Free(LabelImage);
	Image_Free(SegmentImage);
	Coordinate_Free(P);
	Coordinate_Free(SN);
}

	
	












