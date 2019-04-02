#ifndef __OF_H_INCLUDED__
#define __OF_H_INCLUDED__

#define PI 3.14159265359

#define pnum(a) printf("%d\n",a)
#define p(s,a) printf("%s:%d\n",s,a)
#define ps(s,a) printf("%s:%s\n",s,a)
#define pd(s,a) printf("%s:%f\n",s,a)
#define pp(p) printf("x:%d y:%d\n",p.x,p.y)
#define pn printf("\n")
#define Free_ally(ally, w) do {	for(i=0; i < w; i++) { free(ally[i]);}	free(ally);}while(0)

#include "ImageIO/image.h"


typedef struct {
	char descriptor[3];
	int height;
	int width;
	int bright;
	int** data;
}PGM;

typedef struct {
	char descriptor[3];
	int height;
	int width;
	int bright;
	int** dataR;
	int** dataG;
	int** dataB;
}PPM;

typedef struct {
	int x;
	int y;
}Point;


PGM *read_pgm(char *filename);

int write_pgm(char* filename, PGM *pgm) ;

int **create_ally(int width, int height);

double **create_dally(int width, int height);

void copy_ally(int **ally, int **ally2, int w, int h);

void format_ally(int **ally, int w, int h, int bright);

void display_ally(int *ally, int num);

void Free_dally(double **ally, int w);

PGM *copy_pgm(PGM *pgm);

PGM *gaucian_filter(PGM *pgm);

PGM *sobel_filter(PGM *pgm);

PGM *cannyedge_detector(PGM *pgm, double maxValue, double minValue, int thick_min);

void FreePGM(PGM *img);
void FreePPM(PPM *img);

double gause_func(int x, int y, double sigma);

PGM *gaussian_filter(PGM *pgm, double sigma);

PGM *bilateral_filter(PGM *pgm, double sigma, int loopc);

PPM *read_ppm(char *filename);

int write_ppm(char* filename, PPM *pgm);

void devide_ppm(PPM *ppm, PGM* ppmR, PGM* ppmG, PGM* ppmB);

double diffsum_clr(PGM* cmpr, PGM* nimg, Point p, int t, int bright);

PGM *color_gray_conversion(PPM* in);

PPM *copy_ppm(PPM *ppm, int bright);

void normalize_pgm(PGM* in);

PPM *image_to_PPM(image_t *img);
image_t *PPM_to_image(PPM *ppm);
int log_print(char* filename, char *sentence, char *mode);
void Add_dictionary_to_sentence(char* sentence, char *name, int value);

char *get_extension(char *name);

PGM *create_pgm(int width, int height, int bright);
PPM *create_ppm(int width, int height, int bright);
#endif