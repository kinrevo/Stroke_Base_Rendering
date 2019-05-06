#ifndef __SBR_H_INCLUDED__
#define __SBR_H_INCLUDED__

#define PI 3.14159265359
#define Search_DIRECTRY 1
#define Search_FILE 2
#define Reversal_ON 1
#define Reversal_OFF 0
#define NotYetSeach

#define pnum(a) printf("%d\n",a)
#define p(s,a) printf("%s:%d\n",s,a)
#define ps(s,a) printf("%s:%s\n",s,a)
#define pd(s,a) printf("%s:%f\n",s,a)
#define pp(p) printf("x:%f y:%f\n",p.x,p.y)
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
	double x;
	double y;
}Point;

// typedef struct {
	// int x;
	// int y;
// }Point;

typedef struct {
	double a;
	double b;
	double c;
	double d;
}S_para;

typedef struct {
	int R;
	int G;
	int B;
}RGB;

typedef struct {
	int pnum;
	Point* p;
	RGB color;
}Stroke;



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

PGM *CompareImg(PGM *in1, PGM *in2);

PGM *RotateImg(PGM *in, int n);

PGM *SIFT(PGM *pgm, int S);

PGM *bilateral_filter(PGM *pgm, double sigma, int loopc);

PPM *read_ppm(char *filename);

int write_ppm(char* filename, PPM *pgm);

PGM *Illustification(PGM *in, double sigma,int n);

int inner_point(int p1,int p2,float m,float n);

int BezierCurve(int p1, int p2, int p3, int p4, float t);

Point BezierCurve_P(Point p1, Point p2, Point p3, Point p4, float t);

void Paint_Bezier(Point p1, Point s, Point p4, PGM* in_img, int thick, int bright);

void Circle_fill(int xc, int yc, int r, PGM* in_img, int bright, double ratio);

void light_dot(int x, int y, PGM* in_img, int bright);

double sobel_slope(PGM* in_img, int x, int y, int t);

PGM *Illust_brush(PGM *in, char *filename);

void Paint_line(Point p1, Point p4, PGM* in_img, int thick, int bright, double ratio);

int Spline_Curve(double a, double b, double c, double d, int xj, double x);

void Spline_paramater(Point p1, Point p2, Point p3, S_para S[]);

void Paint_Spline(Point p1, Point p2, Point p3, PGM* in_img, int thick, int bright);

Point Rotate_dot(Point p, double theta, Point center);

void sobel_calcu(PGM *pgm, double **sobel_abs, double **sobel_angle);

void Paint_Bezier_ex(Point p[], int pnum, PGM* in_img, int thick, int bright, double ratio);

double calcu_histogram(PGM* cmpr, double **sobel_abs, double **sobel_angle, int partition, 
		double **gauce_filter, int x, int y, int t, int *break_flag);
		
Point calcu_point(PGM *in, Point a, int t, double theta);

void devide_ppm(PPM *ppm, PGM* ppmR, PGM* ppmG, PGM* ppmB);

PPM *c_Illust_brush(PPM *in, char *filename);

double diffsum_clr(PGM* cmpr, PGM* nimg, Point p, int t, int bright);

PGM *color_gray_conversion(PPM* in);

PPM *copy_ppm(PPM *ppm, int bright);

int calcu_color(int** data, int width, int height, int x, int y, int t);
int calcu_color_bi(int** data, int width, int height, int x, int y, int t, double sigma, double** gause_filter);

PGM *calcu_EdgeMap(PGM *canny, int thick_min, double **sobel_angle);
int edge_detect(PGM *canny, PGM *map, double **sobel_angle, int x, int y, int t, Point s, int *flag, int count_init);
PGM *expand_Edge(PGM *canny, int thick_min);

double** color_sobel_abs(PPM *in, int num);
PGM *canny_detect_in(PGM *pgm, double maxValue, double minValue, double** sobel);
void normalize_pgm(PGM* in);

PPM *image_to_PPM(image_t *img);
image_t *PPM_to_image(PPM *ppm);
int log_print(char* filename, char *sentence, char *mode);
void Add_dictionary_to_sentence(char* sentence, char *name, int value);

char** make_file_list(int* file_num, char* search_path, int Search_TYPE);
void multi_file_CIB(char** list, int list_num);
void multi_dir_CIB(char* in_path);

char *get_extension(char *name);

PGM *create_pgm(int width, int height, int bright);
PPM *create_ppm(int width, int height, int bright);

Point *scaling_point(Point p[], int pnum, double canvas_scaling_ratio);

int vec_print(char* filename, Point *p, int pnum, int brightR, int brightG, int brightB, int width, int height);

int test_stroke(PGM* test_Canvas, PGM* cmprR, PGM* cmprG, PGM* cmprB, PGM* nimgR, PGM* nimgG, PGM* nimgB, Point p[], int pnum, int t, int brightR, int brightG, int brightB, double ratio);

int calcu_Stroke_Point(PGM* cmprR, PGM* cmprG, PGM* cmprB, PGM* nimgR, PGM* nimgG, PGM* nimgB,
		double **sobel_abs, double **sobel_angle, double **gauce_filter, 
		int t, int max_stroke, Point Start_P, int brightR, int brightG, int brightB,
		Point* Stroke_P, int direct_reversal_flag);
		
Point search_max_Point(int **ally, int w, int h);

Stroke *create_Stroke(int pnum);
Stroke ***create_Stroke_ally(int width, int height, int max_stroke);

#endif