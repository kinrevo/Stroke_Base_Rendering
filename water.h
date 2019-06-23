#ifndef __WATER_H_INCLUDED__
#define __WATER_H_INCLUDED__

#include "sbr.h"


// UpdateVelocities：Bの係数
#define opt_mhu 0.1
// UpdateVelocities：粘性抵抗
#define opt_kappa 0.01
// RelaxDivergence：減衰の繰り返し回数
#define opt_N 50
// RelaxDivergence：減衰の閾値
#define opt_tau 0.01
// RelaxDivergence：減衰の係数
#define opt_epsilon 0.1
// FlowOutward：ガウスフィルタの直径
#define opt_K 10
// FlowOutward：減衰係数
#define opt_eta 0.5
// TransferPigment：粒状化度、これが大きいほど沈着しづらく、浮上しづらい
#define opt_gamma 0.5
// TransferPigment：顔料の比重、これが大きいほど沈着しやすく浮上しやすい
#define opt_rho 0.20
// TransferPigment：染色力、これが大きいほど浮上しづらい
#define opt_omega 1.5

#define opt_SoakTme 3

double uf(double** u, double x, double y);
double vf(double** v, double x, double y);
void UpdateVelocities(int** M,  double** u, double** v, double** p, double var_t, int width, int height);
void RelaxDivergence(int** M, double** u, double** v, double** p, double var_t, int widht, int height);
void write_Vector_img(PPM* img, double** u, int width, int height);
void FlowOutward(int** M, double** p, double var_t, int width, int height);
void TransferPigment(int** M, double** h, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, double var_t, int width, int height);
void MovePigment(int** M,  double** u, double** v, double** gR, double** gG, double** gB, double var_t, int width, int height);

void calcu_grad_h(double** h, double** grad_hx, double** grad_hy, int width, int height);
void Paint_Water_Stroke(Point StrokeP[], int pnum, int thick, RGB color, int** CanR, int** CanG, int** CanB, 
                            double** h, double** grad_hx, double** grad_hy, double** gauce_filter, int width, int height);
void set_WetStroke(int** M, double** p, double** gR, double** gG, double** gB, Point SP[], int pnum, int thick, RGB color, double** gauce_filter, int width, int height);
void Circle_fill_Water(int** M, double** p, double** gR, double** gG, double** gB, Point SP, int r, RGB color, double** gauce_filter, int width, int height);
void Paint_Water(int** M, double** u, double** v, double** p, double** h, double** grad_hx, double** grad_hy, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, int width, int height);
PPM *c_Illust_brush_Water(PPM *in, char *filename);


void trans_Vector_img(PPM* img, double** u, int width, int height);

double** perlin_img(int width, int height, double freq, int depth);



#endif