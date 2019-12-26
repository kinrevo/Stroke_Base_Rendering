#ifndef __WATER_H_INCLUDED__
#define __WATER_H_INCLUDED__

#include "sbr.h"


// UpdateVelocities：Bの係数(def:0.1)(diffuse:0.)   //ここで水速度の暴走をある程度おさえられる
#define opt_mhu 0.1
// UpdateVelocities：粘性抵抗(def:0.01)
#define opt_kappa 0.01
// RelaxDivergence：減衰の繰り返し回数(def:50)
#define opt_N 50
// RelaxDivergence：減衰の閾値(def:0.01)(diffuse:0.03)
#define opt_tau 0.01
// RelaxDivergence：減衰の係数(def:0.1)
#define opt_xi 0.1
// FlowOutward：ガウスフィルタの直径(def:10)
#define opt_K 5
// FlowOutward：減衰係数(def:0.01-0.05)
#define opt_eta 0.01
// TransferPigment：粒状化度、これが大きいほど沈着しづらく、浮上しづらい(def:0.05)
#define opt_gamma 0.05
// TransferPigment：顔料の比重、これが大きいほど沈着しやすく浮上しやすい(def:0.05)
#define opt_rho 0.05
// TransferPigment：染色力、これが大きいほど浮上しづらい(def:1.0)
#define opt_omega 1.0

// SimulateCapillaryFlow：浅水層から毛細層に水の浸透する割合(def:0.5)(paper:0.5)
#define opt_alpha 0.5
// SimulateCapillaryFlow：毛細層の水がこれより多いと隣に溢れる(def:0.5)(paper:0.62)
#define opt_epsilon 0.5
// SimulateCapillaryFlow：毛細層の水がこれより多いと隣から水が来ない(def:0.65)(paper:0.65)
#define opt_delta 0.65
// SimulateCapillaryFlow：毛細層の水がこれより多いとウェットエリアに加えられる(def:0.001)(paper:0.001)
#define opt_sigma 0.001

// set_WetStroke：水量の初期値をどれだけ中心に寄せるか(def:1.5)(diffuse:10.0)
#define opt_variance_ratio 10.0

// 水と顔料を定着させる時間(diffuse:50,1.0)
#define opt_SoakTime 5
#define opt_SoakTimeStep 1.0
// perlinノイズのパラメータ(def:0.1,6)
#define opt_perlin_freq 0.1
#define opt_perlin_depth 6
// Option機能
#define opt_USE_MoveWater 0
#define opt_USE_Backrun 0
    #define exp_DiffuseNum 1
#define opt_USE_DETAIL_TP 0     //paperでの変数設定法 (diffuse:0.05,0.05,0.8)
    #define opt_deposit 0.05
    #define opt_lift 0.05
    #define opt_exposure 0.8


//////　実験動作  /////
#define exp_RaiseWater 0  // FlowOutward:負の水量が現れたときそれに合わせて全体を底上げする
#define opt_StopSoakVero 2  // UpdateVelocities：水速度が暴走したときに停止する速度条件(diffuse:2.0)
#define opt_RemovePigmentInWater 1
#define opt_FloatPigmentOnPaper 1
// #define _DEBUG_PaintWater



double uf(double** u, double x, double y);
double vf(double** v, double x, double y);
double UpdateVelocities(int** M,  double** u, double** v, double** p, double var_t, int width, int height, Point rectangleP[]);
void RelaxDivergence(int** M, double** u, double** v, double** p, double var_t, int width, int height, Point rectangleP[]);
void FlowOutward(int** M, double** p, int c, double** gause_filter, double var_t, int width, int height, Point rectangleP[]);
void MovePigment(int** M,  double** u, double** v, double** gR, double** gG, double** gB, double var_t, int width, int height, Point rectangleP[]);
void TransferPigment(int** M, double** h, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, double var_t, int width, int height);
void SimulateCapillaryFlow(int** M, double** p, double** c, double** s, double var_t, int width, int height);
void MoveWater(int** M,  double** u, double** v, double** p, double var_t, int width, int height);
double** perlin_img(int width, int height, double freq, int depth);
void calcu_grad_h(double** h, double** grad_hx, double** grad_hy, int width, int height);


void Paint_Water_Stroke(Point StrokeP[], int pnum, int thick, RGB color, int** CanR, int** CanG, int** CanB, double** h, double** grad_hx, double** grad_hy, double** gauce_filter, int width, int height);
void set_WetStroke(int** M, double** p, double** gR, double** gG, double** gB, Point SP[], int pnum, int thick, RGB color, double** gauce_filter, int width, int height);
void Circle_fill_Water(int** M, double** p, double** gR, double** gG, double** gB, Point SP, int r, RGB color, double** gauce_filter, int width, int height);
void Paint_Water(int** M, double** u, double** v, double** p, double** h, double** grad_hx, double** grad_hy, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, int width, int height, Point rectangleP[]);
int test_water_stroke(PPM* test_Canvas, PPM* cmpr, PPM* nimgC, Stroke* stroke, int t, double** h, double** grad_hx, double** grad_hy, double** gauce_filter);

void write_Vector_img(PPM* img, double** u, int width, int height);
void trans_Vector_img(PPM* img, double** u, int width, int height);

void SINGLE_Paint_Water_Stroke(Point StrokeP[], int pnum, int thick, RGB color, int** CanR, int** CanG, int** CanB, double** h, double** grad_hx, double** grad_hy, double** gauce_filter, int width, int height);

__global__ void gpu_calculate_best_stroke(int *dev_GLOBAL_improved_value_map, float *dev_PerlinNoise, int *dev_cmprR,	int *dev_cmprG,	int *dev_cmprB,	int *dev_nimgR, int *dev_nimgG,	int *dev_nimgB,
                                int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y, int *dev_best_stroke_map_R, int *dev_best_stroke_map_G,
                                int *dev_best_stroke_map_B, float *dev_sobel_abs, float *dev_sobel_angle, float *dev_grad_hx, float *dev_grad_hy, float *dev_in_Lab_L, float *dev_in_Lab_a,
                                float *dev_in_Lab_b, char *dev_M, float *dev_u, float *dev_new_u, float *dev_v, float *dev_new_v, float *dev_p, float *dev_gR, float *dev_gG, float *dev_gB,
                                float *dev_dR, float *dev_dG, float *dev_dB, float *dev_new_gR,	float *dev_new_gG, float *dev_new_gB, float *dev_new_dR, float *dev_new_dG, float *dev_new_dB,
                                float *dev_gauss_filter, float *dev_gauss_M, float *dev_h, int width, int height, int t);

__global__ void gpu_select_best_stroke(int *dev_GLOBAL_improved_value_map, int *dev_best_x, int *dev_best_y, int *dev_diff_stroke_max, int width, int height);

__global__ void gpu_draw_best_stroke(float *dev_PerlinNoise, int *dev_nimgR, int *dev_nimgG, int *dev_nimgB, int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y,
                                int *dev_best_stroke_map_R, int *dev_best_stroke_map_G, int *dev_best_stroke_map_B,float *dev_grad_hx, float *dev_grad_hy, char *dev_M, float *dev_u, float *dev_new_u,
                                float *dev_v, float *dev_new_v, float *dev_p, float *dev_gR, float *dev_gG, float *dev_gB, float *dev_dR, float *dev_dG, float *dev_dB, float *dev_new_gR, float *dev_new_gG,
                                float *dev_new_gB, float *dev_new_dR, float *dev_new_dG, float *dev_new_dB, float *dev_gauss_filter, float *dev_gauss_M, float *dev_h, int *dev_best_x, int *dev_best_y, 
                                int width, int height, int t);

__global__ void gpu_reset_improved_value_map(int *dev_GLOBAL_improved_value_map, int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y, int *dev_best_x, int *dev_best_y, int width, int height, int t);


#endif