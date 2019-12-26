#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
//#include <png.h>
//#include <jpeglib.h>
#include <omp.h>
#include "ImageIO/image.h"
#include "sbr.h"
#include "sbr_opt.h"
#include "water.h"

__device__ static float atomicMaxFloat(float* address, float val){
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed, __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

/*
if((threadIdx.x==0)&&(threadIdx.y==0)){
    printf("AAAAAAAAAAAAAAAAAA\n");
    for(int j=0;j<stroke_length_max;j++){
        for(int i=0;i<stroke_length_max;i++){
            printf("%3d ", s_test_Canvas_B[i+j*stroke_length_max]);
        }
        printf("\n");
    }
}
__syncthreads();
*/

//ストロークを沢山描いて各ストロークの改善値を計算する関数
__global__ void gpu_calculate_best_stroke(int *dev_GLOBAL_improved_value_map, float *dev_PerlinNoise, int *dev_cmprR,	int *dev_cmprG,	int *dev_cmprB,	int *dev_nimgR, int *dev_nimgG,	int *dev_nimgB,
                                int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y, int *dev_best_stroke_map_R, int *dev_best_stroke_map_G,
                                int *dev_best_stroke_map_B, float *dev_sobel_abs, float *dev_sobel_angle, float *dev_grad_hx, float *dev_grad_hy, float *dev_in_Lab_L, float *dev_in_Lab_a,
                                float *dev_in_Lab_b, char *dev_M, float *dev_u, float *dev_new_u, float *dev_v, float *dev_new_v, float *dev_p, float *dev_gR, float *dev_gG, float *dev_gB,
                                float *dev_dR, float *dev_dG, float *dev_dB, float *dev_new_gR,	float *dev_new_gG, float *dev_new_gB, float *dev_new_dR, float *dev_new_dG, float *dev_new_dB,
                                float *dev_gauss_filter, float *dev_gauss_M, float *dev_h, int width, int height, int t){

    int global_blockID = blockIdx.x + blockIdx.y * gridDim.x;
    int stroke_length_max = opt_thick_max*(opt_max_stroke+2);
    int SubImageIndex = global_blockID * stroke_length_max * stroke_length_max;
    int SubImageIndex_uv = global_blockID * stroke_length_max * (stroke_length_max+1);

    int x,y,i,j,k,l,yc,xc,sy,sx,ly,lx;
    int pnum, peak, stroke_partition, left_end, right_end, upper_end, lower_end, stroke_length_x, stroke_length_y;
    float temp_R, temp_G, temp_B, Can_L, Can_a, Can_b, testCan_L, testCan_a, testCan_b, Lab_x, Lab_y, Lab_z, error_before, error_after;
    float format_theta, scale, temp_x, temp_y, tmpSP_start_x, tmpSP_start_y, tmpSP_end_x, tmpSP_end_y, tmpSP0_x, tmpSP0_y, tmpSP1_x, tmpSP1_y;
    float UV_var_t, A, B, p_grad, delta, sum, filter_sum, down, up, down_ratio, up_ratio, pigment_density;
    
    int tmp_density_R = 255 - 255 * opt_ratio;//描画色の濃度
    int tmp_density_G = 255 - 255 * opt_ratio;//描画色の濃度
    int tmp_density_B = 255 - 255 * opt_ratio;//描画色の濃度

    int w = (int)(ceil(3.0*opt_K/6.0+0.5)*2-1); //ガウシアンフィルタ用
    int c = (w-1)/2;                            //ガウシアンフィルタ用

	__shared__ int s_improved_value[1];
	__shared__ float s_diff_sum[1];
	__shared__ int s_offscrn_count[1];
    __shared__ float s_theta[1];
    __shared__ float s_histogram[opt_histogram_partition];
    __shared__ float s_error_sum[1];
    __shared__ short s_test_Canvas_R[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ short s_test_Canvas_G[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ short s_test_Canvas_B[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ float s_max_velocity[1];
    __shared__ float s_delta_MAX[1];

    for(y=blockIdx.y; y<height; y=y+gridDim.y) {
		for(x=blockIdx.x; x<width; x=x+gridDim.x) {

			// 改善値が計算済みならSkip
			if(dev_GLOBAL_improved_value_map[x+y*width] != UNCALCULATED) continue;

			pnum = peak = 0;
            if((threadIdx.x==0)&&(threadIdx.y==0)){
                s_diff_sum[0] = 0;
                s_offscrn_count[0] = 0;
            }
            __syncthreads();

			//ウィンドウの中の差分の合計を並列に計算
			for(yc=threadIdx.y-t; yc<=t; yc=yc+blockDim.y) {
				for(xc=threadIdx.x-t; xc<=t; xc=xc+blockDim.x) {
					if((x+xc)<0 || (x+xc)>(width-1) || (y+yc)<0 || (y+yc)>(height-1)){
                        atomicAdd(&s_offscrn_count[0], 1);
					}else{

                        //描画色の平均を求める用
						atomicAdd(&dev_best_stroke_map_R[x+y*width], dev_cmprR[(x+xc)+(y+yc)*width]);
						atomicAdd(&dev_best_stroke_map_G[x+y*width], dev_cmprG[(x+xc)+(y+yc)*width]);
						atomicAdd(&dev_best_stroke_map_B[x+y*width], dev_cmprB[(x+xc)+(y+yc)*width]);

                        //誤差の平均を求める用
                        atomicAdd(&s_diff_sum[0], fabsf((float)(dev_nimgR[(x+xc)+(y+yc)*width] - dev_cmprR[(x+xc)+(y+yc)*width])));
						atomicAdd(&s_diff_sum[0], fabsf((float)(dev_nimgG[(x+xc)+(y+yc)*width] - dev_cmprG[(x+xc)+(y+yc)*width])));
						atomicAdd(&s_diff_sum[0], fabsf((float)(dev_nimgB[(x+xc)+(y+yc)*width] - dev_cmprB[(x+xc)+(y+yc)*width])));
					}
				}
			}
            __syncthreads();

			//1つの代表スレッドを起動
			if((threadIdx.x==0)&&(threadIdx.y==0)){

                //描画色の平均を計算
				dev_best_stroke_map_R[x+y*width] = dev_best_stroke_map_R[x+y*width] / ((2*t+1)*(2*t+1)-s_offscrn_count[0]);
				dev_best_stroke_map_G[x+y*width] = dev_best_stroke_map_G[x+y*width] / ((2*t+1)*(2*t+1)-s_offscrn_count[0]);
				dev_best_stroke_map_B[x+y*width] = dev_best_stroke_map_B[x+y*width] / ((2*t+1)*(2*t+1)-s_offscrn_count[0]);

                //誤差の平均を計算
                s_diff_sum[0] = s_diff_sum[0] / ((2*t+1)*(2*t+1)-s_offscrn_count[0]) / 3;
            }
            __syncthreads();

            //差分の合計平均(画素当たりの差分)が一定以上ならストローク開始位置とする
            if(s_diff_sum[0] < opt_window_diff_border) {
                dev_GLOBAL_improved_value_map[x+y*width] = SMALL_DIFF;
                continue;
            }

            //1つの代表スレッドが第一制御点の座標を格納
			if((threadIdx.x==0)&&(threadIdx.y==0)){
                dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = x+0.5;
				dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = y+0.5;
			}
			__syncthreads();

			pnum = 1;		//第1制御点確定

            //SharedMemory初期化
            for(i=threadIdx.x+threadIdx.y*blockDim.x; i<opt_histogram_partition ; i+= blockDim.x*blockDim.y) {s_histogram[i]=0;}
            __syncthreads();

            //sobelからヒストグラムを作成
            for(sy=threadIdx.y+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]-t-1; sy<=dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]+t+1; sy+=blockDim.y) {
                for(sx=threadIdx.x+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]-t-1; sx<=dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]+t+1; sx+=blockDim.x) {
                    if(!(sx<1 || sx>width-2 || sy<1 || sy>height-2)){
                        atomicAdd(&s_histogram[(int)((dev_sobel_angle[sx+sy*width]/PI)*opt_histogram_partition)], dev_sobel_abs[sx+sy*width]);
                    }
                }
            }
            __syncthreads();

            //1つの代表スレッドを起動
            if((threadIdx.x==0)&&(threadIdx.y==0)){

                //ヒストグラムの中で最も大きい値を探索
                for(i=0; i<opt_histogram_partition; i++){
                    if(s_histogram[peak] < s_histogram[i]) peak = i;
                }

                //thetaを計算
                s_theta[0] = ((float)peak/opt_histogram_partition)*PI+((PI/opt_histogram_partition)/2);

                //次の制御点を方向から計算
                dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = t*cos(s_theta[0])+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
                dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = t*sin(s_theta[0])+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];
                if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] < 0)
                    dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = 0;
                if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] < 0)
                    dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = 0;
                if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] >= width)
                    dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = width-1;
                if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] >= height)
                    dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = height-1;
            }
            __syncthreads();

            //二つ目の描画点周りの色が描画色と一致するか確認する
            error_before = error_after = 0;
            if((threadIdx.x==0)&&(threadIdx.y==0)){
                s_offscrn_count[0] = 0;
                s_error_sum[0] = 0;
            }
            __syncthreads();

            for(ly=threadIdx.y+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]-t; ly<=dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]+t; ly+=blockDim.y) {
                for(lx=threadIdx.x+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]-t; lx<=dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]+t; lx+=blockDim.x) {
                    if(lx<0 || lx>width-1 || ly<0 || ly>height-1) {
                        atomicAdd(&s_offscrn_count[0], 1);
                    }else if( (lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])*(lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])+(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum])*(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]) > t*t ){
                        atomicAdd(&s_offscrn_count[0], 1);
                    }else{

                        //塗る前のLabを計算
                        temp_R = dev_nimgR[lx+ly*width];
                        temp_G = dev_nimgG[lx+ly*width];
                        temp_B = dev_nimgB[lx+ly*width];
                        temp_R = temp_R / 255.0;
                        temp_G = temp_R / 255.0;
                        temp_B = temp_R / 255.0;
                        temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                        temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                        temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                        Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                        Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                        Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                        Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                        Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                        Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                        Can_L = (116 * Lab_y) - 16;
                        Can_a = 500 * (Lab_x - Lab_y);
                        Can_b = 200 * (Lab_y - Lab_z);

                        //塗った後のLabを計算
                        temp_R = dev_nimgR[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_R[x+y*width] * opt_ratio;
                        temp_G = dev_nimgG[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_G[x+y*width] * opt_ratio;
                        temp_B = dev_nimgB[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_B[x+y*width] * opt_ratio;
                        temp_R = temp_R / 255.0;
                        temp_G = temp_R / 255.0;
                        temp_B = temp_R / 255.0;
                        temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                        temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                        temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                        Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                        Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                        Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                        Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                        Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                        Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                        testCan_L = (116 * Lab_y) - 16;
                        testCan_a = 500 * (Lab_x - Lab_y);
                        testCan_b = 200 * (Lab_y - Lab_z);

                        //描画前後のLab空間におけるユークリッド距離誤算の変化を計算
                        error_before = sqrt( powf(dev_in_Lab_L[lx+ly*width]-Can_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-Can_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-Can_b, 2) );
                        error_after = sqrt( powf(dev_in_Lab_L[lx+ly*width]-testCan_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-testCan_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-testCan_b, 2) );
                        atomicAdd(&s_error_sum[0], error_before-error_after);
                    }
                }
            }
            __syncthreads();

            //二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
			if(s_error_sum[0] < opt_color_diff_border){

                //1つの代表スレッドを起動してもう1つの勾配垂直の点を代入
                if((threadIdx.x==0)&&(threadIdx.y==0)){
                    s_theta[0] += PI;
                    dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = t*cos(s_theta[0])+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
                    dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = t*sin(s_theta[0])+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];
                    if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] < 0)
                        dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = 0;
                    if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] < 0)
                        dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = 0;
                    if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] >= width)
                        dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = width-1;
                    if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] >= height)
                        dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = height-1;
                }
                __syncthreads();

                //反対方向の第2点の描画点周りの色が描画色と一致するか確認する
                error_before = error_after = 0;
                if((threadIdx.x==0)&&(threadIdx.y==0)){
                    s_offscrn_count[0] = 0;
                    s_error_sum[0] = 0;
                }
                __syncthreads();

                for(ly=threadIdx.y+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]-t; ly<=dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]+t; ly+=blockDim.y) {
                    for(lx=threadIdx.x+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]-t; lx<=dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]+t; lx+=blockDim.x) {
                        if(lx<0 || lx>width-1 || ly<0 || ly>height-1) {
                            atomicAdd(&s_offscrn_count[0], 1);
                        }else if( (lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])*(lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])+(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum])*(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]) > t*t ){
                            atomicAdd(&s_offscrn_count[0], 1);
                        }else{

                            //塗る前のLabを計算
                            temp_R = dev_nimgR[lx+ly*width];
                            temp_G = dev_nimgG[lx+ly*width];
                            temp_B = dev_nimgB[lx+ly*width];
                            temp_R = temp_R / 255.0;
                            temp_G = temp_R / 255.0;
                            temp_B = temp_R / 255.0;
                            temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                            temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                            temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                            Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                            Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                            Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                            Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                            Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                            Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                            Can_L = (116 * Lab_y) - 16;
                            Can_a = 500 * (Lab_x - Lab_y);
                            Can_b = 200 * (Lab_y - Lab_z);

                            //塗った後のLabを計算
                            temp_R = dev_nimgR[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_R[x+y*width] * opt_ratio;
                            temp_G = dev_nimgG[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_G[x+y*width] * opt_ratio;
                            temp_B = dev_nimgB[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_B[x+y*width] * opt_ratio;
                            temp_R = temp_R / 255.0;
                            temp_G = temp_R / 255.0;
                            temp_B = temp_R / 255.0;
                            temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                            temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                            temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                            Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                            Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                            Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                            Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                            Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                            Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                            testCan_L = (116 * Lab_y) - 16;
                            testCan_a = 500 * (Lab_x - Lab_y);
                            testCan_b = 200 * (Lab_y - Lab_z);

                            //描画前後のLab空間におけるユークリッド距離誤算の変化を計算
                            error_before = sqrt( powf(dev_in_Lab_L[lx+ly*width]-Can_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-Can_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-Can_b, 2) );
                            error_after = sqrt( powf(dev_in_Lab_L[lx+ly*width]-testCan_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-testCan_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-testCan_b, 2) );
                            atomicAdd(&s_error_sum[0], error_before-error_after);
                        }
                    }
                }
                __syncthreads();

                //どちらの第2点も不適切なら描画をせず次のループへ
                if(s_error_sum[0] < opt_color_diff_border){
                    dev_GLOBAL_improved_value_map[x+y*width] = MIN_STROKE;
                    continue;
                }
            }

            pnum=2;		//第2制御点確定

            //ストロークが最大opt_max_strokeになるまでは伸ばし続ける
            while(pnum < opt_max_stroke){

                //thetaの値を保存
                format_theta = s_theta[0];

                //第pnum点周りにおいて、sobelからヒストグラムを作成
                peak = 0;//初期化
                for(i=threadIdx.x+threadIdx.y*blockDim.x; i<opt_histogram_partition ; i+= blockDim.x*blockDim.y) {s_histogram[i]=0;}//ヒストグラムを初期化
                __syncthreads();

                for(sy=threadIdx.y+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]-t-1; sy<=dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]+t+1; sy=sy+blockDim.y) {
                    for(sx=threadIdx.x+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]-t-1; sx<=dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]+t+1; sx=sx+blockDim.x) {
                        if(!(sx<1 || sx>width-2 || sy<1 || sy>height-2))
                            atomicAdd(&s_histogram[(int)((dev_sobel_angle[sx+sy*width]/PI)*opt_histogram_partition)], dev_sobel_abs[sx+sy*width]);
                    }
                }
                __syncthreads();

                //1つの代表スレッドを起動
                if((threadIdx.x==0)&&(threadIdx.y==0)){

                    //ヒストグラムの中で最も大きい値を探索
                    for(i=0; i<opt_histogram_partition; i++){
                        if(s_histogram[peak] < s_histogram[i]) peak = i;
                    }

                    //thetaを計算
                    s_theta[0] = ((float)peak/opt_histogram_partition)*PI+((PI/opt_histogram_partition)/2);

                    //制御点の為す角が急峻になるようなら逆方向に角度を取る
					if( (s_theta[0] < format_theta-PI/2) || (s_theta[0] > format_theta+PI/2) ) {s_theta[0] += PI;}

                    //次の制御点を方向から計算
                    dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = t*cos(s_theta[0])+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
                    dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = t*sin(s_theta[0])+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];
                    if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] < 0)
                        dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = 0;
                    if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] < 0)
                        dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = 0;
                    if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] >= width)
                        dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum] = width-1;
                    if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] >= height)
                        dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum] = height-1;
                }
                __syncthreads();

                //pnum+1目の描画点周りの色が描画色と一致するか確認する
                error_before = error_after = 0;
                if((threadIdx.x==0)&&(threadIdx.y==0)){
                    s_offscrn_count[0] = 0;
                    s_error_sum[0] = 0;
                }
                __syncthreads();

                for(ly=threadIdx.y+dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]-t; ly<=dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]+t; ly=ly+blockDim.y) {
                    for(lx=threadIdx.x+dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]-t; lx<=dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum]+t; lx=lx+blockDim.x) {
                        if(lx<0 || lx>width-1 || ly<0 || ly>height-1) {
                            atomicAdd(&s_offscrn_count[0], 1);
                        }else if( (lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])*(lx-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum])+(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum])*(ly-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum]) > t*t ){
                            atomicAdd(&s_offscrn_count[0], 1);
                        }else{

                            //塗る前のLabを計算
                            temp_R = dev_nimgR[lx+ly*width];
                            temp_G = dev_nimgG[lx+ly*width];
                            temp_B = dev_nimgB[lx+ly*width];
                            temp_R = temp_R / 255.0;
                            temp_G = temp_R / 255.0;
                            temp_B = temp_R / 255.0;
                            temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                            temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                            temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                            Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                            Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                            Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                            Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                            Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                            Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                            Can_L = (116 * Lab_y) - 16;
                            Can_a = 500 * (Lab_x - Lab_y);
                            Can_b = 200 * (Lab_y - Lab_z);

                            //塗った後のLabを計算
                            temp_R = dev_nimgR[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_R[x+y*width] * opt_ratio;
                            temp_G = dev_nimgG[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_G[x+y*width] * opt_ratio;
                            temp_B = dev_nimgB[lx+ly*width] * (1-opt_ratio) + dev_best_stroke_map_B[x+y*width] * opt_ratio;
                            temp_R = temp_R / 255.0;
                            temp_G = temp_R / 255.0;
                            temp_B = temp_R / 255.0;
                            temp_R = temp_R > 0.04045 ? powf(((temp_R + 0.055) / 1.055), 2.4) : (temp_R / 12.92);
                            temp_G = temp_G > 0.04045 ? powf(((temp_G + 0.055) / 1.055), 2.4) : (temp_G / 12.92);
                            temp_B = temp_B > 0.04045 ? powf(((temp_B + 0.055) / 1.055), 2.4) : (temp_B / 12.92);
                            Lab_x = ((temp_R * 0.4124) + (temp_G * 0.3576) + (temp_B * 0.1805)) * 100 / 95.047;
                            Lab_y = (temp_R * 0.2126) + (temp_G * 0.7152) + (temp_B * 0.0722);
                            Lab_z = ((temp_R * 0.0193) + (temp_G * 0.1192) + (temp_B * 0.9505)) * 100 / 108.883;
                            Lab_x = Lab_x > 0.008856 ? powf(Lab_x, 1 / 3.0) : (7.787 * Lab_x) + (4 / 29.0);
                            Lab_y = Lab_y > 0.008856 ? powf(Lab_y, 1 / 3.0) : (7.787 * Lab_y) + (4 / 29.0);
                            Lab_z = Lab_z > 0.008856 ? powf(Lab_z, 1 / 3.0) : (7.787 * Lab_z) + (4 / 29.0);
                            testCan_L = (116 * Lab_y) - 16;
                            testCan_a = 500 * (Lab_x - Lab_y);
                            testCan_b = 200 * (Lab_y - Lab_z);

                            //描画前後のLab空間におけるユークリッド距離誤算の変化を計算
                            error_before = sqrt( powf(dev_in_Lab_L[lx+ly*width]-Can_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-Can_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-Can_b, 2) );
                            error_after = sqrt( powf(dev_in_Lab_L[lx+ly*width]-testCan_L, 2) + powf(dev_in_Lab_a[lx+ly*width]-testCan_a, 2) + powf(dev_in_Lab_b[lx+ly*width]-testCan_b, 2) );
                            atomicAdd(&s_error_sum[0], error_before-error_after);
                        }
                    }
                }
                __syncthreads();

                if(s_error_sum[0] < opt_color_diff_border) break;
                else pnum++;
            }

            //1つの代表スレッドが制御点の個数を記録
            if((threadIdx.x==0)&&(threadIdx.y==0)){
                dev_best_stroke_map_pnum[x+y*width] = pnum;
            }
            __syncthreads();
            
            //////////ストローク形状の計算終了//////////

            //if((threadIdx.x==0)&&(threadIdx.y==0))    printf("(%d,%d)Stroke_form_end\n",x,y);
            //__syncthreads();

            //制御点の個数が足りない場合
            if(pnum < opt_min_stroke){
                dev_GLOBAL_improved_value_map[x+y*width] = MIN_STROKE;
                continue;
            }

            //////////Paint_Water_Stroke(試しに描いてみて誤差を確認)//////////

            //ストローク点を囲む端の座標を特定
            left_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0];       //切り捨て
            right_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]+1;    //切り上げ
            upper_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0];      //切り捨て
            lower_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]+1;    //切り上げ
            for(i=1; i<pnum; i++){
                if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i] < left_end) left_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i];
                if(right_end < dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]) right_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]+1;
                if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i] < upper_end) upper_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i];
                if(lower_end < dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]) lower_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]+1;
            }

            //ストローク半径分、端座標を膨張
            left_end-=t; right_end+=t; upper_end-=t; lower_end+=t;
            if(left_end < 0) left_end = 0;
            if(width < right_end) right_end = width-1;
            if(upper_end < 0) upper_end=0;
            if(height < lower_end) lower_end = height-1;

            stroke_length_x = right_end - left_end; //ストロークの横の長さ
            stroke_length_y = lower_end - upper_end; //ストロークの縦の長さ

            //各パラメータの初期化
            for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                    dev_M[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_p[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_gR[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_gG[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_gB[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_gR[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_gG[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_gB[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_dR[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_dG[SubImageIndex+j+i*stroke_length_max] = 0;
                    dev_new_dB[SubImageIndex+j+i*stroke_length_max] = 0;
                }
            }

            //uとvを初期化
            for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                    dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = 0;
                    dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = 0;
                }
            }
            for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                    dev_v[SubImageIndex_uv+j+i*stroke_length_max] = 0;
                    dev_new_v[SubImageIndex_uv+j+i*stroke_length_max] = 0;
                }
            }

            //影響箇所をコピー
            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    dev_dR[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgR[(left_end+i)+(upper_end+j)*width] / 255.0;
                    dev_dG[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgG[(left_end+i)+(upper_end+j)*width] / 255.0;
                    dev_dB[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgB[(left_end+i)+(upper_end+j)*width] / 255.0;
                    s_test_Canvas_R[i+j*stroke_length_max] = dev_nimgR[(left_end+i)+(upper_end+j)*width];
                    s_test_Canvas_G[i+j*stroke_length_max] = dev_nimgG[(left_end+i)+(upper_end+j)*width];
                    s_test_Canvas_B[i+j*stroke_length_max] = dev_nimgB[(left_end+i)+(upper_end+j)*width];
                    dev_h[SubImageIndex+i+j*stroke_length_max] = dev_PerlinNoise[(left_end+i)+(upper_end+j)*width];
                }
            }
            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] = (-1)*dev_grad_hx[(left_end+i)+(upper_end+j)*(width+1)];
                }
            }
            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    dev_v[SubImageIndex_uv+i+j*stroke_length_max] = (-1)*dev_grad_hy[(left_end+i)+(upper_end+j)*width];
                }
            }
            __syncthreads();

            //////////set_WetStroke(ストローク点に従いウェットエリアと水量を計算)//////////
            if(pnum == 2){//もし制御点が2つの時は直線を引く
                stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1])
                                    + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]); //線分の分割数                    
                
                for(k=0; k<=stroke_partition; k++){
                    scale = (float)k / stroke_partition;
                    temp_x = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] + (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0])*scale;
                    temp_y = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] + (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0])*scale;

                    for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                        for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                            if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                                if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                                    dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                                    dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                                    dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                                    dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                                    dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                                }
                            }
                        }
                    }
                }
            
            //制御点が3つ以上の場合
            }else{

                //1.最初の制御点について
                //端の一つ外の制御点を適当に決める
                tmpSP_start_x = 2*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
                tmpSP_start_y = 2*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];
                
                //ベジエ曲線用の中間点を2点決める
                tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1]-tmpSP_start_x)/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0];
                tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]-tmpSP_start_y)/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0];
                tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+2])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
                tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+2])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];

                //線分の分割数
                stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1])
                                + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]);
                
                for(k=0; k<=stroke_partition; k++){
                    scale = (float)k / stroke_partition;
                    temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
                    temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];

                    for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                        for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                            if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                                if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                                    dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                                    dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                                    dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                                    dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                                    dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                                }
                            }
                        }
                    }
                }
                __syncthreads();

                //2.中間の制御点について
                for(i=1; i<pnum-2; i++){

                    //ベジエ曲線用の中間点を2点決める
                    tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i-1])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i];
                    tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i-1])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i];
                    tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+2])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1];
                    tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+2])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1];

                    //線分の分割数
                    stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1])
                                    + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1]);

                    for(j=0; j<=stroke_partition; j++){
                        scale = (float)j / stroke_partition;
                        temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1];
                        temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1];

                        for(k=threadIdx.y+temp_y-t-upper_end; k<=temp_y+t-upper_end; k+=blockDim.y){
                            for(l=threadIdx.x+temp_x-t-left_end; l<=temp_x+t-left_end; l+=blockDim.x){
                                if((l>=0) && (l<=stroke_length_x) && (k>=0) && (k<=stroke_length_y)){//画像内かどうか
                                    if((l-temp_x+left_end)*(l-temp_x+left_end)+(k-temp_y+upper_end)*(k-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                                        dev_M[SubImageIndex+l+k*stroke_length_max] = 1;
                                        dev_gR[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                                        dev_gG[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                                        dev_gB[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                                        dev_p[SubImageIndex+l+k*stroke_length_max] = 1.0;
                                    }
                                }
                            }
                        }
                    }
                }

                //3.最後の制御点について
                //端の一つ外の制御点を適当に決める
                tmpSP_end_x = 2*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2];
                tmpSP_end_y = 2*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2];

                //ベジエ曲線用の中間点を2点決める
                tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-3])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2];
                tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-3])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2];
                tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2]-tmpSP_end_x)/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
                tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2]-tmpSP_end_y)/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];

                //線分の分割数
                stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1])
                                + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]);

                for(k=0; k<=stroke_partition; k++){
                    scale = (float)k / stroke_partition;
                    temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
                    temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];

                    for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                        for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                            if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                                if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                                    dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                                    dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                                    dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                                    dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                                    dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                                }
                            }
                        }
                    }
                }
            }
            __syncthreads();

            //if((threadIdx.x==0)&&(threadIdx.y==0))    printf("(%d,%d)Set_WetStroke_end\n",x,y);
            //__syncthreads();

            //時間の経過を表すループ
            for (float time=0; time<opt_SoakTime; time+=opt_SoakTimeStep){

                //////////UpdateVelocities(一定時間経過後の速度変化を計算)//////////
                if((threadIdx.x==0)&&(threadIdx.y==0)) s_max_velocity[0] = 0;//最大初速度を初期化
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max]==1){
                            atomicMaxFloat(&s_max_velocity[0], fabsf(dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]));
                            atomicMaxFloat(&s_max_velocity[0], fabsf(dev_v[SubImageIndex_uv+i+j*stroke_length_max]));
                        }
                    }
                }
                __syncthreads();

                if(s_max_velocity[0] > opt_StopSoakVero) break; //水速度が暴走したときは停止
                UV_var_t = fminf(opt_SoakTimeStep, opt_SoakTimeStep/s_max_velocity[0]); //最大初速度が大きいほど細かく更新を行う

                for (float a = 0; a < opt_SoakTimeStep; a+=UV_var_t){

                    for(j=1+threadIdx.y; j<=stroke_length_y-1; j+=blockDim.y){
                        for(i=1+threadIdx.x; i<=stroke_length_x-1; i+=blockDim.x){
                            if(dev_M[SubImageIndex+i+j*stroke_length_max]==1 && dev_M[SubImageIndex+(i+1)+j*stroke_length_max]==1){
                                
                                // paper_DETAIL
                                A = powf((dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)])/2, 2)
                                    - powf((dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+2)+j*(stroke_length_max+1)])/2, 2)
                                    + (dev_u[SubImageIndex_uv+(i+1)+(j-1)*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+j*stroke_length_max]+dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max])/4
                                    - (dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max])/4;
                                B = dev_u[SubImageIndex_uv+(i+2)+j*(stroke_length_max+1)]
                                    + dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]
                                    + dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)]
                                    + dev_u[SubImageIndex_uv+(i+1)+(j-1)*(stroke_length_max+1)]
                                    - 4*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)];
                                p_grad = dev_p[SubImageIndex+i+j*stroke_length_max] - dev_p[SubImageIndex+(i+1)+j*stroke_length_max];

                                dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] + UV_var_t*(A - opt_mhu*B + p_grad - opt_kappa * exp(-0.1*(dev_p[SubImageIndex+i+j*stroke_length_max]+dev_p[SubImageIndex+(i+1)+j*stroke_length_max])/2)*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]);
                            }
                            else{
                                dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = 0; //ウェットエリア外を速度０に
                            }
                        }
                    }
                    __syncthreads();

                    for(j=1+threadIdx.y; j<=stroke_length_y-1; j+=blockDim.y){
                        for(i=1+threadIdx.x; i<=stroke_length_x-1; i+=blockDim.x){
                            if(dev_M[SubImageIndex+i+j*stroke_length_max]==1 && dev_M[SubImageIndex+i+(j+1)*stroke_length_max]==1){
                                
                                // paper_DETAIL
                                A = powf((dev_v[SubImageIndex_uv+i+j*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max])/2, 2)
                                    - powf((dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+2)*stroke_length_max])/2, 2)
                                    + (dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+i+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max])/4
                                    - (dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max])/4;
                                B = dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max]
                                    + dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max]
                                    + dev_v[SubImageIndex_uv+i+(j+2)*stroke_length_max]
                                    + dev_v[SubImageIndex_uv+i+j*stroke_length_max]
                                    - 4*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max];
                                p_grad = dev_p[SubImageIndex+i+j*stroke_length_max] - dev_p[SubImageIndex+i+(j+1)*stroke_length_max];

                                dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] + UV_var_t*(A - opt_mhu*B + p_grad - opt_kappa*exp(-0.1*(dev_p[SubImageIndex_uv+i+j*stroke_length_max]+dev_p[SubImageIndex_uv+i+(j+1)*stroke_length_max])/2)*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]);
                            }
                            else{
                                dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = 0;
                            }
                        }
                    }
                    __syncthreads();

                    //uにnew_uをコピー
                    for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                            dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                        }
                    }
                    //vにnew_vをコピー
                    for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                            dev_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_new_v[SubImageIndex_uv+j+i*stroke_length_max];
                        }
                    }
                    __syncthreads();
                }

                //////////RelaxDivergence(速度ベクトルの発散をある許容範囲τ未満になるまで緩和)//////////
                for(int n = 0; n < opt_N; n++){

                    //new_uにuをコピー
                    for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                            dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                        }
                    }
                    //new_vにvをコピー
                    for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                            dev_new_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_v[SubImageIndex_uv+j+i*stroke_length_max];
                        }
                    }
                    __syncthreads();

                    //初期化
                    if((threadIdx.x==0)&&(threadIdx.y==0)) s_delta_MAX[0] = 0;
                    __syncthreads();

                    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                            if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                                delta = opt_xi * (dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] - dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] + dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] - dev_v[SubImageIndex_uv+i+j*stroke_length_max]);
                                dev_p[SubImageIndex+i+j*stroke_length_max] =  dev_p[SubImageIndex+i+j*stroke_length_max] - delta;
                                if(i!=stroke_length_x && dev_M[SubImageIndex+(i+1)+j*stroke_length_max]==1) dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] - delta;
                                if(j!=stroke_length_y && dev_M[SubImageIndex+i+(j+1)*stroke_length_max]==1) dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] - delta;
                                atomicMaxFloat(&s_delta_MAX[0], fabsf(delta));
                            }
                        }
                    }
                    __syncthreads();

                    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                            if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                                if(i!=0 && dev_M[SubImageIndex+(i-1)+j*stroke_length_max]==1) dev_new_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] + delta;
                                if(j!=0 && dev_M[SubImageIndex+i+(j-1)*stroke_length_max]==1) dev_new_v[SubImageIndex_uv+i+j*stroke_length_max] = dev_new_v[SubImageIndex_uv+i+j*stroke_length_max] + delta;
                            }
                        }
                    }
                    __syncthreads();

                    if(s_delta_MAX[0] < opt_tau) break;

                    //uにnew_uをコピー
                    for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                            dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                        }
                    }
                    //vにnew_vをコピー
                    for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                        for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                            dev_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_new_v[SubImageIndex_uv+j+i*stroke_length_max];
                        }
                    }
                }

                //////////FlowOutward/////////
                for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
                    for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                        dev_gauss_M[SubImageIndex+i+j*stroke_length_max] = 0; //dev_gauss_Mを0で初期化
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        //注目しているピクセルがウェットエリアならガウスフィルタによる拡散された値を周囲に足し込む
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            sum = 0;
                            filter_sum = 1;
                            for(l=-c; l<=c; l++){
                                for(k=-c; k<=c; k++){
                                    //フィルタの端ピクセルがない場合、分子には加算せず分母から減算
                                    if( (i+k)<0 || (i+k)>stroke_length_x || (j+l)<0 || (j+l)>stroke_length_y){
                                        filter_sum -= dev_gauss_filter[(k+c)+(l+c)*w];
                                    }else{
                                        sum += dev_M[SubImageIndex+(i+k)+(j+l)*stroke_length_max] * dev_gauss_filter[(k+c)+(l+c)*w];
                                    }
                                }
                            }
                            dev_gauss_M[SubImageIndex+i+j*stroke_length_max] = sum / filter_sum;
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_p[SubImageIndex+i+j*stroke_length_max] = dev_p[SubImageIndex+i+j*stroke_length_max] - opt_eta * opt_SoakTimeStep * (1-dev_gauss_M[SubImageIndex+i+j*stroke_length_max])*dev_M[SubImageIndex+i+j*stroke_length_max];
                        }
                    }
                }
                __syncthreads();

                //////////MovePigment//////////
                for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
                    for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                        dev_new_gR[SubImageIndex+i+j*stroke_length_max] = dev_gR[SubImageIndex+i+j*stroke_length_max]; //コピー
                        dev_new_gG[SubImageIndex+i+j*stroke_length_max] = dev_gG[SubImageIndex+i+j*stroke_length_max]; //コピー
                        dev_new_gB[SubImageIndex+i+j*stroke_length_max] = dev_gB[SubImageIndex+i+j*stroke_length_max]; //コピー
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_new_gR[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gR[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gG[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gG[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gB[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gB[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_new_gR[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gR[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gG[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gG[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gB[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gB[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_new_gR[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gR[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gG[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gG[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gB[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gB[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_new_gR[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gR[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gG[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gG[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gB[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gB[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                            dev_new_gR[SubImageIndex+i+j*stroke_length_max] = dev_new_gR[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gG[SubImageIndex+i+j*stroke_length_max] = dev_new_gG[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_new_gB[SubImageIndex+i+j*stroke_length_max] = dev_new_gB[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                        }
                    }
                }
                __syncthreads();

                for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
                    for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                        dev_gR[SubImageIndex+i+j*stroke_length_max] = dev_new_gR[SubImageIndex+i+j*stroke_length_max]; //コピー
                        dev_gG[SubImageIndex+i+j*stroke_length_max] = dev_new_gG[SubImageIndex+i+j*stroke_length_max]; //コピー
                        dev_gB[SubImageIndex+i+j*stroke_length_max] = dev_new_gB[SubImageIndex+i+j*stroke_length_max]; //コピー
                    }
                }
                __syncthreads();

                //////////TransferPigment//////////
                for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                    for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                        if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){

                            if(opt_USE_DETAIL_TP){
                                down_ratio = opt_SoakTimeStep * (1-dev_h[SubImageIndex+i+j*stroke_length_max]*opt_exposure) * opt_deposit;
                                up_ratio = opt_SoakTimeStep * (1+(dev_h[SubImageIndex+i+j*stroke_length_max]-1)*opt_exposure) * opt_lift;
                            }else{
                                down_ratio = opt_SoakTimeStep * (1-dev_h[SubImageIndex+i+j*stroke_length_max]*opt_gamma) * opt_rho;
                                up_ratio = opt_SoakTimeStep * (1+(dev_h[SubImageIndex+i+j*stroke_length_max]-1)*opt_gamma) * opt_rho / opt_omega;
                            }

                            //R
                            down = dev_gR[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                            up   = dev_dR[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                            if(dev_dR[SubImageIndex+i+j*stroke_length_max]+down > 1)
                                down = fmaxf(0,1-dev_dR[SubImageIndex+i+j*stroke_length_max]);
                            if(dev_gR[SubImageIndex+i+j*stroke_length_max]+up > 1)
                                up = fmaxf(0,1-dev_gR[SubImageIndex+i+j*stroke_length_max]);
                            dev_dR[SubImageIndex+i+j*stroke_length_max] = dev_dR[SubImageIndex+i+j*stroke_length_max]+down-up;
                            dev_gR[SubImageIndex+i+j*stroke_length_max] = dev_gR[SubImageIndex+i+j*stroke_length_max]+up-down;

                            //G
                            down = dev_gG[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                            up   = dev_dG[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                            if(dev_dG[SubImageIndex+i+j*stroke_length_max]+down > 1)
                                down = fmaxf(0,1-dev_dG[SubImageIndex+i+j*stroke_length_max]);
                            if(dev_gG[SubImageIndex+i+j*stroke_length_max]+up > 1)
                                up = fmaxf(0,1-dev_gG[SubImageIndex+i+j*stroke_length_max]);
                            dev_dG[SubImageIndex+i+j*stroke_length_max] = dev_dG[SubImageIndex+i+j*stroke_length_max]+down-up;
                            dev_gG[SubImageIndex+i+j*stroke_length_max] = dev_gG[SubImageIndex+i+j*stroke_length_max]+up-down;

                            //B
                            down = dev_gB[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                            up   = dev_dB[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                            if(dev_dB[SubImageIndex+i+j*stroke_length_max]+down > 1)
                                down = fmaxf(0,1-dev_dB[SubImageIndex+i+j*stroke_length_max]);
                            if(dev_gB[SubImageIndex+i+j*stroke_length_max]+up > 1)
                                up = fmaxf(0,1-dev_gB[SubImageIndex+i+j*stroke_length_max]);
                            dev_dB[SubImageIndex+i+j*stroke_length_max] = dev_dB[SubImageIndex+i+j*stroke_length_max]+down-up;
                            dev_gB[SubImageIndex+i+j*stroke_length_max] = dev_gB[SubImageIndex+i+j*stroke_length_max]+up-down;
                        }
                    }
                }
            }
            __syncthreads();

            //if((threadIdx.x==0)&&(threadIdx.y==0))    printf("(%d,%d)Time_Roop_end\n",x,y);
            //__syncthreads();

            //堆積顔料をRGBに変換しキャンバスに描画
            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){          

                    // シミュレーション終了時に水中の顔料を落とすか拭き取るか
                    if(opt_FloatPigmentOnPaper){
                        s_test_Canvas_R[i+j*stroke_length_max] = (1 - dev_dR[SubImageIndex+i+j*stroke_length_max]) * 255;    //CMY[0,1]->RGB[0,255]
                        s_test_Canvas_G[i+j*stroke_length_max] = (1 - dev_dG[SubImageIndex+i+j*stroke_length_max]) * 255;
                        s_test_Canvas_B[i+j*stroke_length_max] = (1 - dev_dB[SubImageIndex+i+j*stroke_length_max]) * 255;
                    }
                    else{
                        //R
                        pigment_density = dev_gR[SubImageIndex+i+j*stroke_length_max] + dev_dR[SubImageIndex+i+j*stroke_length_max];
                        LIMIT_RANGE(pigment_density, 0, 1);
                        s_test_Canvas_R[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_R[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_R[x+y*width];

                        //G
                        pigment_density = dev_gG[SubImageIndex+i+j*stroke_length_max] + dev_dG[SubImageIndex+i+j*stroke_length_max];
                        LIMIT_RANGE(pigment_density, 0, 1);
                        s_test_Canvas_G[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_G[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_G[x+y*width];

                        //B
                        pigment_density = dev_gB[SubImageIndex+i+j*stroke_length_max] + dev_dB[SubImageIndex+i+j*stroke_length_max];
                        LIMIT_RANGE(pigment_density, 0, 1);
                        s_test_Canvas_B[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_B[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_B[x+y*width];
                    }
                }
            }
            __syncthreads();

            //if((threadIdx.x==0)&&(threadIdx.y==0))    printf("(%d,%d)Paint_Stroke_end\n",x,y);
            //__syncthreads();

            //実際にストロークを描いた時の誤差を計算
            if((threadIdx.x==0)&&(threadIdx.y==0)){
                s_improved_value[0] = 0; //初期化
            }
            __syncthreads();
            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                        atomicAdd(&s_improved_value[0], abs(dev_cmprR[(left_end+i)+(upper_end+j)*width]-dev_nimgR[(left_end+i)+(upper_end+j)*width]) - abs(dev_cmprR[(left_end+i)+(upper_end+j)*width]-s_test_Canvas_R[i+j*stroke_length_max]));
                        atomicAdd(&s_improved_value[0], abs(dev_cmprG[(left_end+i)+(upper_end+j)*width]-dev_nimgG[(left_end+i)+(upper_end+j)*width]) - abs(dev_cmprG[(left_end+i)+(upper_end+j)*width]-s_test_Canvas_G[i+j*stroke_length_max]));
                        atomicAdd(&s_improved_value[0], abs(dev_cmprB[(left_end+i)+(upper_end+j)*width]-dev_nimgB[(left_end+i)+(upper_end+j)*width]) - abs(dev_cmprB[(left_end+i)+(upper_end+j)*width]-s_test_Canvas_B[i+j*stroke_length_max]));
                    }
                }
            }
            __syncthreads();

            //if((threadIdx.x==0)&&(threadIdx.y==0))    printf("(%d,%d)Error_cul_end\n",x,y);
            //__syncthreads();

            if((threadIdx.x==0)&&(threadIdx.y==0)){
                dev_GLOBAL_improved_value_map[x+y*width] = s_improved_value[0];
                //printf("(%d,%d) = %d\n",x,y,s_improved_value[0]);
            }
            __syncthreads();
        }
    }
}

//改善値マップ中の最大値を探索する関数
__global__ void gpu_select_best_stroke(int *dev_GLOBAL_improved_value_map, int *dev_best_x, int *dev_best_y, int *dev_diff_stroke_max, int width, int height){

    dev_best_x[0] = 0;
    dev_best_y[0] = 0;
    dev_diff_stroke_max[0] = -99999999;

    for(int y=0; y<height; y++){
        for(int x=0; x<width; x++){
            if(dev_GLOBAL_improved_value_map[x+y*width] > dev_diff_stroke_max[0]){
                dev_best_x[0] = x;
                dev_best_y[0] = y;
                dev_diff_stroke_max[0] = dev_GLOBAL_improved_value_map[x+y*width];
            }
        }
    }
}

//実際にベストストロークを描画する関数
__global__ void gpu_draw_best_stroke(float *dev_PerlinNoise, int *dev_nimgR, int *dev_nimgG, int *dev_nimgB, int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y,
                                int *dev_best_stroke_map_R, int *dev_best_stroke_map_G, int *dev_best_stroke_map_B,float *dev_grad_hx, float *dev_grad_hy, char *dev_M, float *dev_u, float *dev_new_u,
                                float *dev_v, float *dev_new_v, float *dev_p, float *dev_gR, float *dev_gG, float *dev_gB, float *dev_dR, float *dev_dG, float *dev_dB, float *dev_new_gR, float *dev_new_gG,
                                float *dev_new_gB, float *dev_new_dR, float *dev_new_dG, float *dev_new_dB, float *dev_gauss_filter, float *dev_gauss_M, float *dev_h, int *dev_best_x, int *dev_best_y, 
                                int width, int height, int t){

    int global_blockID = blockIdx.x + blockIdx.y * gridDim.x; //1ブロックしか起動しないのでglobal_blockID=0
    int stroke_length_max = opt_thick_max*(opt_max_stroke+2);
    int SubImageIndex = global_blockID * stroke_length_max * stroke_length_max;
    int SubImageIndex_uv = global_blockID * stroke_length_max * (stroke_length_max+1);

    int x,y,i,j,k,l;
    int pnum, stroke_partition, left_end, right_end, upper_end, lower_end, stroke_length_x, stroke_length_y;
    float scale, temp_x, temp_y, tmpSP_start_x, tmpSP_start_y, tmpSP_end_x, tmpSP_end_y, tmpSP0_x, tmpSP0_y, tmpSP1_x, tmpSP1_y;
    float UV_var_t, A, B, p_grad, delta, sum, filter_sum, down, up, down_ratio, up_ratio, pigment_density;
    
    int tmp_density_R = 255 - 255 * opt_ratio;//描画色の濃度
    int tmp_density_G = 255 - 255 * opt_ratio;//描画色の濃度
    int tmp_density_B = 255 - 255 * opt_ratio;//描画色の濃度

    int w = (int)(ceil(3.0*opt_K/6.0+0.5)*2-1); //とりあえず動く計算
    int c = (w-1)/2;

    __shared__ short s_test_Canvas_R[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ short s_test_Canvas_G[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ short s_test_Canvas_B[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ float s_max_velocity[1];
    __shared__ float s_delta_MAX[1];

    x = dev_best_x[0];
    y = dev_best_y[0];
    pnum = dev_best_stroke_map_pnum[x+y*width];

    //////////Paint_Water_Stroke(試しに描いてみて誤差を確認)//////////

    //ストローク点を囲む端の座標を特定
    left_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0];       //切り捨て
    right_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]+1;    //切り上げ
    upper_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0];      //切り捨て
    lower_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]+1;    //切り上げ
    for(i=1; i<pnum; i++){
        if(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i] < left_end) left_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i];
        if(right_end < dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]) right_end = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]+1;
        if(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i] < upper_end) upper_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i];
        if(lower_end < dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]) lower_end = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]+1;
    }

    //ストローク半径分、端座標を膨張
    left_end-=t; right_end+=t; upper_end-=t; lower_end+=t;
    if(left_end < 0) left_end = 0;
    if(width < right_end) right_end = width-1;
    if(upper_end < 0) upper_end=0;
    if(height < lower_end) lower_end = height-1;

    stroke_length_x = right_end - left_end; //ストロークの横の長さ
    stroke_length_y = lower_end - upper_end; //ストロークの縦の長さ

    //各パラメータの初期化
    for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
        for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
            dev_M[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_p[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_gR[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_gG[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_gB[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_gR[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_gG[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_gB[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_dR[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_dG[SubImageIndex+j+i*stroke_length_max] = 0;
            dev_new_dB[SubImageIndex+j+i*stroke_length_max] = 0;
        }
    }

    //uとvを初期化
    for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
        for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
            dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = 0;
            dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = 0;
        }
    }
    for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
        for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
            dev_v[SubImageIndex_uv+j+i*stroke_length_max] = 0;
            dev_new_v[SubImageIndex_uv+j+i*stroke_length_max] = 0;
        }
    }

    //影響箇所をコピー
    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
            dev_dR[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgR[(left_end+i)+(upper_end+j)*width] / 255.0;
            dev_dG[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgG[(left_end+i)+(upper_end+j)*width] / 255.0;
            dev_dB[SubImageIndex+i+j*stroke_length_max] = 1 - dev_nimgB[(left_end+i)+(upper_end+j)*width] / 255.0;
            s_test_Canvas_R[i+j*stroke_length_max] = dev_nimgR[(left_end+i)+(upper_end+j)*width];
            s_test_Canvas_G[i+j*stroke_length_max] = dev_nimgG[(left_end+i)+(upper_end+j)*width];
            s_test_Canvas_B[i+j*stroke_length_max] = dev_nimgB[(left_end+i)+(upper_end+j)*width];
            dev_h[SubImageIndex+i+j*stroke_length_max] = dev_PerlinNoise[(left_end+i)+(upper_end+j)*width];
        }
    }
    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
            dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] = (-1)*dev_grad_hx[(left_end+i)+(upper_end+j)*(width+1)];
        }
    }
    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
            dev_v[SubImageIndex_uv+i+j*stroke_length_max] = (-1)*dev_grad_hy[(left_end+i)+(upper_end+j)*width];
        }
    }
    __syncthreads();

    //////////set_WetStroke(ストローク点に従いウェットエリアと水量を計算)//////////
    if(pnum == 2){//もし制御点が2つの時は直線を引く
        stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1])
                            + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]); //線分の分割数                    
        
        for(k=0; k<=stroke_partition; k++){
            scale = (float)k / stroke_partition;
            temp_x = dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] + (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0])*scale;
            temp_y = dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] + (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0])*scale;

            for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                    if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                        if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                            dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                            dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                            dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                            dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                            dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                        }
                    }
                }
            }
        }
    
    //制御点が3つ以上の場合
    }else{

        //1.最初の制御点について
        //端の一つ外の制御点を適当に決める
        tmpSP_start_x = 2*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
        tmpSP_start_y = 2*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];
        
        //ベジエ曲線用の中間点を2点決める
        tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1]-tmpSP_start_x)/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0];
        tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]-tmpSP_start_y)/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0];
        tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+2])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
        tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+2])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];

        //線分の分割数
        stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1])
                        + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1]);
        
        for(k=0; k<=stroke_partition; k++){
            scale = (float)k / stroke_partition;
            temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+0] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+1];
            temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+0] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+1];

            for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                    if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                        if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                            dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                            dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                            dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                            dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                            dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                        }
                    }
                }
            }
        }
        __syncthreads();

        //2.中間の制御点について
        for(i=1; i<pnum-2; i++){

            //ベジエ曲線用の中間点を2点決める
            tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i-1])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i];
            tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i-1])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i];
            tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+2])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1];
            tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+2])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1];

            //線分の分割数
            stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1])
                            + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1]);

            for(j=0; j<=stroke_partition; j++){
                scale = (float)j / stroke_partition;
                temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+i+1];
                temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+i+1];

                for(k=threadIdx.y+temp_y-t-upper_end; k<=temp_y+t-upper_end; k+=blockDim.y){
                    for(l=threadIdx.x+temp_x-t-left_end; l<=temp_x+t-left_end; l+=blockDim.x){
                        if((l>=0) && (l<=stroke_length_x) && (k>=0) && (k<=stroke_length_y)){//画像内かどうか
                            if((l-temp_x+left_end)*(l-temp_x+left_end)+(k-temp_y+upper_end)*(k-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                                dev_M[SubImageIndex+l+k*stroke_length_max] = 1;
                                dev_gR[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                                dev_gG[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                                dev_gB[SubImageIndex+l+k*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                                dev_p[SubImageIndex+l+k*stroke_length_max] = 1.0;
                            }
                        }
                    }
                }
            }
        }

        //3.最後の制御点について
        //端の一つ外の制御点を適当に決める
        tmpSP_end_x = 2*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2];
        tmpSP_end_y = 2*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2];

        //ベジエ曲線用の中間点を2点決める
        tmpSP0_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-3])/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2];
        tmpSP0_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-3])/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2];
        tmpSP1_x = (dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2]-tmpSP_end_x)/6.0 + dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
        tmpSP1_y = (dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2]-tmpSP_end_y)/6.0 + dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];

        //線分の分割数
        stroke_partition = fabsf(dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2]-dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1])
                        + fabsf(dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2]-dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1]);

        for(k=0; k<=stroke_partition; k++){
            scale = (float)k / stroke_partition;
            temp_x = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-2] + 3*scale*(1-scale)*(1-scale)*tmpSP0_x + 3*scale*scale*(1-scale)*tmpSP1_x + scale*scale*scale*dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+pnum-1];
            temp_y = (1-scale)*(1-scale)*(1-scale)*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-2] + 3*scale*(1-scale)*(1-scale)*tmpSP0_y + 3*scale*scale*(1-scale)*tmpSP1_y + scale*scale*scale*dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+pnum-1];

            for(i=threadIdx.y+temp_y-t-upper_end; i<=temp_y+t-upper_end; i+=blockDim.y){
                for(j=threadIdx.x+temp_x-t-left_end; j<=temp_x+t-left_end; j+=blockDim.x){
                    if((j>=0) && (j<=stroke_length_x) && (i>=0) && (i<=stroke_length_y)){//画像内かどうか
                        if((j-temp_x+left_end)*(j-temp_x+left_end)+(i-temp_y+upper_end)*(i-temp_y+upper_end) <= t*t){//円の範囲内かどうか
                            dev_M[SubImageIndex+j+i*stroke_length_max] = 1;
                            dev_gR[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_R[x+y*width]/255.0;
                            dev_gG[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_G[x+y*width]/255.0;
                            dev_gB[SubImageIndex+j+i*stroke_length_max] = 1 - dev_best_stroke_map_B[x+y*width]/255.0;
                            dev_p[SubImageIndex+j+i*stroke_length_max] = 1.0;
                        }
                    }
                }
            }
        }
    }
    __syncthreads();

    //時間の経過を表すループ
    for (float time=0; time<opt_SoakTime; time+=opt_SoakTimeStep){

        //////////UpdateVelocities(一定時間経過後の速度変化を計算)//////////
        if((threadIdx.x==0)&&(threadIdx.y==0)) s_max_velocity[0] = 0;//最大初速度を初期化
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max]==1){
                    atomicMaxFloat(&s_max_velocity[0], fabsf(dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]));
                    atomicMaxFloat(&s_max_velocity[0], fabsf(dev_v[SubImageIndex_uv+i+j*stroke_length_max]));
                }
            }
        }
        __syncthreads();

        if(s_max_velocity[0] > opt_StopSoakVero) break; //水速度が暴走したときは停止
        UV_var_t = fminf(opt_SoakTimeStep, opt_SoakTimeStep/s_max_velocity[0]); //最大初速度が大きいほど細かく更新を行う

        for (float a = 0; a < opt_SoakTimeStep; a+=UV_var_t){

            for(j=1+threadIdx.y; j<=stroke_length_y-1; j+=blockDim.y){
                for(i=1+threadIdx.x; i<=stroke_length_x-1; i+=blockDim.x){
                    if(dev_M[SubImageIndex+i+j*stroke_length_max]==1 && dev_M[SubImageIndex+(i+1)+j*stroke_length_max]==1){
                        
                        // paper_DETAIL
                        A = powf((dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)])/2, 2)
                            - powf((dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+2)+j*(stroke_length_max+1)])/2, 2)
                            + (dev_u[SubImageIndex_uv+(i+1)+(j-1)*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+j*stroke_length_max]+dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max])/4
                            - (dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max])/4;
                        B = dev_u[SubImageIndex_uv+(i+2)+j*(stroke_length_max+1)]
                            + dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]
                            + dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)]
                            + dev_u[SubImageIndex_uv+(i+1)+(j-1)*(stroke_length_max+1)]
                            - 4*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)];
                        p_grad = dev_p[SubImageIndex+i+j*stroke_length_max] - dev_p[SubImageIndex+(i+1)+j*stroke_length_max];

                        dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] + UV_var_t*(A - opt_mhu*B + p_grad - opt_kappa * exp(-0.1*(dev_p[SubImageIndex+i+j*stroke_length_max]+dev_p[SubImageIndex+(i+1)+j*stroke_length_max])/2)*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]);
                    }
                    else{
                        dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = 0; //ウェットエリア外を速度０に
                    }
                }
            }
            __syncthreads();

            for(j=1+threadIdx.y; j<=stroke_length_y-1; j+=blockDim.y){
                for(i=1+threadIdx.x; i<=stroke_length_x-1; i+=blockDim.x){
                    if(dev_M[SubImageIndex+i+j*stroke_length_max]==1 && dev_M[SubImageIndex+i+(j+1)*stroke_length_max]==1){
                        
                        // paper_DETAIL
                        A = powf((dev_v[SubImageIndex_uv+i+j*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max])/2, 2)
                            - powf((dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+2)*stroke_length_max])/2, 2)
                            + (dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+i+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max])/4
                            - (dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]+dev_u[SubImageIndex_uv+(i+1)+(j+1)*(stroke_length_max+1)])*(dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]+dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max])/4;
                        B = dev_v[SubImageIndex_uv+(i+1)+(j+1)*stroke_length_max]
                            + dev_v[SubImageIndex_uv+(i-1)+(j+1)*stroke_length_max]
                            + dev_v[SubImageIndex_uv+i+(j+2)*stroke_length_max]
                            + dev_v[SubImageIndex_uv+i+j*stroke_length_max]
                            - 4*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max];
                        p_grad = dev_p[SubImageIndex+i+j*stroke_length_max] - dev_p[SubImageIndex+i+(j+1)*stroke_length_max];

                        dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] + UV_var_t*(A - opt_mhu*B + p_grad - opt_kappa*exp(-0.1*(dev_p[SubImageIndex_uv+i+j*stroke_length_max]+dev_p[SubImageIndex_uv+i+(j+1)*stroke_length_max])/2)*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]);
                    }
                    else{
                        dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = 0;
                    }
                }
            }
            __syncthreads();

            //uにnew_uをコピー
            for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                    dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                }
            }
            //vにnew_vをコピー
            for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                    dev_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_new_v[SubImageIndex_uv+j+i*stroke_length_max];
                }
            }
            __syncthreads();
        }

        //////////RelaxDivergence(速度ベクトルの発散をある許容範囲τ未満になるまで緩和)//////////
        for(int n = 0; n < opt_N; n++){

            //new_uにuをコピー
            for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                    dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                }
            }
            //new_vにvをコピー
            for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                    dev_new_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_v[SubImageIndex_uv+j+i*stroke_length_max];
                }
            }
            __syncthreads();

            //初期化
            if((threadIdx.x==0)&&(threadIdx.y==0)) s_delta_MAX[0] = 0;
            __syncthreads();

            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                        delta = opt_xi * (dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] - dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] + dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] - dev_v[SubImageIndex_uv+i+j*stroke_length_max]);
                        dev_p[SubImageIndex+i+j*stroke_length_max] =  dev_p[SubImageIndex+i+j*stroke_length_max] - delta;
                        if(i!=stroke_length_x && dev_M[SubImageIndex+(i+1)+j*stroke_length_max]==1) dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)] - delta;
                        if(j!=stroke_length_y && dev_M[SubImageIndex+i+(j+1)*stroke_length_max]==1) dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] = dev_new_v[SubImageIndex_uv+i+(j+1)*stroke_length_max] - delta;
                        atomicMaxFloat(&s_delta_MAX[0], fabsf(delta));
                    }
                }
            }
            __syncthreads();

            for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
                for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                    if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                        if(i!=0 && dev_M[SubImageIndex+(i-1)+j*stroke_length_max]==1) dev_new_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+i+j*(stroke_length_max+1)] + delta;
                        if(j!=0 && dev_M[SubImageIndex+i+(j-1)*stroke_length_max]==1) dev_new_v[SubImageIndex_uv+i+j*stroke_length_max] = dev_new_v[SubImageIndex_uv+i+j*stroke_length_max] + delta;
                    }
                }
            }
            __syncthreads();

            if(s_delta_MAX[0] < opt_tau) break;

            //uにnew_uをコピー
            for(i=threadIdx.y; i<stroke_length_max; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max+1; j+=blockDim.x){
                    dev_u[SubImageIndex_uv+j+i*(stroke_length_max+1)] = dev_new_u[SubImageIndex_uv+j+i*(stroke_length_max+1)];
                }
            }
            //vにnew_vをコピー
            for(i=threadIdx.y; i<stroke_length_max+1; i+=blockDim.y){
                for(j=threadIdx.x; j<stroke_length_max; j+=blockDim.x){
                    dev_v[SubImageIndex_uv+j+i*stroke_length_max] = dev_new_v[SubImageIndex_uv+j+i*stroke_length_max];
                }
            }
        }

        //////////FlowOutward/////////
        for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
            for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                dev_gauss_M[SubImageIndex+i+j*stroke_length_max] = 0; //dev_gauss_Mを0で初期化
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                //注目しているピクセルがウェットエリアならガウスフィルタによる拡散された値を周囲に足し込む
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    sum = 0;
                    filter_sum = 1;
                    for(l=-c; l<=c; l++){
                        for(k=-c; k<=c; k++){
                            //フィルタの端ピクセルがない場合、分子には加算せず分母から減算
                            if( (i+k)<0 || (i+k)>stroke_length_x || (j+l)<0 || (j+l)>stroke_length_y){
                                filter_sum -= dev_gauss_filter[(k+c)+(l+c)*w];
                            }else{
                                sum += dev_M[SubImageIndex+(i+k)+(j+l)*stroke_length_max] * dev_gauss_filter[(k+c)+(l+c)*w];
                            }
                        }
                    }
                    dev_gauss_M[SubImageIndex+i+j*stroke_length_max] = sum / filter_sum;
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_p[SubImageIndex+i+j*stroke_length_max] = dev_p[SubImageIndex+i+j*stroke_length_max] - opt_eta * opt_SoakTimeStep * (1-dev_gauss_M[SubImageIndex+i+j*stroke_length_max])*dev_M[SubImageIndex+i+j*stroke_length_max];
                }
            }
        }
        __syncthreads();

        //////////MovePigment//////////
        for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
            for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                dev_new_gR[SubImageIndex+i+j*stroke_length_max] = dev_gR[SubImageIndex+i+j*stroke_length_max]; //コピー
                dev_new_gG[SubImageIndex+i+j*stroke_length_max] = dev_gG[SubImageIndex+i+j*stroke_length_max]; //コピー
                dev_new_gB[SubImageIndex+i+j*stroke_length_max] = dev_gB[SubImageIndex+i+j*stroke_length_max]; //コピー
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_new_gR[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gR[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gG[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gG[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gB[SubImageIndex+(i+1)+j*stroke_length_max] = dev_new_gB[SubImageIndex+(i+1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_new_gR[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gR[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gG[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gG[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gB[SubImageIndex+(i-1)+j*stroke_length_max] = dev_new_gB[SubImageIndex+(i-1)+j*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_new_gR[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gR[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gG[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gG[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gB[SubImageIndex+i+(j+1)*stroke_length_max] = dev_new_gB[SubImageIndex+i+(j+1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_new_gR[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gR[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gG[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gG[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gB[SubImageIndex+i+(j-1)*stroke_length_max] = dev_new_gB[SubImageIndex+i+(j-1)*stroke_length_max] + fmaxf(0, opt_SoakTimeStep*(-1)*dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                    dev_new_gR[SubImageIndex+i+j*stroke_length_max] = dev_new_gR[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gG[SubImageIndex+i+j*stroke_length_max] = dev_new_gG[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_new_gB[SubImageIndex+i+j*stroke_length_max] = dev_new_gB[SubImageIndex+i+j*stroke_length_max] - fmaxf(0, opt_SoakTimeStep*dev_u[SubImageIndex_uv+(i+1)+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_u[SubImageIndex_uv+i+j*(stroke_length_max+1)]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*dev_v[SubImageIndex_uv+i+(j+1)*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]) - fmaxf(0, opt_SoakTimeStep*-dev_v[SubImageIndex_uv+i+j*stroke_length_max]*dev_gB[SubImageIndex+i+j*stroke_length_max]);
                }
            }
        }
        __syncthreads();

        for(j=threadIdx.y; j<stroke_length_max; j+=blockDim.y){
            for(i=threadIdx.x; i<stroke_length_max; i+=blockDim.x){
                dev_gR[SubImageIndex+i+j*stroke_length_max] = dev_new_gR[SubImageIndex+i+j*stroke_length_max]; //コピー
                dev_gG[SubImageIndex+i+j*stroke_length_max] = dev_new_gG[SubImageIndex+i+j*stroke_length_max]; //コピー
                dev_gB[SubImageIndex+i+j*stroke_length_max] = dev_new_gB[SubImageIndex+i+j*stroke_length_max]; //コピー
            }
        }
        __syncthreads();

        //////////TransferPigment//////////
        for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
            for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
                if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){

                    if(opt_USE_DETAIL_TP){
                        down_ratio = opt_SoakTimeStep * (1-dev_h[SubImageIndex+i+j*stroke_length_max]*opt_exposure) * opt_deposit;
                        up_ratio = opt_SoakTimeStep * (1+(dev_h[SubImageIndex+i+j*stroke_length_max]-1)*opt_exposure) * opt_lift;
                    }else{
                        down_ratio = opt_SoakTimeStep * (1-dev_h[SubImageIndex+i+j*stroke_length_max]*opt_gamma) * opt_rho;
                        up_ratio = opt_SoakTimeStep * (1+(dev_h[SubImageIndex+i+j*stroke_length_max]-1)*opt_gamma) * opt_rho / opt_omega;
                    }

                    //R
                    down = dev_gR[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                    up   = dev_dR[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                    if(dev_dR[SubImageIndex+i+j*stroke_length_max]+down > 1)
                        down = fmaxf(0,1-dev_dR[SubImageIndex+i+j*stroke_length_max]);
                    if(dev_gR[SubImageIndex+i+j*stroke_length_max]+up > 1)
                        up = fmaxf(0,1-dev_gR[SubImageIndex+i+j*stroke_length_max]);
                    dev_dR[SubImageIndex+i+j*stroke_length_max] = dev_dR[SubImageIndex+i+j*stroke_length_max]+down-up;
                    dev_gR[SubImageIndex+i+j*stroke_length_max] = dev_gR[SubImageIndex+i+j*stroke_length_max]+up-down;

                    //G
                    down = dev_gG[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                    up   = dev_dG[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                    if(dev_dG[SubImageIndex+i+j*stroke_length_max]+down > 1)
                        down = fmaxf(0,1-dev_dG[SubImageIndex+i+j*stroke_length_max]);
                    if(dev_gG[SubImageIndex+i+j*stroke_length_max]+up > 1)
                        up = fmaxf(0,1-dev_gG[SubImageIndex+i+j*stroke_length_max]);
                    dev_dG[SubImageIndex+i+j*stroke_length_max] = dev_dG[SubImageIndex+i+j*stroke_length_max]+down-up;
                    dev_gG[SubImageIndex+i+j*stroke_length_max] = dev_gG[SubImageIndex+i+j*stroke_length_max]+up-down;

                    //B
                    down = dev_gB[SubImageIndex+i+j*stroke_length_max] * down_ratio;
                    up   = dev_dB[SubImageIndex+i+j*stroke_length_max] * up_ratio;
                    if(dev_dB[SubImageIndex+i+j*stroke_length_max]+down > 1)
                        down = fmaxf(0,1-dev_dB[SubImageIndex+i+j*stroke_length_max]);
                    if(dev_gB[SubImageIndex+i+j*stroke_length_max]+up > 1)
                        up = fmaxf(0,1-dev_gB[SubImageIndex+i+j*stroke_length_max]);
                    dev_dB[SubImageIndex+i+j*stroke_length_max] = dev_dB[SubImageIndex+i+j*stroke_length_max]+down-up;
                    dev_gB[SubImageIndex+i+j*stroke_length_max] = dev_gB[SubImageIndex+i+j*stroke_length_max]+up-down;
                }
            }
        }
    }
    __syncthreads();

    //堆積顔料をRGBに変換しキャンバスに描画
    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){          

            // シミュレーション終了時に水中の顔料を落とすか拭き取るか
            if(opt_FloatPigmentOnPaper){
                s_test_Canvas_R[i+j*stroke_length_max] = (1 - dev_dR[SubImageIndex+i+j*stroke_length_max]) * 255;    //CMY[0,1]->RGB[0,255]
                s_test_Canvas_G[i+j*stroke_length_max] = (1 - dev_dG[SubImageIndex+i+j*stroke_length_max]) * 255;
                s_test_Canvas_B[i+j*stroke_length_max] = (1 - dev_dB[SubImageIndex+i+j*stroke_length_max]) * 255;
            }
            else{
                //R
                pigment_density = dev_gR[SubImageIndex+i+j*stroke_length_max] + dev_dR[SubImageIndex+i+j*stroke_length_max];
                LIMIT_RANGE(pigment_density, 0, 1);
                s_test_Canvas_R[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_R[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_R[x+y*width];

                //G
                pigment_density = dev_gG[SubImageIndex+i+j*stroke_length_max] + dev_dG[SubImageIndex+i+j*stroke_length_max];
                LIMIT_RANGE(pigment_density, 0, 1);
                s_test_Canvas_G[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_G[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_G[x+y*width];

                //B
                pigment_density = dev_gB[SubImageIndex+i+j*stroke_length_max] + dev_dB[SubImageIndex+i+j*stroke_length_max];
                LIMIT_RANGE(pigment_density, 0, 1);
                s_test_Canvas_B[i+j*stroke_length_max] = (1-pigment_density) * s_test_Canvas_B[i+j*stroke_length_max] + pigment_density * dev_best_stroke_map_B[x+y*width];
            }
        }
    }
    __syncthreads();

    for(j=threadIdx.y; j<=stroke_length_y; j+=blockDim.y){
        for(i=threadIdx.x; i<=stroke_length_x; i+=blockDim.x){
            if(dev_M[SubImageIndex+i+j*stroke_length_max] == 1){
                //キャンバスに部分画像をコピー
                dev_nimgR[(left_end+i)+(upper_end+j)*width] = s_test_Canvas_R[i+j*stroke_length_max];
                dev_nimgG[(left_end+i)+(upper_end+j)*width] = s_test_Canvas_G[i+j*stroke_length_max];
                dev_nimgB[(left_end+i)+(upper_end+j)*width] = s_test_Canvas_B[i+j*stroke_length_max];
            }
        }
    }
    __syncthreads();
}

__global__ void gpu_reset_improved_value_map(int *dev_GLOBAL_improved_value_map, int *dev_best_stroke_map_pnum, float *dev_best_stroke_map_point_x, float *dev_best_stroke_map_point_y, int *dev_best_x, int *dev_best_y, int width, int height, int t){

    int distance;
    int best_x = dev_best_x[0];
    int best_y = dev_best_y[0];
    int pnum = dev_best_stroke_map_pnum[best_x+best_y*width];

    //ストローク点を囲む端の座標を特定
    int left_end = dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+0];
    int right_end = dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+0];
    int upper_end = dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+0];
    int lower_end = dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+0];
    for(int i=1; i<pnum; i++){
        if(dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i] < left_end) left_end = dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i];
        if(right_end < dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i]) right_end = dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i];
        if(dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i] < upper_end) upper_end = dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i];
        if(lower_end < dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i]) lower_end = dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i];
    }

    //端座標を膨張
    left_end  -= t*opt_max_stroke + t+2;
	right_end += t*opt_max_stroke + t+2;
	upper_end -= t*opt_max_stroke + t+2;
	lower_end += t*opt_max_stroke + t+2;
    if(left_end<0) left_end = 0;
	if(width <= right_end) right_end = width-1;
	if(upper_end<0) upper_end = 0;
	if(height <= lower_end) lower_end = height-1;

    for(int y=upper_end; y<=lower_end; y++) {
		for(int x=left_end; x<=right_end; x++) {
			if(dev_GLOBAL_improved_value_map[x+y*width] == UNCALCULATED) continue;
			// 描画したストローク制御点と重なるストロークのみを再計算
			for(int i=0; i<pnum; i++){
				for(int j=0; j<dev_best_stroke_map_pnum[x+y*width]; j++) {
					distance = (dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+j])
                                * (dev_best_stroke_map_point_x[(best_x+best_y*width)*opt_max_stroke+i] - dev_best_stroke_map_point_x[(x+y*width)*opt_max_stroke+j])
                                + (dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+j])
                                * (dev_best_stroke_map_point_y[(best_x+best_y*width)*opt_max_stroke+i] - dev_best_stroke_map_point_y[(x+y*width)*opt_max_stroke+j]);	// 点と点のユークリッド距離
					if(distance < 2*t*2*t){
						dev_GLOBAL_improved_value_map[x+y*width] = UNCALCULATED;
						goto RIM_loopend;
					}
				}
			}
			RIM_loopend: ;
		}
	}
}