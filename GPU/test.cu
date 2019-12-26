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

void printCudaLastError(){
	cudaError_t err = cudaGetLastError();
	printf("cudaGetLastError::%s(code:%d)\n",cudaGetErrorString(err),err);
	if(err)	exit(0);
}

__global__ void test(){
	__shared__ float s_diff_sum[1];
    __shared__ float s_diff_R[1];
    __shared__ float s_diff_G[1];
    __shared__ float s_diff_B[1];
	__shared__ int s_offscrn_count[1];
    __shared__ float s_theta[1];
    __shared__ float s_histogram[opt_histogram_partition];
    __shared__ float s_error_sum[1];
    __shared__ char s_test_Canvas_R[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ char s_test_Canvas_G[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ char s_test_Canvas_B[(opt_thick_max*(opt_max_stroke+2))*(opt_thick_max*(opt_max_stroke+2))];
    __shared__ float s_max_velocity[1];
    __shared__ float s_delta_MAX[1];

	for(int y=blockIdx.y; y<1; y=y+gridDim.y) {
		for(int x=blockIdx.x; x<1; x=x+gridDim.x) {

			if((threadIdx.x==0)&&(threadIdx.y==0)){
				s_diff_sum[0] = 0;
				s_diff_R[0] = 0;
				s_diff_G[0] = 0;
				s_diff_B[0] = 0;
				s_offscrn_count[0] = 0;
			}
			__syncthreads();

			if((threadIdx.x==0)&&(threadIdx.y==0)){
				printf("s_offscrn_count = %d\n",s_offscrn_count[0]);
				printf("s_diff_R = %f\n",s_diff_R[0]);
				printf("s_diff_G = %f\n",s_diff_G[0]);
				printf("s_diff_B = %f\n",s_diff_B[0]);
				printf("s_diff_sum = %f\n",s_diff_sum[0]);
			}
		}
	}
}

int main(){

	dim3 block_num(1, 1);
	dim3 thread_num(32, 32);

	test<<<block_num,thread_num>>>();
	cudaDeviceSynchronize();
	printCudaLastError();
}