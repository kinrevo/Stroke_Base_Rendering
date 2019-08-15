#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include "sbr.h"
#include "water.h"

// #include <mcheck.h>


void ufvf_Test(int argc, char *argv[]);
double UpdateVelocity_Test(int argc, char *argv[]);
double RelaxDivergence_Test(int argc, char *argv[]);
void FlowOutward_Test(int argc, char *argv[]);
double MovePigment_Test(int argc, char *argv[]);
double TransferPigment_Test(int argc, char *argv[]);
void calcu_grad_h_Test(int argc, char *argv[]);
double Paint_Water_Test(int argc, char *argv[]);
double Circle_fill_Water_Test(int argc, char *argv[]);
double set_WetStroke_Test(int argc, char *argv[]);
double Paint_Water_Stroke_Test(int argc, char *argv[]);
void SimulateCapillaryFlow_Test(int argc, char *argv[]);
void SimulateCapillaryFlow_PaintTest(int argc, char *argv[]);



#define p(s,a) printf("%s:%d\n",s,a)


int main(int argc, char *argv[]){
    int i, trials=10;
    double tmp, Ave_TIME=0, min=9999, max=0;
    my_clock();

    #ifdef _OPENMP
        omp_set_num_threads(atoi(argv[2]));
    #endif

    for (i = 0; i < trials; i++) {
        tmp = Paint_Water_Stroke_Test(argc, argv);
        Ave_TIME += tmp;
        if(max<tmp) max=tmp;
        if(min>tmp) min=tmp;
    }
    pn;

	pd("All_Execution_TIME", my_clock());
    printf("Ave_TIME:%f, max:%f, min:%f\n",Ave_TIME/trials, max, min);
    return 0;
}



// int main(int argc, char *argv[])
// {
// 	clock_t start = clock();
	
// 	PPM *in_ppm, *trans_ppm;
// 	char *name;
// 	char *ext;
// 	name = argv[1];
	
//     in_ppm = create_ppm(width, height, bright);

// 	//水彩描画
// 	trans_ppm = (in_ppm);
	
// 	//出力ファイル名に従って画像を出力
// 	ext = get_extension(argv[1]);
// 	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
// 		if(write_ppm(argv[1], trans_ppm)){ printf("WRITE_PPM_ERROR (main)\n");}
// 	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
// 		if(write_jpeg_file(argv[1], PPM_to_image(trans_ppm))){ printf("WRITE JPG ERROR.");}
// 	} else if (strcmp("png", ext) == 0) {
// 		if(write_png_file(argv[1], PPM_to_image(trans_ppm))){ printf("WRITE PNG ERROR.");}
// 	}
	
// 	FreePPM(trans_ppm);
	
// 	pd("TOTAL_TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
// 	return 0;
// }


/*
void ufvf_Test(int argc, char *argv[])
{
    int i,j;
    double tmp;
    int width, height;
    width=height=512;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    
	// uf
    for (i = 0; i < width+1; i++) {  //doubleによる添字の実験
        for (j = 0; j < height; j++) {
            u[i][j] = i+j;
        }
    }
    for (j = 0; j < height; j++) {  //uf(u, i+0.5, j)
        for (i = 0; i < width+1; i++){
            tmp = uf(u, i-0.5, j);
            // printf("%5.1f  ", tmp);
        }
        // pn;
    }
    // pn;
    for (j = 0; j < height-1; j++) {  //uf(u, i+0.5, j+0.5)
        for (i = 0; i < width+1; i++){
            tmp = uf(u, i-0.5, j+0.5);
            // printf("%5.2f  ", tmp);
        }
        // pn;
    }
    // pn;

	// vf
    for (i = 0; i < width; i++) {
        for (j = 0; j < height+1; j++) {
            v[i][j] = i+j;
        }
    }
    for (j = 0; j < height+1; j++) {
        for (i = 0; i < width; i++){
            tmp = vf(v, i, j-0.5);
            // printf("%5.1f  ", tmp);
        }
        // pn;
    }
    // pn;
    for (j = 0; j < height; j++) {
        for (i = 0; i < width; i++){
            tmp = vf(v, i, j);
            // printf("%5.2f  ", tmp);
        }
        // pn;
    }
    
    Free_dally(u, width+1);
    Free_dally(v, width);
}
*/

double UpdateVelocity_Test(int argc, char *argv[])
{
    struct timespec start,end;
    int i,j;
    int width=128, height=128;
    double t;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    PPM* u_img = create_ppm(width+1, height, 255);
    char filename[64]={};
    // char tmp_name[16]={};

    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            v[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            M[i][j] = 1;
            // if(i<width/2 && j<height/2){
            //     M[i][j] = 1;
            // }  else M[i][j]=0;
            p[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
        }
    }

    strcpy(filename, "0_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);

	// UpdateVelocitieの変化の推移を出力
    double var_t = 0.5;
    clock_gettime(CLOCK_MONOTONIC, &start);
    for (t = 0; t < 5; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);
        // if((int)(t/var_t) % 100 == 0 ){
        //     printf("%03d\n", (int)(t/var_t));
        //     snprintf(tmp_name, 16, "%03d", (int)(t/var_t));
        //     strcpy(filename, tmp_name);
        //     strcat(filename, argv[1]);
        //     trans_Vector_img(u_img, u, width+1, height);
        //     write_ppm(filename, u_img);
        // }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    double Ex_TIME = (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
	pd("UV_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}


double RelaxDivergence_Test(int argc, char *argv[])
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int width=128, height=128;
    Point rectangleP[2] = {
        {0,128},
        {0,128}
    };
    double t;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    PPM* u_img = create_ppm(width+1, height, 255);
    PPM* p_img = create_ppm(width, height, 255);
    char filename[64]={};
    char tmp_name[16]={};


    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = -(128-abs(64-i)-abs(64-j))/128.0;
            v[i][j] = -(128-abs(64-i)-abs(64-j))/128.0;
            M[i][j] = 1;
            if(i<width/2 && j<height/2){
                M[i][j] = 1;
            }  else M[i][j]=0;
            p[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0;
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "0p_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);

    double var_t = 0.1;
    for (t = 0; t < 10; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        // if((int)(t/var_t) % 100 == 0 ){
        //     printf("%03d\n", (int)(t/var_t));
        //     snprintf(tmp_name, 16, "%03dUV", (int)(t/var_t));
        //     strcpy(filename, tmp_name);
        //     strcat(filename, argv[1]);
        //     trans_Vector_img(u_img, u, width-1, height);
        //     write_ppm(filename, u_img);
        // }

        clock_gettime(CLOCK_MONOTONIC, &start);
        RelaxDivergence(M, u, v, p, var_t, width, height, rectangleP);		// RelaxDivergenceの変化の推移を出力
        clock_gettime(CLOCK_MONOTONIC, &end);
        Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dRD", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            trans_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);

            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "p%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            trans_Vector_img(p_img, p, width, height);
            write_ppm(filename, p_img);
        }
    }
     
	pd("RD_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}


void FlowOutward_Test(int argc, char *argv[]){
    // struct timespec start,end;
    // double Ex_TIME=0;
    int i,j;
    int width=128, height=128;
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);

    int K=opt_K;
    double eta=opt_eta;
    double** gauss_M = create_dally(width, height);
    double** revers_gauss_M = create_dally(width, height);
    double** remove_p = create_dally(width, height);
    PPM* fig_img = create_ppm(width, height, 256);
    
    for (i = 0; i < width; i++) {	//pの初期化
        for (j = 0; j < height; j++) {
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            if(i<width/2 && j<height/2){
                M[i][j] = 1;
            }  else M[i][j]=0;
        }
    }

    double** dM = create_dally(width, height);  //関数に渡すためにdouble型に変換
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    
    int w = (int)( ceil(3.0*K/6.0+0.5)*2-1 );
	int c=(w-1)/2;
	double** filter = create_dally(w, w);
	
	#ifdef _OPENMP
		#pragma omp parallel for
    #endif
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, K/6.0);
        }
    }
    gauss_M = gaussian_filter_d(dM, c, filter, width, height);

    trans_Vector_img(fig_img, gauss_M, width, height);
    write_ppm("gaussM.ppm", fig_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            revers_gauss_M[i][j] = 1-gauss_M[i][j];
        }
    }

    trans_Vector_img(fig_img, revers_gauss_M, width, height);
    write_ppm("revers_gauss_M.ppm", fig_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            remove_p[i][j] = eta*revers_gauss_M[i][j]*M[i][j];
        }
    }

    trans_Vector_img(fig_img, remove_p, width, height);
    write_ppm("remove_p.ppm", fig_img);
    trans_Vector_img(fig_img, p, width, height);
    write_ppm("orig_p.ppm", fig_img);

    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            p[i][j] = p[i][j] - remove_p[i][j]/eta;
        }
    }

    trans_Vector_img(fig_img, p, width, height);
    write_ppm("res_p.ppm", fig_img);
}


double MovePigment_Test(int argc, char *argv[]){
    struct timespec start,end;
    double UV_TIME=0, RD_TIME=0, FO_TIME=0, MP_TIME=0;
    int i,j;
    int width=128, height=128;
    Point rectangleP[2] = {
        {0,128},
        {0,128}
    };
    double t;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    PPM* u_img = create_ppm(width+1, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    char filename[64]={};
    char tmp_name[16]={};

	int w = (int)( ceil(3.0*opt_K/6+0.5)*2-1 ); //とりあえず動く計算
	int c=(w-1)/2;
	double** filter = create_dally(w, w);	
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, opt_K/6.0);
        }
    }

    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            v[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>0 && j>0 && i<width/1.5 && j<height/1.5){
                M[i][j] = 1;
            }  else M[i][j]=0;
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            if(i!=0 && j!=0) gR[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0; 
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "MPO_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            paint_img->dataR[i][j] = fmax(255-gR[i][j]*255, 0);
            paint_img->dataG[i][j] = fmax(255-gG[i][j]*255, 0);
            paint_img->dataB[i][j] = fmax(255-gB[i][j]*255, 0);
        }
    }
    write_ppm(filename, paint_img);

    double var_t = 0.1;
    for (t = 0; t < 10; t=t+var_t) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        clock_gettime(CLOCK_MONOTONIC, &end);
        UV_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

        clock_gettime(CLOCK_MONOTONIC, &start);
        RelaxDivergence(M, u, v, p, var_t, width, height, rectangleP);		// RelaxDivergenceの変化の推移を出力
        clock_gettime(CLOCK_MONOTONIC, &end);
        RD_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

        clock_gettime(CLOCK_MONOTONIC, &start);
        FlowOutward(M, p, c, filter, var_t, width, height, rectangleP);
        clock_gettime(CLOCK_MONOTONIC, &end);
        FO_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dFO", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            trans_Vector_img(u_img, p, width-1, height);
            write_ppm(filename, u_img);
        }
 
        clock_gettime(CLOCK_MONOTONIC, &start);
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height, rectangleP);	// MovePigmentの変化の推移を出力
        clock_gettime(CLOCK_MONOTONIC, &end);
        MP_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "MP%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            for(i=0; i<width; i++){
                for(j=0; j<height; j++){
                    paint_img->dataR[i][j] = fmax(255-gR[i][j]*255, 0);
                    paint_img->dataG[i][j] = fmax(255-gG[i][j]*255, 0);
                    paint_img->dataB[i][j] = fmax(255-gB[i][j]*255, 0);
                }
            }
            write_ppm(filename, paint_img);
        }
    }
         
	pd("UV_Execution_TIME", UV_TIME);
	pd("RD_Execution_TIME", RD_TIME);
	pd("FO_Execution_TIME", FO_TIME);
	pd("MP_Execution_TIME", MP_TIME);
    return MP_TIME;
}



double TransferPigment_Test(int argc, char *argv[])
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int width=128, height=128;
    Point rectangleP[2] = {
        {0,128},
        {0,128}
    };
    double t;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_dally(p, width, height, 0);
    double** h = create_dally(width, height);
    h = perlin_img(width, height, 0.1, 4);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    double** dR = create_dally(width, height);
    double** dG = create_dally(width, height);
    double** dB = create_dally(width, height);
    format_dally(dR, width, height, 0);
    format_dally(dG, width, height, 0);
    format_dally(dB, width, height, 0);
    PPM* u_img = create_ppm(width+1, height, 255);
    PPM* p_img = create_ppm(width, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    char filename[64]={};
    char tmp_name[16]={};

	int w = (int)( ceil(3.0*opt_K/6+0.5)*2-1 ); //とりあえず動く計算
	int c=(w-1)/2;
	double** filter = create_dally(w, w);	
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, opt_K/6.0);
        }
    }

    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            v[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>30 && j>30 && i<90 && j<90){
                M[i][j] = 1;
                p[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            }  else M[i][j]=0;
            if(i!=0 && j!=0) gR[i][j] = 0.7;//(128-abs(64-i)-abs(64-j))/128.0; 
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "MPO_");		//顔料ｇをCMYで出力
    strcat(filename, argv[1]);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            paint_img->dataR[i][j] = fmax(255-gR[i][j]*255, 0);
            paint_img->dataG[i][j] = fmax(255-gG[i][j]*255, 0);
            paint_img->dataB[i][j] = fmax(255-gB[i][j]*255, 0);
        }
    }
    write_ppm(filename, paint_img);

    double var_t = 0.1;
    for (t = 0; t < 10; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        RelaxDivergence(M, u, v, p, var_t, width, height, rectangleP);		// RelaxDivergenceの変化の推移を出力
        FlowOutward(M, p, c, filter, var_t, width, height, rectangleP);
        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dUV", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            trans_Vector_img(u_img, u, width+1, height);
            write_ppm(filename, u_img);
        }
        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "p%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            trans_Vector_img(p_img, p, width, height);
            write_ppm(filename, p_img);
        }
 
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height, rectangleP);	// MovePigmentの変化の推移を出力
        // if((int)(t/var_t) % 10 == 0 ){
        //     printf("%03d\n", (int)(t/var_t));
        //     snprintf(tmp_name, 16, "MP%03d", (int)(t/var_t));
        //     strcpy(filename, tmp_name);
        //     strcat(filename, argv[1]);
        //     for(i=0; i<width; i++){
        //         for(j=0; j<height; j++){
        //             paint_img->dataR[i][j] = fmax(255-gR[i][j]*255, 0);
        //             paint_img->dataG[i][j] = fmax(255-gG[i][j]*255, 0);
        //             paint_img->dataB[i][j] = fmax(255-gB[i][j]*255, 0);
        //         }
        //     }
        //     write_ppm(filename, paint_img);
        // }
        clock_gettime(CLOCK_MONOTONIC, &start);
        TransferPigment(M, h, gR, gG, gB, dR, dG, dB, var_t, width, height); // TransferPigmentの変化の推移を出力
        clock_gettime(CLOCK_MONOTONIC, &end);
        Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
        if((int)(t/var_t) % 10 == 0 ){
            // printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "TP%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            for(i=0; i<width; i++){
                for(j=0; j<height; j++){
                    paint_img->dataR[i][j] = fmax(255-dR[i][j]*255, 0);
                    paint_img->dataG[i][j] = fmax(255-dG[i][j]*255, 0);
                    paint_img->dataB[i][j] = fmax(255-dB[i][j]*255, 0);
                }
            }
            write_ppm(filename, paint_img);
        }
    }
         
	pd("TP_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}



void calcu_grad_h_Test(int argc, char *argv[])
{
    int i,j;
    int width=8, height=8;
    double** h = create_dally(width, height); 
    double** grad_hx = create_dally(width+1, height); 
    double** grad_hy = create_dally(width, height+1); 

    for (i = 0; i < width; i++) {  //doubleによる添字の実験
        for (j = 0; j < height; j++) {
            h[i][j] = (i+1)*(j+1);
        }
    }
    calcu_grad_h(h, grad_hx, grad_hy, width, height);

    for (j = 0; j < height; j++) {  //uf(u, i+0.5, j)
        for (i = 0; i < width; i++){
            printf("%5.1f  ", h[i][j]);
        }
        pn;
    }pn;
    for (j = 0; j < height; j++) {  //uf(u, i+0.5, j+0.5)
        for (i = 0; i < width+1; i++){
            printf("%5.2f  ", grad_hx[i][j]);
        }
        pn;
    }
    pn;
    for (j = 0; j < height+1; j++) {  //uf(u, i+0.5, j+0.5)
        for (i = 0; i < width; i++){
            printf("%5.2f  ", grad_hy[i][j]);
        }
        pn;
    }
    pn;
}



double Paint_Water_Test(int argc, char *argv[])
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int width=640, height=430;
    Point rectangleP[2] = {
        {0,640},
        {0,430}
    };
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_dally(p, width, height, 0);
    double** h = create_dally(width, height);
    h = perlin_img(width, height, 0.1, 4);
    double** grad_hx = create_dally(width+1, height); 
    double** grad_hy = create_dally(width, height+1); 
    calcu_grad_h(h, grad_hx, grad_hy, width, height);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    double** dR = create_dally(width, height);
    double** dG = create_dally(width, height);
    double** dB = create_dally(width, height);
    format_dally(dR, width, height, 0);
    format_dally(dG, width, height, 0);
    format_dally(dB, width, height, 0);
    PPM* u_img = create_ppm(width+1, height, 255);
    PPM* p_img = create_ppm(width, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    char filename[64]={};


    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            v[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>30 && j>30 && i<90 && j<90){
                M[i][j] = 1;
                p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            }  else M[i][j]=0;
            if(i!=0 && j!=0) gR[i][j] = (128-abs(64-i)-abs(64-j))/128.0; 
        }
    }

    strcpy(filename, "u0_");		//uの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "p0_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);
    strcpy(filename, "dO_");		//顔料ｇをCMYで出力
    strcat(filename, argv[1]);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            paint_img->dataR[i][j] = fmax(255-dR[i][j]*255, 0);
            paint_img->dataG[i][j] = fmax(255-dG[i][j]*255, 0);
            paint_img->dataB[i][j] = fmax(255-dB[i][j]*255, 0);
        }
    }
    write_ppm(filename, paint_img);

    clock_gettime(CLOCK_MONOTONIC, &start);
    Paint_Water(M, u, v, p, h, grad_hx, grad_hy, gR, gG, gB, dR, dG, dB, width, height, rectangleP);
    clock_gettime(CLOCK_MONOTONIC, &end);
    Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

    strcpy(filename, "u1_");		//uの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(u_img, u, width+1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "p1_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);
    strcpy(filename, "d1_");		//顔料ｇをCMYで出力
    strcat(filename, argv[1]);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            paint_img->dataR[i][j] = fmax(255-dR[i][j]*255, 0);
            paint_img->dataG[i][j] = fmax(255-dG[i][j]*255, 0);
            paint_img->dataB[i][j] = fmax(255-dB[i][j]*255, 0);
        }
    }
    write_ppm(filename, paint_img);
             
	pd("PW_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}


double Circle_fill_Water_Test(int argc, char *argv[]) 
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int t=20;
    RGB color = {200, 0, 100};
    double sigma;
    int width=128, height=128;
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_ally(M, width, height, 0);
    format_dally(p, width, height, 0);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    PPM* p_img = create_ppm(width, height, 255);
    Point SP = {0,0};//{width/2, height/2};
    char filename[64]={};

    double** gauce_filter = create_dally(2*t+1, 2*t+1);
    sigma = t/3.0;
    for(i=0; i<2*t+1; i++){
        for(j=0; j<2*t+1; j++){
            gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
        }
    }
    
    double** dM = create_dally(width, height);  //関数に渡すためにdouble型に変換

    strcpy(filename, "p0_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);
    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    strcpy(filename, "M0_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, dM, width, height);
    write_ppm(filename, p_img);

    clock_gettime(CLOCK_MONOTONIC, &start);
    Circle_fill_Water(M, p, gR, gG, gB, SP, t, color, gauce_filter, width, height);
    clock_gettime(CLOCK_MONOTONIC, &end);
    Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            p[i][j] = p[i][j] * 2*t;    //ブラシを動かして塗り重ねないと水量が少なすぎてintだと消える
        }
    }
    strcpy(filename, "p1_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    strcpy(filename, "M1_");		//pの初期速度を画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, dM, width, height);
    write_ppm(filename, p_img);          

	pd("CfW_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}



double set_WetStroke_Test(int argc, char *argv[])
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int t=20;
    int pnum=5;
    RGB color = {200, 0, 100};
    double sigma;
    int width=128, height=128;
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_ally(M, width, height, 0);
    format_dally(p, width, height, 0);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    PPM* p_img = create_ppm(width, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    Point SP[5]= {
        {10,100},
        {20,50},
        {50,20},
        {100,40},
        {120,80}
    };
    char filename[64]={};

    double** gauce_filter = create_dally(2*t+1, 2*t+1);
    sigma = t/3.0;
    for(i=0; i<2*t+1; i++){
        for(j=0; j<2*t+1; j++){
            gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
        }
    }
    
    double** dM = create_dally(width, height);  //関数に渡すためにdouble型に変換

    strcpy(filename, "p0_");		//pを画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);
    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    strcpy(filename, "M0_");		//Mを画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, dM, width, height);
    write_ppm(filename, p_img);

    clock_gettime(CLOCK_MONOTONIC, &start);
    set_WetStroke(M, p, gR, gG, gB, SP, pnum, t, color, gauce_filter, width, height);
    clock_gettime(CLOCK_MONOTONIC, &end);
    Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            p[i][j] = p[i][j];
        }
    }
    strcpy(filename, "p1_");		//pを画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    strcpy(filename, "M1_");		//Mを画像として出力
    strcat(filename, argv[1]);
    trans_Vector_img(p_img, dM, width, height);
    write_ppm(filename, p_img);

    strcpy(filename, "g1_");		//顔料ｇをCMYで出力
    strcat(filename, argv[1]);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            paint_img->dataR[i][j] = fmax(255-gR[i][j]*255, 0);
            paint_img->dataG[i][j] = fmax(255-gG[i][j]*255, 0);
            paint_img->dataB[i][j] = fmax(255-gB[i][j]*255, 0);
        }
    }
    write_ppm(filename, paint_img);
                 
	pd("sWS_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}



double Paint_Water_Stroke_Test(int argc, char *argv[]) 
{
    struct timespec start,end;
    double Ex_TIME=0;
    int i,j;
    int t=10;
    int pnum=6;
    RGB color = {91, 15, 12}, color2={250,20,70};
    // RGB color = {200, 0, 100}, color2={0,0,0};
    double sigma;
    // int width=460, height=360;
    int width=640, height=480;
    double** h = perlin_img(width, height, 0.1, 4);
    double** grad_hx = create_dally(width+1, height); 
    double** grad_hy = create_dally(width, height+1); 
    calcu_grad_h(h, grad_hx, grad_hy, width, height);
    PPM* Canvas_img = create_ppm(width, height, 255);
    Point SP1[6]= {
        {391.5,346.5},
        {400.9,350},
        {410.9,349.5},
        {420.9,349.0},
        {430.8,348.5},
        {440.8,347.9}
    };
    // Point SP1[5]= {
    //     {10,10},
    //     {40,40},
    //     {60,60},
    //     {100,100},
    //     {120,120}
    // };
    Point SP2[6]= {
        {110,10},
        {80,40},
        {60,60},
        {40,100},
        {10,120},
        {0,150}
    };
    char filename[64]={};

    double** gauce_filter = create_dally(2*t+1, 2*t+1);
    sigma = t/3.0;
    for(i=0; i<2*t+1; i++){
        for(j=0; j<2*t+1; j++){
            gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
        }
    }

    strcpy(filename, "C0_");    // ペイント前のキャンバスを出力
    strcat(filename, argv[1]);
    write_ppm(filename, Canvas_img);

    clock_gettime(CLOCK_MONOTONIC, &start);
    Paint_Water_Stroke(SP1, pnum, t, color, Canvas_img->dataR, Canvas_img->dataG, Canvas_img->dataB, h, grad_hx, grad_hy, gauce_filter, width, height);
    clock_gettime(CLOCK_MONOTONIC, &end);
    Ex_TIME += (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;

    strcpy(filename, "C1_");    // ペイント1後のキャンバスを出力
    strcat(filename, argv[1]);
    write_ppm(filename, Canvas_img);

    Paint_Water_Stroke(SP2, pnum, t, color2, Canvas_img->dataR, Canvas_img->dataG, Canvas_img->dataB, h, grad_hx, grad_hy, gauce_filter, width, height);

    strcpy(filename, "C2_");    // ペイント2後のキャンバスを出力
    strcat(filename, argv[1]);
    write_ppm(filename, Canvas_img);

    Free_dally(h,width);
    Free_dally(grad_hx,width+1);
    Free_dally(grad_hy,width);
    Free_dally(gauce_filter,2*t+1);
    FreePPM(Canvas_img);
                 
	pd("PWS_Execution_TIME", Ex_TIME);
    return Ex_TIME;
}



void SimulateCapillaryFlow_Test(int argc, char *argv[])
{
    int i,j;
    int width=128, height=128;
    double t;
    double var_t=0.5;
    int** M = create_ally(width, height);
    double** dM = create_dally(width, height);
    double** p = create_dally(width, height);
    double** h = perlin_img(width, height, 0.1, 4);
    double** s = create_dally(width, height);
    format_dally(s, width, height, 0);

    PPM* fig_img = create_ppm(width, height, 255);
    char count_name[8];
    char out_name[32];
    
    for (i = 0; i < width; i++) {	//pの初期化
        for (j = 0; j < height; j++) {
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            if(i<width/2 && j<height/2){
                M[i][j] = 1;
            }  else M[i][j]=0;
        }
    }

    for (i = 0; i < width; i++) {	
        for (j = 0; j < height; j++) {
            dM[i][j] = M[i][j];
        }
    }

    trans_Vector_img(fig_img, dM, width, height);
    write_ppm("WetAreaO.ppm", fig_img);

    for (t = 0; t < 10; t+=var_t){
        SimulateCapillaryFlow(M, p, h, s, var_t, width, height);
        
        snprintf(count_name, 16, "%02d", (int)(t*2));

        strcpy(out_name, "capillaryWater");
        strcat(out_name, count_name);
        strcat(out_name, ".ppm");
        trans_Vector_img(fig_img, s, width, height);
        write_ppm(out_name, fig_img);

        for (i = 0; i < width; i++) {	
            for (j = 0; j < height; j++) {
                dM[i][j] = M[i][j];
            }
        }
        strcpy(out_name, "WetArea");
        strcat(out_name, count_name);
        strcat(out_name, ".ppm");
        trans_Vector_img(fig_img, dM, width, height);
        write_ppm(out_name, fig_img);
    }
    

}


void SimulateCapillaryFlow_PaintTest(int argc, char *argv[])
{
    int i,j;
    double t;
    double max=0,var_t;
    int width=128, height=128;
    Point rectangleP[2] = {
        {0,128},
        {0,128}
    };
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_dally(p, width, height, 0);
    double** h = perlin_img(width, height, 0.1, 4);
    double** grad_hx = create_dally(width+1, height); 
    double** grad_hy = create_dally(width, height+1); 
    calcu_grad_h(h, grad_hx, grad_hy, width, height);
    PPM* fig_img = create_ppm(width, height, 255);
    PPM* Canvas_img = create_ppm(width, height, 255);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);
    double** dR = create_dally(width, height);
    double** dG = create_dally(width, height);
    double** dB = create_dally(width, height);
    format_dally(dR, width, height, 0);
    format_dally(dG, width, height, 0);
    format_dally(dB, width, height, 0);

	int w = (int)( ceil(3.0*opt_K/6+0.5)*2-1 ); //とりあえず動く計算
	int c=(w-1)/2;
	double** filter = create_dally(w, w);	
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, opt_K/6.0);
        }
    }

    char count_name[8];
    char out_name[32];
    double** dM = create_dally(width, height);
    int paint_count=0;

    double** s = create_dally(width, height);
    format_dally(s, width, height, 0);

    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            u[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            v[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>30 && j>30 && i<90 && j<90){
                M[i][j] = 1;
                p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
                gR[i][j] = 1;
            }  else M[i][j]=0;
            // if(i!=0 && j!=0) (128-abs(64-i)-abs(64-j))/128.0; 
        }
    }


    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            u[i][j] = u[i][j] - grad_hx[i][j];
            max = fmax( max, fabs(u[i][j]) );
            v[i][j] = v[i][j] - grad_hy[i][j];
            max = fmax( max, fabs(v[i][j]) );
        }
    }

    var_t = fmin(1/max, opt_SoakTimeStep);  // maxは１以下なのでおかしい・・・おかしくない？
    pd("var_t", var_t);

    for ( t = 0; t < opt_SoakTime; t=t+var_t)
    {   
        UpdateVelocities(M, u, v, p, var_t, width, height);	
        RelaxDivergence(M, u, v, p, var_t, width, height, rectangleP);
        FlowOutward(M, p, c, filter, var_t, width, height, rectangleP);
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height, rectangleP);
        TransferPigment(M, h, gR, gG, gB, dR, dG, dB, var_t, width, height);
        if(opt_USE_Backrun) SimulateCapillaryFlow(M, p, h, s, var_t, width, height);

        for (i = 0; i < width; i++) {	
            for (j = 0; j < height; j++) {
                dM[i][j] = M[i][j];
                Canvas_img->dataR[i][j] = (1 - dR[i][j]) * 255;    //CMY[0,1]->RGB[0,255]
                Canvas_img->dataG[i][j] = (1 - dG[i][j]) * 255;
                Canvas_img->dataB[i][j] = (1 - dB[i][j]) * 255;
            }
        }
        paint_count++;
        snprintf(count_name, 16, "%02d", paint_count);
        strcpy(out_name, "WetArea");
        strcat(out_name, count_name);
        strcat(out_name, ".ppm");
        trans_Vector_img(fig_img, dM, width, height);
        write_ppm(out_name, fig_img);
        strcpy(out_name, "Can");
        strcat(out_name, count_name);
        strcat(out_name, ".ppm");
        write_ppm(out_name, Canvas_img);
    }

}