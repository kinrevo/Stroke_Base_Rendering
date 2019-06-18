#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "sbr.h"

#define opt_K 20
#define opt_eta 0.25

double uf(double** u, double x, double y);
double vf(double** v, double x, double y);

double** perlin_img(int width, int height, double freq, int depth);
void UpdateVelocities(int** M,  double** u, double** v, double** p, double var_t, int width, int height);
void RelaxDivergence(int** M, double** u, double** v, double** p, double var_t, int widht, int height);
void write_Vector_img(PPM* img, double** u, int width, int height);
void FlowOutward(int** M, double** p, double var_t, int width, int height);
void TransferPigment(int** M, double** h, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, double var_t, int width, int height);
void MovePigment(int** M,  double** u, double** v, double** gR, double** gG, double** gB, double var_t, int width, int height);

void UpdateVelocity_Test(int argc, char *argv[]);
void RelaxDivergence_Test(int argc, char *argv[]);
void FlowOutward_Test(int argc, char *argv[]);
void TransferPigment_Test(int argc, char *argv[]);

#define p(s,a) printf("%s:%d\n",s,a)


int main(int argc, char *argv[]){
    TransferPigment_Test(argc, argv);
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



void ufvf_Test(){
    double i,j;
    int width=8, height=8;
    double** u = create_dally(width-1, height);
    double** v = create_dally(width, height-1);
    
	// uf
    for (i = 0; i < width-1; i++) {  //doubleによる添字の実験
        for (j = 0; j < height; j++) {
            u[(int)i][(int)j] = i+j;
        }
    }
    for (j = 0; j < height; j++) {  //uf(u, i+0.5, j)
        for (i = 0; i < width-1; i++){
            printf("%5.1f  ", uf(u, i+0.5, j));
        }
        pn;
    }
    for (j = 0; j < height-1; j++) {  //uf(u, i+0.5, j+0.5)
        for (i = 1; i < width-2; i++){
            printf("%5.2f  ", uf(u, i+0.5, j+0.5));
        }
        pn;
    }
    pn;

	// vf
    for (i = 0; i < width; i++) {
        for (j = 0; j < height; j++) {
            v[(int)i][(int)j] = i+j;
        }
    }
    for (j = 0; j < height-1; j++) {
        for (i = 0; i < width; i++){
            printf("%5.1f  ", vf(v, i, j+0.5));
        }
        pn;
    }

    for (j = 1; j < height-1; j++) {
        for (i = 0; i < width; i++){
            printf("%5.2f  ", vf(v, i, j));
        }
        pn;
    }
}


void UpdateVelocity_Test(int argc, char *argv[]){
    int i,j;
    int width=128, height=128;
    double t;
    double** u = create_dally(width-1, height);
    double** v = create_dally(width, height-1);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
   	// PGM* u_img = create_pgm(width, height, 256);
    PPM* u_img = create_ppm(width, height, 256);
    char filename[64]={};
    char tmp_name[16]={};

    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            if(i!=width-1) u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(j!=height-1) v[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            M[i][j] = 1;
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
        }
    }

    strcpy(filename, "0_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    write_Vector_img(u_img, u, width-1, height);
    write_ppm(filename, u_img);

	// UpdateVelocitieの変化の推移を出力
    double var_t = 0.01;
    for (t = 0; t < 10; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);
        if((int)(t/var_t) % 100 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);
        }
    }
}


void RelaxDivergence_Test(int argc, char *argv[]){
    int i,j;
    int width=128, height=128;
    double t;
    double** u = create_dally(width-1, height);
    double** v = create_dally(width, height-1);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    PPM* u_img = create_ppm(width-1, height, 256);
    PPM* p_img = create_ppm(width, height, 256);
    char filename[64]={};
    char tmp_name[16]={};


    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            if(i!=width-1) u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(j!=height-1) v[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            M[i][j] = 1;
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    write_Vector_img(u_img, u, width-1, height);
    write_ppm(filename, u_img);
    strcpy(filename, "0p_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    write_Vector_img(p_img, p, width, height);
    write_ppm(filename, p_img);

    double var_t = 0.01;
    for (t = 0; t < 10; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        if((int)(t/var_t) % 100 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dUV", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);
        }
        RelaxDivergence(M, u, v, p, var_t, width, height);		// RelaxDivergenceの変化の推移を出力
        if((int)(t/var_t) % 100 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dRD", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);

            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "p%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(p_img, p, width, height);
            write_ppm(filename, p_img);
        }
    }
}


void FlowOutward_Test(int argc, char *argv[]){
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
    
    gauss_M = gaussian_filter_d(dM, K/6.0, width, height);

    write_Vector_img(fig_img, gauss_M, width, height);
    write_ppm("gaussM.ppm", fig_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            revers_gauss_M[i][j] = 1-gauss_M[i][j];
        }
    }

    write_Vector_img(fig_img, revers_gauss_M, width, height);
    write_ppm("revers_gauss_M.ppm", fig_img);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            remove_p[i][j] = eta*revers_gauss_M[i][j]*M[i][j];
        }
    }

    write_Vector_img(fig_img, remove_p, width, height);
    write_ppm("remove_p.ppm", fig_img);
    write_Vector_img(fig_img, p, width, height);
    write_ppm("orig_p.ppm", fig_img);

    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            p[i][j] = p[i][j] - remove_p[i][j]/eta;
        }
    }

    write_Vector_img(fig_img, p, width, height);
    write_ppm("res_p.ppm", fig_img);
}


void MovePigment_Test(int argc, char *argv[]){
    int i,j;
    int width=128, height=128;
    double t;
    double** u = create_dally(width-1, height);
    double** v = create_dally(width, height-1);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, 0, width, height);
    format_dally(gG, 0, width, height);
    format_dally(gB, 0, width, height);
    PPM* u_img = create_ppm(width-1, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    char filename[64]={};
    char tmp_name[16]={};


    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            if(i!=width-1) u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(j!=height-1) v[i][j] = 0;//(128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>0 && j>0 && i<width/1.5 && j<height/1.5){
                M[i][j] = 1;
            }  else M[i][j]=0;
            p[i][j] = (128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            if(i!=0 && j!=0) gR[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0; 
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として表現
    strcat(filename, argv[1]);
    write_Vector_img(u_img, u, width-1, height);
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

    double var_t = 0.01;
    for (t = 0; t < 10; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        RelaxDivergence(M, u, v, p, var_t, width, height);		// RelaxDivergenceの変化の推移を出力
        FlowOutward(M, p, var_t, width, height);
        if((int)(t/var_t) % 100 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dUV", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);
        }
 
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height);	// MovePigmentの変化の推移を出力
        if((int)(t/var_t) % 100 == 0 ){
            printf("%03d\n", (int)(t/var_t));
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
}



void TransferPigment_Test(int argc, char *argv[])
{
    int i,j;
    int width=128, height=128;
    double t;
    double** u = create_dally(width-1, height);
    double** v = create_dally(width, height-1);
    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_dally(p, 0, width, height);
    double** h = create_dally(width, height);
    h = perlin_img(width, height, 0.1, 4);
    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, 0, width, height);
    format_dally(gG, 0, width, height);
    format_dally(gB, 0, width, height);
    double** dR = create_dally(width, height);
    double** dG = create_dally(width, height);
    double** dB = create_dally(width, height);
    format_dally(dR, 0, width, height);
    format_dally(dG, 0, width, height);
    format_dally(dB, 0, width, height);
    PPM* u_img = create_ppm(width-1, height, 255);
    PPM* p_img = create_ppm(width, height, 255);
    PPM* paint_img = create_ppm(width, height, 255);
    char filename[64]={};
    char tmp_name[16]={};


    for (i = 0; i < width; i++) {	//u,v,pの初期化
        for (j = 0; j < height; j++) {
            if(i!=width-1) u[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(j!=height-1) v[i][j] = (128-abs(64-i)-abs(64-j))/128.0;//(128-i-j)/128.0;
            if(i>30 && j>30 && i<90 && j<90){
                M[i][j] = 1;
                p[i][j] = 0.5;//(128-abs(64-i)-abs(64-j))/128.0; //(256-abs(256-i-j))/256;
            }  else M[i][j]=0;
            if(i!=0 && j!=0) gR[i][j] = 0.7;//(128-abs(64-i)-abs(64-j))/128.0; 
        }
    }

    strcpy(filename, "0u_");		//uの初期速度を画像として出力
    strcat(filename, argv[1]);
    write_Vector_img(u_img, u, width-1, height);
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
    for (t = 0; t < 1000; t=t+var_t) {
        UpdateVelocities(M, u, v, p, var_t, width, height);		// UpdateVelocitieの変化の推移を出力
        RelaxDivergence(M, u, v, p, var_t, width, height);		// RelaxDivergenceの変化の推移を出力
        FlowOutward(M, p, var_t, width, height);
        if((int)(t/var_t) % 10 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "%03dUV", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(u_img, u, width-1, height);
            write_ppm(filename, u_img);
        }
        if((int)(t/var_t) % 10 == 0 ){
            printf("%03d\n", (int)(t/var_t));
            snprintf(tmp_name, 16, "p%03d", (int)(t/var_t));
            strcpy(filename, tmp_name);
            strcat(filename, argv[1]);
            write_Vector_img(p_img, p, width, height);
            write_ppm(filename, p_img);
        }
 
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height);	// MovePigmentの変化の推移を出力
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

        TransferPigment(M, h, gR, gG, gB, dR, dG, dB, var_t, width, height); // TransferPigmentの変化の推移を出力
        if((int)(t/var_t) % 10 == 0 ){
            printf("%03d\n", (int)(t/var_t));
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
}
