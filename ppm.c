#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"ppm.h"


#define p(s,a) printf("%s:%d\n",s,a)
#define pd(s,a) printf("%s:%f\n",s,a)
#define pn printf("\n")
#define t1_8pi 0.4142
#define t3_8pi 2.4142
#define maxValue 0.35
#define minValue 0.25
#define minWindow 0.25
#define simga0 1.6
#define BORDER 12.1



PPM *read_ppm(char *filename)
{
	int i,j;
	FILE *fp;
	PPM *ppm = (PPM *)malloc(sizeof(PPM));
	
	//コマンドラインの名前からファイルをオープン
	if((fp = fopen(filename, "r")) == NULL){
		fprintf(stderr, "Error");
		return NULL;
	}
	fscanf(fp, "%s\n", ppm->descriptor);
	
	
	//最初の行の文字列からファイルがPPMか確認する
	if(strncmp(ppm->descriptor, "P3", 2)) {
		printf("file is not P3\n");
		return NULL;
	}
	
	printf("%s\n",ppm->descriptor);
	
	//コメント行があればスキップする
	if(fgetc(fp)=='#') {
		while(1){
			if(fgetc(fp)=='\n') break;
		}
	} else { fseek(fp,-1,SEEK_CUR); }
	

	
	//ヘッダ情報を取得
	fscanf(fp, "%d %d %d\n", &ppm->width, &ppm->height, &ppm->bright);
	printf("w:%d h:%d b:%d\n",ppm->width, ppm->height, ppm->bright);

	
	int **bufR = (int**)malloc(sizeof(int*)*(ppm->height));
	for(i=0; i<ppm->height; i++) { bufR[i]=(int*)malloc(sizeof(int)*(ppm->width));}
	
	int **bufG = (int**)malloc(sizeof(int*)*(ppm->height));
	for(i=0; i<ppm->height; i++) { bufG[i]=(int*)malloc(sizeof(int)*(ppm->width));}
	
	int **bufB = (int**)malloc(sizeof(int*)*(ppm->height));
	for(i=0; i<ppm->height; i++) { bufB[i]=(int*)malloc(sizeof(int)*(ppm->width));}
	
	for(i=0; i<ppm->height; i++) {
		for(j=0; j<ppm->width; j++) { 
				if(fscanf(fp, "%d ", &bufR[i][j]) == EOF){ printf("ER");}
				if(fscanf(fp, "%d ", &bufG[i][j]) == EOF){ printf("ER");}
				if(fscanf(fp, "%d ", &bufB[i][j]) == EOF){ printf("ER");}
		}
	}	
	
	ppm->dataR = bufR;
	ppm->dataG = bufG;
	ppm->dataB = bufB;

	fclose(fp);


	return ppm;
}



//与えられたイメージをファイルに書き込む関数
int write_ppm(char* filename, PPM *ppm) 
{
	int i,j;
	FILE* fp;

	
	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "W_Open_Error");
		return -1;
	}
	
	//ヘッダ情報を書き込む
	fprintf(fp, "%s\n%d %d\n%d\n", ppm->descriptor, ppm->width,  ppm->height, ppm->bright);

		//ピクセルを一つづつ書き込んでいく
		for(i=0; i<ppm->height; i++) {
			p("i",i);
			for(j=0; j<ppm->width; j++) { 
				if(fprintf(fp, "%d ", ppm->dataR[i][j]) < 0){ printf("ER");}
				if(fprintf(fp, "%d ", ppm->dataG[i][j]) < 0){ printf("ER");}
				if(fprintf(fp, "%d ", ppm->dataB[i][j]) < 0){ printf("ER");}
			}
			fprintf(fp, "\n");
		}
	
	fclose(fp);
	
	return 0;
}

//配列を動的に確保する関数
int **create_ally(int width, int height) {
	int i;
	int **buf = (int**)malloc(sizeof(int*)*(height));
	for(i=0; i<height; i++) { buf[i]=(int*)malloc(sizeof(int)*(width)); }
	return buf;
}



//double配列を動的に確保する関数
double **create_dally(int width, int height) {
	int i;
	double **buf = (double**)malloc(sizeof(double*)*(height));
	for(i=0; i<height; i++) { buf[i]=(double*)malloc(sizeof(double)*(width)); }
	return buf;
}


//領域を開放する関数
void Free_dally(double **ally, int w) {
	int i;
	for(i=0; i < w; i++) { free(ally[i]);}
	free(ally);;
}


//領域を開放する関数
void Freeimg(PPM *img) {
	int i;
	for(i=0; i < img->height; i++) { free(img->dataR[i]);}
	for(i=0; i < img->height; i++) { free(img->dataG[i]);}
	for(i=0; i < img->height; i++) { free(img->dataB[i]);}
	free(img->dataR);
	free(img->dataG);
	free(img->dataB);
	free(img);
}



//画像イメージを複製する関数
PPM *copy_ppm(PPM *ppm){
	PPM *img = (PPM *)malloc(sizeof(PPM));
	memcpy(img->descriptor, ppm->descriptor, 3);
	img->height = ppm->height;
	img->width = ppm->width;
	img->bright = ppm->bright;
	img->dataR = create_ally(ppm->width, ppm->height);
	img->dataG = create_ally(ppm->width, ppm->height);
	img->dataB = create_ally(ppm->width, ppm->height);
	
	return img;
}





//ガウス関数
double gause_func(int x, int y, double sigma){
	double f;
	f = 1/(2*PI*sigma*sigma)*exp(-(x*x+y*y)/(2*sigma*sigma));
	return f;
}





//バイラテラルフィルタ
PPM *bilateral_filter(PPM *ppm, double sigma)
{
	int i,j,k,l;
	double sum=0;
	//新しい画像データをコピーして作っておく
	PPM *nimg = copy_ppm(ppm);
	
	int w = (int)( ceil(3.0*sigma+0.5)*2-1 );
	int c=(w-1)/2;
	
	double filter[w][w];
	double fil_sum;
	double brightdiff;  //バイラテル追加　：カーネル中心と注目画素との輝度差
	double bi_fil;
	double sigma_2 = 2*sigma*sigma;
	
	
    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, sigma);
        }
    }
	
	for(i=0; i<ppm->height; i++) {
		for(j=0; j<ppm->width; j++) {
			//注目している周囲一ピクセルの値を合計し平均する
			//R
			for(k=-c; k<=c; k++) { 
				for(l=-c; l<=c; l++)  {
					//フィルタの端ピクセルがない場合分母にも分子にも加算しない
					if( (i+k)<=0 || (i+k)>=ppm->height-1 || (j+l)<=0 || (j+l)>=ppm->width-1) {}
					else {
						brightdiff = ppm->dataR[i][j] - ppm->dataR[i+k][j+l];  //bi追加
						//bai追加：輝度差で重み付けしたガウスフィルタ
						bi_fil = filter[k+c][l+c] * exp( -(brightdiff*brightdiff)/sigma_2 );
						sum += ppm->dataR[i+k][j+l] * bi_fil;
						fil_sum += bi_fil;
					}
				}
			}
			nimg->dataR[i][j] = sum/fil_sum+0.5;
			sum = 0;
			fil_sum = 0;
			
			//G
			for(k=-c; k<=c; k++) { 
				for(l=-c; l<=c; l++)  {
					//フィルタの端ピクセルがない場合分母にも分子にも加算しない
					if( (i+k)<=0 || (i+k)>=ppm->height-1 || (j+l)<=0 || (j+l)>=ppm->width-1) {}
					else {
						brightdiff = ppm->dataG[i][j] - ppm->dataG[i+k][j+l];  //bi追加
						//bai追加：輝度差で重み付けしたガウスフィルタ
						bi_fil = filter[k+c][l+c] * exp( -(brightdiff*brightdiff)/sigma_2 );
						sum += ppm->dataG[i+k][j+l] * bi_fil;
						fil_sum += bi_fil;
					}
				}
			}
			nimg->dataG[i][j] = sum/fil_sum+0.5;
			sum = 0;
			fil_sum = 0;
			
			//B
			for(k=-c; k<=c; k++) { 
				for(l=-c; l<=c; l++)  {
					//フィルタの端ピクセルがない場合分母にも分子にも加算しない
					if( (i+k)<=0 || (i+k)>=ppm->height-1 || (j+l)<=0 || (j+l)>=ppm->width-1) {}
					else {
						brightdiff = ppm->dataB[i][j] - ppm->dataB[i+k][j+l];  //bi追加
						//bai追加：輝度差で重み付けしたガウスフィルタ
						bi_fil = filter[k+c][l+c] * exp( -(brightdiff*brightdiff)/sigma_2 );
						sum += ppm->dataB[i+k][j+l] * bi_fil;
						fil_sum += bi_fil;
					}
				}
			}
			nimg->dataB[i][j] = sum/fil_sum+0.5;
			sum = 0;
			fil_sum = 0;
			
		}
	}
	
	
	return nimg;
}


/*
//与えられたイメージを絵画調にする
PPM *Illustification(PPM *in, int n) {
	
	int i,j,k;
	//新しい画像データをコピーして作っておく
	PPM *nimg = copy_ppm(in);
	
	for(k=0; k<n; k++) {
		for(i=0; i<in->height; i++) {
			for(j=0; j<in->width; j++) {
				nimg->data[j][in->height-1-i] = in->data[i][j];
			}
		}
	}
	
	

	return nimg;
}
*/

