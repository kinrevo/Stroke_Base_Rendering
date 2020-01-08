#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <png.h>
#include <jpeglib.h>
#include "ImageIO/image.h"
#include <stdlib.h>
#include "sbr.h"
// #include "other_func.h"



/////////////////////////////////////////////
//			配列メモリを処理する関数群
/////////////////////////////////////////////


//配列を動的に確保する関数
int **create_ally(int width, int height) {
	int i;
	int **buf = (int**)malloc(sizeof(int*)*(width));
	for(i=0; i<width; i++) { buf[i]=(int*)malloc(sizeof(int)*(height)); }
	return buf;
}


//double配列を動的に確保する関数
double **create_dally(int width, int height) {
	int i;
	double **buf = (double**)malloc(sizeof(double*)*(width));
	for(i=0; i<width; i++) { buf[i]=(double*)malloc(sizeof(double)*(height)); }
	return buf;
}


//int配列をコピーする関数
void copy_ally(int **ally, int **ally2, int w, int h) {
	int i, j;
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally2[i][j] = ally[i][j];
		}
	}
}


//double配列をコピーする関数
void copy_dally(double **ally, double **ally2, int w, int h) {
	int i, j;
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally2[i][j] = ally[i][j];
		}
	}
}


//double配列を限定的にコピーする関数
void rect_copy_dally(double **ally, double **ally2, int x_min, int x_max, int y_min, int y_max) {
	int i, j;
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
	for(i=x_min; i<=x_max; i++) {
		for(j=y_min; j<=y_max; j++) {
			ally2[i][j] = ally[i][j];
		}
	}
}


//配列を端末表示する関数
void display_ally(int *ally, int num) {
	int i;
	printf("(%d",ally[0]);
	for(i=1; i<num; i++) {
		printf(" %d", ally[i]);
	}
	printf(")\n");
}


//配列を初期化する関数
void format_ally(int **ally, int w, int h, int bright) {
	int i, j;
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally[i][j] = bright;
		}
	}
}

//配列を初期化する関数
void format_dally(double **ally, int w, int h, double bright) {
	int i, j;
	#ifdef _OPENMP
		#pragma omp parallel for private(i,j)
    #endif
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally[i][j] = bright;
		}
	}
}


//double領域を開放する関数
void Free_dally(double **ally, int w) {
	int i;
	#ifdef _OPENMP
		#pragma omp parallel for private(i)
    #endif
	for(i=0; i < w; i++) { free(ally[i]);}
	free(ally);
}


// １次元int配列を２次元に変形
int** ReshapeInt_1to2(int* x, int width, int height)
{
    int i,j;
    int** y = create_ally(width, height);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            y[i][j] = x[j + i*height];
        }
    }

    return y;
}


// 配列中の最大値を探索しPointを返す
Point search_max_Point(int **ally, int w, int h) {
	int y, x, max_value=-99999999;
	Point max;
	max.x=0; max.y=0;
	for(y=0; y<h; y++) {
		for(x=0; x<w; x++) {
			if(ally[x][y]>max_value) {
				max.x=x;
				max.y=y;
				max_value=ally[x][y];
			}
		}
	}

	return max;
}


//　ストローク構造体メモリを確保
Stroke *create_Stroke(int pnum){
	Stroke* sp = (Stroke *)malloc(sizeof(Stroke));
	sp->pnum = pnum;
	sp->p = (Point*)malloc(sizeof(Point)*(sp->pnum));

	return sp;
}


// Stroke構造体のアドレス変数の二重配列を生成
Stroke ***create_Stroke_ally(int width, int height, int max_stroke) {
	int i,j;

	Stroke*** stroke_map = (Stroke***)malloc(sizeof(Stroke**)*(width));
	for(i=0; i<width; i++){
		stroke_map[i] = (Stroke**)malloc(sizeof(Stroke*)*(height));
	}
	for(i=0; i<width; i++){
		for(j=0; j<height; j++){
			stroke_map[i][j]=create_Stroke(max_stroke);
		}
	}

	return stroke_map;
}


/////////////////////////////////////////////
//			PGM,PPM構造体を処理する関数群
/////////////////////////////////////////////


//PPM画像ファイルを読み込みPPM構造体を返す
PPM *read_ppm(char *filename)
{
	int i,j;
	FILE *fp;
	PPM *ppm = (PPM *)malloc(sizeof(PPM));

	//コマンドラインの名前からファイルをオープン
	if((fp = fopen(filename, "r")) == NULL){
		fprintf(stderr, "Error");
		free(ppm);
		return NULL;
	}
	fscanf(fp, "%s\n", ppm->descriptor);


	//最初の行の文字列からファイルがPPMか確認する
	if(strncmp(ppm->descriptor, "P3", 2)) {
		printf("file is not P3\n");
		free(ppm);	fclose(fp);
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


	int **bufR = (int**)malloc(sizeof(int*)*(ppm->width));
	for(i=0; i<ppm->width; i++) { bufR[i]=(int*)malloc(sizeof(int)*(ppm->height));}

	int **bufG = (int**)malloc(sizeof(int*)*(ppm->width));
	for(i=0; i<ppm->width; i++) { bufG[i]=(int*)malloc(sizeof(int)*(ppm->height));}

	int **bufB = (int**)malloc(sizeof(int*)*(ppm->width));
	for(i=0; i<ppm->width; i++) { bufB[i]=(int*)malloc(sizeof(int)*(ppm->height));}

	for(i=0; i<ppm->height; i++) {
		for(j=0; j<ppm->width; j++) {
				if(fscanf(fp, "%d ", &bufR[j][i]) == EOF){ printf("ER");}
				if(fscanf(fp, "%d ", &bufG[j][i]) == EOF){ printf("ER");}
				if(fscanf(fp, "%d ", &bufB[j][i]) == EOF){ printf("ER");}
		}
	}

	ppm->dataR = bufR;
	ppm->dataG = bufG;
	ppm->dataB = bufB;

	fclose(fp);


	return ppm;
}


//PGM画像ファイルを読み込みPGM構造体を返す
PGM *read_pgm(char *filename)
{
	int i,j;
	FILE *fp;
	PGM *pgm = (PGM *)malloc(sizeof(PGM));

	//コマンドラインの名前からファイルをオープン
	if((fp = fopen(filename, "r")) == NULL){
		fprintf(stderr, "Error");
		free(pgm);
		return NULL;
	}
	fscanf(fp, "%s\n", pgm->descriptor);


	//最初の行の文字列からファイルがPGMか確認する
	if(strncmp(pgm->descriptor, "P2", 2)) {
		free(pgm);  fclose(fp);
		printf("not P2\n");
		return NULL;
	}
	printf("%s\n",pgm->descriptor);

	//コメント行があればスキップする
	if(fgetc(fp)=='#') {
		while(1){
			if(fgetc(fp)=='\n') break;
		}
	} else { fseek(fp,-1,SEEK_CUR); }



	//ヘッダ情報を取得
	fscanf(fp, "%d %d %d\n", &pgm->width, &pgm->height, &pgm->bright);
	printf("w:%d h:%d b:%d\n",pgm->width, pgm->height, pgm->bright);


	int **buf = (int**)malloc(sizeof(int*)*(pgm->width));
	for(i=0; i<pgm->width; i++) { buf[i]=(int*)malloc(sizeof(int)*(pgm->height));}


	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
				if(fscanf(fp, "%d ", &buf[j][i]) == EOF){ printf("ER");}
		}
	}

	pgm->data = buf;



	fclose(fp);

	return pgm;
}


//与えられたPGMイメージをファイルに書き込む関数
int write_pgm(char* filename, PGM *pgm)
{
	int i,j;
	FILE* fp;

	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "W_Open_Error(write_pgm)\n");
		return -1;
	}

	//ヘッダ情報を書き込む
	fprintf(fp, "%s\n%d %d\n%d\n", pgm->descriptor, pgm->width,  pgm->height, pgm->bright);

	if(!(strncmp(pgm->descriptor, "P2", 2))) {
		//ピクセルを一つづつ書き込んでいく
		for(i=0; i<pgm->height; i++) {
			for(j=0; j<pgm->width; j++) {
				if(fprintf(fp, "%d ", pgm->data[j][i]) < 0){ printf("ER");}
			}
			fprintf(fp, "\n");
		}
	}
	else{fprintf(stderr, "Write_PGM_data_Error");}
	fclose(fp);

	return 0;
}


//与えられたPPMイメージをファイルに書き込む関数
int write_ppm(char* filename, PPM *pgm)
{
	int i,j;
	FILE* fp;

	if((fp = fopen(filename, "w")) == NULL){
		fprintf(stderr, "W_Open_Error(write_ppm)\n");
		return -1;
	}

	//ヘッダ情報を書き込む
	fprintf(fp, "%s\n%d %d\n%d\n", pgm->descriptor, pgm->width,  pgm->height, pgm->bright);

	if(!(strncmp(pgm->descriptor, "P3", 2))) {
		//ピクセルを一つづつ書き込んでいく
		for(i=0; i<pgm->height; i++) {
			for(j=0; j<pgm->width; j++) {
				if(fprintf(fp, "%d ", pgm->dataR[j][i]) < 0){ printf("ER");}
				if(fprintf(fp, "%d ", pgm->dataG[j][i]) < 0){ printf("ER");}
				if(fprintf(fp, "%d ", pgm->dataB[j][i]) < 0){ printf("ER");}
			}
			fprintf(fp, "\n");
		}
	}
	else{fprintf(stderr, "Write_PPM_data_Error");}
	fclose(fp);

	return 0;
}


//PGM領域を開放する関数
void FreePGM(PGM *img) {
	int i;
	for(i=0; i < img->width; i++) { free(img->data[i]);}
	free(img->data);
	free(img);
}


//PPM領域を開放する関数
void FreePPM(PPM *img) {
	int i;
	for(i=0; i < img->width; i++) {
		free(img->dataR[i]);
		free(img->dataG[i]);
		free(img->dataB[i]);
	}
	free(img->dataR);
	free(img->dataG);
	free(img->dataB);
	free(img);
}


//image形式データ（png用）をPPM形式データ（自定義）に変換
PPM *image_to_PPM(image_t *img){
	int i,j;

	PPM *ppm = (PPM *)malloc(sizeof(PPM));
	ppm->dataR = create_ally(img->width, img->height);
	ppm->dataG = create_ally(img->width, img->height);
	ppm->dataB = create_ally(img->width, img->height);

	memcpy(ppm->descriptor, "P3", 3);
	ppm->height = img->height;
	ppm->width = img->width;
	ppm->bright = 255;

	for(i=0; i<img->width; i++) {
		for(j=0; j<img->height; j++) {
			ppm->dataR[i][j] = img->map[j][i].c.r;
			ppm->dataG[i][j] = img->map[j][i].c.g;
			ppm->dataB[i][j] = img->map[j][i].c.b;
		}
	}

	return ppm;
}


//PPM形式データ（自定義）をimage形式データ（汎用）に変換
image_t *PPM_to_image(PPM *ppm){
	int i,j;
	//image_t *nimg = clone_image(img);
	image_t *nimg = allocate_image(ppm->width, ppm->height, COLOR_TYPE_RGB);
	if (nimg == NULL) {
	  return NULL;
	}


	for(i=0; i<ppm->width; i++) {
		for(j=0; j<ppm->height; j++) {
			nimg->map[j][i].c.r = ppm->dataR[i][j];
			nimg->map[j][i].c.g = ppm->dataG[i][j];
			nimg->map[j][i].c.b = ppm->dataB[i][j];
		}
	}

	return nimg;
}


// PPMカラー画像データをPGMグレイ画像に変換して返す
PGM *color_gray_conversion(PPM* in){
	int i,j;
	PGM *gray = (PGM *)malloc(sizeof(PGM));
	memcpy(gray->descriptor, "P2", 3);
	gray->height = in->height;
	gray->width = in->width;
	gray->bright = in->bright;
	gray->data = create_ally(in->width, in->height);

	for(i=0; i<in->width; i++) {
		for(j=0; j<in->height; j++) {
			gray->data[i][j] = 0.299*in->dataR[i][j] + 0.587*in->dataG[i][j] + 0.114*in->dataB[i][j] + 0.5;
		}
	}

	return gray;
}


// PPMカラーデータを三つのPGMグレイデータに分割
void devide_ppm(PPM *ppm, PGM* ppmR, PGM* ppmG, PGM* ppmB){
	int i,j;
	for(i=0; i<ppm->width; i++) {
		for(j=0; j<ppm->height; j++) {
			ppmR->data[i][j] = ppm->dataR[i][j];
			ppmG->data[i][j] = ppm->dataG[i][j];
			ppmB->data[i][j] = ppm->dataB[i][j];
		}
	}
}


//　PGMデータを画素値を除いて複製する関数
PGM *copy_pgm(PGM *pgm){
	PGM *img = (PGM *)malloc(sizeof(PGM));
	memcpy(img->descriptor, pgm->descriptor, 3);
	img->height = pgm->height;
	img->width = pgm->width;
	img->bright = pgm->bright;
	img->data = create_ally(pgm->width, pgm->height);

	return img;
}


//　PGM画像領域を確保（最大輝度で初期化）
PGM *create_pgm(int width, int height, int bright){
	PGM *img = (PGM *)malloc(sizeof(PGM));
	memcpy(img->descriptor, "P2", 3);
	img->height = height;
	img->width = width;
	img->bright = bright;
	img->data = create_ally(width, height);
	format_ally(img->data, img->width, img->height, bright);

	return img;
}


//　PPMデータを画素値を除いて複製する関数
PPM *copy_ppm(PPM *ppm, int bright){
	PPM *img = (PPM *)malloc(sizeof(PPM));
	memcpy(img->descriptor, ppm->descriptor, 3);
	img->height = ppm->height;
	img->width = ppm->width;
	img->bright = ppm->bright;
	img->dataR = create_ally(ppm->width, ppm->height);
	img->dataG = create_ally(ppm->width, ppm->height);
	img->dataB = create_ally(ppm->width, ppm->height);
	format_ally(img->dataR, ppm->width, ppm->height, bright);
	format_ally(img->dataG, ppm->width, ppm->height, bright);
	format_ally(img->dataB, ppm->width, ppm->height, bright);

	return img;
}


// PPM画像領域を確保（最大輝度で初期化）
PPM *create_ppm(int width, int height, int bright){
	PPM *img = (PPM *)malloc(sizeof(PPM));
	memcpy(img->descriptor, "P3", 3);
	img->height = height;
	img->width = width;
	img->bright = bright;
	img->dataR = create_ally(width, height);
	img->dataG = create_ally(width, height);
	img->dataB = create_ally(width, height);
	format_ally(img->dataR, img->width, img->height, bright);
	format_ally(img->dataG, img->width, img->height, bright);
	format_ally(img->dataB, img->width, img->height, bright);

	return img;
}


// PPMデータを任意の画像形式でファイルに書き込む
void MY__write_img(PPM* out_ppm, char dir_path[], char in_filename[], char last_chara[], char ext[])
{
	char out_filename[128];
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, last_chara);

	image_t *out_img;
	out_img = PPM_to_image(out_ppm);

	if(strncmp(ext, "png", 3) == 0){
		strcat(out_filename, ".png");
		if(write_png_file(out_filename, out_img)) { printf("WRITE PNG ERROR."); }
		free_image(out_img);
	}else if(strncmp(ext, "jpg", 3) == 0){
		strcat(out_filename, ".jpg");
		if(write_jpeg_file(out_filename, out_img)) { printf("WRITE JPG ERROR."); }
		free_image(out_img);
	}else if(strncmp(ext, "ppm", 3) == 0){
		strcat(out_filename, ".ppm");
		if(write_ppm(out_filename, out_ppm)) { printf("WRITE PPM ERROR."); }
	}else {
		printf("ERROR: please enter default extension !\n");
	}

	strcpy(out_filename, in_filename);
	strcat(out_filename, last_chara);
	strcat(out_filename, ".");
	strcat(out_filename, ext);
	printf("WRITE: %s\n",out_filename);
}




/////////////////////////////////////////////
//			画像処理を行う関数群
/////////////////////////////////////////////

//ガウス関数
double gause_func(int x, int y, double sigma){
	double f;
	f = 1/(2*PI*sigma*sigma)*exp(-(x*x+y*y)/(2*sigma*sigma));
	return f;
}


//ガウシアンフィルタ可変半径
PGM *gaussian_filter(PGM *pgm, double sigma)
{
	int i,j,k,l;
	double sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);

	int w = (int)( ceil(3.0*sigma+0.5)*2-1 );
	int c=(w-1)/2;

	double filter[w][w];
	double fil_sum;

    for(i=0;i<w;i++){
        for(j=0;j<w;j++){
        	filter[i][j] = gause_func(i-c, j-c, sigma);
        }
    }

	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			//注目している周囲cピクセルの値を合計し平均する
			for(k=-c; k<=c; k++) {
				for(l=-c; l<=c; l++)  {
					//フィルタの端ピクセルがない場合分母にも分子にも加算しない
					if( (i+k)<0 || (i+k)>=pgm->width-1 || (j+l)<0 || (j+l)>=pgm->height-1) {}
					else {
						sum += pgm->data[i+k][j+l] * filter[k+c][l+c];
						fil_sum += filter[k+c][l+c];
					}
				}
			}
			nimg->data[i][j] = sum/fil_sum+0.5;
			sum = 0;
			fil_sum = 0;
		}
	}

	return nimg;
}


//ガウシアンフィルタ可変半径(double)
double** gaussian_filter_d(double** in, int c, double** filter, int width, int height)
{
	int i,j,k,l;
	double sum=0;
	double fil_sum;
	double** new = create_dally(width, height);
	//入力データをコピー
	copy_dally(in, new, width, height);


	#ifdef _OPENMP
		#pragma omp parallel for private(i,j,sum,fil_sum)
    #endif
	for(i=0; i<width; i++) {
		for(j=0; j<height; j++) {
			sum = 0;
			fil_sum = 0;
			//注目している周囲cピクセルの値を合計し平均する
			for(k=-c; k<=c; k++) {
				for(l=-c; l<=c; l++)  {
					//フィルタの端ピクセルがない場合分母にも分子にも加算しない
					if( (i+k)<0 || (i+k)>=width-1 || (j+l)<0 || (j+l)>=height-1) {}
					else {
						sum += in[i+k][j+l] * filter[k+c][l+c];
						fil_sum += filter[k+c][l+c];
					}
				}
			}
			new[i][j] = sum/fil_sum;
		}
	}

	return new;
}


//バイラテラルフィルタ
PGM *bilateral_filter(PGM *pgm, double sigma, int loopc)
{
	int i,j,k,l,n;
	double sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);  //変換後データ
	PGM *nimg2 = copy_pgm(pgm);  //変換前データ

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

	copy_ally(pgm->data, nimg2->data, pgm->width, pgm->height);

	for(n=0; n<loopc; n++) {
		for(i=0; i<pgm->width; i++) {
			for(j=0; j<pgm->height; j++) {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-c; k<=c; k++) {
					for(l=-c; l<=c; l++)  {
						//フィルタの端ピクセルがない場合分母にも分子にも加算しない
						if( (i+k)<=0 || (i+k)>=pgm->width-1 || (j+l)<=0 || (j+l)>=pgm->height-1) {}
						else {
							brightdiff = pgm->data[i][j] - pgm->data[i+k][j+l];  //bi追加
							//bai追加：輝度差で重み付けしたガウスフィルタ
							bi_fil = filter[k+c][l+c] * exp( -(brightdiff*brightdiff)/sigma_2 );
							sum += nimg2->data[i+k][j+l] * bi_fil;
							fil_sum += bi_fil;
						}
					}
				}
				nimg->data[i][j] = sum/fil_sum+0.5;
				sum = 0;
				fil_sum = 0;
			}
		}
		//入力ループ毎に変換データを保存して更新して行く
		copy_ally(nimg->data, nimg2->data, pgm->width, pgm->height);
	}

	return nimg;
}


//内分点を返す関数
int inner_point(int p1,int p2,float m,float n){
	int r = ( (n*p1+m*p2)/(m+n) + 0.5);
	return r;
}


//キャニーエッジ検出器
PGM *cannyedge_detector(PGM *pgm, double maxValue, double minValue, int thick_min)
{
	int i,j,k,l,flag=0,sum=0,sumx=0,sumy=0,min=pgm->bright,max=0;
	//double maxValue=0.30, minValue=0.05;
	double t1_8pi=PI/8, t3_8pi=PI*3/8;
	//新しい画像データを作りガウシアンフィルタを掛けたものを入れる(nimg1)
	PGM *nimg1 = gaussian_filter(pgm,thick_min/3.0);
	printf("Canny_step1 finish\n");
	//画像データを作りsobelフィルタをかける(nimg2)
	PGM *nimg2 = copy_pgm(nimg1);

	double **theta = (double**)malloc(sizeof(double*)*(pgm->width));
	if(theta==NULL) {fprintf(stderr, "fail memory\n");}
	for(i=0; i<pgm->width; i++) { theta[i]=(double*)malloc(sizeof(double)*(pgm->height));}
	//double theta[pgm->height][pgm->height];


	printf("Canny_step2-1 finish\n");
	int filterx[3][3] = {
		{-1,0,1},
		{-2,0,2},
		{-1,0,1},
	};
	int filtery[3][3] = {
		{-1,-2,-1},
		{0,0,0},
		{1,2,1},
	};


	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			if(i==0 || i==pgm->width-1 || j==0 || j==pgm->height-1) {
				nimg2->data[i][j] = 0;
			}else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) {
					for(l=-1; l<=1; l++)  {
						sumx += nimg1->data[i+k][j+l] * filterx[l+1][k+1];
						sumy += nimg1->data[i+k][j+l] * filtery[l+1][k+1];
					}
				}

				sum = (int)sqrt((double)(sumx*sumx+sumy*sumy));
				/*
				if(sum > nimg1->bright){ nimg2->data[i][j]=nimg1->bright;
				}else if(sum<0){ nimg2->data[i][j]=0;
				}else { nimg2->data[i][j]=sum; }
				*/ nimg2->data[i][j]=sum;
				if(sumx==0) {theta[i][j] = (double)(sumy*1000);} //分母がゼロにならないようにする
				else { theta[i][j] = (double)(sumy/sumx);}  //次の処理のため角度を計算し格納しておく
				theta[i][j] = atan(theta[i][j]);		//ラジアン変換と鉛直方向に変換

				if(sum>max) {max=sum;}  //最大小値をスケーリングのため保存しておく
				if(sum<min) {min=sum;}

				sum = sumx = sumy = 0;
			}
		}
	}

	//スケーリング処理を行う
	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			nimg2->data[i][j] += (-1)*min;
			nimg2->data[i][j] = (nimg2->data[i][j]*pgm->bright)/(max-min);
			if(nimg2->data[i][j] > 255) printf("error\n");
		}
	}
	printf("max:%d min:%d\n",max,min);

	printf("Canny_step2 finish\n");


	//非極大値抑制処理を行う(nimg3)
	PGM *nimg3 = copy_pgm(nimg2);

	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			if(i==0 || i==pgm->width-1 || j==0 || j==pgm->height-1) {
				nimg3->data[i][j] = nimg2->data[i][j];
			}else {
				//注目しているピクセルのエッジ方向の鉛直方向の角度を取る　(x,y)=(k,l)
				//-1/8pi < 1/8pi
				if( (-1)*t1_8pi <= theta[i][j] && theta[i][j] <= t1_8pi ) { k=1; l=0; }
				//1/8pi < 3/8pi
				else if( t1_8pi <= theta[i][j] && theta[i][j] <= t3_8pi ) { k=1; l=1; }
				//3/8pi < 5/8pi
				else if( (-1)*t3_8pi <= theta[i][j] && theta[i][j] <= (-1)*t1_8pi ) { k=1; l=-1; }
				//5/8pi < 7/8pi
				else  {	k=0; l=1; }

				//注目しているピクセルをエッジ鉛直方向の隣接画素と比較し最大でなければ0とする
				if( nimg2->data[i][j] < nimg2->data[i-k][j-l] || nimg2->data[i][j] < nimg2->data[i+k][j+l] ) {
					nimg3->data[i][j] =0;
				}
				else { nimg3->data[i][j] = nimg2->data[i][j]; }


			}
		}
	}



	printf("Canny_step3 finish\n");
	//ヒステリシス閾処理を行う(nimg4)
	PGM *nimg4 = copy_pgm(nimg3);

	//ハイパスを超える値を最大、ローパスより低い値を最小とする
	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			if(nimg3->data[i][j] > maxValue*nimg3->bright){ nimg4->data[i][j]=nimg3->bright; }
			else if(nimg3->data[i][j] < minValue*nimg3->bright){ nimg4->data[i][j]=0; }
			else { nimg4->data[i][j] = nimg3->data[i][j]; }

		}
	}


	//2つの閾値の間の値はエッジに繋がっていればエッジとする
	int wflag=1, flag2=0;
	//エッジの接続を全て確認するまで画素の全探索を続ける
	while(wflag){
		wflag=0;
		for(i=0; i<pgm->width; i++) {
			for(j=0; j<pgm->height; j++) {
				if(i==0 || i==pgm->width-1 || j==0 || j==pgm->height-1) {
				}else if(nimg4->data[i][j]==0 || nimg4->data[i][j]==pgm->bright){
				}else {
					//注目しているピクセルの周りにエッジがあるか確認する
					for(k=-1; k<=1; k++) {
						for(l=-1; l<=1; l++)  {
							if(nimg4->data[i+k][j+l]==nimg4->bright) { flag=1; }
							if(nimg4->data[i+k][j+l]>minValue*nimg4->bright && nimg4->data[i+k][j+l]<maxValue*nimg4->bright)
								{ flag2=1; }
						}
					}

					if(flag) {
						nimg4->data[i][j]=nimg3->bright;
						wflag=1;
					}else if(flag2){}
					else {
						nimg4->data[i][j] = 0;
						wflag=1;
					}

					flag=0;
					flag2=0;

				}
			}
		}
	}

	//グレー画素を全て最低輝度に
	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			if(nimg4->data[i][j]!=pgm->bright){
				nimg4->data[i][j]=0;
			}
		}
	}

	//return nimg3;

	FreePGM(nimg1);
	FreePGM(nimg2);
	FreePGM(nimg3);
	for(i=0; i < pgm->width; i++) { free(theta[i]);}
	free(theta);

	return nimg4;
}


//sobelフィルタを適応した計算結果を返す
void sobel_calcu(PGM *pgm, double **sobel_abs, double **sobel_angle)
{
	int i,j,k,l,sumx=0,sumy=0;
	double theta;

	int filterx[3][3] = {
		{-1,0,1},
		{-2,0,2},
		{-1,0,1},
	};
	int filtery[3][3] = {
		{-1,-2,-1},
		{0,0,0},
		{1,2,1},
	};


	for(i=0; i<pgm->width; i++) {
		for(j=0; j<pgm->height; j++) {
			if(i==0 || i==pgm->width-1 || j==0 || j==pgm->height-1) {
				sobel_abs[i][j] = sobel_angle[i][j] = 0;
			}else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) {
					for(l=-1; l<=1; l++)  {
						sumx += pgm->data[i+k][j+l] * filterx[l+1][k+1];
						sumy += pgm->data[i+k][j+l] * filtery[l+1][k+1];
					}
				}
				sobel_abs[i][j] = (int)sqrt((double)(sumx*sumx+sumy*sumy));
				if(sumx==0) theta=10000000;
				else theta=(double)sumy/sumx;
				sobel_angle[i][j] = atan(theta)+PI/2;		//0~πに90度変換
				sumx = sumy = 0;
			}
		}
	}
}



//与えられた二つのイメージの差分の平均（RGBの平均）を取る（in2-in1）
double image_MSE(PPM *in1, PPM *in2) {

	int i,j;
	int error_sum=0;
	//新しい画像データをコピーして作っておく

	for(i=0; i<in1->width; i++) {
		for(j=0; j<in1->height; j++) {
			error_sum += abs(in2->dataR[i][j] - in1->dataR[i][j]);
			error_sum += abs(in2->dataG[i][j] - in1->dataG[i][j]);
			error_sum += abs(in2->dataB[i][j] - in1->dataB[i][j]);
		}
	}
	error_sum /= (in1->height * in1->width * 3);

	return error_sum;
}


//与えられた座標を目立たせるように打鋲する
void light_dot(int x, int y, PGM* in_img, int bright) {
	if(x>(in_img->width-1) || x<1 || y>(in_img->height-1) || y<1) return;
	in_img->data[x][y]=bright;
	in_img->data[x+1][y]=bright;
	in_img->data[x][y+1]=bright;
	in_img->data[x-1][y]=bright;
	in_img->data[x][y-1]=bright;
}



//配列を端末表示する関数
void display_Point_ally(Point *ally, int pnum) {
	int i;
	printf("(%.1f,%.1f)",ally[0].x,ally[0].y);
	for(i=1; i<pnum; i++) {
		printf("->(%.1f,%.1f)", ally[i].x,ally[i].y);
	}
	printf("\n");
}



// RGB->lRGB->CIExyz->CIELab  https://qiita.com/hachisukansw/items/09caabe6bec46a2a0858
Lab RGB2Lab(RGB in) {
    Lab out;
    double r = in.R / 255.0;
    double g = in.G / 255.0;
    double b = in.B / 255.0;

    double lr = r > 0.04045 ? pow(((r + 0.055) / 1.055), 2.4) : (r / 12.92);
    double lg = g > 0.04045 ? pow(((g + 0.055) / 1.055), 2.4) : (g / 12.92);
    double lb = b > 0.04045 ? pow(((b + 0.055) / 1.055), 2.4) : (b / 12.92);

    double x = (lr * 0.4124) + (lg * 0.3576) + (lb * 0.1805);
    double y = (lr * 0.2126) + (lg * 0.7152) + (lb * 0.0722);
    double z = (lr * 0.0193) + (lg * 0.1192) + (lb * 0.9505);

    x = x * 100 / 95.047;
    y = y * 100 / 100;
    z = z * 100 / 108.883;

    x = x > 0.008856 ? pow(x, 1 / 3.0) : (7.787 * x) + (4 / 29.0);
    y = y > 0.008856 ? pow(y, 1 / 3.0) : (7.787 * y) + (4 / 29.0);
    z = z > 0.008856 ? pow(z, 1 / 3.0) : (7.787 * z) + (4 / 29.0);

    out.L = (116 * y) - 16;
    out.a = 500 * (x - y);
    out.b = 200 * (y - z);

    return out;
}


// CIELab->CIExyz->lRGB->RGB  http://story-monoroch.blogspot.com/2014/09/2rgbcielab-rgbcielab-rgb-rgb-xyz-cielab.html
RGB Lab2RGB(Lab in) {
    RGB out;
    double yn = (in.L + 16.0) / 116.0;
    double xn = yn + in.a / 500.0;
    double zn = yn - in.b / 200.0;

    double fxn = xn > 6.0/29 ? pow(xn, 3) : 3*pow(6.0/29, 2.0) * (xn - 4.0/29);
    double fyn = yn > 6.0/29 ? pow(yn, 3) : 3*pow(6.0/29, 2.0) * (yn - 4.0/29);
    double fzn = zn > 6.0/29 ? pow(zn, 3) : 3*pow(6.0/29, 2.0) * (zn - 4.0/29);

    double x = fxn * 0.95047;
    double y = fyn * 1.00000;
    double z = fzn * 1.08883;

    double lr =   3.24062547732005470 * x - 1.53720797221031910 * y - 0.49862859869824794 * z;
    double lg = - 0.96893071472931970 * x + 1.87575606088524200 * y + 0.04151752384295397 * z;
    double lb =   0.05571012044551063 * x - 0.20402105059848677 * y + 1.05699594225438860 * z;

    double r = lr <= 0.0031308 ? lr * 12.92 : pow(lr, 1 / 2.4) * 1.055 - 0.055;
    double g = lg <= 0.0031308 ? lg * 12.92 : pow(lg, 1 / 2.4) * 1.055 - 0.055;
    double b = lb <= 0.0031308 ? lb * 12.92 : pow(lb, 1 / 2.4) * 1.055 - 0.055;

    out.R = r * 255;
    out.G = g * 255;
    out.B = b * 255;

    return out;
}


// 画像の代表構成色をKmeansにより抽出
RGB* Kmeans_ImageLab3D(PPM* img, int cluster_num, int maxiter, int* x_centlabel, int* num_cluster)
{
    int i,j,data;
    RGB* color_set = (RGB *)malloc(sizeof(RGB) * cluster_num);
    RGB rgb;
    Lab lab;
    int data_num =img->width*img->height;
    int dimension_num = 3;
    double* x = (double*)malloc(sizeof(double) * data_num*dimension_num);
    double* centroid = (double*)malloc(sizeof(double) * cluster_num*dimension_num);
    double* cluster_weight = (double*)malloc(sizeof(double) * cluster_num);

    // 二次元画像のab２つのデータを一次元に変換
    for(i=0; i<img->width; i++){
        for(j=0; j<img->height; j++){
            data = i*img->height + j;   //データのINDEX
            rgb.R=img->dataR[i][j];
            rgb.G=img->dataG[i][j];
            rgb.B=img->dataB[i][j];
            lab = RGB2Lab(rgb);

            x[data + data_num*0] = lab.L;
            x[data + data_num*1] = lab.a;
            x[data + data_num*2] = lab.b;
        }
    }

    // KmeansPlusによりセントロイド初期位置を決定
    KmeansPlus(x, data_num, dimension_num, cluster_num, centroid);

    Kmeans(x, data_num, dimension_num, centroid, cluster_num, x_centlabel, maxiter, num_cluster, cluster_weight);

    for (i = 0; i<cluster_num; i++) {
        lab.L = centroid[i + cluster_num*0];
        lab.a = centroid[i + cluster_num*1];
        lab.b = centroid[i + cluster_num*2];
        rgb = Lab2RGB(lab);
        color_set[i] = rgb;
        printf("NumberCluster[%d]:%d \n", i, num_cluster[i]);
        printf("l:%f,a:%f,b:%f \n", lab.L, lab.a, lab.b);
        printf("r:%d,g:%d,b:%d \n", rgb.R, rgb.G, rgb.B);
        pn;
    }

    free(x);
    free(centroid);
    free(cluster_weight);

    return color_set;
}


// Kmeansの初期値を決定：セントロイド間の距離が最大になるように置いていく
void KmeansPlus(double* x, int data_num, int dimension_num, int cluster_num, double* centroid)
{
    int cen=0,c=0,d=0,dim,tmp;

    // 第一セントロイドを適当に決定
    for (dim = 0; dim < dimension_num; dim++) {
        centroid[cen + cluster_num*dim] =x[d + data_num*dim];
    }

    //残りのセントロイドを順に決定
    for (cen = 1; cen<cluster_num; cen++) {
        double max_dist = 0;
        /* データn個に対するループ */
        /* それぞれの点で、最も近いセントロイドの距離を測る */
        for(d = 0; d < data_num; d++) {
            double min_cent_dist = 99999999;
            /* クラスタ数に対するループ */
            /* それぞれのセントロイドとの距離を測り最小の距離を見つける */
            for(c = 0; c < cen; c++) {
                double dd = 0.0;
                /* 次元に対するループ。点x_iにおける、各次元でのj番目のクラスタの中心との二乗和を取っている */
                for(dim = 0; dim < dimension_num; dim++) {
                    tmp = x[d+data_num*dim] - centroid[c+cluster_num*dim];
                    dd += tmp * tmp;
                }
                /* x_iと一番近いクラスタとの距離に更新 */
                if(dd < min_cent_dist) {
                    min_cent_dist = dd;
                }
            }

            /* この点における最小セントロイド距離が最大なら更新 */
            if(max_dist < min_cent_dist){
                max_dist = min_cent_dist;
                for (dim = 0; dim < dimension_num; dim++) {
                    centroid[cen + cluster_num*dim] =x[d + data_num*dim];
                }
            }
        }
    }
}


// 与えられたデータ、クラスタ数を用いてKmeansによりクラスタリング
void Kmeans(double *x, int data_num, int dimension_num, double *cen, int cluster_num, int *cl, int maxiter, int *nc, double *wss)
{
    // int data_num = *pn, cluster_num = *pk, dimension_num = *pp, maxiter = *pmaxiter;
    int iter, d, c, dim, it, inew = 0;
    double best, dd, tmp;
    int updated;
    /* cl[d]はd番目のデータが第何クラスタ[1~]に所属しているかを表わす */
    for(d = 0; d < data_num; d++) cl[d] = -1;

    for(iter = 0; iter < maxiter; iter++) {
        // printf("iteration:%d\n",iter+1);
        updated = FALSE;
        /* データn個に対するループ */
        /* それぞれの点で、最も近いクラスタの中心(重心)を見つける */
        for(d = 0; d < data_num; d++) {
            best = 99999999;
            /* クラスタ数に対するループ */
            for(c = 0; c < cluster_num; c++) {
                dd = 0.0;
                /* 次元に対するループ。点x_iにおける、各次元でのj番目のクラスタの中心との二乗和を取っている */
                for(dim = 0; dim < dimension_num; dim++) {
                    tmp = x[d+data_num*dim] - cen[c+cluster_num*dim];
                    dd += tmp * tmp;
                }
                /* x_iと一番近いクラスタとの距離に更新 */
                if(dd < best) {
                    best = dd;
                    inew = c+1;
                }
            }
            if(cl[d] != inew) {
                updated = TRUE;
                cl[d] = inew;
            }
        }
        if(!updated) break;

        /* 各クラスタの中心の点を更新 */
        for(c = 0; c < cluster_num*dimension_num; c++) cen[c] = 0.0;
        for(c = 0; c < cluster_num; c++) nc[c] = 0;
        for(d = 0; d < data_num; d++) {
            it = cl[d] - 1; nc[it]++;
            for(dim = 0; dim < dimension_num; dim++) cen[it+dim*cluster_num] += x[d+dim*data_num];
        }
        for(c = 0; c < cluster_num*dimension_num; c++) cen[c] /= nc[c % cluster_num];
    }

    /* 各クラスタのセントロイドと属する点のの距離の二乗総和wss */
    for(c = 0; c < cluster_num; c++) wss[c] = 0.0;
    for(d = 0; d < data_num; d++) {
        it = cl[d] - 1;
        for(dim = 0; dim < dimension_num; dim++) {
            tmp = x[d+data_num*dim] - cen[it+cluster_num*dim];
            wss[it] += tmp * tmp;
        }
    }
}


// カラーセットを帯画像で割合表現
PPM* Visualize_ColorSet(RGB* color_set, int cluster_num, int* num_cluster) {
	int c,i,j;

	PPM *ppm = (PPM *)malloc(sizeof(PPM));
	memcpy(ppm->descriptor, "P3", 3);
	ppm->height = 64;
	ppm->width = 512;
	ppm->bright = 255;
	ppm->dataR = create_ally(ppm->width, ppm->height);
	ppm->dataG = create_ally(ppm->width, ppm->height);
	ppm->dataB = create_ally(ppm->width, ppm->height);
    format_ally(ppm->dataR, ppm->width, ppm->height, ppm->bright);
    format_ally(ppm->dataG, ppm->width, ppm->height, ppm->bright);
    format_ally(ppm->dataB, ppm->width, ppm->height, ppm->bright);


    int sum_cn = 0;
    for (i = 0; i < cluster_num; i++) {
        sum_cn += num_cluster[i];
    }

    int paint_sum=0;
    int paint_num=0;
    // クラスタの割合分、左から塗りつぶしていく
    for (c = 0; c < cluster_num; c++)
    {
        paint_num = ppm->width * (double)num_cluster[c]/sum_cn + 0.5;
        for(i=paint_sum; i<paint_sum + paint_num; i++) {
            if(i >= ppm->width) break;
            for(j=0; j<ppm->height; j++) {
                ppm->dataR[i][j] = color_set[c].R;
                ppm->dataG[i][j] = color_set[c].G;
                ppm->dataB[i][j] = color_set[c].B;
            }
        }
        paint_sum += paint_num;
    }

	return ppm;
}


// カラーセットを帯画像で割合表現
PPM* Visualize_KmeanImg(PPM* in, RGB* color_set, int** x_centlabel) {
	int i,j;
    int CentLabel;
	PPM *KmeanImg = create_ppm(in->width, in->height, in->bright);

    // データのセントロイドラベルに従いカラーセットを配置していく
    for(i=0; i<in->width; i++) {
        for(j=0; j<in->height; j++) {
            CentLabel = x_centlabel[i][j]-1;
            KmeanImg->dataR[i][j] = color_set[CentLabel].R;
            KmeanImg->dataG[i][j] = color_set[CentLabel].G;
            KmeanImg->dataB[i][j] = color_set[CentLabel].B;
        }
    }

	return KmeanImg;
}


//JIS規格カラーセットを作成
RGB* create_JIS_ColorSet(int cluster_num)
{
    RGB* color_set = (RGB *)malloc(sizeof(RGB) * cluster_num);

	for (int i = 0; i < cluster_num; i++){
		switch (i)
		{
			//rose madder
			case 0:
				color_set[0].R = 241;
				color_set[0].G = 91;
				color_set[0].B = 134;
				break;
			//red
			case 1:
				color_set[1].R = 206;
				color_set[1].G = 34;
				color_set[1].B = 60;
				break;
			//vermilion
			case 2:
				color_set[2].R = 222;
				color_set[2].G = 56;
				color_set[2].B = 56;
			//brown
			case 3:
				color_set[3].R = 148;
				color_set[3].G = 73;
				color_set[3].B = 61;
				break;
			//vandyke brown
			case 4:
				color_set[4].R = 91;
				color_set[4].G = 68;
				color_set[4].B = 49;
				break;
			//raw sienna
			case 5:
				color_set[5].R = 177;
				color_set[5].G = 99;
				color_set[5].B = 42;
				break;
			//orange
			case 6:
				color_set[6].R = 255;
				color_set[6].G = 113;
				color_set[6].B = 68;
				break;
			//burnt umber
			case 7:
				color_set[7].R = 94;
				color_set[7].G = 66;
				color_set[7].B = 52;
				break;
			//yellow ochre
			case 8:
				color_set[8].R = 192;
				color_set[8].G = 132;
				color_set[8].B = 42;
				break;
			//chrome yellow
			case 9:
				color_set[9].R = 250;
				color_set[9].G = 190;
				color_set[9].B = 0;
				break;
			//naples yellow
			case 10:
				color_set[10].R = 254;
				color_set[10].G = 205;
				color_set[10].B = 121;
				break;
			//yellow
			case 11:
				color_set[11].R = 244;
				color_set[11].G = 213;
				color_set[11].B = 0;
				break;
			//olive
			case 12:
				color_set[12].R = 129;
				color_set[12].G = 121;
				color_set[12].B = 35;
				break;
			//lemon yellow
			case 13:
				color_set[13].R = 238;
				color_set[13].G = 234;
				color_set[13].B = 0;
				break;
			//olive green
			case 14:
				color_set[14].R = 85;
				color_set[14].G = 101;
				color_set[14].B = 45;
				break;
			//yellow green
			case 15:
				color_set[15].R = 101;
				color_set[15].G = 192;
				color_set[15].B = 66;
				break;
			//grass green
			case 16:
				color_set[16].R = 96;
				color_set[16].G = 145;
				color_set[16].B = 54;
				break;
			//green
			case 17:
				color_set[17].R = 0;
				color_set[17].G = 140;
				color_set[17].B = 75;
				break;
			//emerald green
			case 18:
				color_set[18].R = 0;
				color_set[18].G = 168;
				color_set[18].B = 105;
				break;
			//cobalt green
			case 19:
				color_set[19].R = 0;
				color_set[19].G = 164;
				color_set[19].B = 129;
				break;
			//deep green
			case 20:
				color_set[20].R = 27;
				color_set[20].G = 81;
				color_set[20].B = 63;
				break;
			//viridian
			case 21:
				color_set[21].R = 0;
				color_set[21].G = 73;
				color_set[21].B = 73;
				break;
			//blue green
			case 22:
				color_set[22].R = 0;
				color_set[22].G = 123;
				color_set[22].B = 122;
				break;
			//sky blue
			case 23:
				color_set[23].R = 46;
				color_set[23].G = 152;
				color_set[23].B = 204;
				break;
			//cerulean blue
			case 24:
				color_set[24].R = 0;
				color_set[24].G = 130;
				color_set[24].B = 181;
				break;
			//cobalt blue
			case 25:
				color_set[25].R = 0;
				color_set[25].G = 96;
				color_set[25].B = 175;
				break;
			//prussian blue
			case 26:
				color_set[26].R = 42;
				color_set[26].G = 61;
				color_set[26].B = 97;
				break;
			//ultramarine light
			case 27:
				color_set[27].R = 82;
				color_set[27].G = 112;
				color_set[27].B = 201;
				break;
			//ultramarine
			case 28:
				color_set[28].R = 14;
				color_set[28].G = 63;
				color_set[28].B = 158;
				break;
			//violet
			case 29:
				color_set[29].R = 71;
				color_set[29].G = 59;
				color_set[29].B = 135;
				break;
			//light purple
			case 30:
				color_set[30].R = 159;
				color_set[30].G = 90;
				color_set[30].B = 174;
				break;
			//purple
			case 31:
				color_set[31].R = 97;
				color_set[31].G = 52;
				color_set[31].B = 116;
				break;
			//red purple
			case 32:
				color_set[32].R = 194;
				color_set[32].G = 42;
				color_set[32].B = 119;
				break;
			//white
			case 33:
				color_set[33].R = 234;
				color_set[33].G = 234;
				color_set[33].B = 234;
				break;
			//gray
			case 34:
				color_set[34].R = 156;
				color_set[34].G = 156;
				color_set[34].B = 156;
				break;
			//black
			case 35:
				color_set[35].R = 52;
				color_set[35].G = 52;
				color_set[35].B = 52;
				break;
		}
	}

    return color_set;
}


// RGBカラーセットをLabの明るい順に整列
void arrange_ColorSet_inLight(RGB* ColorSet, int cluster_num)
{
	// 比較のため入力カラーセットをLabに変換
	Lab* Lab_ColorSet= (Lab*)malloc(sizeof(Lab)*cluster_num);
	for (int i = 0; i < cluster_num; i++) {
		Lab_ColorSet[i] = RGB2Lab(ColorSet[i]);
	}

	// 入力カラーセットをLab明度に従い整列
	for (int i = 0; i < cluster_num; i++) {
		for (int j = 0; j < cluster_num; j++) {
			if(Lab_ColorSet[i].L > Lab_ColorSet[j].L){
				RGB tmpRGB = ColorSet[i];
				ColorSet[i] = ColorSet[j];
				ColorSet[j] = tmpRGB;
				Lab tmpLab = Lab_ColorSet[i];
				Lab_ColorSet[i] = Lab_ColorSet[j];
				Lab_ColorSet[j] = tmpLab;
			}
		}
	}
}




/////////////////////////////////////////////
//			その他の関数群
/////////////////////////////////////////////

//[value] -> "~~value"
void Add_NUM_to_sentence(char* sentence, int letterNum, const char * NumFormat,int value){
	char value_word[16];

	snprintf(value_word, letterNum, NumFormat, value);
	strcat(sentence, value_word);
}


//[name, value] -> "name : value\n"
void Add_dictionary_to_sentence(char* sentence, char *name, int value){
	char value_word[16];

	strcat(sentence, name);
	strcat(sentence, " : ");
	snprintf(value_word, 16, "%d", value);
	strcat(sentence, value_word);
	strcat(sentence, "\r\n");
}


//[name, value] -> "name : value\n"
void Add_dictionary_to_sentence_d(char* sentence, char *name, double value){
	char value_word[16];

	strcat(sentence, name);
	strcat(sentence, " : ");
	snprintf(value_word, 16, "%f", value);
	strcat(sentence, value_word);
	strcat(sentence, "\r\n");
}


//logテキストファイルを生成する
int log_print(char* filename, char *sentence, char *mode){
	FILE* fp;

	if((fp = fopen(filename, mode)) == NULL){
		fprintf(stderr, "W_Open_Error(log_print)\n");
		return -1;
	}

	fprintf(fp, "%s\n", sentence);

	fclose(fp);
	return 0;
}


//ストロークデータをファイルに追加
int vec_print(char* filename, Point *p, int pnum, int brightR, int brightG, int brightB, int width, int height){
	char vec_sentence[256] = "";
	char tmp_sentence[64] = "";
	int i;

	snprintf(tmp_sentence, 32, "%d", brightR);
	strcat(vec_sentence, tmp_sentence);
	strcat(vec_sentence, ",");
	snprintf(tmp_sentence, 32, "%d", brightG);
	strcat(vec_sentence, tmp_sentence);
	strcat(vec_sentence, ",");
	snprintf(tmp_sentence, 32, "%d", brightB);
	strcat(vec_sentence, tmp_sentence);
	strcat(vec_sentence, " ");

	for(i=0; i<pnum; i++){
		snprintf(tmp_sentence, 32, "%f", p[i].x/(double)width);
		strcat(vec_sentence, tmp_sentence);
		strcat(vec_sentence, ",");
		snprintf(tmp_sentence, 32, "%f", p[i].y/(double)height);
		strcat(vec_sentence, tmp_sentence);
		strcat(vec_sentence, " ");
	}

	log_print(filename, vec_sentence, "a");
	return 0;
}


// 文字列から拡張子を取得
char *get_extension(char *name) {
  int i;
  for (i = strlen(name) - 1; i >= 0; i--) {
    if (name[i] == '.') {
      return &name[i + 1];
    }
  }
  return NULL;
}


// 整数をchar文字列の最後尾に文字として付ける
void Add_num(char sentence[], int num)
{
	char num_chara[32];

	snprintf(num_chara, 32, "%d", num);
	strcat(sentence, num_chara);
}


double my_clock() {
	static int start_flag=1;
    static struct timespec start;
	struct timespec end;

	if(start_flag){
    	clock_gettime(CLOCK_MONOTONIC, &start);
		start_flag=0;
		return 0;
	}

	clock_gettime(CLOCK_MONOTONIC, &end);

	return (double)(end.tv_sec-start.tv_sec)+(double)(end.tv_nsec-start.tv_nsec)/1e+9;
}



/////////////////////////////////////////////
//			ストロークレンダリング関数群
/////////////////////////////////////////////


//与えられたパスの中のフォルダの中の画像すべてに対してC_illust_brushを実行
void multi_dir_CIB(char* in_path){
	int i;
	int dir_num, file_num;
	char** file_list;

	char** dir_list = make_file_list(&dir_num, in_path, Search_DIRECTRY); //フォルダリストを取得

	for(i=0; i<dir_num; i++) {
		//i番目フォルダのファイルリストを取得
		file_list = make_file_list(&file_num, dir_list[i], Search_FILE);
		//ファイルリストの画像それぞれに対してC_illust_brushを実行
		multi_file_CIB(file_list, file_num);

		for(i=0; i < file_num; i++) { free(file_list[i]);}
		free(file_list);
	}

	for(i=0; i < dir_num; i++) { free(dir_list[i]);}
	free(dir_list);
}

//与えられたパスリストの画像すべてに対してC_illust_brushを実行
void multi_file_CIB(char** list, int list_num){
	int i;
	PPM *in_data, *trans_data;
	image_t *in_img;
	char out_dir_name[256];

	for(i=0; i<list_num; i++) {
		strcpy(out_dir_name, list[i]);
		strtok(out_dir_name,".");		//拡張子は取る
		printf("out_dir_name:%s\n",out_dir_name);
		printf("list[i]:%s\n",list[i]);

		if((in_img = read_jpeg_file(list[i])) == NULL){  // <<<<<====対象ファイルの拡張子によって変更が必要
			printf("PNG_READ_ERROR\n");		exit(1);
		}
		dump_image_info(in_img);	//PNG情報出力
		in_data = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
		free_image(in_img);

		trans_data = c_Illust_brush(in_data, list[i]);
		FreePPM(in_data);
		FreePPM(trans_data);

	}

}


//与えられたパスのディレクトリ中にあるファイルもしくはディレクトリをリスト化し配列で返す
//引数に従ってどちらかしか検索できない
char** make_file_list(int* file_num, char* search_path, int Search_TYPE){
	int i=0;
	DIR *dir;
	char path[256];
    struct dirent *dp;
    struct stat fi;
	*file_num=0;

	dir = opendir(search_path);
	p("Search_TYPE",Search_TYPE);
    for (dp = readdir(dir); dp != NULL; dp = readdir(dir)) {
        if (dp->d_name[0] == '.') continue;

	    strcpy(path, search_path);
	    strcat(path, "/");
	    strcat(path, dp->d_name);
        stat(path, &fi);	//pathのファイルの情報を取得

    	switch(Search_TYPE){
	    	case Search_FILE:
				if (!S_ISDIR(fi.st_mode)) {
					i++;
		            printf("FILE[%d]:%s\n",i, path);
		        }
    			break;
	    	case Search_DIRECTRY:
		        if (S_ISDIR(fi.st_mode)) {
					i++;
		            printf("DIR[%d]:%s\n",i, path);
		        }
    			break;
    	}
    }
	*file_num = i;
	printf("file_num:%d\n", *file_num);

	char** list = (char**)malloc(sizeof(char*)*(*file_num));
	for(i=0; i<*file_num; i++) { list[i]=(char*)malloc(sizeof(char)*(256));}

	i=0;
	rewinddir(dir);
	for (dp = readdir(dir); dp != NULL; dp = readdir(dir)) {
        if (dp->d_name[0] == '.') continue;

	    strcpy(path, search_path);
	    strcat(path, "/");
	    strcat(path, dp->d_name);
        stat(path, &fi);	//pathのファイルの情報を取得

    	switch(Search_TYPE){
	    	case Search_FILE:
				if (!S_ISDIR(fi.st_mode)) {
					strcpy(list[i], path);
					i++;
		            printf("FILE[%d]:%s\n",i, path);
		        }
    			break;
	    	case Search_DIRECTRY:
		        if (S_ISDIR(fi.st_mode)) {
					strcpy(list[i], path);
					i++;
		            printf("DIR[%d]:%s\n",i, path);
		        }
    			break;
    	}
    }

    closedir(dir);
	return list;//file_num;
}

