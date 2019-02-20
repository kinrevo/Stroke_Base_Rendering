#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"pgm.h"

#define t1_8pi 0.4142
#define t3_8pi 2.4142
#define maxValue 0.35
#define minValue 0.25



//平均化フィルタ
PGM *average_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *img = (PGM *)malloc(sizeof(PGM));
	memcpy(img->descriptor, pgm->descriptor, 3);
	img->height = pgm->height;
	img->width = pgm->width;
	img->bright = pgm->bright;
	img->data = create_ally(pgm->width, pgm->height);
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				img->data[i][j] = pgm->data[i][j];}
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l];}
				img->data[i][j] = sum/9;
				sum = 0;
			}
		}
	}
	
	
	return img;
}



//与えられた閾値でイメージを二値化する関数
PGM *binarization_filter(PGM *pgm, int boundary) 
{
	int i,j;
	PGM *nimg = copy_pgm(pgm);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			//ピクセルを閾値で分けて二値化する
			if(pgm->data[i][j] < boundary) {
				nimg->data[i][j] = 0;				
			}else {
				nimg->data[i][j] = pgm->bright;
			}
		}
	}
	
	
	return nimg;
	
}



//白黒反転フィルタ
PGM *inversion_filter(PGM *pgm) 
{
	int i,j;
	PGM *nimg = copy_pgm(pgm);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			//ピクセルを輝度の最大値から反転する
			nimg->data[i][j] = pgm->bright - pgm->data[i][j];	
		}
	}
	
	
	return nimg;
	
}



//加重平均フィルタ
PGM *weightaverage_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filter[3][3] = {
		{1,2,1},
		{2,6,2},
		{1,2,1},
	};
	int fil_sum = 18;
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = pgm->data[i][j];}
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l] * filter[k+1][l+1];
				}
				nimg->data[i][j] = sum/fil_sum;
				sum = 0;
			}
		}
	}
	
	
	return nimg;
}



//縦平滑フィルタ
PGM *verticalsmooth_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filter[3][3] = {
		{0,1,0},
		{0,1,0},
		{0,1,0},
	};
	int fil_sum = 3;
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = pgm->data[i][j];}
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l] * filter[k+1][l+1];
				}
				nimg->data[i][j] = sum/fil_sum;
				sum = 0;
			}
		}
	}
	
	
	return nimg;
}




//横平滑フィルタ
PGM *horizonsmooth_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filter[3][3] = {
		{0,0,0},
		{1,1,1},
		{0,0,0},
	};
	int fil_sum = 3;
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = pgm->data[i][j];}
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l] * filter[k+1][l+1];
				}
				nimg->data[i][j] = sum/fil_sum;
				sum = 0;
			}
		}
	}
	
	
	return nimg;
}




//prewittフィルタ
PGM *prewitt_filter(PGM *pgm)
{
	int i,j,k,l,sum=0,sumx=0,sumy=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filterx[3][3] = {
		{-1,0,1},
		{-1,0,1},
		{-1,0,1},
	};
	int filtery[3][3] = {
		{-1,-1,-1},
		{0,0,0},
		{1,1,1},
	};
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = 0;
			}else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						sumx += pgm->data[i+k][j+l] * filterx[k+1][l+1];
						sumy += pgm->data[i+k][j+l] * filtery[k+1][l+1];
					}
				}
				sum = (int)sqrt((double)(sumx*sumx*+sumy*sumy));
				if(sum > pgm->bright){ nimg->data[i][j]=pgm->bright;
				}else if(sum<0){ nimg->data[i][j]=0;
				}else { nimg->data[i][j]=sum; }
				sum = sumx = sumy = 0;
			}
		}
	}
	
	
	return nimg;
}




//Laplacianフィルタ
PGM *Laplacian_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filter[3][3] = {
		{0,1,0},
		{1,-4,1},
		{0,1,0},
	};
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = 0;}
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l] * filter[k+1][l+1];
				}
				if(sum > pgm->bright){ nimg->data[i][j]=pgm->bright;
				}else if(sum<0){ nimg->data[i][j]=0;
				}else { nimg->data[i][j]=sum; }
				sum = 0;
			}
		}
	}
	
	
	return nimg;
}



//ガウシアンフィルタ半径１
PGM *gaucian_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
	int filter[3][3] = {
		{1,2,1},
		{2,4,2},
		{1,2,1},
	};
	int fil_sum = 16;
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = pgm->data[i][j]; }
			else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l] * filter[k+1][l+1];
				}
				nimg->data[i][j] = sum/fil_sum;
				sum = 0;
			}
		}
	}
	
	
	return nimg;
}


//sobelフィルタ
PGM *sobel_filter(PGM *pgm)
{
	int i,j,k,l,sum=0,sumx=0,sumy=0;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(pgm);
	
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
	
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg->data[i][j] = 0;
			}else {
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						sumx += pgm->data[i+k][j+l] * filterx[k+1][l+1];
						sumy += pgm->data[i+k][j+l] * filtery[k+1][l+1];
					}
				}
				sum = (int)sqrt((double)(sumx*sumx*+sumy*sumy));
				if(sum > pgm->bright){ nimg->data[i][j]=pgm->bright;
				}else if(sum<0){ nimg->data[i][j]=0;
				}else { nimg->data[i][j]=sum; }
				sum = sumx = sumy = 0;
			}
		}
	}
	
	
	return nimg;
}



//与えられた二つのイメージの差分を取る（in2-in1）
PGM *DifferentImg(PGM *in1, PGM *in2) {
	
	int i,j,min=255;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(in1);
	
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//差分を取った後二値化する
			nimg->data[i][j] = in2->data[i][j] - in1->data[i][j];
			//if(nimg->data[i][j]<0) nimg->data[i][j]=0;
			//nimg->data[i][j] = (nimg->data[i][j]>0 ?nimg->bright:0);
			if(nimg->data[i][j] < min) min=nimg->data[i][j];
		}
	}
	
	min = (-1)*min;
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//底上げ
			nimg->data[i][j] += min; 
			nimg->data[i][j] = (nimg->data[i][j] * nimg->bright)/(min+nimg->bright);
		}
	}
	
	return nimg;
}



//与えられた二つのイメージを比較する（in2-in1）
PGM *CompareImg(PGM *in1, PGM *in2) {
	
	int i,j;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(in1);
	
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//差分を取った後二値化する
			if(in2->data[i][j] != in1->data[i][j]) {
				nimg->data[i][j] = 0;
				printf("[%d][%d]:%d\n",i ,j ,in2->data[i][j] - in1->data[i][j]);
			}
			else{
				nimg->data[i][j] = nimg->bright;
			}
		}
	}
	
	return nimg;
}
	

//与えられたイメージを時計回りにn回転する
PGM *RotateImg(PGM *in, int n) {
	
	int i,j,k;
	//新しい画像データをコピーして作っておく
	PGM *nimg = copy_pgm(in);
	
	for(k=0; k<n; k++) {
		for(i=0; i<in->height; i++) {
			for(j=0; j<in->width; j++) {
				nimg->data[j][in->height-1-i] = in->data[i][j];
			}
		}
	}
	
	

	return nimg;
}



//キャニーエッジ検出器
PGM *cannyedge_detector(PGM *pgm)
{
	int i,j,k,l,flag=0,sum=0,sumx=0,sumy=0,min=pgm->bright,max=0;
	//新しい画像データを作りガウシアンフィルタを掛けたものを入れる(nimg1)
	PGM *nimg1 = gaucian_filter(pgm);
	printf("step1 finish\n");
	//画像データを作りsobelフィルタをかける(nimg2)
	PGM *nimg2 = copy_pgm(nimg1);
	
	double **theta = (double**)malloc(sizeof(double*)*(pgm->height));
	if(theta==NULL) {fprintf(stderr, "fail memory\n");}
	for(i=0; i<pgm->height; i++) { theta[i]=(double*)malloc(sizeof(double)*(pgm->width));}
	//double theta[pgm->height][pgm->height];
	//double theta[100][100];
	
	printf("step2-1 finish\n");
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
	
	printf("ph=%d\n",pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) {
				nimg2->data[i][j] = 0;
			}else { 
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						sumx += nimg1->data[i+k][j+l] * filterx[k+1][l+1];
						sumy += nimg1->data[i+k][j+l] * filtery[k+1][l+1];
					}
				}
				
				sum = (int)sqrt((double)(sumx*sumx+sumy*sumy));
				/*
				if(sum > nimg1->bright){ nimg2->data[i][j]=nimg1->bright;
				}else if(sum<0){ nimg2->data[i][j]=0;
				}else { nimg2->data[i][j]=sum; }
				*/ nimg2->data[i][j]=sum;
				if(sumx==0) {theta[i][j] = (double)(sumy*10000);} //分母がゼロにならないようにする
				else { theta[i][j] = (double)(sumy/sumx);}  //次の処理のため角度を計算し格納しておく
				
				if(sum>max) {max=sum;}  //最大小値をスケーリングのため保存しておく
				if(sum<min) {min=sum;}
				
				sum = sumx = sumy = 0;
			}
		}
	}
	
	//スケーリング処理を行う
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			nimg2->data[i][j] += (-1)*min;
			nimg2->data[i][j] = (nimg2->data[i][j]*pgm->bright)/(max-min);
			if(nimg2->data[i][j] > 255) printf("error\n");
		}
	}
	printf("max:%d min:%d\n",max,min);
	
	printf("step2 finish\n");
	
	
	//非極大値抑制処理を行う(nimg3)
	PGM *nimg3 = copy_pgm(nimg2);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) { 
				nimg3->data[i][j] = nimg2->data[i][j];
			}else {
				//注目しているピクセルのエッジ方向の鉛直方向の角度を取る　(x,y)=(k,l)
				//-1/8pi < 1/8pi
				if( (-1)*t1_8pi <= theta[i][j] && theta[i][j] <= t1_8pi ) { k=1; l=0; }
				//1/8pi < 3/8pi
				else if( t1_8pi <= theta[i][j] && theta[i][j] <= t3_8pi ) { k=1; l=1; }
				//3/8pi < 5/8pi
				else if( t3_8pi <= theta[i][j] && theta[i][j] <= (-1)*t3_8pi ) { k=0; l=1; }
				//5/8pi < 7/8pi
				else  {	k=-1; l=1; }
				
				//注目しているピクセルをエッジ鉛直方向の隣接画素と比較し最大でなければ0とする
				if( nimg2->data[i][j] < nimg2->data[i-l][j-k] || nimg2->data[i][j] < nimg2->data[i+l][j+k] ) {
					nimg3->data[i][j] =0;
				}
				else { nimg3->data[i][j] = nimg2->data[i][j]; }
					
				
			}
		}
	}
	
	
	
	printf("step3 finish\n");
	//ヒステリシス閾処理を行う(nimg4)
	PGM *nimg4 = copy_pgm(nimg3);
	
	//ハイパスを超える値を最大、ローパスより低い値を最小とする
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
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
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) { 
			}else if(nimg4->data[i][j]==0 || nimg4->data[i][j]==pgm->bright){
			}else {				
				//注目しているピクセルの周りにエッジがあるか確認する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						if(nimg4->data[i+k][j+l]==nimg4->bright) { flag=1; }
						if(nimg4->data[i+k][j+l]>minValue*nimg4->bright && nimg4->data[i+k][j+l]<maxValue*nimg4->bright) { flag2=1; }
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
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(nimg4->data[i][j]!=pgm->bright){
				nimg4->data[i][j]=0;
			}
		}
	}
	
	
	Freeimg(nimg1);
	Freeimg(nimg2);
	Freeimg(nimg3);
	for(i=0; i < pgm->height; i++) { free(theta[i]);}
	free(theta);
	
	return nimg4;
}






//SIFT変換
PGM *SIFT(PGM *pgm, int S) {
	
	int i,j,x,y,s,ds,dx,dy,temp;
	double sigma = 1.6;
	PGM *Gimg[S+3];  //平滑化画像を入れる配列
	PGM *DoG[S+2];  //DoG
	int is_max,is_min;
	double size = sigma;
	
	//増加率kを算出
	double k = pow(4, 1/(double)S);
	pd("k",k);
	
	
	//
	for(j=0; j<(S+3); j++){
		Gimg[j] = gaussian_filter(pgm, sigma);
		sigma *= k;
	}
	

	for(j=0; j<(S+2); j++){
		DoG[S+1-j] = DifferentImg(Gimg[S+1-j], Gimg[S+2-j]);
	}
	
	/*
	int **keypt[S];  //キーポイント
	
	for(s=0; s<S; s++){
		keypt[s] = create_ally(pgm->width, pgm->height);
	}
	*/
	
	//keyptの値は平滑化画像の指定子（例えば2ならスケールはk^2σである）
	double **keypt = create_dally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { keypt[i][j]=0; }
	}
	
	int **flag = create_ally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { flag[i][j]=0; }
	}
	
	int kpcount=0;

	//極値の検出
	for(s=1; s<=S; s++){
		//size = size*k;
		p("s",s);
		for(x=0; x<pgm->height; x++){
        	for(y=0; y<pgm->width; y++){
        		is_max = is_min = 1;
        		
        		//既にキーポイントと検出しているか計算できない領域だと0
        		if(x==0 || x==pgm->height-1 || y==0 || y==pgm->width-1) {}
        		else if(keypt[x][y]>0) {}
        		//近傍と比較
        		else {
        			for(ds=-1; ds<=1; ds++) { 
						for(dx=-1; dx<=1; dx++)  {
							for(dy=-1; dy<=1; dy++)  {
								if(ds==0 && dx==0 && dy==0){} //中心は無視
								else { 
									if(DoG[s+ds]->data[x+dx][y+dy] >= DoG[s]->data[x][y]) {
										is_max=0;  //極大値でないと判定
									}
									if(DoG[s+ds]->data[x+dx][y+dy] <= DoG[s]->data[x][y]) {
										is_min=0;  //極小値でないと判定
									}
								}
							}
						}
					}
        			if(is_max==0 && is_min==0) keypt[x][y]=0;  //極大でも極小でもないなら0とする
        			else { 
        				keypt[x][y]=s; //size; 
        				kpcount++;
        			}
				}
        	}
		}
	}
	
	p("kp",kpcount);
	
	
	//-----------------------
	//キーポイント削除
	//-----------------------
	double Dxx, Dyy, Dxy, Tr, Det;
	double Dx, Dy, Ds, Dss, Dxs, Dys, detmD, Dpow;
	double mD[3][3], iD[3][3], X[3], xD[3];
	int sm2, sp2;
	double TH_POW = 0.05*pgm->bright; //しきい値（最大輝度の3％）

	for(s=1; s<=S; s++){
		for(x=2; x<pgm->height-2; x++){
        	for(y=2; y<pgm->width-2; y++){
        		//主曲線によるキーポイント削除
        		if(keypt[x][y]>0){
        			
        			Dxx = DoG[s]->data[x-2][y] + DoG[s]->data[x+2][y] - 2*DoG[s]->data[x][y];
        			Dyy = DoG[s]->data[x][y-2] + DoG[s]->data[x][y+2] - 2*DoG[s]->data[x][y];
        			Dxy = (DoG[s]->data[x-1][y-1] - DoG[s]->data[x+1][y-1]) - (DoG[s]->data[x-1][y+1] - DoG[s]->data[x+1][y+1]);
        			
        			Tr = Dxx+Dyy;
        			Det = Dxx*Dyy-(Dxy*Dxy);
        			
        			if((Tr*Tr/Det) >  BORDER) {
        				keypt[x][y]=0;
        				kpcount--;
        			}
        		}
        		
        		//サブピクセル推定
        		if(keypt[x][y]>0){
        			sm2=(s-2<0) ? 0:s-2;
        			sp2=(s+2>S+1) ? S+1:s+2;
        			
        			Dx=(DoG[s]->data[x-1][y]-DoG[s]->data[x+1][y]);
        			Dy=(DoG[s]->data[x][y-1]-DoG[s]->data[x][y+1]);
        			Ds=(DoG[s-1]->data[x][y]-DoG[s+1]->data[x][y]);
        			
        			Dss=(DoG[sm2]->data[x][y]-DoG[sp2]->data[x][y]+2*DoG[s]->data[x][y]);
        			Dxs=(DoG[s-1]->data[x-1][y]-DoG[s-1]->data[x+1][y]+DoG[s-1]->data[x-1][y]-DoG[s+1]->data[x+1][y]);
        			Dys=(DoG[s-1]->data[x][y-1]-DoG[s-1]->data[x][y+1]+DoG[s+1]->data[x][y-1]-DoG[s+1]->data[x][y+1]);

        			mD[0][0]=Dxx; mD[0][1]=Dxy; mD[0][2]=Dxs;
        			mD[1][0]=Dxy; mD[1][1]=Dyy; mD[1][2]=Dys;
        			mD[2][0]=Dxs; mD[2][1]=Dys; mD[2][2]=Dss;
        			
        			xD[0]=-Dx; xD[1]=-Dy; xD[2]=-Ds;
        			
        			//逆行列計算(mDの逆行列をiDに)
        			detmD = mD[0][0]*mD[1][1]*mD[2][2] + mD[1][0]*mD[2][1]*mD[0][2] + mD[2][0]*mD[0][1]*mD[1][2]
        					- mD[0][0]*mD[2][1]*mD[1][2] - mD[2][0]*mD[1][1]*mD[0][2] - mD[1][0]*mD[0][1]*mD[2][2];
        			if(detmD==0) continue;
        			
        			iD[0][0]=mD[1][1]*mD[2][2]-mD[1][2]*mD[2][1]; iD[0][1]=mD[0][2]*mD[2][1]-mD[0][1]*mD[2][2]; iD[0][2]=mD[0][1]*mD[1][2]-mD[0][2]*mD[1][1];
        			iD[1][0]=mD[1][2]*mD[2][0]-mD[1][0]*mD[2][2]; iD[1][1]=mD[0][0]*mD[2][2]-mD[0][2]*mD[2][0]; iD[1][2]=mD[0][2]*mD[1][0]-mD[0][0]*mD[1][2];
        			iD[2][0]=mD[1][0]*mD[2][1]-mD[1][1]*mD[2][0]; iD[2][1]=mD[0][1]*mD[2][0]-mD[0][0]*mD[2][1]; iD[2][2]=mD[0][0]*mD[1][1]-mD[0][1]*mD[1][0];
        			
        			for(i=0; i<3; i++){
        				for(j=0; j<3; j++){ iD[i][j] = iD[i][j]/detmD;   }
        			}
        			
        			//サブピクセル位置(行列の積)
        			X[0]=iD[0][0]*xD[0]+iD[0][1]*xD[1]+iD[0][2]*xD[2];
        			X[1]=iD[1][0]*xD[0]+iD[1][1]*xD[1]+iD[1][2]*xD[2];
        			X[2]=iD[2][0]*xD[0]+iD[2][1]*xD[1]+iD[2][2]*xD[2];
        			//サブピクセル位置での出力
        			Dpow=fabs(DoG[s]->data[x][y]+(xD[0]*X[0]+xD[1]*X[1]+xD[2]*X[2])/2);
        			
        			if(Dpow<TH_POW) {
        				keypt[x][y]=0;
        				kpcount--;
        			}
        		}
        		
        	}	
        }
	}
	
	
	
	p("kp2",kpcount);
    
	
	
	
	//----------------------
	//オリエンテーションの算出
	//----------------------

	//キーポイントごとのオリエンテーションを保存する配列
	int **orient = create_ally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { orient[i][j]=0; }
	}
	
	double fu, fv;
	double histogram[36]={};
	int w = (int)( ceil(3.0*sigma+0.5)*2-1 ); //ガウスの値が大体0になる距離（直径）
	int c=(w-1)/2;  //ガウスフィルタの半径
	double **filter; //ガウスフィルタを格納する配列のポインタ
	double Fpow;
	int Farg, peak, im1,ip1;
	
	
	//キーポイントから平滑化画像とガウシアン半径を取り出す
	for(s=1; s<=S; s++){
		sigma=1.6; //σを初期化
		sigma = sigma*pow(k, s);
		w = (int)( ceil(3.0*sigma+0.5)*2-1 );  //ガウスの値が大体0になる距離（直径）
		c=(w-1)/2;  //ガウスフィルタの半径
				
		//ガウスフィルタの領域を確保
		filter = create_dally(w, w);
		for(i=0; i<w; i++) {
			for(j=0; j<w; j++) { filter[i][j]=0; }
		}
		
		//ガウスフィルタを構成
    	for(i=0;i<w;i++){
    	    for(j=0;j<w;j++){
     		   	filter[i][j] = gause_func(i-c, j-c, sigma);
       		}
  		}
		
		for(x=0; x<pgm->height; x++){
			for(y=0; y<pgm->width; y++){
				//探索しているσのガウスフィルタを持っているキーポイントならヒストグラムを計算
				if(keypt[x][y]==s){
					//histogramを初期化
    				for(i=0;i<36;i++){ histogram[i] = 0;	}
					
					//注目しているピクセルのガウス窓のピクセル一つ一つの勾配を調べる
					for(i=-c; i<=c; i++) { 
						for(j=-c; j<=c; j++)  {
							//計算出来ない領域はスキップする
							if( (x+i-1)<0 || (x+i+1)>pgm->height-1 || (y+j-1)<0 || (y+j+1)>pgm->width-1) {}
							else {
								fu = Gimg[s]->data[x+i+1][y+j] - Gimg[s]->data[x+i-1][y+j];
								fv = Gimg[s]->data[x+i][y+j+1] - Gimg[s]->data[x+i][y+j-1];
								
								Fpow = sqrt(fu*fu+fv*fv);
								Farg = (atan2(fv,fu)/PI+1)*18;
								//ガウス窓で強度を重み付けしそれぞれの方向に加算
								histogram[Farg] += Fpow*filter[i+c][j+c];
							}
						}
					}
					
					//ヒストグラムの中で最も大きい値を探す
					peak=0;
					for(i=0;i<36;i++){
						if(peak < histogram[i]) peak=i;
					}
					orient[x][y] = peak;
					//最大値の80％以上の勾配で極値ならオリエンテーションに割り当てる
					for(i=0;i<36;i++){
						if(histogram[i] > 0.8*histogram[peak] && i!=peak){
							im1 = (i-1<0) ? 35:i-1;
							ip1 = (i+1>35) ? 0:i+1;
							if(histogram[i]>histogram[im1] && histogram[i]>histogram[ip1]){
								//orient配列には二桁毎に勾配方向の情報が保存される
								//（orient=1209なら12と9がオリエンテーション）
								orient[x][y] = orient[x][y]*100 + i;
							}
						}
					}
				}
			}
        }
		
		Free_dally(filter, w);
	}
	
	
	
	
	
	
	
	
	PGM *kimg = copy_pgm(pgm);
	for(x=0; x<kimg->height; x++){  //画像データをコピー
        for(y=0; y<kimg->width; y++){
        	kimg->data[x][y] = pgm->data[x][y]; 
        }
	}
	/*
	for(x=0; x<kimg->height; x++){  //キーポイント点表示
        for(y=0; y<kimg->width; y++){
        	if(keypt[x][y]>0){ kimg->data[x][y] = kimg->bright;	}
        }
	}
	*/
	
	double fai;
	
	for(x=0; x<kimg->height; x++){  //キーポイント円表示
        for(y=0; y<kimg->width; y++){
        	if(keypt[x][y]>0){ 
        		size = keypt[x][y]*keypt[x][y];
        		for(i=0; i<size*2; i++){
        			fai = (i/size)*2*PI;
        			dx = keypt[x][y]*cos(fai);
        			dy = keypt[x][y]*sin(fai);
        			//画像領域内に円を描く
        			if( (x+dx)>0 && (x+dx)<(pgm->height-1) && (y+dy)>0 && (y+dy)<(pgm->width-1) ){
        				kimg->data[x+dx][y+dy] = kimg->bright;
        			}
        		}
        		//オリエンテーションを描画
        		size= keypt[x][y];
        		temp=orient[x][y];
        		while(temp!=0){
        			fai=((double)(temp%100)/18 - 1)*PI;
        			
        			for(i=0; i<size; i++){
        				dx = i*cos(fai);
        				dy = i*sin(fai);
        				//円の中に線を描く
        				if( (x+dx)>0 && (x+dx)<(pgm->height-1) && (y+dy)>0 && (y+dy)<(pgm->width-1) ){
        					kimg->data[x+dx][y+dy] = kimg->bright;
        				}
        			}
        			temp = temp/100;
        		}
        	}
        }
	}
	/**/
	
	/*
	for(x=0; x<kimg->height; x++){  //キーポイント点画像作成
        for(y=0; y<kimg->width; y++){
        	kimg->data[x][y] = keypt[x][y];
        }
	}
	/**/
	return kimg;
}




//バイラテラルフィルタを用いて与えられたイメージを絵画調にする
PGM *Illustification(PGM *in, double sigma, int n) {
	
	int i,j;
	//新しい画像データをコピーして作っておく
	PGM *nimg = bilateral_filter(in, sigma, n);
	PGM *nimg2 = copy_pgm(in);
	/*for(i=0; i<in->height; i++) {
		for(j=0; j<in->width; j++) {
			//自己商画像の生成
			if( ((double)in->data[i][j] / nimg->data[i][j]) < 0.95){
				nimg->data[i][j] = 0;
			}
			//else nimg2->data[i][j] = 255;
		}
	}
	*/
	
	nimg2 = cannyedge_detector(in);
	for(i=0; i<in->height; i++) {
		for(j=0; j<in->width; j++) {
			nimg->data[i][j] = ((nimg->data[i][j]-nimg2->data[i][j])>0 ? nimg->data[i][j]-nimg2->data[i][j]:0 );
			
			//else nimg2->data[i][j] = 255;
		}
	}
	
	
	return nimg;
}



//sobelを適応した範囲の合計の勾配を返す
double sobel_slope(PGM* in_img, int x, int y, int t) {
	double theta;
	int i,j,k,l, sumx=0, sumy=0;
	
	//sobelフィルタ
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
	
	//sobelフィルタによりストローク方向を決定
	for(i=x-t; i<=x+t; i++) {
		for(j=y-t; j<=y+t; j++) {
			if(i<0 || i>in_img->width-1 || j<0 || j>in_img->height-1) {
			}else { 
				//注目している周囲一ピクセルの値を合計し平均する
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						sumx += in_img->data[i+k][j+l] * filterx[k+1][l+1];
						sumy += in_img->data[i+k][j+l] * filtery[k+1][l+1];
					}
				}
			}
		}
	}
					
	//勾配の偏りからストローク方向を計算
	if(sumx==0) {theta = (double)(sumy*10000);} //分母ゼロの例外処理
	else { theta = (double)(sumy/sumx);}  //次の処理のため角度を計算し格納しておく
	return theta;
}