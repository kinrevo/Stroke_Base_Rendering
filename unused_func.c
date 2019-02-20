#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"pgm.h"

#define t1_8pi 0.4142
#define t3_8pi 2.4142
#define maxValue 0.35
#define minValue 0.25



//���ω��t�B���^
PGM *average_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
				for(k=-1; k<=1; k++) { for(l=-1; l<=1; l++)  sum += pgm->data[i+k][j+l];}
				img->data[i][j] = sum/9;
				sum = 0;
			}
		}
	}
	
	
	return img;
}



//�^����ꂽ臒l�ŃC���[�W���l������֐�
PGM *binarization_filter(PGM *pgm, int boundary) 
{
	int i,j;
	PGM *nimg = copy_pgm(pgm);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			//�s�N�Z����臒l�ŕ����ē�l������
			if(pgm->data[i][j] < boundary) {
				nimg->data[i][j] = 0;				
			}else {
				nimg->data[i][j] = pgm->bright;
			}
		}
	}
	
	
	return nimg;
	
}



//�������]�t�B���^
PGM *inversion_filter(PGM *pgm) 
{
	int i,j;
	PGM *nimg = copy_pgm(pgm);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			//�s�N�Z�����P�x�̍ő�l���甽�]����
			nimg->data[i][j] = pgm->bright - pgm->data[i][j];	
		}
	}
	
	
	return nimg;
	
}



//���d���σt�B���^
PGM *weightaverage_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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



//�c�����t�B���^
PGM *verticalsmooth_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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




//�������t�B���^
PGM *horizonsmooth_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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




//prewitt�t�B���^
PGM *prewitt_filter(PGM *pgm)
{
	int i,j,k,l,sum=0,sumx=0,sumy=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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




//Laplacian�t�B���^
PGM *Laplacian_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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



//�K�E�V�A���t�B���^���a�P
PGM *gaucian_filter(PGM *pgm)
{
	int i,j,k,l,sum=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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


//sobel�t�B���^
PGM *sobel_filter(PGM *pgm)
{
	int i,j,k,l,sum=0,sumx=0,sumy=0;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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



//�^����ꂽ��̃C���[�W�̍��������iin2-in1�j
PGM *DifferentImg(PGM *in1, PGM *in2) {
	
	int i,j,min=255;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
	PGM *nimg = copy_pgm(in1);
	
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//��������������l������
			nimg->data[i][j] = in2->data[i][j] - in1->data[i][j];
			//if(nimg->data[i][j]<0) nimg->data[i][j]=0;
			//nimg->data[i][j] = (nimg->data[i][j]>0 ?nimg->bright:0);
			if(nimg->data[i][j] < min) min=nimg->data[i][j];
		}
	}
	
	min = (-1)*min;
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//��グ
			nimg->data[i][j] += min; 
			nimg->data[i][j] = (nimg->data[i][j] * nimg->bright)/(min+nimg->bright);
		}
	}
	
	return nimg;
}



//�^����ꂽ��̃C���[�W���r����iin2-in1�j
PGM *CompareImg(PGM *in1, PGM *in2) {
	
	int i,j;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
	PGM *nimg = copy_pgm(in1);
	
	for(i=0; i<in1->height; i++) {
		for(j=0; j<in1->width; j++) {
			//��������������l������
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
	

//�^����ꂽ�C���[�W�����v����n��]����
PGM *RotateImg(PGM *in, int n) {
	
	int i,j,k;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
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



//�L���j�[�G�b�W���o��
PGM *cannyedge_detector(PGM *pgm)
{
	int i,j,k,l,flag=0,sum=0,sumx=0,sumy=0,min=pgm->bright,max=0;
	//�V�����摜�f�[�^�����K�E�V�A���t�B���^���|�������̂�����(nimg1)
	PGM *nimg1 = gaucian_filter(pgm);
	printf("step1 finish\n");
	//�摜�f�[�^�����sobel�t�B���^��������(nimg2)
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
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
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
				if(sumx==0) {theta[i][j] = (double)(sumy*10000);} //���ꂪ�[���ɂȂ�Ȃ��悤�ɂ���
				else { theta[i][j] = (double)(sumy/sumx);}  //���̏����̂��ߊp�x���v�Z���i�[���Ă���
				
				if(sum>max) {max=sum;}  //�ő召�l���X�P�[�����O�̂��ߕۑ����Ă���
				if(sum<min) {min=sum;}
				
				sum = sumx = sumy = 0;
			}
		}
	}
	
	//�X�P�[�����O�������s��
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			nimg2->data[i][j] += (-1)*min;
			nimg2->data[i][j] = (nimg2->data[i][j]*pgm->bright)/(max-min);
			if(nimg2->data[i][j] > 255) printf("error\n");
		}
	}
	printf("max:%d min:%d\n",max,min);
	
	printf("step2 finish\n");
	
	
	//��ɑ�l�}���������s��(nimg3)
	PGM *nimg3 = copy_pgm(nimg2);
	
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) { 
				nimg3->data[i][j] = nimg2->data[i][j];
			}else {
				//���ڂ��Ă���s�N�Z���̃G�b�W�����̉��������̊p�x�����@(x,y)=(k,l)
				//-1/8pi < 1/8pi
				if( (-1)*t1_8pi <= theta[i][j] && theta[i][j] <= t1_8pi ) { k=1; l=0; }
				//1/8pi < 3/8pi
				else if( t1_8pi <= theta[i][j] && theta[i][j] <= t3_8pi ) { k=1; l=1; }
				//3/8pi < 5/8pi
				else if( t3_8pi <= theta[i][j] && theta[i][j] <= (-1)*t3_8pi ) { k=0; l=1; }
				//5/8pi < 7/8pi
				else  {	k=-1; l=1; }
				
				//���ڂ��Ă���s�N�Z�����G�b�W���������̗אډ�f�Ɣ�r���ő�łȂ����0�Ƃ���
				if( nimg2->data[i][j] < nimg2->data[i-l][j-k] || nimg2->data[i][j] < nimg2->data[i+l][j+k] ) {
					nimg3->data[i][j] =0;
				}
				else { nimg3->data[i][j] = nimg2->data[i][j]; }
					
				
			}
		}
	}
	
	
	
	printf("step3 finish\n");
	//�q�X�e���V�X臏������s��(nimg4)
	PGM *nimg4 = copy_pgm(nimg3);
	
	//�n�C�p�X�𒴂���l���ő�A���[�p�X���Ⴂ�l���ŏ��Ƃ���
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(nimg3->data[i][j] > maxValue*nimg3->bright){ nimg4->data[i][j]=nimg3->bright; }
			else if(nimg3->data[i][j] < minValue*nimg3->bright){ nimg4->data[i][j]=0; }
			else { nimg4->data[i][j] = nimg3->data[i][j]; }
			
		}
	}
	
	
	//2��臒l�̊Ԃ̒l�̓G�b�W�Ɍq�����Ă���΃G�b�W�Ƃ���
	int wflag=1, flag2=0;
	//�G�b�W�̐ڑ���S�Ċm�F����܂ŉ�f�̑S�T���𑱂���
	while(wflag){
		wflag=0;
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->width; j++) {
			if(i==0 || i==pgm->height-1 || j==0 || j==pgm->width-1) { 
			}else if(nimg4->data[i][j]==0 || nimg4->data[i][j]==pgm->bright){
			}else {				
				//���ڂ��Ă���s�N�Z���̎���ɃG�b�W�����邩�m�F����
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
	
	//�O���[��f��S�čŒ�P�x��
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






//SIFT�ϊ�
PGM *SIFT(PGM *pgm, int S) {
	
	int i,j,x,y,s,ds,dx,dy,temp;
	double sigma = 1.6;
	PGM *Gimg[S+3];  //�������摜������z��
	PGM *DoG[S+2];  //DoG
	int is_max,is_min;
	double size = sigma;
	
	//������k���Z�o
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
	int **keypt[S];  //�L�[�|�C���g
	
	for(s=0; s<S; s++){
		keypt[s] = create_ally(pgm->width, pgm->height);
	}
	*/
	
	//keypt�̒l�͕������摜�̎w��q�i�Ⴆ��2�Ȃ�X�P�[����k^2�Ђł���j
	double **keypt = create_dally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { keypt[i][j]=0; }
	}
	
	int **flag = create_ally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { flag[i][j]=0; }
	}
	
	int kpcount=0;

	//�ɒl�̌��o
	for(s=1; s<=S; s++){
		//size = size*k;
		p("s",s);
		for(x=0; x<pgm->height; x++){
        	for(y=0; y<pgm->width; y++){
        		is_max = is_min = 1;
        		
        		//���ɃL�[�|�C���g�ƌ��o���Ă��邩�v�Z�ł��Ȃ��̈悾��0
        		if(x==0 || x==pgm->height-1 || y==0 || y==pgm->width-1) {}
        		else if(keypt[x][y]>0) {}
        		//�ߖT�Ɣ�r
        		else {
        			for(ds=-1; ds<=1; ds++) { 
						for(dx=-1; dx<=1; dx++)  {
							for(dy=-1; dy<=1; dy++)  {
								if(ds==0 && dx==0 && dy==0){} //���S�͖���
								else { 
									if(DoG[s+ds]->data[x+dx][y+dy] >= DoG[s]->data[x][y]) {
										is_max=0;  //�ɑ�l�łȂ��Ɣ���
									}
									if(DoG[s+ds]->data[x+dx][y+dy] <= DoG[s]->data[x][y]) {
										is_min=0;  //�ɏ��l�łȂ��Ɣ���
									}
								}
							}
						}
					}
        			if(is_max==0 && is_min==0) keypt[x][y]=0;  //�ɑ�ł��ɏ��ł��Ȃ��Ȃ�0�Ƃ���
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
	//�L�[�|�C���g�폜
	//-----------------------
	double Dxx, Dyy, Dxy, Tr, Det;
	double Dx, Dy, Ds, Dss, Dxs, Dys, detmD, Dpow;
	double mD[3][3], iD[3][3], X[3], xD[3];
	int sm2, sp2;
	double TH_POW = 0.05*pgm->bright; //�������l�i�ő�P�x��3���j

	for(s=1; s<=S; s++){
		for(x=2; x<pgm->height-2; x++){
        	for(y=2; y<pgm->width-2; y++){
        		//��Ȑ��ɂ��L�[�|�C���g�폜
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
        		
        		//�T�u�s�N�Z������
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
        			
        			//�t�s��v�Z(mD�̋t�s���iD��)
        			detmD = mD[0][0]*mD[1][1]*mD[2][2] + mD[1][0]*mD[2][1]*mD[0][2] + mD[2][0]*mD[0][1]*mD[1][2]
        					- mD[0][0]*mD[2][1]*mD[1][2] - mD[2][0]*mD[1][1]*mD[0][2] - mD[1][0]*mD[0][1]*mD[2][2];
        			if(detmD==0) continue;
        			
        			iD[0][0]=mD[1][1]*mD[2][2]-mD[1][2]*mD[2][1]; iD[0][1]=mD[0][2]*mD[2][1]-mD[0][1]*mD[2][2]; iD[0][2]=mD[0][1]*mD[1][2]-mD[0][2]*mD[1][1];
        			iD[1][0]=mD[1][2]*mD[2][0]-mD[1][0]*mD[2][2]; iD[1][1]=mD[0][0]*mD[2][2]-mD[0][2]*mD[2][0]; iD[1][2]=mD[0][2]*mD[1][0]-mD[0][0]*mD[1][2];
        			iD[2][0]=mD[1][0]*mD[2][1]-mD[1][1]*mD[2][0]; iD[2][1]=mD[0][1]*mD[2][0]-mD[0][0]*mD[2][1]; iD[2][2]=mD[0][0]*mD[1][1]-mD[0][1]*mD[1][0];
        			
        			for(i=0; i<3; i++){
        				for(j=0; j<3; j++){ iD[i][j] = iD[i][j]/detmD;   }
        			}
        			
        			//�T�u�s�N�Z���ʒu(�s��̐�)
        			X[0]=iD[0][0]*xD[0]+iD[0][1]*xD[1]+iD[0][2]*xD[2];
        			X[1]=iD[1][0]*xD[0]+iD[1][1]*xD[1]+iD[1][2]*xD[2];
        			X[2]=iD[2][0]*xD[0]+iD[2][1]*xD[1]+iD[2][2]*xD[2];
        			//�T�u�s�N�Z���ʒu�ł̏o��
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
	//�I���G���e�[�V�����̎Z�o
	//----------------------

	//�L�[�|�C���g���Ƃ̃I���G���e�[�V������ۑ�����z��
	int **orient = create_ally(pgm->width, pgm->height);
	for(i=0; i<pgm->height; i++) {
		for(j=0; j<pgm->height; j++) { orient[i][j]=0; }
	}
	
	double fu, fv;
	double histogram[36]={};
	int w = (int)( ceil(3.0*sigma+0.5)*2-1 ); //�K�E�X�̒l�����0�ɂȂ鋗���i���a�j
	int c=(w-1)/2;  //�K�E�X�t�B���^�̔��a
	double **filter; //�K�E�X�t�B���^���i�[����z��̃|�C���^
	double Fpow;
	int Farg, peak, im1,ip1;
	
	
	//�L�[�|�C���g���畽�����摜�ƃK�E�V�A�����a�����o��
	for(s=1; s<=S; s++){
		sigma=1.6; //�Ђ�������
		sigma = sigma*pow(k, s);
		w = (int)( ceil(3.0*sigma+0.5)*2-1 );  //�K�E�X�̒l�����0�ɂȂ鋗���i���a�j
		c=(w-1)/2;  //�K�E�X�t�B���^�̔��a
				
		//�K�E�X�t�B���^�̗̈���m��
		filter = create_dally(w, w);
		for(i=0; i<w; i++) {
			for(j=0; j<w; j++) { filter[i][j]=0; }
		}
		
		//�K�E�X�t�B���^���\��
    	for(i=0;i<w;i++){
    	    for(j=0;j<w;j++){
     		   	filter[i][j] = gause_func(i-c, j-c, sigma);
       		}
  		}
		
		for(x=0; x<pgm->height; x++){
			for(y=0; y<pgm->width; y++){
				//�T�����Ă���Ђ̃K�E�X�t�B���^�������Ă���L�[�|�C���g�Ȃ�q�X�g�O�������v�Z
				if(keypt[x][y]==s){
					//histogram��������
    				for(i=0;i<36;i++){ histogram[i] = 0;	}
					
					//���ڂ��Ă���s�N�Z���̃K�E�X���̃s�N�Z�����̌��z�𒲂ׂ�
					for(i=-c; i<=c; i++) { 
						for(j=-c; j<=c; j++)  {
							//�v�Z�o���Ȃ��̈�̓X�L�b�v����
							if( (x+i-1)<0 || (x+i+1)>pgm->height-1 || (y+j-1)<0 || (y+j+1)>pgm->width-1) {}
							else {
								fu = Gimg[s]->data[x+i+1][y+j] - Gimg[s]->data[x+i-1][y+j];
								fv = Gimg[s]->data[x+i][y+j+1] - Gimg[s]->data[x+i][y+j-1];
								
								Fpow = sqrt(fu*fu+fv*fv);
								Farg = (atan2(fv,fu)/PI+1)*18;
								//�K�E�X���ŋ��x���d�ݕt�������ꂼ��̕����ɉ��Z
								histogram[Farg] += Fpow*filter[i+c][j+c];
							}
						}
					}
					
					//�q�X�g�O�����̒��ōł��傫���l��T��
					peak=0;
					for(i=0;i<36;i++){
						if(peak < histogram[i]) peak=i;
					}
					orient[x][y] = peak;
					//�ő�l��80���ȏ�̌��z�ŋɒl�Ȃ�I���G���e�[�V�����Ɋ��蓖�Ă�
					for(i=0;i<36;i++){
						if(histogram[i] > 0.8*histogram[peak] && i!=peak){
							im1 = (i-1<0) ? 35:i-1;
							ip1 = (i+1>35) ? 0:i+1;
							if(histogram[i]>histogram[im1] && histogram[i]>histogram[ip1]){
								//orient�z��ɂ͓񌅖��Ɍ��z�����̏�񂪕ۑ������
								//�iorient=1209�Ȃ�12��9���I���G���e�[�V�����j
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
	for(x=0; x<kimg->height; x++){  //�摜�f�[�^���R�s�[
        for(y=0; y<kimg->width; y++){
        	kimg->data[x][y] = pgm->data[x][y]; 
        }
	}
	/*
	for(x=0; x<kimg->height; x++){  //�L�[�|�C���g�_�\��
        for(y=0; y<kimg->width; y++){
        	if(keypt[x][y]>0){ kimg->data[x][y] = kimg->bright;	}
        }
	}
	*/
	
	double fai;
	
	for(x=0; x<kimg->height; x++){  //�L�[�|�C���g�~�\��
        for(y=0; y<kimg->width; y++){
        	if(keypt[x][y]>0){ 
        		size = keypt[x][y]*keypt[x][y];
        		for(i=0; i<size*2; i++){
        			fai = (i/size)*2*PI;
        			dx = keypt[x][y]*cos(fai);
        			dy = keypt[x][y]*sin(fai);
        			//�摜�̈���ɉ~��`��
        			if( (x+dx)>0 && (x+dx)<(pgm->height-1) && (y+dy)>0 && (y+dy)<(pgm->width-1) ){
        				kimg->data[x+dx][y+dy] = kimg->bright;
        			}
        		}
        		//�I���G���e�[�V������`��
        		size= keypt[x][y];
        		temp=orient[x][y];
        		while(temp!=0){
        			fai=((double)(temp%100)/18 - 1)*PI;
        			
        			for(i=0; i<size; i++){
        				dx = i*cos(fai);
        				dy = i*sin(fai);
        				//�~�̒��ɐ���`��
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
	for(x=0; x<kimg->height; x++){  //�L�[�|�C���g�_�摜�쐬
        for(y=0; y<kimg->width; y++){
        	kimg->data[x][y] = keypt[x][y];
        }
	}
	/**/
	return kimg;
}




//�o�C���e�����t�B���^��p���ė^����ꂽ�C���[�W���G�撲�ɂ���
PGM *Illustification(PGM *in, double sigma, int n) {
	
	int i,j;
	//�V�����摜�f�[�^���R�s�[���č���Ă���
	PGM *nimg = bilateral_filter(in, sigma, n);
	PGM *nimg2 = copy_pgm(in);
	/*for(i=0; i<in->height; i++) {
		for(j=0; j<in->width; j++) {
			//���ȏ��摜�̐���
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



//sobel��K�������͈͂̍��v�̌��z��Ԃ�
double sobel_slope(PGM* in_img, int x, int y, int t) {
	double theta;
	int i,j,k,l, sumx=0, sumy=0;
	
	//sobel�t�B���^
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
	
	//sobel�t�B���^�ɂ��X�g���[�N����������
	for(i=x-t; i<=x+t; i++) {
		for(j=y-t; j<=y+t; j++) {
			if(i<0 || i>in_img->width-1 || j<0 || j>in_img->height-1) {
			}else { 
				//���ڂ��Ă�����͈�s�N�Z���̒l�����v�����ς���
				for(k=-1; k<=1; k++) { 
					for(l=-1; l<=1; l++)  {
						sumx += in_img->data[i+k][j+l] * filterx[k+1][l+1];
						sumy += in_img->data[i+k][j+l] * filtery[k+1][l+1];
					}
				}
			}
		}
	}
					
	//���z�̕΂肩��X�g���[�N�������v�Z
	if(sumx==0) {theta = (double)(sumy*10000);} //����[���̗�O����
	else { theta = (double)(sumy/sumx);}  //���̏����̂��ߊp�x���v�Z���i�[���Ă���
	return theta;
}