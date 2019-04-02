#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <png.h>
#include <jpeglib.h>
#include "sbr.h"
#include "sbr_opt.h"
//#define _CRTDBG_MAP_ALLOC #include <stdlib.h> #include <crtdbg.h>  
//#define malloc _malloc_dbg
//#define _DEBUG



#define minWindow 0.25
#define simga0 1.6
#define BORDER 12.1


//�^����ꂽ�p�X�̒��̃t�H���_�̒��̉摜���ׂĂɑ΂���C_illust_brush�����s
void multi_dir_CIB(char* in_path){
	int i;
	int dir_num, file_num;
	char** file_list;
	
	char** dir_list = make_file_list(&dir_num, in_path, Search_DIRECTRY); //�t�H���_���X�g���擾
	
	for(i=0; i<dir_num; i++) {
		//i�Ԗڃt�H���_�̃t�@�C�����X�g���擾
		file_list = make_file_list(&file_num, dir_list[i], Search_FILE); 
		//�t�@�C�����X�g�̉摜���ꂼ��ɑ΂���C_illust_brush�����s
		multi_file_CIB(file_list, file_num);
		
		for(i=0; i < file_num; i++) { free(file_list[i]);}
		free(file_list);
	}
	
	for(i=0; i < dir_num; i++) { free(dir_list[i]);}
	free(dir_list);
}
	
//�^����ꂽ�p�X���X�g�̉摜���ׂĂɑ΂���C_illust_brush�����s
void multi_file_CIB(char** list, int list_num){
	int i;
	PPM *in_data, *trans_data;
	image_t *in_img;
	char out_dir_name[256];
	
	for(i=0; i<list_num; i++) {
		strcpy(out_dir_name, list[i]);	
		strtok(out_dir_name,".");		//�g���q�͎��
		printf("out_dir_name:%s\n",out_dir_name);
		printf("list[i]:%s\n",list[i]);
		
		if((in_img = read_jpeg_file(list[i])) == NULL){  // <<<<<====�Ώۃt�@�C���̊g���q�ɂ���ĕύX���K�v
			printf("PNG_READ_ERROR\n");		exit(1);
		}
		dump_image_info(in_img);	//PNG���o��
		in_data = image_to_PPM(in_img);		//�����₷���f�[�^�\���ɕϊ�
		free_image(in_img);
		
		trans_data = c_Illust_brush(in_data, list[i]);
		FreePPM(in_data);
		FreePPM(trans_data);
		
	}
	
}


//�^����ꂽ�p�X�̃f�B���N�g�����ɂ���t�@�C���������̓f�B���N�g�������X�g�����z��ŕԂ�
//�����ɏ]���Ăǂ��炩���������ł��Ȃ�
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
        stat(path, &fi);	//path�̃t�@�C���̏����擾
    	
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
        stat(path, &fi);	//path�̃t�@�C���̏����擾
    	
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



//�����_��Ԃ��֐�
int inner_point(int p1,int p2,float m,float n){
	int r = ( (n*p1+m*p2)/(m+n) + 0.5);
	return r;
}


//�x�W�G�Ȑ��̕Ѝ��W���v�Z����֐�
int BezierCurve(int p1, int p2, int p3, int p4, float t){
	return (1-t)*(1-t)*(1-t)*p1 + 3*t*(1-t)*(1-t)*p2 + 3*t*t*(1-t)*p3 + t*t*t*p4;
}

//�x�W�G�Ȑ��̍��W���v�Z����֐�
Point BezierCurve_P(Point p1, Point p2, Point p3, Point p4, float t){
	Point r={BezierCurve(p1.x, p2.x, p3.x, p4.x, t), BezierCurve(p1.y, p2.y, p3.y, p4.y, t)};
	return r;
}


//�^����ꂽ�摜�Ƀx�W�G�Ȑ���`�悷��֐��i�O�_�j
void Paint_Bezier(Point p1, Point s, Point p4, PGM* in_img, int thick, int bright) {
	int i,r;
	int partition = abs(p1.x-p4.x) + abs(p1.y-p4.y); //�����̕�����
	p("parti",partition);
	float t;
	Point temp, p2, p3;
	
	p2.x = ( ((s.x-p4.x)/5 + s.x)>0 ? (s.x-p4.x)/5 + s.x : 0);
	if(p2.x >= in_img->width) p2.x=in_img->width-1;
	p2.y = ( ((s.y-p4.y)/5 + s.y)>0 ? (s.y-p4.y)/5 + s.y : 0);
	if(p2.y >= in_img->height) p2.y=in_img->height-1;
	p3.x = ( ((s.x-p1.x)/5 + s.x)>0 ? (s.x-p1.x)/5 + s.x : 0);
	if(p3.x >= in_img->width) p3.x=in_img->width-1;
	p3.y = ( ((s.y-p1.y)/5 + s.y)>0 ? (s.y-p1.y)/5 + s.y : 0);
	if(p3.y >= in_img->height) p3.y=in_img->height-1;
	

	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick*sin((10*t<PI/2 ? 10*t:PI/2));
		temp = BezierCurve_P(p1, p2, p3, p4, t);
		Circle_fill(temp.x, temp.y, r, in_img, bright, 0.03);
		//light_dot(temp.x, temp.y, in_img, 255);
	}
}



//�^����ꂽ�摜�Ƀx�W�G�Ȑ���`�悷��֐��ipnum�_�j
void Paint_Bezier_ex(Point p[], int pnum, PGM* in_img, int thick, int bright, double ratio) {
	//��_�܂ł����^�����Ȃ���Β���������
	if(pnum==2){
		Paint_line(p[0], p[1], in_img, thick, bright, ratio);
		return;
	}
	Point p0={2*p[0].x-p[1].x, 2*p[0].y-p[1].y}, 	//���[�̈�O�̐���_��K���Ɍ��߂�
		p_np1={2*p[pnum-1].x-p[pnum-2].x, 2*p[pnum-1].y-p[pnum-2].y};
	//light_dot(p0.x, p0.y, in_img, 0);
	//light_dot(p_np1.x, p_np1.y, in_img, 0);
	//printf("%d,%d\n",p_np1.x, p_np1.y);
	Point temp;  //�`�����ʂ�Ȃ�����_
	int i,j,r;
	double t;
	int partition = abs(p[0].x-p[1].x)+abs(p[0].y-p[1].y); //�����̕�����
	
	Point bp0 = {(p[1].x-p0.x)/6.0 + p[0].x, (p[1].y-p0.y)/6.0 + p[0].y}
		,bp1 = {(p[0].x-p[2].x)/6.0 + p[1].x, (p[0].y-p[2].y)/6.0 + p[1].y};
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick*sin((10*t<PI/2 ? 10*t:PI/2));
		temp = BezierCurve_P(p[0], bp0, bp1, p[1], t);
		Circle_fill(temp.x, temp.y, r, in_img, bright, 0.03);
		//light_dot(temp.x, temp.y, in_img, 0);
	}
	//light_dot(p[0].x, p[0].y, in_img, 0);
	
	for(i=1; i<pnum-2; i++){
		bp0.x = (p[i+1].x-p[i-1].x)/6.0 + p[i].x;	bp0.y = (p[i+1].y-p[i-1].y)/6.0 + p[i].y;
		bp1.x = (p[i].x-p[i+2].x)/6.0 + p[i+1].x;	bp1.y = (p[i].y-p[i+2].y)/6.0 + p[i+1].y;
		partition = abs(p[i].x-p[i+1].x)+abs(p[i].y-p[i+1].y);
		for(j=0; j <= partition; j++) {
			t = (double)j/partition;
			r = thick;
			temp = BezierCurve_P(p[i], bp0, bp1, p[i+1], t);
			Circle_fill(temp.x, temp.y, r, in_img, bright, 0.03);
			//light_dot(temp.x, temp.y, in_img, 0);
		}
		//light_dot(p[i].x, p[i].y, in_img, 0);
	}
	

	bp0.x = (p[pnum-1].x-p[pnum-3].x)/6.0 + p[pnum-2].x;
	bp0.y = (p[pnum-1].y-p[pnum-3].y)/6.0 + p[pnum-2].y;
	bp1.x = (p[pnum-2].x-p_np1.x)/6.0 + p[pnum-1].x;
	bp1.y = (p[pnum-2].y-p_np1.y)/6.0 + p[pnum-1].y;
	partition = abs(p[i].x-p[i+1].x)+abs(p[i].y-p[i+1].y);
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick;
		temp = BezierCurve_P(p[pnum-2], bp0, bp1, p[pnum-1], t);
		Circle_fill(temp.x, temp.y, r, in_img, bright, 0.03);
		//light_dot(temp.x, temp.y, in_img, 0);
	}
	//light_dot(p[pnum-2].x, p[pnum-2].y, in_img, 0);
	//light_dot(p[pnum-1].x, p[pnum-1].y, in_img, 0);
}



//�^����ꂽ�摜�ɒ�����`�悷��֐�
void Paint_line(Point p1, Point p4, PGM* in_img, int thick, int bright, double ratio) {
	int i,r;
	int partition = abs(p1.x-p4.x) + abs(p1.y-p4.y); //�����̕�����
	p("parti",partition);
	float t;
	Point temp;
	
	p("b",bright);
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick*sin((10*t<PI/2 ? 10*t:PI/2));
		temp.x = p1.x+(p4.x-p1.x)*t;
		temp.y = p1.y+(p4.y-p1.y)*t;
		Circle_fill(temp.x, temp.y, r, in_img, bright, ratio);
	}
}


//�^����ꂽ�摜�̍��W�𒆐S�Ƃ���~��h��ׂ��֐�
void Circle_fill(int xc, int yc, int r, PGM* in_img, int bright, double ratio) {
	int x,y;
	for(x=xc-r; x <= xc+r; x++) {
		for(y=yc-r; y <= yc+r; y++) {
			if(x<0 || x>in_img->width-1 || y<0 || y>in_img->height-1) {}
			else if( (x-xc)*(x-xc)+(y-yc)*(y-yc) <= r*r ){
				in_img->data[x][y] = (1-ratio)*in_img->data[x][y] + ratio*bright;
			}
		}
	}
}


//�^����ꂽ���W��ڗ�������悤�ɑŕe����
void light_dot(int x, int y, PGM* in_img, int bright) {
	if(x>(in_img->width-1) || x<1 || y>(in_img->height-1) || y<1) return;
	in_img->data[x][y]=bright;
	in_img->data[x+1][y]=bright;
	in_img->data[x][y+1]=bright;
	in_img->data[x-1][y]=bright;
	in_img->data[x][y-1]=bright;
}


//������p�xtheta�����ɋ���t�̓_�����߂�
Point calcu_point(PGM *in, Point a, int t, double theta){
	Point b;
	b.x = t*cos(theta)+a.x;  
	if(b.x<0) b.x=0;  if(b.x>=in->width) b.x=in->width-1; 
	b.y = t*sin(theta)+a.y;
	if(b.y<0) b.y=0; if(b.y>=in->height) b.y=in->height-1;
	
	return b;
}


//�^����ꂽ�_����̐F���`��F�Ƃ����ق����߂�
double diffsum_clr(PGM* cmpr, PGM* nimg, Point p, int t, int bright){
	//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
	int i,j, offscrn_count=0;
	double sum=0;
	for(i=p.x-t; i<=p.x+t; i++) {
		for(j=p.y-t; j<=p.y+t; j++) {
			if(i<0 || i>cmpr->width-1 || j<0 || j>cmpr->height-1) {
				offscrn_count++;
			}else { 
				//���ڂ��Ă������t�s�N�Z���̒l�����v�����ς���
				sum += abs(cmpr->data[i][j]-nimg->data[i][j]) - abs(cmpr->data[i][j]-bright);
				//sum += abs(cmpr->data[i][j]-bright) -30
			}
		}
	}
	sum = sum/((2*t+1)*(2*t+1)-offscrn_count);  //�_������̍��ٕ���
	return sum;
}


//���ڃs�N�Z�����肔�̕��ς����`��F���v�Z
int calcu_color(int** data, int width, int height, int x, int y, int t){
	int i,j,bright=0, offscrn_count=0;
	for(i=x-t; i<=x+t; i++) {
		for(j=y-t; j<=y+t; j++) {
			if(i<0 || i>width-1 || j<0 || j>height-1) {
				offscrn_count++;
			}else { 
				//���ڂ��Ă������t�s�N�Z���̒l�����v�����ς���
				bright += data[i][j];
			}
		}
	}
	bright = bright/((2*t+1)*(2*t+1)-offscrn_count);  //�E�B���h�E�T�C�Y�Ŋ���(t�͔��a)
	return bright;
}


//���ڃs�N�Z�����ӂ��Ńo�C���e�����t�B���^�ɂ����d���ς����`��F���v�Z
int calcu_color_bi(int** data, int width, int height, int x, int y, int t, double sigma, double** gause_filter){
	int i,j, bright=0, br_diff;
	double sum=0, co_sum=0, sigma_2=2*sigma*sigma, coefficient;
	for(i=x-t; i<=x+t; i++) {
		if(i<0 || i>width-1) continue;
		for(j=y-t; j<=y+t; j++) {
			if(j<0 || j>height-1) continue;
			//���ڂ��Ă������t�s�N�Z���̒l���K�E�X����
			br_diff = (data[x][y] - data[i][j]);
			coefficient = exp(-(br_diff*br_diff)/sigma_2);
			sum += data[i][j] * coefficient;
			co_sum += coefficient;
		}
	}
	//pd("co",co_sum);
	bright = sum/co_sum + 0.5;  
	return bright;
}


//�@�^����ꂽ�͈͂ɂ�������z�̃q�X�g�O�������Ƃ�
double calcu_histogram(PGM* cmpr, double **sobel_abs, double **sobel_angle, int partition, 
					double **gauce_filter, int x, int y, int t, int *histogram_direct, int *break_flag)
{
	int i,j,Farg, im1,ip1;//,im2,ip2;
	int peak=partition/2+0.5;		//�l�̌ܓ�
	double histogram[partition];
	double theta;
	for(i=0; i<partition; i++) {histogram[i]=0;}
	
	//sobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
	for(i=x-t; i<=x+t; i++) {
		for(j=y-t; j<=y+t; j++) {
			if(i<1 || i>cmpr->width-2 || j<1 || j>cmpr->height-2) continue;
			Farg = (sobel_angle[i][j]/PI)*partition;
			histogram[Farg] += sobel_abs[i][j]*gauce_filter[i-(x-t)][j-(y-t)];
		}
	}
	
	//�q�X�g�O�����̒��ōł��傫���l��T��
	for(i=0;i<partition;i++){
		if(histogram[peak] < histogram[i]) peak=i;
	}
	
	//�ő�l��80���ȏ�̌��z�����ɑ��݂����continue
	for(i=0;i<partition;i++){
		im1 = (i-1<0) ? partition-1:i-1;
		ip1 = (i+1>partition-1) ? 0:i+1;
		//im2 = (im1-1<0) ? partition-1:im1-1;
		//ip2 = (ip1+1>partition-1) ? 0:ip1+1;
		if(histogram[i] > 0.8*histogram[peak] && i!=peak && im1!=peak && ip1!=peak && im1!=peak && ip1!=peak){
			*break_flag=1;
			break;
		}
	}
	//printf("peak:%d  ",peak);
	histogram_direct[peak]++;
	theta = ((double)peak/partition)*PI+((PI/partition)/2);
	
	return theta;
}




//�G�b�W�̓���g�񂾓_��T�����A���ӂɓ���g�񂾃G�b�W�����_�̃}�b�v��Ԃ�
PGM *expand_Edge(PGM *canny, int thick_min){
	PGM *nimg = copy_pgm(canny);
	int x,y ,t=thick_min;
	format_ally(nimg->data, canny->width, canny->height, 0);
	
	for(y=0; y<canny->height; y++) {  
		for(x=0; x<canny->width; x++) {
			//�G�b�W�𔭌������Ȃ�
			if(canny->data[x][y] == 255){
				//���̎��͂���h��ׂ�
				Circle_fill(x, y, t, nimg, 255, 1);
			}
		}
	}
		
	return nimg;
}



//�G�b�W�̓���g�񂾓_��T�����A���ӂɓ���g�񂾃G�b�W�����_�̃}�b�v��Ԃ�
PGM *calcu_EdgeMap(PGM *canny, int thick_min, double **sobel_angle){
	PGM *nimg = copy_pgm(canny);
	PGM *map = copy_pgm(canny);
	Point p;
	int x,y,xc,yc ,t=thick_min, flag1, flag2a=0,flag2b=0;
	int px,py;
	int flag1_c=0,flag2_c=0;
	format_ally(map->data, canny->width, canny->height, 0);
	format_ally(nimg->data, canny->width, canny->height, 0);
	
	for(y=0; y<canny->height; y++) {  //�E�B���h�E�̑傫���ɍ��킹��
		for(x=0; x<canny->width; x++) {  //�E�B���h�E�����炷������ς��Ƃ�
			flag1=flag2a=flag2b=0;
			for(xc=-t; xc<=t; xc++) {
				if((x+xc)<0 || (x+xc)>canny->width-1) continue;
				for(yc=-t; yc<=t; yc++) {
					if( xc*xc+yc*yc > t*t ) continue;
					if((y+yc)<0 || (y+yc)>canny->height-1) continue;
					//�G�b�W�Ȃ玟�̃X�e�b�v��
					if(canny->data[x+xc][y+yc] != canny->bright) continue;		
					//�������G�b�W�̓_����G�b�W��H��G�b�W�̌`��T��
					p.x=x+xc, p.y=y+yc;
					//�T���ς݂łȂ���ڂ̃G�b�W�𔭌������Ȃ�
					if(flag2a==1 && map->data[x+xc][y+yc]!=255) {
						flag2b=1;	//��̃G�b�W�𔭌�
						flag2_c++;printf("flag2\n");
						break;
					}
					px=p.x; py=p.y;
					map->data[px][py]=255;
					if(flag2a!=1){
						//�������G�b�W�̏ꏊ���L�^
						pp(p); //getchar();
						//�������G�b�W�̉ӏ�����Ȃ����Ă���G�b�W��S�ĒT��
						edge_detect(canny, map, sobel_angle, x, y, t, p, &flag1, 0);
						printf("detect_done\n");
						
					}
					if(flag1) {flag1_c++; //break;
					}
					flag2a=1; //��̃G�b�W�𔭌�
				}
				if(flag2b==1) {
				//if(flag1==1 || flag2b==1) {
					nimg->data[x][y]=255;
					break;
				}
			}
			format_ally(map->data, canny->width, canny->height, 0);
		}
	}
	p("F1",flag1_c);
	p("F2",flag2_c);
	//getchar();
	FreePGM(map); 
	return nimg;
}

// �אڂ���G�b�W�̐�����
int edge_detect(PGM *canny, PGM *map, double **sobel_angle, int x, int y, int t, Point s, int *flag, int count_init){
	int i, j, ec=0, count=count_init;	//�אڂ���G�b�W�̐�
	Point p;
	if(s.x==0 || s.x==canny->width-1 || s.y==0 || s.y==canny->height-1) return 0; //��ʒ[�ɓ��B
	if(s.x<(x-t) || s.x>(x+t) || s.y<(y-t) || s.y>(y+t)) return 0;		//�E�B���h�E�[�ɓ��B
	map->data[(int)s.x][(int)s.y] = 255;
	
	//�^����ꂽ�G�b�W�̎��ӂ̃G�b�W��T��
	for(i=-1; i<=1; i++) {
		if((s.x+i)<0 || (s.x+i)>canny->width-1) continue;
		for(j=-1; j<=1; j++) {
			if((s.y+j)<0 || (s.y+j)>canny->height-1) continue;
			if(canny->data[(int)s.x+i][(int)s.y+j] == canny->bright){
				if(map->data[(int)s.x+i][(int)s.y+j]!=255)	{
					ec=1;	//��ȏ�G�b�W���Ȃ����Ă���
					map->data[(int)s.x+i][(int)s.y+j]=255;
					p.x=s.x+i;	p.y=s.y+j;	
					count += edge_detect(canny, map, sobel_angle, x, y, t, p, flag, 1);
				}
			}
		}
	}
	
	if(count > 2) *flag=1;		//��ȏ�A�������G�b�W���O�Ȃ����Ă���Ȃ�
	
	return ec;
}


//�^����ꂽ�{���œ_�ʒu���X�P�[�����O
Point *scaling_point(Point p[], int pnum, double canvas_scaling_ratio){
	int i;
	Point *new_p = (Point *)malloc(sizeof(Point)*(pnum));
	
	for(i=0; i<pnum; i++){
		new_p[i].x = p[i].x * canvas_scaling_ratio;
		new_p[i].y = p[i].y * canvas_scaling_ratio;
	}
	
	return new_p;
}

//�^����ꂽ�C���[�W����u���b�V���O�ɂ��G������
PPM *c_Illust_brush(PPM *in, char *filename) {
	clock_t start = clock();
	int i,j,x,y,xc,yc,t,bright,brightR,brightG,brightB,break_flag,pnum, offscrn_count;
	int P2a_count, P2b_count, inapp_count, direct_count;
	P2a_count = P2b_count = inapp_count = direct_count = 0;
	int window_diff_border = opt_window_diff_border; 	//�X�g���[�N�ʒu�T���̂������l
	int color_diff_border = opt_color_diff_border;  	//�`��F�̍��ق̂������l
	int max_stroke = opt_max_stroke;
	int min_stroke = opt_min_stroke;
	Point p[max_stroke];
//	Point *scaling_p;
	int stroke_histogram[max_stroke+1];
	for(i=0; i<max_stroke+1; i++){stroke_histogram[i]=0;}
	double ratio=opt_ratio;		//�X�g���[�N�̔Z�x
	double theta, former_theta;	
	double **gauce_filter;
	double sigma, diff_sum, sum;
	int histogram_partition=opt_histogram_partition;
	int histogram_direct[31]={};  int histogram_direct2[31]={};
	int paint_count=0, nc=0, tc=-1;
	//char input_char='i';
	int loop_cont=opt_loop_cont, lc=0, x_defo=0, y_defo=0;
	double maxValue, minValue;
//	double canvas_scaling_ratio = 1.0;
	
	//�ő召�X�g���[�N���a�i�������F��ʂ�1/10,�ő��1/10�j
	int thick_max = opt_thick_max;//(in->height < in->width ? in->height : in->width)/10;
	int thick_min = opt_thick_min;//(thick_max/15 > 3 ? thick_max/15 : 3);
	
	
	//�o�̓t�@�C�����̃T�C�Y���擾
	int namesize = strlen(filename)*2 + 16;
	char out_filename[namesize], log_filename[namesize], dir_path[namesize], vec_filename[namesize];
	char in_filename[namesize-16];
	char count_name[16];
	char log_sentence[2028] = "";
	char vec_sentence[128] = "";
	char tmp_sentence[32] = "";
	image_t *out_png;
	//�p�X������̓t�@�C�����i�g���q�܂܂Ȃ��j���擾
	char tmp[namesize];
	char* tp;
	strcpy(in_filename, filename);
	strcpy(tmp, filename);
	strtok(tmp, "/\\");
	while((tp = strtok(NULL, "/\\")) != NULL ) { 
		ps("in_filename",in_filename);
		ps("tmp",tmp);
		ps("tp",tp);
		strcpy(in_filename, tp);
	} 
	strtok(in_filename, ".");
	ps("in_filename",in_filename);
	
	//�o�͂���t�H���_�𐶐����t�H���_�ւ̃p�X���i�[�A���O�t�@�C�����쐬
	strcpy(tmp, filename);
	strcpy(dir_path, strtok(tmp, "."));
	if(mkdir(tmp, 0775)){ printf("FAIL TO CREATE DIRECTRY\n"); }
	strcat(dir_path, "/");
	strcpy(log_filename, dir_path);
	strcat(log_filename, in_filename);
	strcat(log_filename, ".log");
	
	//log�f�[�^�i�[
	strcat(log_sentence, "<");
	strcat(log_sentence, in_filename);
	strcat(log_sentence, ">\r\n");
	Add_dictionary_to_sentence(log_sentence, "width", in->width);
	Add_dictionary_to_sentence(log_sentence, "height", in->height);
	Add_dictionary_to_sentence(log_sentence, "thick_max", thick_max);
	Add_dictionary_to_sentence(log_sentence, "thick_min", thick_min);
	Add_dictionary_to_sentence(log_sentence, "max_stroke", max_stroke);
	Add_dictionary_to_sentence(log_sentence, "min_stroke", min_stroke);
	Add_dictionary_to_sentence(log_sentence, "window_diff_border", window_diff_border);
	Add_dictionary_to_sentence(log_sentence, "color_diff_border", color_diff_border);
	Add_dictionary_to_sentence(log_sentence, "ratio", (int)(ratio*100));
	Add_dictionary_to_sentence(log_sentence, "histogram_partition", histogram_partition);
	Add_dictionary_to_sentence(log_sentence, "loop_cont", loop_cont);
	strcat(log_sentence, "CompareImage : Origin\r\n");
	
	//vector�f�[�^�̃w�b�_���i�[
	// snprintf(tmp_sentence, 32, "%d", in->width);
	// strcat(vec_sentence, tmp_sentence);
	// strcat(vec_sentence, " ");
	// snprintf(tmp_sentence, 32, "%d", in->height);
	// strcat(vec_sentence, tmp_sentence);
	// strcpy(vec_filename, dir_path);
	// strcat(vec_filename, in_filename);
	// strcat(vec_filename, ".vec");
	// log_print(vec_filename, vec_sentence, "w");
	
	
	//�J���[�摜����
	PGM *gray = color_gray_conversion(in);
	PGM *nimgR = create_pgm(gray->width, gray->height, gray->bright); 
	PGM *nimgG = create_pgm(gray->width, gray->height, gray->bright); 
	PGM *nimgB = create_pgm(gray->width, gray->height, gray->bright); 
	// PGM *nimgR_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
	// PGM *nimgG_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
	// PGM *nimgB_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
	PGM *inR = create_pgm(gray->width, gray->height, gray->bright); 
	PGM *inG = create_pgm(gray->width, gray->height, gray->bright); 
	PGM *inB = create_pgm(gray->width, gray->height, gray->bright); 
	devide_ppm(in, inR, inG, inB);
	//��r�p�̃C���[�W����
	//PGM *gauss = gaussian_filter(gray, thick_min);
	//PGM *cmpr = gauss;
//	PGM *cmpr = bilateral_filter(gray, thick_min, 2);  
//	PGM *cmprR = bilateral_filter(inR, thick_min, 2);  
//	PGM *cmprG = bilateral_filter(inG, thick_min, 2);  
//	PGM *cmprB = bilateral_filter(inB, thick_min, 2);  
//	PGM *cmpr  = gaussian_filter(gray,thick_min);  
//	PGM *cmprR = gaussian_filter(inR, thick_min);  
//	PGM *cmprG = gaussian_filter(inG, thick_min);  
//	PGM *cmprB = gaussian_filter(inB, thick_min);  
	PGM *cmpr  = gray;  
	PGM *cmprR = inR;  
	PGM *cmprG = inG;  
	PGM *cmprB = inB;  
	//printf("Birateral Filtar done.\n");
	//�L�����o�X�C���[�W����
	PGM *nimgV = create_pgm(gray->width, gray->height, gray->bright); //���x�݂̂̃L�����o�X�i��r�p�j  
	PPM *nimgC = create_ppm(in->width, in->height, in->bright); //���ۂɕ`�悷��L�����o�X
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;
	// PPM *nimgC_Scaling = create_ppm(in->width*canvas_scaling_ratio, in->height*canvas_scaling_ratio, in->bright); //���ۂɕ`�悷��L�����o�X�i�g�k�`��p�j
	// nimgC_Scaling->dataR = nimgR_Scaling->data;
	// nimgC_Scaling->dataG = nimgG_Scaling->data;
	// nimgC_Scaling->dataB = nimgB_Scaling->data;
	
	
	//sobel�t�B���^��K�������v�Z���ʂ�\�ߊi�[���Ă���
	double **sobel_abs = create_dally(in->width, in->height);
	double **sobel_angle = create_dally(in->width, in->height);
	sobel_calcu(gray, sobel_abs, sobel_angle);
	
	
	//�����X�g���[�N���珇�ԂɃX�g���[�N�����������Ă����܂��ɊG�̌`������Ă���
	for(t=thick_max; t>=thick_min; t--){
		//if(t!=12 && t!=8 && t!=4) continue;
			
		//�X�g���[�N�T�C�Y�̃K�E�X�t�B���^�𐶐�
		gauce_filter = create_dally(2*t+1, 2*t+1);
		sigma = t/3.0;
		for(i=0; i<2*t+1; i++){
	        for(j=0; j<2*t+1; j++){
	        	gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
	        }
	    }
		
		//vec�f�[�^��������
		// strcpy(vec_sentence, "t");
		// snprintf(tmp_sentence, 32, "%d", t);
		// strcat(vec_sentence, tmp_sentence);
		// log_print(vec_filename, vec_sentence, "a");
		
		
		for(y=y_defo; y<in->height; y=y+t) {  //�E�B���h�E�̑傫���ɍ��킹��
			for(x=x_defo; x<in->width; x=x+t) {  //�E�B���h�E�����炷������ς��Ƃ�
				//�E�B���h�E�̒��̍����̍��v�����
				diff_sum = break_flag = pnum = 0;
				//printf("x:%d y:%d   ",x,y);
				
				offscrn_count = 0;
				for(xc=-t; xc<=t; xc++) {
					if((x+xc)<0 || (x+xc)>in->width-1) {offscrn_count += 2*t+1;		continue;	}
					for(yc=-t; yc<=t; yc++) {
						if((y+yc)<0 || (y+yc)>in->height-1) {	offscrn_count++;
						}else{
							//diff_sum += abs(nimgV->data[x+xc][y+yc] - cmpr->data[x+xc][y+yc]);
							diff_sum += abs(nimgR->data[x+xc][y+yc] - cmprR->data[x+xc][y+yc]);
							diff_sum += abs(nimgG->data[x+xc][y+yc] - cmprG->data[x+xc][y+yc]);
							diff_sum += abs(nimgB->data[x+xc][y+yc] - cmprB->data[x+xc][y+yc]);
						}
					}
				}
				diff_sum = diff_sum/((2*t+1)*(2*t+1)-offscrn_count)/3;
				//printf("diffsum:%f \n", diff_sum);
				
				//�����̍��v����(��f������̍���)�����ȏ�Ȃ�X�g���[�N�J�n�ʒu�Ƃ���
				if(diff_sum < window_diff_border) {
					stroke_histogram[pnum]++;
					continue;
				}
				pnum=1;		//���_�m��
				p[0].x=x+0.5; p[0].y=y+0.5;
				
				//��ڂ̕`��̈悩��`��F�𕽋ς�����Ď擾
			//	bright  = calcu_color(cmpr->data, cmpr->width, cmpr->height, x, y, t);
			//	brightR = calcu_color(in->dataR, in->width, in->height, x, y, t);
			//	brightG = calcu_color(in->dataG, in->width, in->height, x, y, t);
			//	brightB = calcu_color(in->dataB, in->width, in->height, x, y, t);
				bright  = calcu_color_bi(cmpr->data, cmpr->width, cmpr->height, x, y, t, 50, gauce_filter);
				brightR = calcu_color_bi(cmprR->data, in->width, in->height, x, y, t, 50, gauce_filter);
				brightG = calcu_color_bi(cmprG->data, in->width, in->height, x, y, t, 50, gauce_filter);
				brightB = calcu_color_bi(cmprB->data, in->width, in->height, x, y, t, 50, gauce_filter);
				//p("bright",bright);
				
				theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
						gauce_filter, p[0].x, p[0].y, t, histogram_direct, &break_flag);
				
				// if(0){//if(break_flag) {}
					// stroke_histogram[pnum]++;
					// direct_count++;
					// printf("DIRECTION IS UNSTABLE.\n");
					// continue;
				// }
				
				//pd("theta1",theta* 180.0 / PI);
				
				//����_���������v�Z�����
				p[1] = calcu_point(cmpr, p[0], t, theta);
				
				
				//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, brightR);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, brightG);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, brightB);
				//pd("sum1",sum);
				
				//��ڂ̐���_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�`�悹�����Ε����̐���_������
				if( sum < color_diff_border){
					//������̌��z�����̓_����
					theta += PI; 	//printf("+PI\n");
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//���Ε����̑��_�̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, brightR);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, brightG);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, brightB);
					//pd("sum1.2",sum);
					
					//�ǂ���̑��_���s�K�؂Ȃ�`����������̃��[�v��
					if( sum < color_diff_border) {
						//printf("BOTH POINT2 IS INAPRROPRIATE\n");
						inapp_count++;
						stroke_histogram[pnum]++;
						continue;
					}
					else{ P2b_count++;}
				}
				//�K�؂ȑ��_��������Ύ���
				else{  P2a_count++; }
				pnum=2;		//���_�m��
				
				
				/*
					POINT2������z�ɂ�莟�̐���_��T���Ă���
				*/
				
				while(pnum!=(max_stroke)){
					former_theta=theta;
					
					//��pnum�_����ɂ����āAsobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
					theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
											gauce_filter, p[pnum-1].x, p[pnum-1].y, t, histogram_direct2, &break_flag);
					//if(break_flag==1) break;
					
					//����_�ׂ̈��p���}�s�ɂȂ�悤�Ȃ�t�����Ɋp�x�����
					if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;} 
					p[pnum] = calcu_point(cmpr, p[pnum-1], t, theta);
					//printf("theta%d:%f\n",pnum,theta* 180.0 / PI);
					
					//pnum+1�ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[pnum], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, brightR);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, brightG);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, brightB);
					//printf("sum%d:%f\n",pnum,sum);
					
					/*
						pnum+1�ڂ�(����)����_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�
						����܂ł̐���_��p���Đ���`��
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}
					
				}
				
				//�Z�o����pnum�̐���_��p���ăX�g���[�N��`��
				if(pnum>=min_stroke) { 
					// printf("C%f: (%f:%f)",pnum,p[0].x,p[0].y);
					// for(i=1; i<pnum; i++){
						// printf("->(%f:%f)", p[i].x, p[i].y);
					// }
					
					Paint_Bezier_ex(p, pnum, nimgV, t, bright, ratio);	
					Paint_Bezier_ex(p, pnum, nimgR, t, brightR, ratio);	
					Paint_Bezier_ex(p, pnum, nimgG, t, brightG, ratio);	
					Paint_Bezier_ex(p, pnum, nimgB, t, brightB, ratio);	
					stroke_histogram[pnum]++;  
					
					// ����_��S�Ăƃu���V�T�C�Y���g�嗦�ɏ]���X�P�[�����O���A�g��L�����o�X�ɕ`��
					// scaling_p = scaling_point(p, pnum, canvas_scaling_ratio);

					// Paint_Bezier_ex(scaling_p, pnum, nimgR_Scaling, t*canvas_scaling_ratio, brightR, ratio);	
					// Paint_Bezier_ex(scaling_p, pnum, nimgG_Scaling, t*canvas_scaling_ratio, brightG, ratio);	
					// Paint_Bezier_ex(scaling_p, pnum, nimgB_Scaling, t*canvas_scaling_ratio, brightB, ratio);
					
					//�X�g���[�N�f�[�^���t�@�C���ɒǉ�
					// vec_print(vec_filename, p, pnum, brightR, brightG, brightB, nimgV->width, nimgV->height);
					
					nc++;
				}
				
				
				
				//���X�g���[�N���Ƃɓr���o�߉摜�������o��
				//if(nc%1000==0)
				//if(0)
				if(t!=tc){
					if(t%2==0){ //�u���V���a��2�ς�邲��
						paint_count++;
						tc=t;	
						snprintf(count_name, 16, "%03d", tc+1);
						strcpy(out_filename, dir_path);
						strcat(out_filename, in_filename);
						strcat(out_filename, "_t");
						strcat(out_filename, count_name);
						strcat(out_filename, ".jpg");
						out_png = PPM_to_image(nimgC);
						// out_png = PPM_to_image(nimgC_Scaling);
						if(write_jpeg_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
						free_image(out_png);
						printf("%s\n",out_filename);
						printf("%d:",t);
						pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
					}
				}
				//�f�o�b�O�p
				/*
				if(input_char!='c'){
					do{
						if((input_char = getchar())=='e')exit(0);
						if(input_char == 'p'){
							for(i=x-t; i<=x+t; i++) {
								for(j=y-t; j<=y+t; j++) {
									if(i<0 || i>cmpr->width-1 || j<0 || j>cmpr->height-1) continue;
									printf("%d ",cmpr->data[i][j]);
								}
								pn;
							}
						}
					}while(input_char != '\n' && input_char != 'c');
				}
				pn;
				*/
			}
		}
		
		printf("t:%d , paint_num:%d\n", t,nc);
		printf("////////////////////\nt%d done.\n////////////////////\n\n",t);
		Free_dally(gauce_filter, 2*t+1);
		
		lc++;		//�������a�ł̃��[�v��cont�񂷂�
		if((lc%loop_cont) != 0){
			p("lc",lc);
			t++;
			x_defo += t/loop_cont;
			y_defo += t/loop_cont;
			if(t/loop_cont==0){
				x_defo++;	y_defo++;
			}
		}
		
	}
	
	//���i�K�`���̒��ԉ摜���o��
	paint_count++;
	snprintf(count_name, 16, "%03d", tc);
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "_t");
	strcat(out_filename, count_name);
	strcat(out_filename, ".jpg");
	out_png = PPM_to_image(nimgC);
	// out_png = PPM_to_image(nimgC_Scaling);
	if(write_jpeg_file(out_filename, out_png)){ printf("WRITE JPG ERROR.");}
	free_image(out_png);
	printf("%s\n",out_filename);
	tc=t;
	printf("%d",t);
	pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	
	
	
	
	//---------------------------
	//�G�b�W�}�b�v���v�Z���A�G�b�W�̕��G�Ȏ��ӂ����ɕ`����s��
	//---------------------------
	//color_diff_border = 0;
	//min_stroke = 4;
	//window_diff_border = 2;
	//ratio = 0.30;
	loop_cont = 2;
	tc=-1;
	maxValue=0.30, minValue=0.10;
	PGM *canny;
	PGM *EdgeMap;
	
	thick_max = opt2_thick_max;
	if(thick_max){
		canny = cannyedge_detector(gray, maxValue, minValue, thick_min);
		EdgeMap = calcu_EdgeMap(canny, thick_min, sobel_angle);
		//EdgeMap = expand_Edge(canny, thick_min);
	}
	thick_min = opt2_thick_min;
	
	for(t=thick_max; t>=thick_min; t--){
		color_diff_border = color_diff_border*(2*t+1)*(2*t+1);
		//�X�g���[�N�T�C�Y�̃K�E�X�t�B���^�𐶐�
		gauce_filter = create_dally(2*t+1, 2*t+1);
		sigma = t/3.0;
		for(i=0; i<2*t+1; i++){
	        for(j=0; j<2*t+1; j++){
	        	gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
	        }
	    }
		
		for(y=y_defo; y<in->height; y=y+t) {  //�E�B���h�E�̑傫���ɍ��킹��
			for(x=x_defo; x<in->width; x=x+t) {  //�E�B���h�E�����炷������ς��Ƃ�
				//�G�b�W�}�b�v�ɂȂ��ӏ��Ȃ玟�̃��[�v��
				if(EdgeMap->data[x][y]!=255) continue;
				
				//�E�B���h�E�̒��̍����̍��v�����
				diff_sum = break_flag = pnum = 0;
				printf("x:%d y:%d   ",x,y);
				
				offscrn_count = 0;
				for(xc=-t; xc<=t; xc++) {
					if((x+xc)<0 || (x+xc)>in->width-1) {offscrn_count += 2*t+1;		continue;	}
					for(yc=-t; yc<=t; yc++) {
						if((y+yc)<0 || (y+yc)>in->height-1) {	offscrn_count++;
						}else{
							//diff_sum += abs(nimgV->data[x+xc][y+yc] - cmpr->data[x+xc][y+yc]);
							diff_sum += abs(nimgR->data[x+xc][y+yc] - cmprR->data[x+xc][y+yc]);
							diff_sum += abs(nimgG->data[x+xc][y+yc] - cmprG->data[x+xc][y+yc]);
							diff_sum += abs(nimgB->data[x+xc][y+yc] - cmprB->data[x+xc][y+yc]);
						}
					}
				}
				diff_sum = diff_sum/((2*t+1)*(2*t+1)-offscrn_count)/3;
				printf("diffsum:%f \n", diff_sum);
				
				//�����̍��v����(��f������̍���)�����ȏ�Ȃ�X�g���[�N�J�n�ʒu�Ƃ���
				if(diff_sum < window_diff_border) {
					stroke_histogram[pnum]++;
					continue;
				}
				pnum=1;		//���_�m��
				
				//��ڂ̕`��̈悩��`��F�𕽋ς�����Ď擾
			//	bright  = calcu_color(cmpr->data, cmpr->width, cmpr->height, x, y, t);
			//	brightR = calcu_color(in->dataR, in->width, in->height, x, y, t);
			//	brightG = calcu_color(in->dataG, in->width, in->height, x, y, t);
			//	brightB = calcu_color(in->dataB, in->width, in->height, x, y, t);
				bright  = calcu_color_bi(cmpr->data, cmpr->width, cmpr->height, x, y, t, 50, gauce_filter);
				brightR = calcu_color_bi(cmprR->data, in->width, in->height, x, y, t, 50, gauce_filter);
				brightG = calcu_color_bi(cmprG->data, in->width, in->height, x, y, t, 50, gauce_filter);
				brightB = calcu_color_bi(cmprB->data, in->width, in->height, x, y, t, 50, gauce_filter);
				p("bright",bright);
				
				
				theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
						gauce_filter, p[0].x, p[0].y, t, histogram_direct, &break_flag);
				
				if(0){//if(break_flag) {}
					stroke_histogram[pnum]++;
					direct_count++;
					printf("DIRECTION IS UNSTABLE.\n");
					continue;
				}
				
				pd("theta1",theta* 180.0 / PI);
				
				//����_���������v�Z�����
				p[0].x=x; p[0].y=y;
				p[1] = calcu_point(cmpr, p[0], t, theta);
				
				
				//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, brightR);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, brightG);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, brightB);
				pd("sum1",sum);
				
				//��ڂ̐���_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�`�悹�����Ε����̐���_������
				if( sum < color_diff_border){
					//������̌��z�����̓_����
					theta += PI; 	printf("+PI\n");
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//���Ε����̑��_�̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, brightR);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, brightG);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, brightB);
					pd("sum1.2",sum);
					
					//�ǂ���̑��_���s�K�؂Ȃ�`����������̃��[�v��
					if( sum < color_diff_border) {
						printf("BOTH POINT2 IS INAPRROPRIATE\n");
						inapp_count++;
						stroke_histogram[pnum]++;
						continue;
					}
					else{ P2b_count++;}
				}
				//�K�؂ȑ��_��������Ύ���
				else{  P2a_count++; }
				pnum=2;		//���_�m��
				
				
				/*
					POINT2������z�ɂ�莟�̐���_��T���Ă���
				*/
				
				while(pnum!=(max_stroke)){
					former_theta=theta;
					
					//��pnum�_����ɂ����āAsobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
					theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
											gauce_filter, p[pnum-1].x, p[pnum-1].y, t, histogram_direct2, &break_flag);
					//if(break_flag==1) break;
					
					//����_�ׂ̈��p���}�s�ɂȂ�悤�Ȃ�t�����Ɋp�x�����
					if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;} 
					p[pnum] = calcu_point(cmpr, p[pnum-1], t, theta);
					printf("theta%d:%f\n",pnum,theta* 180.0 / PI);
					
					//pnum+1�ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[pnum], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, brightR);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, brightG);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, brightB);
					printf("sum%d:%f\n",pnum,sum);
					
					/*
						pnum+1�ڂ�(����)����_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�
						����܂ł̐���_��p���Đ���`��
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}
					
				}
				
				if(pnum>=min_stroke) { 
					//�Z�o����pnum�̐���_��p���ăX�g���[�N��`��
					printf("C%d: (%d:%d)",pnum,p[0].x,p[0].y);
					for(i=1; i<pnum; i++){
						printf("->(%d:%d)", p[i].x, p[i].y);
					}
					pn;
					Paint_Bezier_ex(p, pnum, nimgV, t, bright, ratio);	
					Paint_Bezier_ex(p, pnum, nimgR, t, brightR, ratio);	
					Paint_Bezier_ex(p, pnum, nimgG, t, brightG, ratio);	
					Paint_Bezier_ex(p, pnum, nimgB, t, brightB, ratio);	
					stroke_histogram[pnum]++;  
				}
				
				
				
				
				//�X�g���[�N���Ƃɉ摜�t�@�C���ɏ����o��
				nc++;
				//if(nc%1000==0)
				if(t!=tc)
				if(1){
					paint_count++;
					snprintf(count_name, 16, "%08d", paint_count);
					strcpy(in_filename, filename);
					strcpy(out_filename, strtok(in_filename,".") );
					strcat(out_filename, count_name);
					strcat(out_filename, ".pnm");
					if(write_ppm(out_filename, nimgC)){	exit(1);}
					printf("%s\n",out_filename);
					tc=t;
				}
				
				/*
				if(input_char!='c'){
					do{
						if((input_char = getchar())=='e')exit(0);
						if(input_char == 'p'){
							for(i=x-t; i<=x+t; i++) {
								for(j=y-t; j<=y+t; j++) {
									if(i<0 || i>cmpr->width-1 || j<0 || j>cmpr->height-1) continue;
									printf("%d ",cmpr->data[i][j]);
								}
								pn;
							}
						}
					}while(input_char != '\n' && input_char != 'c');
				}
				pn;
				*/
			}
		}
		
	
		printf("////////////////////\nt%d done.\n////////////////////\n\n",t);
		Free_dally(gauce_filter, 2*t+1);
		
		lc++;		//�������a�ł̃��[�v��cont�񂷂�
		if((lc%loop_cont) != 0){
			p("lc",lc);
			t++;
			x_defo += t/loop_cont;
			y_defo += t/loop_cont;
			if(t/loop_cont==0){
				x_defo++;	y_defo++;
			}
		}
		
	}
	
	/*
	strcpy(in_filename, filename);
	strcpy(out_filename, strtok(in_filename,".") );
	strcat(out_filename, "EdgeMap");
	strcat(out_filename, ".pnm");
	if(write_pgm(out_filename, EdgeMap)){	exit(1);}
	printf("%s\n",out_filename);
	*/
//	for(i=0; i<in->width; i++) {
//		for(j=0; j<in->height; j++) {
//			nimgV->data[i][j] = sobel_abs[i][j];
//		}
//	}
	/*
	strcpy(in_filename, filename);
	strcpy(out_filename, strtok(in_filename,".") );
	strcat(out_filename, "Canny");
	
	strcat(out_filename, ".pnm");
	if(write_pgm(out_filename, canny)){	exit(1);}
	printf("%s\n",out_filename);
	*/	
	
	
	printf("%s\n", log_filename);
	if(log_print(log_filename, log_sentence, "w") ){ printf("LOG_PRINTING_FAIL\n"); }
	
	//p("diffsum_border",window_diff_border);
	//p("color_border",color_diff_border);
	//p("thick_max",thick_max);
	//p("thick_min",thick_min);
	//p("max_stroke",max_stroke);
	//p("min_stroke",min_stroke);
	//p("loop_cont",loop_cont);
	//pd("ratio",ratio);
	//pn;
	
	//p("direction_unstable",direct_count);
	//printf("histogram_direction");  display_ally(histogram_direct, histogram_partition);
	//printf("histogram_direction2");  display_ally(histogram_direct2, histogram_partition);
	//p("Point alpha",P2a_count);
   //	p("Point beta",P2b_count);
	//p("BOTH INAPRROPRIATE", inapp_count);
	//printf("stroke_histogram");  display_ally(stroke_histogram, max_stroke+1);
	Free_dally(sobel_abs, in->width);
	Free_dally(sobel_angle, in->width);
	FreePGM(gray);
	FreePGM(inR);
	FreePGM(inG);
	FreePGM(inB);
	FreePGM(nimgV);
	
	//_CrtDumpMemoryLeaks();  
	
	pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);	
	return nimgC;
}
