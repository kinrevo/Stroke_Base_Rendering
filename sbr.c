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

	Point temp;  //�`�����ʂ�Ȃ�����_
	int i,j,r;
	double t;
	int partition = abs(p[0].x-p[1].x)+abs(p[0].y-p[1].y); //�����̕�����
	
	Point bp0 = {(p[1].x-p0.x)/6.0 + p[0].x, (p[1].y-p0.y)/6.0 + p[0].y}
		,bp1 = {(p[0].x-p[2].x)/6.0 + p[1].x, (p[0].y-p[2].y)/6.0 + p[1].y};
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick*sin((10*t<PI/2 ? 10*t:PI/2));     //�������x�݂̂����ɏ]������
		temp = BezierCurve_P(p[0], bp0, bp1, p[1], t);
		Circle_fill(temp.x, temp.y, r, in_img, bright, ratio);
	}
	
	for(i=1; i<pnum-2; i++){
		bp0.x = (p[i+1].x-p[i-1].x)/6.0 + p[i].x;	bp0.y = (p[i+1].y-p[i-1].y)/6.0 + p[i].y;
		bp1.x = (p[i].x-p[i+2].x)/6.0 + p[i+1].x;	bp1.y = (p[i].y-p[i+2].y)/6.0 + p[i+1].y;
		partition = abs(p[i].x-p[i+1].x)+abs(p[i].y-p[i+1].y);
		for(j=0; j <= partition; j++) {
			t = (double)j/partition;
			r = thick;
			temp = BezierCurve_P(p[i], bp0, bp1, p[i+1], t);
			Circle_fill(temp.x, temp.y, r, in_img, bright, ratio);
		}
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
		Circle_fill(temp.x, temp.y, r, in_img, bright, ratio);
	}
}



//�^����ꂽ�摜�ɒ�����`�悷��֐�
void Paint_line(Point p1, Point p4, PGM* in_img, int thick, int bright, double ratio) {
	int i,r;
	int partition = abs(p1.x-p4.x) + abs(p1.y-p4.y); //�����̕�����
	float t;
	Point temp;
	
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



//������p�xtheta�����ɋ���t�̓_�����߂�
Point calcu_point(PGM *in, Point a, int t, double theta){
	Point b;
	b.x = t*cos(theta)+a.x;  
	if(b.x<0) b.x=0;  
	if(b.x>=in->width) b.x=in->width-1; 
	b.y = t*sin(theta)+a.y;
	if(b.y<0) b.y=0; 
	if(b.y>=in->height) b.y=in->height-1;
	
	return b;
}


//�^����ꂽ�_����ŕ`��F��h�����ۂ̉��P�l�����߂�
double diffsum_clr(PGM* cmpr, PGM* nimg, Point p, int t, int bright){
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
	sum = sum / ((2*t+1)*(2*t+1)-offscrn_count);  //�_������̍��ٕ���
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
					double **gauce_filter, int x, int y, int t, int *break_flag)
{
	int i,j,Farg, im1,ip1;//,im2,ip2;
	int peak=partition/2+0.5;		//�l�̌ܓ�
	double histogram[partition];
	double theta;
	for(i=0; i<partition; i++) {histogram[i]=0;}
	
	//sobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
	for(i=x-t-1; i<=x+t+1; i++) {
		for(j=y-t-1; j<=y+t+1; j++) {
			if(i<1 || i>cmpr->width-2 || j<1 || j>cmpr->height-2) continue;
			Farg = (sobel_angle[i][j]/PI)*partition;
			histogram[Farg] += sobel_abs[i][j];
			if(opt_USE_gause_histogram) histogram[Farg] *= gauce_filter[i-(x-t)][j-(y-t)];
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
						flag2_c++;
						// printf("flag2\n");
						break;
					}
					px=p.x; py=p.y;
					map->data[px][py]=255;
					if(flag2a!=1){
						//�������G�b�W�̏ꏊ���L�^
						// pp(p); 
						//�������G�b�W�̉ӏ�����Ȃ����Ă���G�b�W��S�ĒT��
						edge_detect(canny, map, sobel_angle, x, y, t, p, &flag1, 0);
						// printf("detect_done\n");
						
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

//�@�����ɕ`���Ă݂Č덷���m�F
int test_stroke(PGM* test_Canvas, PGM* cmprR, PGM* cmprG, PGM* cmprB, PGM* nimgR, PGM* nimgG, PGM* nimgB, Point p[], int pnum, int t, int brightR, int brightG, int brightB, double ratio)
{
	int i,j,diffsum=0;
	double left_end,right_end,upper_end,lower_end;//upper,lower�͉摜�̌������̏㉺�i�l�̏㉺�ł͂Ȃ��j
	
	//�X�g���[�N�_���͂ޒ[�̍��W�����
	left_end=right_end=p[0].x;
	upper_end=lower_end=p[0].y;
	for(i=1; i<pnum; i++){
		if(p[i].x < left_end) left_end=p[i].x;
		if(right_end < p[i].x) right_end=p[i].x;
		if(p[i].y < upper_end) upper_end=p[i].y;
		if(lower_end < p[i].y) lower_end=p[i].y;
	}
	//�X�g���[�N���a���A�[���W��c��
	left_end-=2*t; right_end+=2*t; upper_end-=2*t; lower_end+=2*t;
	if(left_end<0)left_end=0; 
	if(cmprR->width <= right_end)right_end=cmprR->width-1; 
	if(upper_end<0)upper_end=0; 
	if(cmprR->height <= lower_end)lower_end=cmprR->height-1; 
	
	// printf("Test: (%3.0f:%3.0f)", p[0].x, p[0].y);
	// for(i=1; i<pnum; i++){
		// printf("->(%3.0f:%3.0f)", p[i].x, p[i].y);
	// }
	// printf("l:%3.0f,r:%3.0f,u:%3.0f,l:%3.0f\n", left_end,right_end,upper_end,lower_end);
	
	//�e�X�g�L�����o�X�ɂ�����X�g���[�N�ɂ��e���ӏ����R�s�[
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			//if(i<0 || i>cmprR->width-1 || j<0 || j>cmprR->height-1)continue;
			test_Canvas->data[i][j] = nimgR->data[i][j];
		}
	}
	//�����`��
	Paint_Bezier_ex(p, pnum, test_Canvas, t, brightR, ratio);
	//�����`���ɂ��ω����m�F
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			if(i<0 || i>cmprR->width-1 || j<0 || j>cmprR->height-1)continue;
			diffsum += abs(cmprR->data[i][j]-nimgR->data[i][j]) - abs(cmprR->data[i][j]-test_Canvas->data[i][j]);
		}
	}
	
	//�e�X�g�L�����o�X�ɂ�����X�g���[�N�ɂ��e���ӏ����R�s�[
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			test_Canvas->data[i][j] = nimgG->data[i][j];
		}
	}
	//�����`��
	Paint_Bezier_ex(p, pnum, test_Canvas, t, brightG, ratio);
	//�����`���ɂ��ω����m�F
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			diffsum += abs(cmprG->data[i][j]-nimgG->data[i][j]) - abs(cmprG->data[i][j]-test_Canvas->data[i][j]);
		}
	}
	
	//�e�X�g�L�����o�X�ɂ�����X�g���[�N�ɂ��e���ӏ����R�s�[
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			test_Canvas->data[i][j] = nimgB->data[i][j];
		}
	}
	//�����`��
	Paint_Bezier_ex(p, pnum, test_Canvas, t, brightB, ratio);
	//�����`���ɂ��ω����m�F
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			diffsum += abs(cmprB->data[i][j]-nimgB->data[i][j]) - abs(cmprB->data[i][j]-test_Canvas->data[i][j]);
		}
	}
	
	return diffsum;
}



// ����Point���J�n�_�Ƃ����Ƃ��̐���_�̑g���v�Z
// �t���O�ŏ����X�g���[�N�������������������������
int calcu_Stroke_Point(PGM* cmprR, PGM* cmprG, PGM* cmprB, PGM* nimgR, PGM* nimgG, PGM* nimgB,
		double **sobel_abs, double **sobel_angle, double **gauce_filter, 
		int t, int max_stroke, Point Start_P, int brightR, int brightG, int brightB,
		Point* Stroke_P, int direct_reversal_flag){
	int sum,pnum;
	double theta,former_theta;
	int dummy_flag;
	
	Stroke_P[0].x=Start_P.x; Stroke_P[0].y=Start_P.y;
	pnum=1;		//�@���͂ɗ^����ꂽ�_�͐������ƍl����
	
	// �����_�ɂ�����������v�Z
	theta =  calcu_histogram(cmprR, sobel_abs, sobel_angle, opt_histogram_partition, 
					gauce_filter, Start_P.x, Start_P.y, t, &dummy_flag);
	
	//�@���͐ݒ�ɏ]�������������t�]
	// if(direct_reversal_flag){
		// theta += PI;		
	// }
	
	//����_���������v�Z�����
	Stroke_P[1] = calcu_point(cmprR, Stroke_P[0], t, theta);
	
	
	//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
	sum = 0;
	sum += diffsum_clr(cmprR, nimgR, Stroke_P[1], t, brightR);
	sum += diffsum_clr(cmprG, nimgG, Stroke_P[1], t, brightG);
	sum += diffsum_clr(cmprB, nimgB, Stroke_P[1], t, brightB);
	
	
	//��ڂ̐���_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�`�悹�����Ε����̐���_������
	if(sum < opt_color_diff_border){
		//������̌��z�����̓_����
		theta += PI;
		Stroke_P[1] = calcu_point(cmprR, Stroke_P[0], t, theta);
		
		//���Ε����̑��_�̕`��_����̐F���`��F�ƈ�v���邩�m�F����
		//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //�_������̍��ٕ���
		sum = 0;
		sum += diffsum_clr(cmprR, nimgR, Stroke_P[1], t, brightR);
		sum += diffsum_clr(cmprG, nimgG, Stroke_P[1], t, brightG);
		sum += diffsum_clr(cmprB, nimgB, Stroke_P[1], t, brightB);
		
		//�ǂ���̑��_���s�K�؂Ȃ�`����������̃��[�v��
		if( sum < opt_color_diff_border) {
			return pnum;
		}
		else{ }
	}

	pnum=2;		//���P�l�����ȏ�ł���Α��_�m��
	
	/*
		POINT2������z�ɂ�莟�̐���_��T���Ă���
	*/
	while(pnum<max_stroke){
		former_theta=theta;
		
		//��pnum�_����ɂ����āAsobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
		theta =  calcu_histogram(cmprR, sobel_abs, sobel_angle, opt_histogram_partition, 
								gauce_filter, Stroke_P[pnum-1].x, Stroke_P[pnum-1].y, t, &dummy_flag);
		
		//����_�ׂ̈��p���}�s�ɂȂ�悤�Ȃ�t�����Ɋp�x�����
		if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;} 
		Stroke_P[pnum] = calcu_point(cmprR, Stroke_P[pnum-1], t, theta);
		
		//pnum+1�ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
		sum = 0;
		sum += diffsum_clr(cmprR, nimgR, Stroke_P[pnum], t, brightR);
		sum += diffsum_clr(cmprG, nimgG, Stroke_P[pnum], t, brightG);
		sum += diffsum_clr(cmprB, nimgB, Stroke_P[pnum], t, brightB);
		
		/*
			pnum+1�ڂ�(����)����_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�
			����܂ł̐���_��p���Đ���`��
		*/
		if(sum < opt_color_diff_border) { break;
		}else { pnum++; }
	}
	
	return pnum;
}


// �`�悵���X�g���[�N�̎��͂̂݉��P�l�}�b�v�����Z�b�g
void reset_improved_value_map(PGM* map, Point* p, int pnum, int t, int max_stroke){
	int i,x,y;
	int left_end,right_end,upper_end,lower_end;
	
	// �X�g���[�N�̊e����_����ő�X�g���[�N���������ꂽ���W�܂Ń��Z�b�g
	for(i=0; i<pnum; i++){
		// �X�g���[�N�̍ő咷�{�`�悳�ꂽ�X�g���[�N�̔��a�{�덷�z���V��
		left_end  = p[i].x - t*max_stroke - t -2;
		right_end = p[i].x + t*max_stroke + t +2;
		upper_end = p[i].y - t*max_stroke - t -2;
		lower_end = p[i].y + t*max_stroke + t +2; 
		if(left_end<0)left_end=0; 
		if(map->width <= right_end)right_end=map->width-1; 
		if(upper_end<0)upper_end=0; 
		if(map->height <= lower_end)lower_end=map->height-1; 
		
		for(y=upper_end; y<=lower_end; y++) { 
			for(x=left_end; x<=right_end; x++) {  
				map->data[x][y] = UNCALCULATED;
			}
		}
	}
	
}




//�^����ꂽ�C���[�W����u���b�V���O�ɂ��G������
PPM *c_Illust_brush(PPM *in, char *filename) {
	clock_t start = clock();
	int i,j,x,y,xc,yc,t,break_flag,pnum, offscrn_count;
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
	// int histogram_direct[31]={}; 
	int paint_count=0, nc=0, tc=-1;
	int loop_cont=opt_loop_cont, lc=0, x_defo=0, y_defo=0;
	double maxValue, minValue;
//	double canvas_scaling_ratio = 1.0;
	RGB bright;
	
	
	//�ő召�X�g���[�N���a�i�������F��ʂ�1/10,�ő��1/10�j
	int thick_max = opt_thick_max;//(in->height < in->width ? in->height : in->width)/10;
	int thick_min = opt_thick_min;//(thick_max/15 > 3 ? thick_max/15 : 3);
	
	
	//�o�̓t�@�C�����̃T�C�Y���擾
	int namesize = strlen(filename)*2 + 16;
	char out_filename[namesize], log_filename[namesize], dir_path[namesize];
	char in_filename[namesize-16];
	char count_name[16];
	char log_sentence[2028] = "";
	// char vec_sentence[128] = "", tmp_sentence[32] = "", vec_filename[namesize];
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
	Add_dictionary_to_sentence(log_sentence, "USE_calcu_color_bi", opt_USE_calcu_color_bi);
	Add_dictionary_to_sentence(log_sentence, "USE_gause_histogram", opt_USE_gause_histogram);
	Add_dictionary_to_sentence(log_sentence, "optimal_improved_value_border", opt_optimal_improved_value_border);
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
	PGM *cmpr  = gray;  
	PGM *cmprR = inR;  
	PGM *cmprG = inG;  
	PGM *cmprB = inB;  
	//�L�����o�X�C���[�W����
	PGM *nimgV = create_pgm(gray->width, gray->height, gray->bright); //���x�݂̂̃L�����o�X�i��r�p�j  
	PPM *nimgC = create_ppm(in->width, in->height, in->bright); //���ۂɕ`�悷��L�����o�X
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;
	PGM *test_Canvas = create_pgm(in->width, in->height, in->bright); //�덷�\���Ɏg���e�X�g�L�����o�X
	// PGM *improved_value_map = create_pgm(in->width, in->height, in->bright); //���P�l�}�b�v
	// PPM *nimgC_Scaling = create_ppm(in->width*canvas_scaling_ratio, in->height*canvas_scaling_ratio, in->bright); //���ۂɕ`�悷��L�����o�X�i�g�k�`��p�j
	// nimgC_Scaling->dataR = nimgR_Scaling->data;
	// nimgC_Scaling->dataG = nimgG_Scaling->data;
	// nimgC_Scaling->dataB = nimgB_Scaling->data;
	
	// sobel�t�B���^��K�������v�Z���ʂ�\�ߊi�[���Ă���
	double **sobel_abs = create_dally(in->width, in->height);
	double **sobel_angle = create_dally(in->width, in->height);
	sobel_calcu(gray, sobel_abs, sobel_angle);
	
	//Greedy�A�v���[�`
	// int s_count,stroke_num, best_pnum;
	// int diff_stroke_max=0;
	// int tmp_num;
	// int diff_stroke_max_ave[50]={};
	// int **reversal_map = create_ally(in->width, in->height);
	// int best_x=0, best_y=0;
	// int miss_stroke_count=0;
	// Stroke*** best_stroke_map = create_Stroke_ally(in->width, in->height, max_stroke);
	// Point best_P;
	// int optimal_improved_value_border = opt_optimal_improved_value_border;
	int diff_stroke=0;

	
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

		
		// stroke_num=9999;
		// stroke_num = in->width/t * in->height/t / (max_stroke+min_stroke)/2;
		// p("stroke_num", stroke_num);
		
		// format_ally(improved_value_map->data, improved_value_map->width, improved_value_map->height, UNCALCULATED);
		// format_ally(reversal_map, improved_value_map->width, improved_value_map->height, Reversal_OFF);
		
		// for(s_count=0; s_count<stroke_num; s_count++) {  //���锼�a�ɂ�����X�g���[�N��
			// �덷�T���̕ϐ���������
			// diff_stroke_max = UNCALCULATED;
			// format_ally(improved_value_map->data, improved_value_map->width, improved_value_map->height, UNCALCULATED);

			
		// for(y=0; y<in->height; y++) {  
		// 	for(x=0; x<in->width; x++) {  
		for(y=y_defo; y<in->height; y=y+t) {  //�E�B���h�E�̑傫���ɍ��킹��
			for(x=x_defo; x<in->width; x=x+t) {  //�E�B���h�E�����炷������ς��Ƃ�
				// ���P�l���v�Z�ς݂Ȃ�Skip
				// if(improved_value_map->data[x][y] != UNCALCULATED) continue;
			
				diff_sum = break_flag = pnum = 0;
				// diff_stroke=0;	//Greedy�A�v���[�`
				
				//�E�B���h�E�̒��̍����̍��v�����
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
				
				//�����̍��v����(��f������̍���)�����ȏ�Ȃ�X�g���[�N�J�n�ʒu�Ƃ���
				if(diff_sum < window_diff_border) {
					stroke_histogram[pnum]++;
					continue;
				}
				pnum=1;		//���_�m��
				p[0].x=x+0.5; p[0].y=y+0.5;
				
				//��ڂ̕`��̈悩��`��F���擾
				if(opt_USE_calcu_color_bi){
					bright.R = calcu_color_bi(cmprR->data, cmprR->width, cmprR->height, x, y, t, 50, gauce_filter);
					bright.G = calcu_color_bi(cmprG->data, cmprG->width, cmprG->height, x, y, t, 50, gauce_filter);
					bright.B = calcu_color_bi(cmprB->data, cmprB->width, cmprB->height, x, y, t, 50, gauce_filter);
				} else{
					bright.R = calcu_color(cmprR->data, cmprR->width, cmprR->height, x, y, t);
					bright.G = calcu_color(cmprG->data, cmprG->width, cmprG->height, x, y, t);
					bright.B = calcu_color(cmprB->data, cmprB->width, cmprB->height, x, y, t);
				}

	
				theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
						gauce_filter, p[0].x, p[0].y, t, &break_flag);
				
				
				//����_���������v�Z�����
				p[1] = calcu_point(cmpr, p[0], t, theta);
				
				
				//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
				
				
				//��ڂ̐���_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�`�悹�����Ε����̐���_������
				if(sum < color_diff_border){
					//������̌��z�����̓_����
					theta += PI;
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//���Ε����̑��_�̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					
					//�ǂ���̑��_���s�K�؂Ȃ�`����������̃��[�v��
					if( sum < color_diff_border) {
						stroke_histogram[pnum]++;
						// improved_value_map->data[x][y] = MIN_STROKE;
						continue;
					}
					else{ 
						// reversal_map[x][y]=Reversal_ON;
					}
				}

				pnum=2;		//���_�m��
				
				
				/*
					POINT2������z�ɂ�莟�̐���_��T���Ă���
				*/
				
				while(pnum<max_stroke){
					former_theta=theta;
					
					//��pnum�_����ɂ����āAsobel����q�X�g�O�������쐬���ő�̂��̂����z�Ƃ���
					theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
											gauce_filter, p[pnum-1].x, p[pnum-1].y, t, &break_flag);
					
					//����_�ׂ̈��p���}�s�ɂȂ�悤�Ȃ�t�����Ɋp�x�����
					if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;} 
					p[pnum] = calcu_point(cmpr, p[pnum-1], t, theta);
					
					//pnum+1�ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
					
					/*
						pnum+1�ڂ�(����)����_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�
						����܂ł̐���_��p���Đ���`��
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}
				}
				
				//�@�����ɕ`���Ă݂Č덷���m�F
				if(pnum>=min_stroke){
					diff_stroke = test_stroke(test_Canvas, cmprR, cmprG, cmprB, nimgR, nimgG, nimgB, p, pnum, t, bright.R, bright.G, bright.B, ratio);
					// improved_value_map->data[x][y] = diff_stroke;
				}else if(pnum<min_stroke){
					// improved_value_map->data[x][y] = MIN_STROKE;
				}

				//�Z�o����pnum�̐���_��p���ăX�g���[�N��`��
				if(diff_stroke>0)
				if(pnum>=min_stroke) { 
					Paint_Bezier_ex(p, pnum, nimgR, t, bright.R, ratio);	
					Paint_Bezier_ex(p, pnum, nimgG, t, bright.G, ratio);	
					Paint_Bezier_ex(p, pnum, nimgB, t, bright.B, ratio);	
					stroke_histogram[pnum]++;  
					paint_count++;
				}
				//���݂̃X�g���[�N���ۑ����Ă���J�n�_�̂��̂��ǂ��Ȃ�X�V
				// if(pnum>=min_stroke){
				// 	for(i=0; i<pnum; i++){
				// 		best_stroke_map[x][y]->p[i] = p[i];
				// 	}
				// 	best_stroke_map[x][y]->color=bright;
				// 	best_stroke_map[x][y]->pnum=pnum;
				// }
				// if( (diff_stroke>diff_stroke_max) && (pnum>=min_stroke) ){
					// best_x=x; best_y=y;
					// diff_stroke_max = diff_stroke;
				// }
				
			}
		}
			
			// ���P�l�}�b�v���̍ő�l��T��
			// best_P = search_max_Point(improved_value_map->data, improved_value_map->width, improved_value_map->height);
			// best_x=best_P.x;  best_y=best_P.y;
			
			
			//���߂�50�̃X�g���[�N�ł̕ω����قƂ�ǂȂ���Ύ��̔��a�X�e�b�v�Ɉڍs
			// if(t==1)optimal_improved_value_border = 1;
			// diff_stroke_max_ave[s_count%50] = diff_stroke_max;
			// // if(1)
			// if(s_count%50 == 49)
			// {
			// 	tmp_num=0;
			// 	for(i=0; i<50; i++) tmp_num+=diff_stroke_max_ave[i];
			// 	if(tmp_num/50 < optimal_improved_value_border) {
			// 		strcat(log_sentence, "\r\n");
			// 		Add_dictionary_to_sentence(log_sentence, "t", t);
			// 		Add_dictionary_to_sentence(log_sentence, "s_count", s_count);
			// 		snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
			// 		strcat(log_sentence, count_name);
			// 		strcat(log_sentence, "\r\n");
			// 		break;
			// 	}
			// }
			
			
			//�Z�o����pnum�̐���_��p���ăX�g���[�N��`��
			// Paint_Bezier_ex(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, nimgR, t, best_stroke_map[best_x][best_y]->color.R, ratio);	
			// Paint_Bezier_ex(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, nimgG, t, best_stroke_map[best_x][best_y]->color.G, ratio);	
			// Paint_Bezier_ex(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, nimgB, t, best_stroke_map[best_x][best_y]->color.B, ratio);	
			// stroke_histogram[best_stroke_map[best_x][best_y]->pnum]++;  
			
			// �`�悵���X�g���[�N�̎��͂̂݉��P�l�}�b�v�����Z�b�g
			// reset_improved_value_map(improved_value_map, best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, t, max_stroke);
			
					// ����_��S�Ăƃu���V�T�C�Y���g�嗦�ɏ]���X�P�[�����O���A�g��L�����o�X�ɕ`��
					// scaling_p = scaling_point(p, pnum, canvas_scaling_ratio);
					// Paint_Bezier_ex(scaling_p, pnum, nimgR_Scaling, t*canvas_scaling_ratio, brightR, ratio);	
					// Paint_Bezier_ex(scaling_p, pnum, nimgG_Scaling, t*canvas_scaling_ratio, brightG, ratio);	
					// Paint_Bezier_ex(scaling_p, pnum, nimgB_Scaling, t*canvas_scaling_ratio, brightB, ratio);
					//�X�g���[�N�f�[�^���t�@�C���ɒǉ�
					// vec_print(vec_filename, p, pnum, brightR, brightG, brightB, nimgV->width, nimgV->height);
					
			//���X�g���[�N���Ƃɓr���o�߉摜�������o��
			nc++;
			// if(nc%100==0)
			//if(t!=tc)
			// {
			// 		tc=t;	
			// 		strcpy(out_filename, dir_path);
			// 		strcat(out_filename, in_filename);
			// 		snprintf(count_name, 16, "%02d", t);
			// 		strcat(out_filename, "_t");
			// 		strcat(out_filename, count_name);
			// 		snprintf(count_name, 16, "%d", paint_count);
			// 		strcat(out_filename, "_s");
			// 		strcat(out_filename, count_name);
			// 		strcat(out_filename, ".jpg");
			// 		out_png = PPM_to_image(nimgC);
			// 		// out_png = PPM_to_image(nimgC_Scaling);
			// 		if(write_jpeg_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
			// 		free_image(out_png);
			// 		printf("%s\n",out_filename);
			// 		printf("%d:",t);
			// 		pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
			// }
						
		if(t!=tc)
		{
			tc=t;	
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			snprintf(count_name, 16, "%02d", t+1);
			strcat(out_filename, "__t");
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
			
			// Greedy�A�v���[�`�F�œK�X�g���[�N��T�����邲�ƂɒT���ʒu�����炷
			// x_defo += 1;
			// if(x_defo>t) {
			// 	x_defo -= t;
			// 	y_defo += 1;
			// 	if(y_defo>t) y_defo -= t;
			// }
		//}
		
		// printf("miss_stroke/t%d:%d\n",t,miss_stroke_count);
	
		printf("\n////////////////////\nt%d done.\n////////////////////\n",t);
		p("Paint_num",paint_count); paint_count=0;
		Free_dally(gauce_filter, 2*t+1);
		x_defo=y_defo=paint_count=0;		
	}
	
	//���i�K�`���̒��ԉ摜���o��
	snprintf(count_name, 16, "%03d", tc);
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__t");
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
				if(opt_USE_calcu_color_bi){
					bright.R = calcu_color_bi(cmprR->data, cmprR->width, cmprR->height, x, y, t, 50, gauce_filter);
					bright.G = calcu_color_bi(cmprG->data, cmprG->width, cmprG->height, x, y, t, 50, gauce_filter);
					bright.B = calcu_color_bi(cmprB->data, cmprB->width, cmprB->height, x, y, t, 50, gauce_filter);
				} else{
					bright.R = calcu_color(cmprR->data, cmprR->width, cmprR->height, x, y, t);
					bright.G = calcu_color(cmprG->data, cmprG->width, cmprG->height, x, y, t);
					bright.B = calcu_color(cmprB->data, cmprB->width, cmprB->height, x, y, t);
				}
				
				
				theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition, 
						gauce_filter, p[0].x, p[0].y, t, &break_flag);
				
				
				pd("theta1",theta* 180.0 / PI);
				
				//����_���������v�Z�����
				p[0].x=x; p[0].y=y;
				p[1] = calcu_point(cmpr, p[0], t, theta);
				
				
				//��ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
				
				//��ڂ̐���_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�`�悹�����Ε����̐���_������
				if( sum < color_diff_border){
					//������̌��z�����̓_����
					theta += PI; 	printf("+PI\n");
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//���Ε����̑��_�̕`��_����̐F���`��F�ƈ�v���邩�m�F����
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					
					//�ǂ���̑��_���s�K�؂Ȃ�`����������̃��[�v��
					if( sum < color_diff_border) {
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
											gauce_filter, p[pnum-1].x, p[pnum-1].y, t, &break_flag);
					
					//����_�ׂ̈��p���}�s�ɂȂ�悤�Ȃ�t�����Ɋp�x�����
					if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;} 
					p[pnum] = calcu_point(cmpr, p[pnum-1], t, theta);
					
					//pnum+1�ڂ̕`��_����̐F���`��F�ƈ�v���邩�m�F����  //�_������̍��ٕ���
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
					
					/*
						pnum+1�ڂ�(����)����_����̐F���`��F�Ƃ������l�ȏ�̍������Ȃ�
						����܂ł̐���_��p���Đ���`��
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}
					
				}
				
				if(pnum>=min_stroke) { 
					//�Z�o����pnum�̐���_��p���ăX�g���[�N��`��
					pn;
					Paint_Bezier_ex(p, pnum, nimgR, t, bright.R, ratio);	
					Paint_Bezier_ex(p, pnum, nimgG, t, bright.G, ratio);	
					Paint_Bezier_ex(p, pnum, nimgB, t, bright.B, ratio);	
					stroke_histogram[pnum]++;  
				}
				
				
				
				
				//�X�g���[�N���Ƃɉ摜�t�@�C���ɏ����o��
				nc++;
				//if(nc%1000==0)
				if(t!=tc){
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
	
	double MSE = image_MSE(nimgC, in);
	Add_dictionary_to_sentence(log_sentence, "MSE", MSE);
	
	Add_dictionary_to_sentence(log_sentence, "All_Execution_TIME", (double)(clock()-start)/CLOCKS_PER_SEC);
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
