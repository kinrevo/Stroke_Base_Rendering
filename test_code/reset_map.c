#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../sbr.h"


//配列を動的に確保する関数
int **create_ally(int width, int height) {
	int i;
	int **buf = (int**)malloc(sizeof(int*)*(width));
	for(i=0; i<width; i++) { buf[i]=(int*)malloc(sizeof(int)*(height)); }
	return buf;
}


//　PGM画像領域を確保（最大輝度で初期化）
PGM *create_pgm(int width, int height, int bright){
	PGM *img = (PGM *)malloc(sizeof(PGM));
	memcpy(img->descriptor, "P3", 3);
	img->height = height;
	img->width = width;
	img->bright = bright;
	img->data = create_ally(width, height);
	format_ally(img->data, img->width, img->height, bright);
	
	return img;
}

void display_ally(int *ally, int num) {
	int i;
	printf("(%d",ally[0]);
	for(i=1; i<num; i++) {
		printf(" %d", ally[i]);
	}
	printf(")\n");
}


void reset_improved_value_map(PGM* map, Point* p, int pnum, int t, int max_stroke){
	int i,j,x,y;
	int left_end,right_end,upper_end,lower_end;
	
	// ストロークの各制御点から最大ストローク長だけ離れた座標までリセット
	for(i=0; i<pnum; i++){
		// ストロークの最大長＋描画されたストロークの半径＋誤差吸収遊び
		left_end  = p[i].x - t*max_stroke - t -2;
		right_end = p[i].x + t*max_stroke + t +2;
		upper_end = p[i].y - t*max_stroke - t -2;
		lower_end = p[i].y + t*max_stroke + t +2; 
		if(left_end<0)left_end=0; if(map->width <= right_end)right_end=map->width-1; 
		if(upper_end<0)upper_end=0; if(map->height <= lower_end)lower_end=map->height-1; 
		
		for(y=upper_end; y<=lower_end; y++) { 
			for(x=left_end; x<=right_end; x++) {  
				map->data[x][y] = UNCALCULATED;
			}
		}
	}
}
void format_ally(int **ally, int w, int h, int bright) {
	int i, j;
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally[i][j] = bright;
		}
	}
}

int main(){
	int i,j;
	int pnum =2, t=1;
	Point p[4];
	p[0].x=1; p[0].y=1;
	p[1].x=1; p[1].y=3;
	
	PGM* map = create_pgm(16, 16, 0);

	for(i=0; i<16; i++){
		display_ally(map->data[i], 16);
	}
	
	reset_improved_value_map(map, p, pnum, t, 2);
	
	p("t",t);
	
	for(i=0; i<16; i++){
		display_ally(map->data[i], 16);
	}
	
	
	return 0;
}