#include <stdio.h>
#include <stdlib.h>
#include "../sbr.h"


//配列を動的に確保する関数
int **create_ally(int width, int height) {
	int i;
	int **buf = (int**)malloc(sizeof(int*)*(width));
	for(i=0; i<width; i++) { buf[i]=(int*)malloc(sizeof(int)*(height)); }
	return buf;
}

//配列を初期化する関数
void format_ally(int **ally, int w, int h, int bright) {
	int i, j;
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			ally[i][j] = bright;
		}
	}
}

// 配列中の最大値を探索しPointを返す
Point search_max_Point(int **ally, int w, int h) {
	int i, j, max_value=0;
	Point max;
	max.x=0; max.y=0;
	for(i=0; i<w; i++) {
		for(j=0; j<h; j++) {
			if(ally[i][j]>max_value) {
				max.x=i;
				max.y=j;
				max_value=ally[i][j];
			}
		}
	}
	
	return max;
}

int main(){

	int** img=create_ally(10,10);
	format_ally(img, 10, 10, 0);
	
	
	img[2][2]=-999;
	img[9][9]=999;
	
	Point max=search_max_Point(img,10,10);
	
	printf("%d:%d\n",(int)max.x,(int)max.y);	
	printf("%d\n",img[(int)max.x][(int)max.y]);
	
	return 0;
}