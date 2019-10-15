#include <stdio.h>
#include <stdlib.h>
#include "../sbr.h"




//　ストローク構造体メモリを確保
Stroke *create_Stroke(int pnum){	
	Stroke* sp = (Stroke *)malloc(sizeof(Stroke));
	sp->pnum = pnum;
	sp->p = (Point*)malloc(sizeof(Point)*(sp->pnum));
	
	return sp;
}

int main(){
	int i,j;

	Stroke* sp = (Stroke *)malloc(sizeof(Stroke));
	sp->pnum = 10;
	sp->p = (Point*)malloc(sizeof(Point)*(sp->pnum));
	//sp->color = (RGB *)malloc(sizeof(RGB));
	sp->color.R = 1; 
	sp->color.G = 2;
	sp->color.B = 3;
	sp->p[9].x=29;
	sp->p[9].y=229;
	
	printf("%d",sp->color.R);
	printf("%d",sp->color.G);
	printf("%d",sp->color.B);
	printf("%f",sp->p[9].x);
	printf("%f",sp->p[9].y);
	pn;
	
	Stroke*** sp_map = (Stroke***)malloc(sizeof(Stroke**)*(512));
	for(i=0; i<512; i++){
		sp_map[i] = (Stroke**)malloc(sizeof(Stroke*)*(512));
	}
	
	for(i=0; i<512; i++){
		for(j=0; j<512; j++){
			sp_map[i][j]=create_Stroke(10);
		}
	}
	
	// sp_map[1][1]->color.R = 1; 
	// sp_map[1][1]->color.G = 2;
	// sp_map[1][1]->color.B = 3;
	sp_map[1][1]->color = sp->color;
	// sp_map[1][1]->p[9].x=29;
	// sp_map[1][1]->p[9].y=229;
	sp_map[1][1]->p[9] = sp->p[9];
	
	printf("%d",sp_map[1][1]->color.R);
	printf("%d",sp_map[1][1]->color.G);
	printf("%d",sp_map[1][1]->color.B);
	printf("%f",sp_map[1][1]->p[9].x);
	printf("%f",sp_map[1][1]->p[9].y);
	pn;
	
	sp->color.R = -11111; 
	sp->p[9].x=-999999;
	
	printf("%d",sp_map[1][1]->color.R);
	printf("%d",sp_map[1][1]->color.G);
	printf("%d",sp_map[1][1]->color.B);
	printf("%f",sp_map[1][1]->p[9].x);
	printf("%f",sp_map[1][1]->p[9].y);
	
	
	return 0;
}