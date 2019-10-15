#include <stdio.h>
#include "sbr.h"



int func(Point* P){
	Point tmp = {3,2};
	P[2] = tmp;
	P[3] = tmp;
	P[4] = tmp;
	return 0;
}

int main(){
	Point P[10];
	int i;
	for(i=0; i<10; i++) {
		P[i].x = i;
	}
	for(i=0; i<10; i++) {
		printf("Point[%d]:%f\n",i,P[i].x);
		pn;
	}
	
	func(P);
	
	for(i=0; i<10; i++) {
		printf("Point[%d]:%f\n",i,P[i].x);
		pn;
	}
	
	return 0;
}