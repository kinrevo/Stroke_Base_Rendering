#include"pgm.h"
#include<stdlib.h>
#include<stdio.h>

int main(int argc, char *argv[])
{
	if(argc>4 || argc<2){
		fprintf(stderr, "Usage: program <inputfile> <outputfile>\n");
		exit(1);
	}

	PPM *in_img, *trans_img;	//, *in_img2, *trans_img2, *trans_img3;

	int i;
	/*
	if((in_img = read_pgm(argv[1])) == NULL){
		exit(1);
	}
	*/
	if((in_img = read_ppm(argv[1])) == NULL){
		printf("read_error\n");
		exit(1);
	}
	if(argc ==3) {
		//PGM *inR,*inG,*inB,*inL;
		//devide_ppm(in_img, inR, inG, inB);
		//trans_img = inR;
		trans_img = c_Illust_brush(in_img, argv[2]);
	//	trans_img = bilateral_filter(in_img, 4);
	//	for(i=0; i<3; i++){			
	//		trans_img2 = bilateral_filter(trans_img, 4);
	//		Freeimg(trans_img);
	//		trans_img = bilateral_filter(trans_img2, 4);
	//		Freeimg(trans_img2);
	//	}
		
		if(write_ppm(argv[2], trans_img)){
			exit(1);
		}
	}
	/*
	if(argc ==4) {
		if((in_img2 = read_pgm(argv[2])) == NULL){ exit(1);}
		
		trans_img = SIFT(in_img, 2);
		trans_img2 = SIFT(in_img2, 2);
		
		trans_img = RotateImg(trans_img,1);
		trans_img3 = CompareImg(trans_img, trans_img2);
		
		if(write_pgm(argv[3], trans_img3)){
			printf("WError\n");
			exit(1);
		}		
		
		Freeimg(trans_img2);
		Freeimg(trans_img3);
		Freeimg(in_img2);
	}
	*/
		
	


	//Freeimg(trans_img);
	//Freeimg(in_img);
	/*
	for(i=0; i < in_img->height; i++) { free(in_img->data[i]);}
	free(in_img->data);
	for(i=0; i < trans_img->height; i++) { free(trans_img->data[i]);}
	free(trans_img->data);
	*/
	return 0;
}

