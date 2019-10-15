#include "sbr.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <png.h>
#include "ImageIO/image.h"
#include <dirent.h>
#include <sys/stat.h>

#define p(s,a) printf("%s:%d\n",s,a)



int main(int argc, char *argv[])
{
	clock_t start = clock();

	if(argc>4 || argc<2){
		fprintf(stderr, "Usage: program <inputfile> <outputfile>\n");
		exit(1);
	}
	
	image_t *in_img;
	PPM *in_ppm, *in_ppm2;

	
	char *name;
	char *ext;
	
	//入力画像を読み込む
	name = argv[1];
	ext = get_extension(name);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		in_ppm = read_ppm(name);
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		in_img = read_jpeg_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	} else if (strcmp("png", ext) == 0) {
		in_img = read_png_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	} else {
		printf("Plese use JPEG,PNG or PPM!\n");
		exit(1);
	}
	free_image(in_img);
	
	//入力画像を読み込む
	name = argv[2];
	ext = get_extension(name);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		in_ppm2 = read_ppm(name);
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		in_img = read_jpeg_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm2 = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	} else if (strcmp("png", ext) == 0) {
		in_img = read_png_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm2 = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	} else {
		printf("Plese use JPEG,PNG or PPM!\n");
		exit(1);
	}
	free_image(in_img);
	
	double error = image_MSE(in_ppm, in_ppm2);
	printf("MSE:%f\n", error);
	
	FreePPM(in_ppm);
		
	return 0;
}