#include "sbr.h"
//#include "other_func.h"
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
	PPM *in_ppm, *trans_ppm;

	
	char *name;
	char *ext;
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
	
	//入力画像の絵画化
	trans_ppm = c_Illust_brush(in_ppm, argv[2]);
	
	//if(write_jpeg_file(argv[2], trans_ppm)){ printf("WRITE PNG ERROR.");}
	//if(write_png_file(argv[2], trans_ppm)){ printf("WRITE PNG ERROR."); exit(1); }
	if(write_ppm(argv[2], trans_ppm)){ printf("WRITE_PPM_ERROR (main)\n");};
	
	FreePPM(in_ppm);
	FreePPM(trans_ppm);
		
	// 複数ファイルをまとめて処理するコード
	/*
	int i;
	int file_num;
	char** file_list = make_file_list(&file_num, argv[1], Search_FILE);
	multi_file_CIB(file_list, file_num);

	for(i=0; i < file_num; i++) { free(file_list[i]);}
	free(file_list);
	*/
	
	pd("TOTAL_TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	return 0;
}