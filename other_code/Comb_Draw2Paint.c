#include "../sbr.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <png.h>
#include "../ImageIO/image.h"
#include <dirent.h>
#include <sys/stat.h>

#define p(s,a) printf("%s:%d\n",s,a)



int main(int argc, char *argv[])
{
	my_clock();

	if(argc>=6 || argc<=2){
		fprintf(stderr, "Usage: program <PaintFile> <DrawFile> <OutputFile> (NumThread)\n");
		exit(1);
	}

	#ifdef _OPENMP
        omp_set_num_threads(atoi(argv[4]));
    #endif

	image_t *in_img;
	PPM *Paint_ppm, *Draw_ppm;

	
	char *name;
	char *ext;
	
	//塗り画像を読み込む
	name = argv[1];
	ext = get_extension(name);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		Paint_ppm = read_ppm(name);
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		in_img = read_jpeg_file(name);
		dump_image_info(in_img);	//画像情報出力
		Paint_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	    free_image(in_img);
	} else if (strcmp("png", ext) == 0) {
		in_img = read_png_file(name);
		dump_image_info(in_img);	//画像情報出力
		Paint_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	    free_image(in_img);
	} else {
		printf("Plese use JPEG,PNG or PPM!\n");
		exit(1);
	}
	
	//線画画像を読み込む
	name = argv[2];
	ext = get_extension(name);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		Draw_ppm = read_ppm(name);
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		in_img = read_jpeg_file(name);
		dump_image_info(in_img);	//画像情報出力
		Draw_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	    free_image(in_img);
	} else if (strcmp("png", ext) == 0) {
		in_img = read_png_file(name);
		dump_image_info(in_img);	//画像情報出力
		Draw_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
	    free_image(in_img);
	} else {
		printf("Plese use JPEG,PNG or PPM!\n");
		exit(1);
	}
	

    // 塗り画像に線画を重ねる
	combine_LineDrawing2PaintIllust(Paint_ppm, Draw_ppm);
	
	//出力ファイル名に従って画像を出力
	ext = get_extension(argv[3]);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		if(write_ppm(argv[3], Paint_ppm)){ printf("WRITE_PPM_ERROR (main)\n");}
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		if(write_jpeg_file(argv[3], PPM_to_image(Paint_ppm))){ printf("WRITE JPG ERROR.");}
	} else if (strcmp("png", ext) == 0) {
		if(write_png_file(argv[3], PPM_to_image(Paint_ppm))){ printf("WRITE PNG ERROR.");}
	}

	FreePPM(Paint_ppm);
	FreePPM(Draw_ppm);

	pd("TOTAL_TIME[s]",my_clock());
	return 0;
}