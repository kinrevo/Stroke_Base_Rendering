#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <png.h>
#include <jpeglib.h>
#include <omp.h>
#include "ImageIO/image.h"
#include "sbr.h"


#define p(s,a) printf("%s:%d\n",s,a)


int main(int argc, char *argv[])
{
	my_clock();

	if(argc<2){
		fprintf(stderr, "Usage: program <inputfile> <outputfile>\n");
		exit(1);
	}

	#ifdef _OPENMP
        omp_set_num_threads(atoi(argv[3]));
    #endif

	image_t *in_img;
	PPM *in_ppm, *trans_ppm;


	char *name;
	char *ext;
	name = argv[1];

	//入力画像を読み込む
	ext = get_extension(name);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		in_ppm = read_ppm(name);
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		in_img = read_jpeg_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
		free_image(in_img);
	} else if (strcmp("png", ext) == 0) {
		in_img = read_png_file(name);
		dump_image_info(in_img);	//画像情報出力
		in_ppm = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
		free_image(in_img);
	} else {
		printf("Plese use JPEG,PNG or PPM!\n");
		exit(1);
	}

	//入力画像にキャニー
	PGM *gray = color_gray_conversion(in_ppm);
	double maxValue=0.30, minValue=0.10;
    int thick_min = 5;
    PGM * canny_pgm = cannyedge_detector(gray, maxValue, minValue, thick_min);

	//出力ファイル名に従って画像を出力
    write_pgm(argv[2],canny_pgm);

	FreePPM(in_ppm);
	FreePGM(gray);
	FreePGM(canny_pgm);


	pd("TOTAL_TIME[s]",my_clock());
	return 0;
}
