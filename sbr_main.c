#include"sbr.h"
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <png.h>
#include "ImageIO/image.h"
#include <dirent.h>
#include <sys/stat.h>

#define p(s,a) printf("%s:%d\n",s,a)

int main(int argc, char *argv[])
{
	if(argc>4 || argc<2){
		fprintf(stderr, "Usage: program <inputfile> <outputfile>\n");
		exit(1);
	}

	clock_t start = clock();
	
	//image_t *in_img;, *out_img;
	
	//コピペソース（ImageIO）
	/*
	char *name;
	char *ext;
	name = argv[1];
    ext = get_extension(name);
    if (strcmp("bmp", ext) == 0) {
      in_img = read_bmp_file(name);
    } else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
      in_img = read_jpeg_file(name);
    } else if (strcmp("png", ext) == 0) {
      in_img = read_png_file(name);
    } else if (strcmp("ppm", ext) == 0 || strcmp("pbm", ext) == 0
        || strcmp("pgm", ext) == 0) {
      in_img = read_pnm_file(name);
    } */
	
	// 複数ファイルをまとめて処理するコード
	int i;
	int file_num;
	char** file_list = make_file_list(&file_num, argv[1], Search_FILE);
	multi_file_CIB(file_list, file_num);

	for(i=0; i < file_num; i++) { free(file_list[i]);}
	free(file_list);
	
	//multi_dir_CIB(argv[1]);
	
	
	//if(write_png_file(argv[2], out_img)){ printf("WRITE PNG ERROR."); exit(1); }
	pd("TOTAL_TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	return 0;
}



// 文字列から拡張子を取得
static char *get_extension(char *name) {
  int i;
  for (i = strlen(name) - 1; i >= 0; i--) {
    if (name[i] == '.') {
      return &name[i + 1];
    }
  }
  return NULL;
}
