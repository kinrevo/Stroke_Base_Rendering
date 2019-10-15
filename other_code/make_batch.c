#include<stdio.h>
#include<string.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdlib.h>

#define Search_DIRECTRY 0
#define Search_FILE 1

char** make_file_list(int* file_num, char* search_path, int Search_TYPE);


int main(int argc, char *argv[]){
	
	int i;
	int dir_num;
	char** dir_list = make_file_list(&dir_num, argv[1], Search_DIRECTRY);
	
	//前に生成したバッチファイルをあらかじめ削除
	if( remove( "make_dataset_batch.sh" )){
		fprintf(stderr, "%sRemove_Error", argv[1]);
	}
	
	FILE *fp;
	
	if((fp = fopen("make_dataset_batch.sh", "w")) == NULL){
		fprintf(stderr, "W_Open_Error");
		return -1;
	}
	
	//読み込んだファイル名から順に実行コマンドを記述しバッチファイルを生成
	for(i=0; i<dir_num; i++) {
		fprintf(fp, "./multi_file_CIB ");
		fprintf(fp, "%s\n", dir_list[i]);
	}
		
		
	fclose(fp);
	return 0;
	
}



//与えられたパスのディレクトリ中にあるファイルもしくはディレクトリをリスト化し配列で返す
//引数に従ってどちらかしか検索できない
char** make_file_list(int* file_num, char* search_path, int Search_TYPE){
	int i=0;
	DIR *dir;
	char path[256];
    struct dirent *dp;
    struct stat fi;
	*file_num=0;
	
	dir = opendir(search_path);
	printf("Search_TYPE:%d\n",Search_TYPE);
    for (dp = readdir(dir); dp != NULL; dp = readdir(dir)) {
        if (dp->d_name[0] == '.') continue;              
        
	    strcpy(path, search_path);
	    strcat(path, "/");
	    strcat(path, dp->d_name);
        stat(path, &fi);	//pathのファイルの情報を取得
    	
    	switch(Search_TYPE){
	    	case Search_FILE:
				if (!S_ISDIR(fi.st_mode)) {
					i++;
		            printf("FILE[%d]:%s\n",i, path);
		        }
    			break;
	    	case Search_DIRECTRY:
		        if (S_ISDIR(fi.st_mode)) {
					i++;
		            printf("DIR[%d]:%s\n",i, path);
		        }
    			break;
    	}
    }
	*file_num = i;
	printf("file_num:%d\n", *file_num);
	
	char** list = (char**)malloc(sizeof(char*)*(*file_num));
	for(i=0; i<*file_num; i++) { list[i]=(char*)malloc(sizeof(char)*(256));}
	
	i=0;
	rewinddir(dir);
	for (dp = readdir(dir); dp != NULL; dp = readdir(dir)) {
        if (dp->d_name[0] == '.') continue;              
        
	    strcpy(path, search_path);
	    strcat(path, "/");
	    strcat(path, dp->d_name);
        stat(path, &fi);	//pathのファイルの情報を取得
    	
    	switch(Search_TYPE){
	    	case Search_FILE:
				if (!S_ISDIR(fi.st_mode)) {
					strcpy(list[i], path);
					i++;
		            printf("FILE[%d]:%s\n",i, path);
		        }
    			break;
	    	case Search_DIRECTRY:
		        if (S_ISDIR(fi.st_mode)) {
					strcpy(list[i], path);
					i++;
		            printf("DIR[%d]:%s\n",i, path);
		        }
    			break;
    	}
    }
	
    closedir(dir);
	return list;//file_num;
}