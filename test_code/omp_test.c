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


void main(int argc, char *argv[])
{
    int i,j;

    omp_set_num_threads(atoi(argv[1]));

    // #pragma omp parallel for private(i,j)
    for(i=0; i<3; i++){
        printf("使用中のスレッド数「%d」", omp_get_num_threads());
        printf("i:%d \n", i);
        #pragma omp parallel for private(j)
        for(j=0; j<5; j++){
            printf("使用中のスレッド数「%d」", omp_get_num_threads());
            printf("[%d,%d]",i,j);
        }
        printf("\n");
    }

}