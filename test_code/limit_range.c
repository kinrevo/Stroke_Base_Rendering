#include <stdio.h>

#define HELLO "HELLO"
#define arr_num 3
#define arr {1,3,6}

#define LIMIT_RANGE(x,min,max) ((x= (x<min  ? min : x<max ? x : max)))

void limit_range_TEST();
void define_allay_TEST();


void main(){
    limit_range_TEST();
}



void limit_range_TEST(){
    int i,j;

    int int_small=-5;
    int int_big=555;
    double double_small=-2.7;
    double double_big=1.9;

    LIMIT_RANGE(int_small, 0, 255);
    printf("int_small:%d \n", int_small);
    LIMIT_RANGE(int_big, 0, 255);
    printf("int_big:%d \n", int_big);
    LIMIT_RANGE(double_small, 0, 1);
    printf("double_small:%f \n", double_small);
    LIMIT_RANGE(double_big, 0, 1);
    printf("double_big:%f \n", double_big);
}




void define_allay_TEST(){
    int i,j;
    int array[arr_num] = arr;

    int u=-5;
    printf("%d,%d\n", u%10, (-1*u)%10);
    printf(HELLO);
    printf("\n");

    for(i=0; i<10; i++){
        for(j=0; j<3; j++){
                if(i==array[j]) printf("%d \n", i);
        }
    }
}