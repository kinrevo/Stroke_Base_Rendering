#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include "sbr.h"


double Test_RGB2Lab(int argc, char *argv[]);
double Test_Kmeans_ImageLab3D(int argc, char *argv[]);
void Kmeans(double *x, int data_num, int dimension_num, double *cen, int cluster_num, int *cl, int maxiter, int *nc, double *wss);


// RGB->lRGB->CIExyz->CIELab  https://qiita.com/hachisukansw/items/09caabe6bec46a2a0858
Lab RGB2Lab(RGB in) {
    Lab out;
    double r = in.R / 255.0;
    double g = in.G / 255.0;
    double b = in.B / 255.0;

    double lr = r > 0.04045 ? pow(((r + 0.055) / 1.055), 2.4) : (r / 12.92);
    double lg = g > 0.04045 ? pow(((g + 0.055) / 1.055), 2.4) : (g / 12.92);
    double lb = b > 0.04045 ? pow(((b + 0.055) / 1.055), 2.4) : (b / 12.92);

    double x = (lr * 0.4124) + (lg * 0.3576) + (lb * 0.1805);
    double y = (lr * 0.2126) + (lg * 0.7152) + (lb * 0.0722);
    double z = (lr * 0.0193) + (lg * 0.1192) + (lb * 0.9505);

    x = x * 100 / 95.047;
    y = y * 100 / 100;
    z = z * 100 / 108.883;

    x = x > 0.008856 ? pow(x, 1 / 3.0) : (7.787 * x) + (4 / 29.0);
    y = y > 0.008856 ? pow(y, 1 / 3.0) : (7.787 * y) + (4 / 29.0);
    z = z > 0.008856 ? pow(z, 1 / 3.0) : (7.787 * z) + (4 / 29.0);

    out.L = (116 * y) - 16;
    out.a = 500 * (x - y);
    out.b = 200 * (y - z);

    return out;
}


// CIELab->CIExyz->lRGB->RGB  http://story-monoroch.blogspot.com/2014/09/2rgbcielab-rgbcielab-rgb-rgb-xyz-cielab.html
RGB Lab2RGB(Lab in) {
    RGB out;
    double yn = (in.L + 16.0) / 116.0;
    double xn = yn + in.a / 500.0;
    double zn = yn - in.b / 200.0;
    
    double fxn = xn > 6.0/29 ? pow(xn, 3) : 3*pow(6.0/29, 2.0) * (xn - 4.0/29);
    double fyn = yn > 6.0/29 ? pow(yn, 3) : 3*pow(6.0/29, 2.0) * (yn - 4.0/29);
    double fzn = zn > 6.0/29 ? pow(zn, 3) : 3*pow(6.0/29, 2.0) * (zn - 4.0/29);
    
    double x = fxn * 0.95047;
    double y = fyn * 1.00000;
    double z = fzn * 1.08883;
    
    double lr =   3.24062547732005470 * x - 1.53720797221031910 * y - 0.49862859869824794 * z;
    double lg = - 0.96893071472931970 * x + 1.87575606088524200 * y + 0.04151752384295397 * z;
    double lb =   0.05571012044551063 * x - 0.20402105059848677 * y + 1.05699594225438860 * z;
    
    double r = lr <= 0.0031308 ? lr * 12.92 : pow(lr, 1 / 2.4) * 1.055 - 0.055;
    double g = lg <= 0.0031308 ? lg * 12.92 : pow(lg, 1 / 2.4) * 1.055 - 0.055;
    double b = lb <= 0.0031308 ? lb * 12.92 : pow(lb, 1 / 2.4) * 1.055 - 0.055;
    
    out.R = r * 255;
    out.G = g * 255;
    out.B = b * 255;

    return out;
}



RGB* Kmeans_ImageLab3D(PPM* img, int cluster_num, int maxiter, int* x_centlabel, int* num_cluster)
{
    int i,j,data;
    double maxL=-999999, minL=999999,maxa=-999999, mina=999999, maxb=-999999,minb=999999;
    RGB* color_set = (RGB *)malloc(sizeof(RGB) * cluster_num);
    RGB rgb;
    Lab lab;
    int data_num =img->width*img->height;
    int dimension_num = 3;
    double* x = (double*)malloc(sizeof(double) * data_num*dimension_num);
    double* centroid = (double*)malloc(sizeof(double) * cluster_num*dimension_num);
    // x_centlabel = (int*)malloc(sizeof(int) * data_num);
    // num_cluster = (int*)malloc(sizeof(int) * cluster_num);
    double* cluster_weight = (double*)malloc(sizeof(double) * cluster_num);

    // 二次元画像のab２つのデータを一次元に変換
    for(i=0; i<img->width; i++){
        for(j=0; j<img->height; j++){
            data = i*img->height + j;   //データのINDEX
            rgb.R=img->dataR[i][j];
            rgb.G=img->dataG[i][j];
            rgb.B=img->dataB[i][j];
            lab = RGB2Lab(rgb);

            x[data + data_num*0] = lab.L;
            x[data + data_num*1] = lab.a;
            x[data + data_num*2] = lab.b;

            maxL = fmax(maxa, lab.L);
            maxa = fmax(maxa, lab.a);
            maxb = fmax(maxb, lab.b);
            minL = fmin(mina, lab.L);
            mina = fmin(mina, lab.a);
            minb = fmin(minb, lab.b);
        }
    }

    // 適当に初期セントロイド位置を決定
    for (i = 0; i < cluster_num; i++){
        centroid[i + cluster_num*0] = (maxL-minL)/(cluster_num+1)*(i+1) + minL;
        centroid[i + cluster_num*1] = (maxa-mina)/(cluster_num+1)*(i+1) + mina;
        centroid[i + cluster_num*2] = (maxb-minb)/(cluster_num+1)*(i+1) + minb;
    }

    Kmeans(x, data_num, dimension_num, centroid, cluster_num, x_centlabel, maxiter, num_cluster, cluster_weight);

    for (i = 0; i<cluster_num; i++) {
        lab.L = centroid[i + cluster_num*0];
        lab.a = centroid[i + cluster_num*1];
        lab.b = centroid[i + cluster_num*2];
        rgb = Lab2RGB(lab);
        color_set[i] = rgb;
        printf("ClusterWeight[%d]:%f \n", i, cluster_weight[i]);
        printf("NumberCluster[%d]:%d \n", i, num_cluster[i]);
        printf("l:%f,a:%f,b:%f \n", lab.L, lab.a, lab.b);
        printf("r:%d,g:%d,b:%d \n", rgb.R, rgb.G, rgb.B);
        pn;
    }

    free(x);
    free(centroid);
    // free(x_centlabel);
    // free(num_cluster);
    free(cluster_weight);

    return color_set;
}



void Kmeans(double *x, int data_num, int dimension_num, double *cen, int cluster_num, int *cl, int maxiter, int *nc, double *wss)
{
    // int data_num = *pn, cluster_num = *pk, dimension_num = *pp, maxiter = *pmaxiter;
    int iter, d, c, dim, it, inew = 0;
    double best, dd, tmp;
    int updated;
    /* cl[d]はd番目のデータが第何クラスタに所属しているかを表わす */
    for(d = 0; d < data_num; d++) cl[d] = -1;

    for(iter = 0; iter < maxiter; iter++) {
        // printf("iteration:%d\n",iter+1);
        updated = FALSE;
        /* データn個に対するループ */
        /* それぞれの点で、最も近いクラスタの中心(重心)を見つける */
        for(d = 0; d < data_num; d++) {
            best = 99999999;
            /* クラスタ数に対するループ */
            for(c = 0; c < cluster_num; c++) {
                dd = 0.0;
                /* 次元に対するループ。点x_iにおける、各次元でのj番目のクラスタの中心との二乗和を取っている */
                for(dim = 0; dim < dimension_num; dim++) {
                    tmp = x[d+data_num*dim] - cen[c+cluster_num*dim];
                    dd += tmp * tmp;
                }
                /* x_iと一番近いクラスタとの距離に更新 */
                if(dd < best) {
                    best = dd;
                    inew = c+1;
                    // printf(" %3d番目のデータ(",d+1);
                    // for(dim = 0; dim < dimension_num; dim++) {
                    //     if(dim < dimension_num - 1){
                    //         printf("% .2f, ",x[d+data_num*dim]);
                    //     }else{
                    //         printf("% .2f",x[d+data_num*dim]);
                    //     }
                    // }
                    // printf(")がクラスタ%dに移動しました。\n",c+1);
                }
            }
            if(cl[d] != inew) {
                updated = TRUE;
                cl[d] = inew;
            }
        }
        if(!updated) break;

        /* 各クラスタの中心の点を更新 */
        for(c = 0; c < cluster_num*dimension_num; c++) cen[c] = 0.0;
        for(c = 0; c < cluster_num; c++) nc[c] = 0;
        for(d = 0; d < data_num; d++) {
            it = cl[d] - 1; nc[it]++;
            for(dim = 0; dim < dimension_num; dim++) cen[it+dim*cluster_num] += x[d+dim*data_num];
        }
        for(c = 0; c < cluster_num*dimension_num; c++) cen[c] /= nc[c % cluster_num];
    }

    /* 各クラスタのセントロイドと属する点のの距離の二乗総和wss */
    for(c = 0; c < cluster_num; c++) wss[c] = 0.0;
    for(d = 0; d < data_num; d++) {
        it = cl[d] - 1;
        for(dim = 0; dim < dimension_num; dim++) {
            tmp = x[d+data_num*dim] - cen[it+cluster_num*dim];
            wss[it] += tmp * tmp;
        }
    }
}


//image形式データ（png用）をPPM形式データ（自定義）に変換
PPM* Visualize_ColorSet(RGB* color_set, int cluster_num, int* num_cluster) {
	int c,i,j;

	PPM *ppm = (PPM *)malloc(sizeof(PPM));
	memcpy(ppm->descriptor, "P3", 3);
	ppm->height = 64;
	ppm->width = 512;
	ppm->bright = 255;
	ppm->dataR = create_ally(ppm->width, ppm->height);
	ppm->dataG = create_ally(ppm->width, ppm->height);
	ppm->dataB = create_ally(ppm->width, ppm->height);
    format_ally(ppm->dataR, ppm->width, ppm->height, ppm->bright);
    format_ally(ppm->dataG, ppm->width, ppm->height, ppm->bright);
    format_ally(ppm->dataB, ppm->width, ppm->height, ppm->bright);


    int sum_cn = 0;
    for (i = 0; i < cluster_num; i++) {
        sum_cn += num_cluster[i];
    }

    int paint_sum=0;
    int paint_num=0;
    // クラスタの割合分、左から塗りつぶしていく
    for (c = 0; c < cluster_num; c++) 
    {
        paint_num = ppm->width * (double)num_cluster[c]/sum_cn + 0.5;
        for(i=paint_sum; i<paint_sum + paint_num; i++) {
            if(i >= ppm->width) break;
            for(j=0; j<ppm->height; j++) {
                ppm->dataR[i][j] = color_set[c].R;
                ppm->dataG[i][j] = color_set[c].G;
                ppm->dataB[i][j] = color_set[c].B;
            }
        }
        paint_sum += paint_num;
    }

	return ppm;
}



int main(int argc, char *argv[])
{
    my_clock();

    Test_Kmeans_ImageLab3D(argc, argv);

    return 0;
}


double Test_Kmeans_ImageLab3D(int argc, char *argv[])
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

    int cluster_num = 28;
    int *x_centlabel = (int*)malloc(sizeof(int) * in_ppm->width*in_ppm->height);
    int *num_cluster = (int*)malloc(sizeof(int) * cluster_num);
    RGB* color_set = Kmeans_ImageLab3D(in_ppm, cluster_num, 50, x_centlabel, num_cluster);
    for (int i = 0; i < 4; i++)
    {
        printf("NumberCluster[%d]:%d \n", i, num_cluster[i]);
        printf("r:%d,g:%d,b:%d \n", color_set[i].R, color_set[i].G, color_set[i].B);
    }
    

	//入力画像の絵画化
    trans_ppm = Visualize_ColorSet(color_set, cluster_num, num_cluster);


	//出力ファイル名に従って画像を出力
	ext = get_extension(argv[2]);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		if(write_ppm(argv[2], trans_ppm)){ printf("WRITE_PPM_ERROR (main)\n");}
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		if(write_jpeg_file(argv[2], PPM_to_image(trans_ppm))){ printf("WRITE JPG ERROR.");}
	} else if (strcmp("png", ext) == 0) {
		if(write_png_file(argv[2], PPM_to_image(trans_ppm))){ printf("WRITE PNG ERROR.");}
	}

	FreePPM(in_ppm);
	FreePPM(trans_ppm);

	pd("TOTAL_TIME[s]",my_clock());
	return my_clock();
}



double Test_RGB2Lab(int argc, char *argv[])
{
	double start = my_clock();

    RGB in = {68,89,59};
    printf("%d,%d,%d \n",in.R, in.G, in.B);

    Lab lab = RGB2Lab(in);
    printf("%f,%f,%f \n",lab.L, lab.a, lab.b);

    RGB rgb = Lab2RGB(lab);
    printf("%d,%d,%d \n",rgb.R, rgb.G, rgb.B);
    
    return my_clock()-start;
}

