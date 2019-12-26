#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <omp.h>
#include "sbr.h"
#include "sbr_opt.h"
#include "water.h"

#define p(s,a) printf("%s:%d\n",s,a)

PPM *gpu_c_Illust_brush_Water_best(PPM *in, char *filename);

PGM* GLOBAL_improved_value_map;

//GPUのエラーを検出する関数
void printCudaLastError(){
	cudaError_t err = cudaGetLastError();
	printf("cudaGetLastError::%s(code:%d)\n",cudaGetErrorString(err),err);
	if(err)	exit(0);
}

int main(int argc, char *argv[])
{

	//実行時間計測用
	float total_time = 0;
	cudaEvent_t total_time_start, total_time_end;
	cudaEventCreate(&total_time_start);
	cudaEventCreate(&total_time_end);

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
	/*
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
	*/
	in_ppm = read_ppm(name);

	cudaEventRecord(total_time_start, 0);//計測スタート

	//入力画像の絵画化
	trans_ppm = gpu_c_Illust_brush_Water_best(in_ppm, argv[2]);

	cudaEventRecord(total_time_end, 0);//計測ストップ

	//実行時間表示
	cudaEventElapsedTime(&total_time, total_time_start, total_time_end);
	printf("Total_Time:%f[ms]\n", total_time);

	//出力ファイル名に従って画像を出力
	/*
	ext = get_extension(argv[2]);
	if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
		if(write_ppm(argv[2], trans_ppm)){ printf("WRITE_PPM_ERROR (main)\n");}
	} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
		if(write_jpeg_file(argv[2], PPM_to_image(trans_ppm))){ printf("WRITE JPG ERROR.");}
	} else if (strcmp("png", ext) == 0) {
		if(write_png_file(argv[2], PPM_to_image(trans_ppm))){ printf("WRITE PNG ERROR.");}
	}
	*/
	write_ppm(argv[2], trans_ppm);

	FreePPM(in_ppm);
	FreePPM(trans_ppm);

	pd("TOTAL_TIME[s]",my_clock());
	return 0;
}

//与えられたイメージからブラッシングによる絵画を作る
PPM *gpu_c_Illust_brush_Water_best(PPM *in, char *filename)
{

	//実行時間計測用
	float ftimer = 0;
	cudaEvent_t start, end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);

	my_clock();
	int i,j,x,y,xc,yc,t=1000,break_flag,pnum, offscrn_count;
	int window_diff_border = opt_window_diff_border; 	//ストローク位置探索のしきい値
	int color_diff_border = opt_color_diff_border;  	//描画色の差異のしきい値
	int max_stroke = opt_max_stroke;
	int min_stroke = opt_min_stroke;
	Point p[max_stroke];
	int stroke_histogram[max_stroke+1];
	for(i=0; i<max_stroke+1; i++){stroke_histogram[i]=0;}
	double ratio=opt_ratio;		//ストロークの濃度
	double theta, former_theta;
	double sigma, diff_sum, sum;
	int histogram_partition=opt_histogram_partition;
	int paint_count=0, nc=0, tc=-1;
	int loop_cont=opt_loop_cont, x_defo=0, y_defo=0;
	// int lc=0;
	// double maxValue, minValue;
	RGB bright;


	//最大小ストローク半径（自動化：画面の1/10,最大の1/10）
	int thick_max = opt_thick_max;//(in->height < in->width ? in->height : in->width)/10;
	int thick_min = opt_thick_min;//(thick_max/15 > 3 ? thick_max/15 : 3);


	//出力ファイル名のサイズを取得
	int namesize = strlen(filename)*2 + 16;
	char out_filename[namesize], log_filename[namesize], dir_path[namesize];
	char in_filename[namesize-16];
	char count_name[16];
	char log_sentence[2028] = "";
	image_t *out_png;
	//パスから入力ファイル名（拡張子含まない）を取得
	char tmp[namesize];
	char* tp;
	strcpy(in_filename, filename);
	strcpy(tmp, filename);
	strtok(tmp, "/\\");
	while((tp = strtok(NULL, "/\\")) != NULL ) {
		ps("in_filename",in_filename);
		ps("tmp",tmp);
		ps("tp",tp);
		strcpy(in_filename, tp);
	}
	strtok(in_filename, ".");
	ps("in_filename",in_filename);

	//出力するフォルダを生成しフォルダへのパスを格納、ログファイルを作成
	strcpy(tmp, filename);
	strcpy(dir_path, strtok(tmp, "."));
	if(mkdir(tmp, 0775)){ printf("FAIL TO CREATE DIRECTRY\n"); }
	strcat(dir_path, "/");
	strcpy(log_filename, dir_path);
	strcat(log_filename, in_filename);
	strcat(log_filename, ".log");

	//logデータ格納
	strcat(log_sentence, "<");
	strcat(log_sentence, in_filename);
	strcat(log_sentence, ">\r\n");
	strcat(log_sentence, "Stroke_Method:BestStroke\r\n");
	Add_dictionary_to_sentence(log_sentence, "width", in->width);
	Add_dictionary_to_sentence(log_sentence, "height", in->height);
	Add_dictionary_to_sentence(log_sentence, "thick_max", thick_max);
	Add_dictionary_to_sentence(log_sentence, "thick_min", thick_min);
	Add_dictionary_to_sentence(log_sentence, "max_stroke", max_stroke);
	Add_dictionary_to_sentence(log_sentence, "min_stroke", min_stroke);
	Add_dictionary_to_sentence(log_sentence, "window_diff_border", window_diff_border);
	Add_dictionary_to_sentence(log_sentence, "color_diff_border", color_diff_border);
	Add_dictionary_to_sentence_d(log_sentence, "ratio", ratio);
	Add_dictionary_to_sentence(log_sentence, "histogram_partition", histogram_partition);
	Add_dictionary_to_sentence(log_sentence, "loop_cont", loop_cont);
	Add_dictionary_to_sentence(log_sentence, "USE_Lab_ColorDiff", opt_USE_Lab_ColorDiff);
	Add_dictionary_to_sentence(log_sentence, "USE_calcu_color_bi", opt_USE_calcu_color_bi);
	Add_dictionary_to_sentence(log_sentence, "USE_gause_histogram", opt_USE_gause_histogram);
	Add_dictionary_to_sentence(log_sentence, "optimal_improved_value_border", opt_optimal_improved_value_border);
	Add_dictionary_to_sentence_d(log_sentence, "StrokeWindowStep", opt_StrokeWindowStep/t);
	if(opt_USE_calcu_Kmean_ColorSet){
		Add_dictionary_to_sentence(log_sentence, "Kmean_ClusterNum", opt_Kmean_ClusterNum);
	}else if(opt_USE_calcu_JIS_ColorSet){
		Add_dictionary_to_sentence(log_sentence, "JIS_ClusterNum", opt_JIS_ClusterNum);
	}
	if(opt2_thick_max){
		Add_dictionary_to_sentence(log_sentence, "thick_max[2]", opt2_thick_max);
		Add_dictionary_to_sentence(log_sentence, "thick_min[2]", opt2_thick_min);
		Add_dictionary_to_sentence(log_sentence, "min_stroke[2]", opt2_min_stroke);
		Add_dictionary_to_sentence_d(log_sentence, "ratio[2]", opt2_ratio);
		Add_dictionary_to_sentence(log_sentence, "loop_cont[2]", opt2_loop_cont);
	}
	strcat(log_sentence, "\r\n[Water Option]\r\n");
	Add_dictionary_to_sentence_d(log_sentence, "mhu", opt_mhu);
	Add_dictionary_to_sentence_d(log_sentence, "kappa", opt_kappa);
	Add_dictionary_to_sentence(log_sentence, "N", opt_N);
	Add_dictionary_to_sentence_d(log_sentence, "tau", opt_tau);
	Add_dictionary_to_sentence_d(log_sentence, "xi", opt_xi);
	Add_dictionary_to_sentence(log_sentence, "K", opt_K);
	Add_dictionary_to_sentence_d(log_sentence, "eta", opt_eta);
	if(opt_USE_DETAIL_TP){
		Add_dictionary_to_sentence_d(log_sentence, "deposit", opt_deposit);
		Add_dictionary_to_sentence_d(log_sentence, "lift", opt_lift);
		Add_dictionary_to_sentence_d(log_sentence, "exposure", opt_exposure);
	}else{
		Add_dictionary_to_sentence_d(log_sentence, "gamma", opt_gamma);
		Add_dictionary_to_sentence_d(log_sentence, "rho", opt_rho);
		Add_dictionary_to_sentence_d(log_sentence, "omega", opt_omega);
	}
	Add_dictionary_to_sentence(log_sentence, "SoakTme", opt_SoakTime);
	Add_dictionary_to_sentence_d(log_sentence, "SoakTimeStep", opt_SoakTimeStep);
	Add_dictionary_to_sentence_d(log_sentence, "perlin_freq", opt_perlin_freq);
	Add_dictionary_to_sentence(log_sentence, "perlin_depth", opt_perlin_depth);
	Add_dictionary_to_sentence_d(log_sentence, "variance_ratio", opt_variance_ratio);
	if(opt_USE_Backrun){
		Add_dictionary_to_sentence_d(log_sentence, "alpha", opt_alpha);
		Add_dictionary_to_sentence_d(log_sentence, "epsilon", opt_epsilon);
		Add_dictionary_to_sentence_d(log_sentence, "delta", opt_delta);
		Add_dictionary_to_sentence_d(log_sentence, "sigma", opt_sigma);
	}
	Add_dictionary_to_sentence(log_sentence, "RemovePigmentInWater", opt_RemovePigmentInWater);
	Add_dictionary_to_sentence(log_sentence, "FloatPigmentOnPaper", opt_FloatPigmentOnPaper);
	#ifdef _OPENMP
		Add_dictionary_to_sentence(log_sentence, "OMP_NUM_THREADS", omp_get_max_threads());
		printf("OMP_NUM_T:%d \n",(int)omp_get_max_threads());
	#endif
	Add_dictionary_to_sentence(log_sentence, "opt_GPU_Block_Num", opt_GPU_Block_Num_x*opt_GPU_Block_Num_y);
	Add_dictionary_to_sentence(log_sentence, "opt_GPU_Thread_Num", opt_GPU_Thread_Num_x*opt_GPU_Thread_Num_y);

	pn;

	//カラー画像分割
	PGM *gray = color_gray_conversion(in);
	PGM *nimgR = create_pgm(gray->width, gray->height, gray->bright);
	PGM *nimgG = create_pgm(gray->width, gray->height, gray->bright);
	PGM *nimgB = create_pgm(gray->width, gray->height, gray->bright);
	PGM *inR = create_pgm(gray->width, gray->height, gray->bright);
	PGM *inG = create_pgm(gray->width, gray->height, gray->bright);
	PGM *inB = create_pgm(gray->width, gray->height, gray->bright);
	devide_ppm(in, inR, inG, inB);
	//比較用のイメージ生成
	PGM *cmpr  = gray;
	PGM *cmprR = inR;
	PGM *cmprG = inG;
	PGM *cmprB = inB;
	//キャンバスイメージ生成
	PGM *nimgV = create_pgm(gray->width, gray->height, gray->bright); //明度のみのキャンバス（比較用）
	PPM *nimgC;
	if(opt_USE_input_progress_image){
		//描画中キャンバス画像を読み込む
		image_t *in_img;
		char name[128] = opt_progress_image_address;
		char* ext = get_extension(name);

		//if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
			nimgC = read_ppm(name);
		//}
		/* else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
			in_img = read_jpeg_file(name);
			dump_image_info(in_img);	//画像情報出力
			nimgC = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
			free_image(in_img);
		} else if (strcmp("png", ext) == 0) {
			in_img = read_png_file(name);
			dump_image_info(in_img);	//画像情報出力
			nimgC = image_to_PPM(in_img);		//扱いやすいデータ構造に変換
			free_image(in_img);
		} else {
			printf("Plese use JPEG,PNG or PPM!\n");
			exit(1);
		}
		*/

		devide_ppm(nimgC, nimgR, nimgG, nimgB);
	} else{
		nimgC = create_ppm(in->width, in->height, in->bright); //実際に描画するキャンバス
	}
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;

	PPM *test_Canvas;
	// #ifndef _OPENMP
		test_Canvas = create_ppm(in->width, in->height, in->bright);
	// #endif
	// PGM *improved_value_map = create_pgm(in->width, in->height, in->bright); //改善値マップ
	GLOBAL_improved_value_map = create_pgm(in->width, in->height, in->bright); //改善値マップ


	// sobelフィルタを適応した計算結果を予め格納しておく
	double **sobel_abs = create_dally(in->width, in->height);
	double **sobel_angle = create_dally(in->width, in->height);
	sobel_calcu(gray, sobel_abs, sobel_angle);

	//Greedyアプローチ
	int s_count,stroke_num;
	// int best_pnum;
	int diff_stroke_max[1];
	int tmp_num;
	int diff_stroke_max_ave[50]={};
	int best_x=0, best_y=0;
	// int miss_stroke_count=0;
	Stroke*** best_stroke_map = create_Stroke_ally(in->width, in->height, max_stroke);
	Point best_P;
	int optimal_improved_value_border = opt_optimal_improved_value_border;
	// Point before_P={0,0};

    double** h = perlin_img(in->width, in->height, opt_perlin_freq, opt_perlin_depth);
    double** grad_hx = create_dally(in->width+1, in->height);
    double** grad_hy = create_dally(in->width, in->height+1);
    calcu_grad_h(h, grad_hx, grad_hy, in->width, in->height);

	// KmeanカラーセットまたはJISカラーセット
	int CentLabel = 0, JISLabel = 0;
	int *num_cluster, *x_centlabel;
	int** x_centlabel_2D;
	RGB* ColorSet;
	if(opt_USE_calcu_Kmean_ColorSet){
		x_centlabel = (int*)malloc(sizeof(int) * in->width*in->height);
		num_cluster = (int*)malloc(sizeof(int) * opt_Kmean_ClusterNum);
		ColorSet = Kmeans_ImageLab3D(in, opt_Kmean_ClusterNum, 50, x_centlabel, num_cluster);
		x_centlabel_2D = ReshapeInt_1to2(x_centlabel, in->width, in->height);	//クラスタ番号は１から数えている
		free(x_centlabel);
		PPM* Kmean_Img = Visualize_KmeanImg(in, ColorSet, x_centlabel_2D);
		PPM* ColorSet_Img = Visualize_ColorSet(ColorSet, opt_Kmean_ClusterNum, num_cluster);
		strcpy(out_filename, dir_path);
		strcat(out_filename, in_filename);
		strcat(out_filename, "__Kmean");
		strcat(out_filename, ".png");
		if(write_png_file(out_filename, PPM_to_image(Kmean_Img))){ printf("WRITE PNG ERROR.");}
		strcpy(out_filename, dir_path);
		strcat(out_filename, in_filename);
		strcat(out_filename, "__ColorSet");
		strcat(out_filename, ".png");
		if(write_png_file(out_filename, PPM_to_image(ColorSet_Img))){ printf("WRITE PNG ERROR.");}
	}
	else if(opt_USE_calcu_JIS_ColorSet){
		ColorSet = create_JIS_ColorSet(opt_JIS_ClusterNum);
	}

	// Lab誤差
	Lab** in_Lab;
	float *in_Lab_L;
	float *in_Lab_a;
	float *in_Lab_b;
	if(opt_USE_Lab_ColorDiff){
		//　RGB入力画像をLabに変換した配列を用意
		RGB CanRGB;
		in_Lab = (Lab**)malloc(sizeof(Lab*)*(in->width));
		for(i=0; i<in->width; i++){
			in_Lab[i] = (Lab*)malloc(sizeof(Lab)*(in->height));
		}

		for(i=0; i<in->width; i++){
			for(j=0; j<in->height; j++){
				CanRGB.R = in->dataR[i][j];
				CanRGB.G = in->dataG[i][j];
				CanRGB.B = in->dataB[i][j];
				in_Lab[i][j] = RGB2Lab(CanRGB);
			}
		}

		in_Lab_L = (float*)malloc(sizeof(float)*in->width*in->height);
		in_Lab_a = (float*)malloc(sizeof(float)*in->width*in->height);
		in_Lab_b = (float*)malloc(sizeof(float)*in->width*in->height);

		for(i=0; i<in->width; i++){
			for(j=0; j<in->height; j++){
				in_Lab_L[i+j*in->width] = in_Lab[i][j].L;
				in_Lab_a[i+j*in->width] = in_Lab[i][j].a;
				in_Lab_b[i+j*in->width] = in_Lab[i][j].b;
			}
		}
	}

	//ガウスフィルタを生成
	double** gauce_filter;
    int w = (int)(ceil(3.0*opt_K/6.0+0.5)*2-1); //とりあえず動く計算
    int c=(w-1)/2;
    gauce_filter = create_dally(w, w);
	for(i=0;i<w;i++){
		for(j=0;j<w;j++){
			gauce_filter[i][j] = gause_func(i-c, j-c, opt_K/6.0);
			//printf("%.2f ",gauce_filter[i][j]);
		}
		//printf("\n");
	}

	//GPUの領域確保
	int stroke_length_max = opt_thick_max*(opt_max_stroke+2); //ストローク長さの最大値（本来はopt_thick_max*(opt_max_stroke+1)だが余裕を持たせている）
	int *dev_GLOBAL_improved_value_map;
	float *dev_PerlinNoise;
	float *dev_in_Lab_L;
	float *dev_in_Lab_a;
	float *dev_in_Lab_b;
	int *dev_cmpr;
	int *dev_cmprR;
	int *dev_cmprG;
	int *dev_cmprB;
	int *dev_nimgR;
	int *dev_nimgG;
	int *dev_nimgB;
	int *dev_best_stroke_map_pnum;
	float *dev_best_stroke_map_point_x;
	float *dev_best_stroke_map_point_y;
	int *dev_best_stroke_map_R;
	int *dev_best_stroke_map_G;
	int *dev_best_stroke_map_B;
	float *dev_gauss_filter;
	float *dev_sobel_abs;
	float *dev_sobel_angle;
	float *dev_grad_hx;
    float *dev_grad_hy;
	char *dev_M;
    float *dev_u;
    float *dev_new_u;
	float *dev_v;
	float *dev_new_v;
	float *dev_p;
	float *dev_gR;
	float *dev_gG;
	float *dev_gB;
	float *dev_dR;
	float *dev_dG;
	float *dev_dB;
	float *dev_new_gR;
	float *dev_new_gG;
	float *dev_new_gB;
	float *dev_new_dR;
	float *dev_new_dG;
	float *dev_new_dB;
	float *dev_gauss_M;
	float *dev_h;
	int *dev_best_x;
	int *dev_best_y;
	int *dev_diff_stroke_max;

	cudaMalloc(&dev_GLOBAL_improved_value_map , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_PerlinNoise , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_in_Lab_L , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_in_Lab_a , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_in_Lab_b , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_cmpr , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_cmprR , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_cmprG , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_cmprB , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_nimgR , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_nimgG , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_nimgB , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_best_stroke_map_pnum , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_best_stroke_map_point_x , in->width*in->height*max_stroke*sizeof(float));
	cudaMalloc(&dev_best_stroke_map_point_y , in->width*in->height*max_stroke*sizeof(float));
	cudaMalloc(&dev_best_stroke_map_R , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_best_stroke_map_G , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_best_stroke_map_B , in->width*in->height*sizeof(int));
	cudaMalloc(&dev_sobel_abs , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_sobel_angle , in->width*in->height*sizeof(float));
	cudaMalloc(&dev_grad_hx , (in->width+1)*in->height*sizeof(float));
	cudaMalloc(&dev_grad_hy , in->width*(in->height+1)*sizeof(float));
	cudaMalloc(&dev_M , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(char));
	cudaMalloc(&dev_u , (stroke_length_max+1)*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_u , (stroke_length_max+1)*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_v , stroke_length_max*(stroke_length_max+1)*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_v , stroke_length_max*(stroke_length_max+1)*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_p , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_gR , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_gG , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_gB , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_dR , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_dG , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_dB , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_gR , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_gG , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_gB , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_dR , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_dG , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_new_dB , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_gauss_filter , w*w*sizeof(float));
	cudaMalloc(&dev_gauss_M , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_h , stroke_length_max*stroke_length_max*opt_GPU_Block_Num_x*opt_GPU_Block_Num_y*sizeof(float));
	cudaMalloc(&dev_best_x , 1*sizeof(int));
	cudaMalloc(&dev_best_y , 1*sizeof(int));
	cudaMalloc(&dev_diff_stroke_max , 1*sizeof(int));

	float *PerlinNoise_1dim = ReshapeDouble_2to1(h, in->height, in->width);
	int *cmpr_1dim = ReshapeInt_2to1(cmpr->data, in->height, in->width);
	int *cmprR_1dim = ReshapeInt_2to1(cmprR->data, in->height, in->width);
	int *cmprG_1dim = ReshapeInt_2to1(cmprG->data, in->height, in->width);
	int *cmprB_1dim = ReshapeInt_2to1(cmprB->data, in->height, in->width);
	float *sobel_abs_1dim = ReshapeDouble_2to1(sobel_abs, in->height, in->width);
	float *sobel_angle_1dim = ReshapeDouble_2to1(sobel_angle, in->height, in->width);
	float *grad_hx_1dim = ReshapeDouble_2to1(grad_hx, in->height, in->width+1);
	float *grad_hy_1dim = ReshapeDouble_2to1(grad_hy, in->height+1, in->width);
	float *gauss_filter_1dim = ReshapeDouble_2to1(gauce_filter, w, w);
	
	//キャンバスを0で初期化
	int *nimgR_1dim = (int*)malloc(sizeof(int)*in->height*in->width);
	int *nimgG_1dim = (int*)malloc(sizeof(int)*in->height*in->width);
	int *nimgB_1dim = (int*)malloc(sizeof(int)*in->height*in->width);
	for(int i=0; i<in->width*in->height; i++){
		nimgR_1dim[i] = 255;
		nimgG_1dim[i] = 255;
		nimgB_1dim[i] = 255;
	}

	cudaMemcpy(dev_PerlinNoise, PerlinNoise_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_Lab_L, in_Lab_L, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_Lab_a, in_Lab_a, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_in_Lab_b, in_Lab_b, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_cmpr, cmpr_1dim, in->width*in->height*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_cmprR, cmprR_1dim, in->width*in->height*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_cmprG, cmprG_1dim, in->width*in->height*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_cmprB, cmprB_1dim, in->width*in->height*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sobel_abs, sobel_abs_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_sobel_angle, sobel_angle_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_grad_hx, grad_hx_1dim, (in->width+1)*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_grad_hy, grad_hy_1dim, in->width*(in->height+1)*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nimgR, nimgR_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nimgG, nimgG_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_nimgB, nimgB_1dim, in->width*in->height*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_gauss_filter, gauss_filter_1dim, w*w*sizeof(float), cudaMemcpyHostToDevice);

	dim3 block_num(opt_GPU_Block_Num_x, opt_GPU_Block_Num_y);	//カーネルの起動ブロック数を定義
	dim3 thread_num(opt_GPU_Thread_Num_x, opt_GPU_Thread_Num_y);	//カーネルの起動スレッド数を定義

	///////////////////preprocess終了/////////////////
	Add_dictionary_to_sentence_d(log_sentence, "PreProsessTIME[s]", my_clock());
	pd("PreProsessTIME[s]",my_clock());

	//太いストロークから順番にストロークを小さくしておおまかに絵の形を取っていく
	for(t=thick_max; t>=thick_min; t--){
		if(opt_num_thick){
			int thick_arr[opt_num_thick] = opt_thick_assignment;
			int thick_flag=1;
			for (i = 0; i < opt_num_thick; i++){
				if(t==thick_arr[i]) thick_flag=0;
			}
			if(thick_flag) continue;
		}

		printf("t = %d\n",t);

		// 最適なストロークに関するデータを初期化
		format_ally(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height, UNCALCULATED);
		int *GLOBAL_improved_value_map_1dim = ReshapeInt_2to1(GLOBAL_improved_value_map->data, in->height, in->width);
		cudaMemcpy(dev_GLOBAL_improved_value_map, GLOBAL_improved_value_map_1dim, in->width*in->height*sizeof(int), cudaMemcpyHostToDevice);

		stroke_num = 99999;

		for(s_count=0; s_count<stroke_num; s_count++) {  //ある半径におけるストローク回数

			//cudaEventRecord(start, 0);//計測スタート

			//各ストロークの改善値を並列に計算するカーネル
			gpu_calculate_best_stroke<<<block_num,thread_num>>>(dev_GLOBAL_improved_value_map, dev_PerlinNoise, dev_cmprR, dev_cmprG, dev_cmprB, dev_nimgR,
														dev_nimgG, dev_nimgB, dev_best_stroke_map_pnum, dev_best_stroke_map_point_x,
														dev_best_stroke_map_point_y, dev_best_stroke_map_R, dev_best_stroke_map_G, 
														dev_best_stroke_map_B, dev_sobel_abs, dev_sobel_angle, dev_grad_hx, dev_grad_hy, dev_in_Lab_L,
														dev_in_Lab_a, dev_in_Lab_b, dev_M, dev_u, dev_new_u,
														dev_v, dev_new_v, dev_p, dev_gR, dev_gG, dev_gB, dev_dR, dev_dG, dev_dB,
														dev_new_gR, dev_new_gG, dev_new_gB, dev_new_dR, dev_new_dG, dev_new_dB,
														dev_gauss_filter, dev_gauss_M, dev_h, in->width, in->height, t);
			cudaDeviceSynchronize();
			//printCudaLastError();

			//改善値マップ中の最大値を探索するカーネル
			gpu_select_best_stroke<<<1,1>>>(dev_GLOBAL_improved_value_map, dev_best_x, dev_best_y, dev_diff_stroke_max, in->width, in->height);
			cudaDeviceSynchronize();
			//printCudaLastError();

			//改善値マップ中の最大値をCPUにコピー
			cudaMemcpy(diff_stroke_max, dev_diff_stroke_max, 1*sizeof(int), cudaMemcpyDeviceToHost);
			printf("%d : diff_stroke_max = %d\n" , s_count, diff_stroke_max[0]);

			//改善値マップ中の最大値が閾値以下ならば次の半径の処理へ進む
			if(diff_stroke_max[0] < optimal_improved_value_border) {
				strcat(log_sentence, "\r\n");
				Add_dictionary_to_sentence(log_sentence, "t", t);
				Add_dictionary_to_sentence(log_sentence, "s_count", s_count);
				Add_dictionary_to_sentence_d(log_sentence, "TIME[s]", my_clock());
				{
					strcpy(out_filename, dir_path);
					strcat(out_filename, in_filename);
					snprintf(count_name, 16, "%02d", t);
					strcat(out_filename, "_t");
					strcat(out_filename, count_name);
					strcat(out_filename, "_BSM.pgm");
					if(write_pgm(out_filename, GLOBAL_improved_value_map)){ printf("WRITE PGM ERROR.");}
					printf("%s\n",out_filename);
				}
				break;
			}

			//実際にストロークを描画する関数
			gpu_draw_best_stroke<<<1,thread_num>>>(dev_PerlinNoise, dev_nimgR, dev_nimgG, dev_nimgB, dev_best_stroke_map_pnum, dev_best_stroke_map_point_x,
														dev_best_stroke_map_point_y, dev_best_stroke_map_R, dev_best_stroke_map_G, 
														dev_best_stroke_map_B, dev_grad_hx, dev_grad_hy, dev_M, dev_u, dev_new_u,
														dev_v, dev_new_v, dev_p, dev_gR, dev_gG, dev_gB, dev_dR, dev_dG, dev_dB,
														dev_new_gR, dev_new_gG, dev_new_gB, dev_new_dR, dev_new_dG, dev_new_dB,
														dev_gauss_filter, dev_gauss_M, dev_h, dev_best_x, dev_best_y, in->width, in->height, t);
			
			cudaDeviceSynchronize();
			//printCudaLastError();

			//描画したストロークの周囲のみ改善値マップをリセット
			gpu_reset_improved_value_map<<<1,1>>>(dev_GLOBAL_improved_value_map, dev_best_stroke_map_pnum, dev_best_stroke_map_point_x, dev_best_stroke_map_point_y, dev_best_x, dev_best_y, in->width, in->height, t);
			cudaDeviceSynchronize();
			//printCudaLastError();

			//cudaEventRecord(end, 0);//計測ストップ

			//実行時間表示
			//cudaEventElapsedTime(&ftimer, start, end);
			//printf("%f,ms\n", ftimer);

			//キャンバスをGPUからCPUにコピー
			cudaMemcpy(nimgR_1dim, dev_nimgR, in->width*in->height*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(nimgG_1dim, dev_nimgG, in->width*in->height*sizeof(int), cudaMemcpyDeviceToHost);
			cudaMemcpy(nimgB_1dim, dev_nimgB, in->width*in->height*sizeof(int), cudaMemcpyDeviceToHost);
			
			//1次元配列を2次元配列に変換
			nimgR->data = ReshapeInt_1to2(nimgR_1dim, in->width, in->height);
			nimgG->data = ReshapeInt_1to2(nimgG_1dim, in->width, in->height);
			nimgB->data = ReshapeInt_1to2(nimgB_1dim, in->width, in->height);
			nimgC->dataR = nimgR->data;
			nimgC->dataG = nimgG->data;
			nimgC->dataB = nimgB->data;
			
			//一定ストロークごとに途中経過画像を書き出す
			paint_count++;
			nc++;
			if(nc%100==0 || nc<=100){
				strcpy(out_filename, dir_path);
				strcat(out_filename, in_filename);
				snprintf(count_name, 16, "%d", nc);
				strcat(out_filename, "_s");
				strcat(out_filename, count_name);
				strcat(out_filename, ".ppm");
				write_ppm(out_filename, nimgC);
				//printf("%s\n",out_filename);
				//printf("%d:",t);
				//pd("TIME[s]",my_clock());
			}
		}

		//現在の半径の完成画像を出力
		strcpy(out_filename, dir_path);
		strcat(out_filename, in_filename);
		snprintf(count_name, 16, "%02d", t);
		strcat(out_filename, "__t");
		strcat(out_filename, count_name);
		strcat(out_filename, ".ppm");
		write_ppm(out_filename, nimgC);
		//printf("%s\n",out_filename);
		//printf("%d:",t);
		//pd("TIME[s]",my_clock());

		printf("\n////////////////////\nt%d done.\n////////////////////\n\n\n",t);
		p("Paint_num",paint_count); paint_count=0;
		paint_count=0;
	}

	//第一段階描画後の中間画像を出力
	snprintf(count_name, 16, "%03d", tc);
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__t");
	strcat(out_filename, count_name);
	strcat(out_filename, ".ppm");
	write_ppm(out_filename, nimgC);
	printf("%s\n",out_filename);
	tc=t;
	printf("%d",t);
	pd("TIME[s]",my_clock());

	//---------------------------
	//エッジマップを計算し、エッジの複雑な周辺だけに描画を行う
	//---------------------------
	int lc=0;
	loop_cont = opt2_loop_cont;
	min_stroke = opt2_min_stroke;
	ratio = opt2_ratio;
	tc=-1;
	double maxValue=0.30, minValue=0.10;
	PGM *canny;
	PGM *EdgeMap;

	thick_max = opt2_thick_max;
	if(thick_max){
		strcat(log_sentence, "\r\n\r\n[EdgeStroke]\r\n");
		canny = cannyedge_detector(gray, maxValue, minValue, thick_min);
		EdgeMap = calcu_EdgeMap(canny, thick_min, sobel_angle);
		//EdgeMap = expand_Edge(canny, thick_min);
	}
	thick_min = opt2_thick_min;

	Add_dictionary_to_sentence_d(log_sentence, "CannyTIME[s]", my_clock());
	pd("Canny:TIME[s]",my_clock());

	for(t=thick_max; t>=thick_min; t--){
		paint_count=0;

		//ストロークサイズのガウスフィルタを生成
		gauce_filter = create_dally(2*t+1, 2*t+1);
		sigma = t/3.0;
		for(i=0; i<2*t+1; i++){
			for(j=0; j<2*t+1; j++){
				gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
			}
		}

		for(y=y_defo; y<in->height; y=y+t) {  //ウィンドウの大きさに合わせて
			for(x=x_defo; x<in->width; x=x+t) {  //ウィンドウをずらす距離を変えとく
				//エッジマップにない箇所なら次のループに
				if(EdgeMap->data[x][y]!=255) continue;

				//ウィンドウの中の差分の合計を取る
				diff_sum = break_flag = pnum = 0;

				offscrn_count = 0;
				for(xc=-t; xc<=t; xc++) {
					if((x+xc)<0 || (x+xc)>in->width-1) {offscrn_count += 2*t+1;		continue;	}
					for(yc=-t; yc<=t; yc++) {
						if((y+yc)<0 || (y+yc)>in->height-1) {	offscrn_count++;
						}else{
							diff_sum += abs(nimgR->data[x+xc][y+yc] - cmprR->data[x+xc][y+yc]);
							diff_sum += abs(nimgG->data[x+xc][y+yc] - cmprG->data[x+xc][y+yc]);
							diff_sum += abs(nimgB->data[x+xc][y+yc] - cmprB->data[x+xc][y+yc]);
						}
					}
				}
				diff_sum = diff_sum/((2*t+1)*(2*t+1)-offscrn_count)/3;

				//差分の合計平均(画素当たりの差分)が一定以上ならストローク開始位置とする
				if(diff_sum < window_diff_border) {
					stroke_histogram[pnum]++;
					continue;
				}
				pnum=1;		//第一点確定
				p[0].x=x+0.5; p[0].y=y+0.5;

				//一つ目の描画領域から描画色を平均を取って取得
				if(opt_USE_calcu_color_bi){
					bright.R = calcu_color_bi(cmprR->data, cmprR->width, cmprR->height, x, y, t, 50, gauce_filter);
					bright.G = calcu_color_bi(cmprG->data, cmprG->width, cmprG->height, x, y, t, 50, gauce_filter);
					bright.B = calcu_color_bi(cmprB->data, cmprB->width, cmprB->height, x, y, t, 50, gauce_filter);
				} else{
					bright.R = calcu_color(cmprR->data, cmprR->width, cmprR->height, x, y, t);
					bright.G = calcu_color(cmprG->data, cmprG->width, cmprG->height, x, y, t);
					bright.B = calcu_color(cmprB->data, cmprB->width, cmprB->height, x, y, t);
				}


				theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
						gauce_filter, p[0].x, p[0].y, t, &break_flag);


				//制御点を方向から計算し代入
				p[1] = calcu_point(cmpr, p[0], t, theta);


				//二つ目の描画点周りの色が描画色と一致するか確認する
				if(opt_USE_Lab_ColorDiff){
					sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, 1.0);
				} else{
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
				}

				//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
				if( sum < color_diff_border){
					//もう一つの勾配垂直の点を代入
					theta += PI;
					p[1] = calcu_point(cmpr, p[0], t, theta);

					//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
					if(opt_USE_Lab_ColorDiff){
						sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, 1.0);
					} else{
						sum = 0;
						sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
						sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
						sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					}

					//どちらの第二点も不適切なら描画をせず次のループへ
					if( sum < color_diff_border) {
						continue;
					}
				}
				//適切な第二点が見つかれば次へ
				pnum=2;		//第二点確定


				/*
					POINT2から勾配により次の制御点を探していく
				*/

				while(pnum<max_stroke){
					former_theta=theta;

					//第pnum点周りにおいて、sobelからヒストグラムを作成し最大のものを勾配とする
					theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
											gauce_filter, p[pnum-1].x, p[pnum-1].y, t, &break_flag);

					//制御点の為す角が急峻になるようなら逆方向に角度を取る
					if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;}
					p[pnum] = calcu_point(cmpr, p[pnum-1], t, theta);

					//pnum+1目の描画点周りの色が描画色と一致するか確認する  //点当たりの差異平均
					if(opt_USE_Lab_ColorDiff){
						sum = diffsum_Lab(in_Lab, nimgC, p[pnum], t, bright, 1.0);
					} else{
						sum = 0;
						sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
						sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
						sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
					}

					/*
						pnum+1目の(次の)制御点周りの色が描画色としきい値以上の差を持つなら
						それまでの制御点を用いて線を描画
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}

				}

				if(pnum>=min_stroke) {
					//算出したpnum個の制御点を用いてストロークを描画
					Paint_Bezier_ex(p, pnum, nimgR, t, bright.R, ratio);
					Paint_Bezier_ex(p, pnum, nimgG, t, bright.G, ratio);
					Paint_Bezier_ex(p, pnum, nimgB, t, bright.B, ratio);

				}

				paint_count++;
				nc++;
			}
		}


		printf("////////////////////\nt%d done.\n////////////////////\n\n\n",t);
		Free_dally(gauce_filter, 2*t+1);


		{
			tc=t;
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			snprintf(count_name, 16, "%02d", t);
			strcat(out_filename, "__st");
			strcat(out_filename, count_name);
			snprintf(count_name, 16, "%02d", lc);
			strcat(out_filename, "_lc");
			strcat(out_filename, count_name);
			strcat(out_filename, ".ppm");
			write_ppm(out_filename, nimgC);
			printf("%s\n",out_filename);
			printf("%d:",t);
			strcat(log_sentence, "\r\n");
			Add_dictionary_to_sentence(log_sentence, "t", t);
			Add_dictionary_to_sentence(log_sentence, "s_count", paint_count);
			Add_dictionary_to_sentence_d(log_sentence, "TIME[s]", my_clock());
			pd("TIME[s]",my_clock());
		}

		lc++;		//同じ半径でのループをcont回する
		if((lc%loop_cont) != 0){
			p("lc",lc);
			t++;
			x_defo += t/loop_cont;
			y_defo += t/loop_cont;
			if(t/loop_cont==0){
				x_defo++;	y_defo++;
			}
		}
	}

	if(thick_max){
		strcpy(out_filename, dir_path);
		strcat(out_filename, in_filename);
		strcat(out_filename, "__EdgeMap");
		strcat(out_filename, ".pgm");
		if(write_pgm(out_filename, EdgeMap)){ printf("WRITE PNG ERROR.");}
	}

	double MSE = image_MSE(nimgC, in);
	Add_dictionary_to_sentence(log_sentence, "MSE", (int)MSE);

	Add_dictionary_to_sentence_d(log_sentence, "All_Execution_TIME", my_clock());
	printf("%s\n", log_filename);
	if(log_print(log_filename, log_sentence, "w") ){ printf("LOG_PRINTING_FAIL\n"); }


	Free_dally(sobel_abs, in->width);
	Free_dally(sobel_angle, in->width);
	Free_dally(h, in->width);
	Free_dally(grad_hx, in->width+1);
	Free_dally(grad_hy, in->width);
	FreePGM(gray);
	FreePGM(inR);
	FreePGM(inG);
	FreePGM(inB);
	FreePGM(nimgV);
	free(nimgR);
	free(nimgG);
	free(nimgB);
	if(thick_max){
		FreePGM(canny);
		FreePGM(EdgeMap);
	}

	cudaFree(dev_GLOBAL_improved_value_map);
	cudaFree(dev_PerlinNoise);
	cudaFree(dev_in_Lab_L);
	cudaFree(dev_in_Lab_a);
	cudaFree(dev_in_Lab_b);
	cudaFree(dev_cmpr);
	cudaFree(dev_cmprR);
	cudaFree(dev_cmprG);
	cudaFree(dev_cmprB);
	cudaFree(dev_nimgR);
	cudaFree(dev_nimgG);
	cudaFree(dev_nimgB);
	cudaFree(dev_best_stroke_map_pnum);
	cudaFree(dev_best_stroke_map_point_x);
	cudaFree(dev_best_stroke_map_point_y);
	cudaFree(dev_best_stroke_map_R);
	cudaFree(dev_best_stroke_map_G);
	cudaFree(dev_best_stroke_map_B);
	cudaFree(dev_sobel_abs);
	cudaFree(dev_sobel_angle);
	cudaFree(dev_grad_hx);
	cudaFree(dev_grad_hy);
	cudaFree(dev_M);
	cudaFree(dev_u);
	cudaFree(dev_new_u);
	cudaFree(dev_v);
	cudaFree(dev_new_v);
	cudaFree(dev_p);
	cudaFree(dev_gR);
	cudaFree(dev_gG);
	cudaFree(dev_gB);
	cudaFree(dev_dR);
	cudaFree(dev_dG);
	cudaFree(dev_dB);
	cudaFree(dev_gauss_M);
	cudaFree(dev_h);

    return nimgC;
}