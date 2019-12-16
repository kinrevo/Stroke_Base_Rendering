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
#include "sbr_opt.h"
#include "water.h"


#define p(s,a) printf("%s:%d\n",s,a)

PPM *c_Illust_brush_Water(PPM *in, char *filename);
PPM *c_Illust_brush_Water_best(PPM *in, char *filename);

PGM* GLOBAL_improved_value_map;
// extern const int opt_Stroke_Method;

int main(int argc, char *argv[])
{
	my_clock();

	if(argc<=2){
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

	//入力画像の絵画化
	trans_ppm = c_Illust_brush_Water_INTEGRATED(in_ppm, argv[2]);

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
	return 0;
}


//与えられたイメージからブラッシングによる絵画を作る
PPM *c_Illust_brush_Water_INTEGRATED(PPM *in, char *filename)
{
	my_clock();
	int i,j,x,y,xc,yc,t=100,break_flag,pnum, offscrn_count;
	int window_diff_border = opt_window_diff_border; 	//ストローク位置探索のしきい値
	int color_diff_border = opt_color_diff_border;  	//描画色の差異のしきい値
	int max_stroke = opt_max_stroke;
	int min_stroke = opt_min_stroke;
	Point p[max_stroke];
	// int stroke_histogram[max_stroke+1];
	// for(i=0; i<max_stroke+1; i++){stroke_histogram[i]=0;}
	double ratio=opt_ratio;
	double theta, former_theta;
	double **gauce_filter;
	double sigma, diff_sum, sum;
	int histogram_partition=opt_histogram_partition;
	int paint_count=0, nc=0;
	int loop_cont=opt_loop_cont, x_defo=0, y_defo=0;
	int lc=0;
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
	// char vec_sentence[128] = "", tmp_sentence[32] = "", vec_filename[namesize];
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
	{
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
		if(opt_Stroke_Method==Best_StrokeOrder)	strcat(log_sentence, "Stroke_Method:Best_StrokeOrder\r\n");
		if(opt_Stroke_Method==Raster_StrokeOrder) strcat(log_sentence, "Stroke_Method:Raster_StrokeOrder\r\n");
		if(opt_Stroke_Method==EveryColor_StrokeOrder) strcat(log_sentence, "Stroke_Method:EveryColor_StrokeOrder\r\n");
		Add_dictionary_to_sentence(log_sentence, "width", in->width);
		Add_dictionary_to_sentence(log_sentence, "height", in->height);
		Add_dictionary_to_sentence(log_sentence, "thick_max", thick_max);
		Add_dictionary_to_sentence(log_sentence, "thick_min", thick_min);
		Add_dictionary_to_sentence(log_sentence, "max_stroke", max_stroke);
		Add_dictionary_to_sentence(log_sentence, "min_stroke", min_stroke);
		Add_dictionary_to_sentence(log_sentence, "window_diff_border", window_diff_border);
		Add_dictionary_to_sentence(log_sentence, "color_diff_border", color_diff_border);
		Add_dictionary_to_sentence(log_sentence, "histogram_partition", histogram_partition);
		Add_dictionary_to_sentence(log_sentence, "loop_cont", loop_cont);
		if(opt_USE_best_ratio){
			strcat(log_sentence, "# USE_best_ratio\r\n");
			Add_dictionary_to_sentence_d(log_sentence, "max_ratio", opt_max_ratio);
			Add_dictionary_to_sentence_d(log_sentence, "min_ratio", opt_min_ratio);
			Add_dictionary_to_sentence_d(log_sentence, "ratio_step", opt_ratio_step);
		}else{
			Add_dictionary_to_sentence_d(log_sentence, "ratio", opt_ratio);
		}
		if(opt_USE_Lab_ColorDiff){
			strcat(log_sentence, "# USE_Lab_ColorDiff\r\n");
			if(opt_USE_best_ratio){
				strcat(log_sentence, "  # USE_Lab_ColorDiff_Weight\r\n");
				Add_dictionary_to_sentence(log_sentence, "    Lab_Weight", opt_Lab_Weight);
			}
		}else strcat(log_sentence, "# USE_RGB_ColorDiff\r\n");
		if(opt_USE_Lab_PaintColorDiff){
			strcat(log_sentence, "# USE_Lab_PaintColorDiff\r\n");
			if(opt_USE_best_ratio){
				strcat(log_sentence, "  # USE_Lab_ColorDiff_Weight\r\n");
				Add_dictionary_to_sentence(log_sentence, "    Lab_Weight", opt_Lab_Weight);
			}
		}else strcat(log_sentence, "# USE_RGB_PaintColorDiff\r\n");
		if(opt_USE_calcu_color_bi) strcat(log_sentence, "# USE_calcu_color_bi\r\n");
		else if(opt_USE_calcu_Kmean_ColorSet){
			strcat(log_sentence, "# USE_calcu_Kmean_ColorSet\r\n");
			Add_dictionary_to_sentence(log_sentence, "  Kmean_ClusterNum", opt_Kmean_ClusterNum);
		}else if(opt_USE_calcu_JIS_ColorSet){
			strcat(log_sentence, "# USE_calcu_JIS_ColorSet\r\n");
			Add_dictionary_to_sentence(log_sentence, "  JIS_ClusterNum", opt_JIS_ClusterNum);
		}else strcat(log_sentence, "# USE_average_PaintColor\r\n");
		if(opt_USE_gause_histogram) strcat(log_sentence, "  # USE_gause_histogram\r\n");
		if(opt_Stroke_Method==Best_StrokeOrder) Add_dictionary_to_sentence(log_sentence, "optimal_improved_value_border", opt_optimal_improved_value_border);
		Add_dictionary_to_sentence_d(log_sentence, "StrokeWindowStep", opt_StrokeWindowStep/t);
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
		pn;
	}

	//vectorデータのヘッダを格納
	// snprintf(tmp_sentence, 32, "%d", in->width);
	// strcat(vec_sentence, tmp_sentence);
	// strcat(vec_sentence, " ");
	// snprintf(tmp_sentence, 32, "%d", in->height);
	// strcat(vec_sentence, tmp_sentence);
	// strcpy(vec_filename, dir_path);
	// strcat(vec_filename, in_filename);
	// strcat(vec_filename, ".vec");
	// log_print(vec_filename, vec_sentence, "w");


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
	PPM *cmprC  = in;
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

		if (strcmp("ppm", ext) == 0 || strcmp("pnm", ext) == 0) {
			nimgC = read_ppm(name);
		} else if (strcmp("jpg", ext) == 0 || strcmp("jpeg", ext) == 0) {
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

		devide_ppm(nimgC, nimgR, nimgG, nimgB);
	} else{
		nimgC = create_ppm(in->width, in->height, in->bright); //実際に描画するキャンバス
	}
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;


	// スケーリングキャンバスの生成
	Point *scaling_p;
	double **gauce_filter_Scaling;
	PGM *nimgR_Scaling = create_pgm(gray->width*opt_canvas_scaling_ratio, gray->height*opt_canvas_scaling_ratio, gray->bright);
	PGM *nimgG_Scaling = create_pgm(gray->width*opt_canvas_scaling_ratio, gray->height*opt_canvas_scaling_ratio, gray->bright);
	PGM *nimgB_Scaling = create_pgm(gray->width*opt_canvas_scaling_ratio, gray->height*opt_canvas_scaling_ratio, gray->bright);
	PPM *nimgC_Scaling = create_ppm(in->width*opt_canvas_scaling_ratio, in->height*opt_canvas_scaling_ratio, in->bright); //実際に描画するキャンバス（拡縮描画用）
	nimgC_Scaling->dataR = nimgR_Scaling->data;
	nimgC_Scaling->dataG = nimgG_Scaling->data;
	nimgC_Scaling->dataB = nimgB_Scaling->data;
	double** h_Scaling = perlin_img(nimgC_Scaling->width, nimgC_Scaling->height, opt_perlin_freq, opt_perlin_depth);
    double** grad_hx_Scaling = create_dally(nimgC_Scaling->width+1, nimgC_Scaling->height);
    double** grad_hy_Scaling = create_dally(nimgC_Scaling->width, nimgC_Scaling->height+1);
    calcu_grad_h(h_Scaling, grad_hx_Scaling, grad_hy_Scaling, nimgC_Scaling->width, nimgC_Scaling->height);


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
	int diff_stroke_max=0;
	int tmp_num;
	int diff_stroke_max_ave[50]={};
	int best_x=0, best_y=0;
	Stroke*** best_stroke_map = create_Stroke_ally(in->width, in->height, max_stroke);
	Point best_P;
	int optimal_improved_value_border = opt_optimal_improved_value_border;

    double** h = perlin_img(in->width, in->height, opt_perlin_freq, opt_perlin_depth);
    double** grad_hx = create_dally(in->width+1, in->height);
    double** grad_hy = create_dally(in->width, in->height+1);
    calcu_grad_h(h, grad_hx, grad_hy, in->width, in->height);

	// KmeanカラーセットまたはJISカラーセット
	int *num_cluster, *x_centlabel;
	int** x_centlabel_2D;
	RGB* ColorSet;
	int ColorNum;
	if(opt_USE_calcu_Kmean_ColorSet){
		ColorNum = opt_Kmean_ClusterNum;
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
		FreePPM(Kmean_Img);
		FreePPM(ColorSet_Img);
	}
	else if(opt_USE_calcu_JIS_ColorSet){
		ColorNum = opt_JIS_ClusterNum;
		ColorSet = create_JIS_ColorSet(opt_JIS_ClusterNum);
	}

	//　RGB入力画像をLabに変換した配列を用意
	RGB CanRGB;
	Lab** in_Lab = (Lab**)malloc(sizeof(Lab*)*(in->width));
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

	///////////////////preprocess終了/////////////////
	Add_dictionary_to_sentence_d(log_sentence, "\r\nPreProsessTIME[s]", my_clock());
	pd("PreProsessTIME[s]",my_clock());



	// 最も誤差の減らせるストロークから順に配置して行く
	if(opt_Stroke_Method==Best_StrokeOrder)
	{
		printf("Best_Stroke \n");
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


			//ストロークサイズのガウスフィルタを生成
			gauce_filter = create_dally(2*t+1, 2*t+1);
			sigma = t/3.0*opt_variance_ratio;
			sum = 0;
			for(i=0; i<2*t+1; i++){
				for(j=0; j<2*t+1; j++){
					gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
					sum += gauce_filter[i][j];
				}
			}
			for(i=0; i<2*t+1; i++){
				for(j=0; j<2*t+1; j++){
					gauce_filter[i][j] = gauce_filter[i][j] / sum;
				}
			}

			//ストロークサイズのガウスフィルタを生成(スケーリングキャンバス))
			int t_Scaling = t*opt_canvas_scaling_ratio;
			gauce_filter_Scaling = create_dally(2*t_Scaling+1, 2*t_Scaling+1);
			sigma = t_Scaling/3.0*opt_variance_ratio;
			sum = 0;
			for(i=0; i<2*t_Scaling+1; i++){
				for(j=0; j<2*t_Scaling+1; j++){
					gauce_filter_Scaling[i][j] = gause_func(i-t_Scaling, j-t_Scaling, sigma);
					sum += gauce_filter_Scaling[i][j];
				}
			}
			for(i=0; i<2*t+1; i++){
				for(j=0; j<2*t+1; j++){
					gauce_filter_Scaling[i][j] = gauce_filter_Scaling[i][j] / sum;
				}
			}


			stroke_num=99999;

			// 最適なストロークに関するデータを初期化
			format_ally(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height, UNCALCULATED);

			// #ifdef _OPENMP
			// 	Format_SINGLE_Paint_Water();
			// #endif
			for(s_count=0; s_count<stroke_num; s_count++) {  //ある半径におけるストローク回数
				// 誤差探索の変数を初期化
				diff_stroke_max = UNCALCULATED;

				int x_step=opt_StrokeWindowStep, y_step=opt_StrokeWindowStep;
				// #pragma omp parallel for private(x,xc,yc,test_Canvas,diff_sum,break_flag,pnum,offscrn_count,theta,sum,former_theta) schedule(static, 1)
				for(y=0; y<in->height; y=y+y_step) {
					// #ifdef _OPENMP
					// 	test_Canvas = create_ppm(in->width, in->height, in->bright);
					// #endif
					for(x=0; x<in->width; x=x+x_step) {
					//for(y=y_defo; y<in->height; y=y+t) {  //ウィンドウの大きさに合わせて
					//for(x=x_defo; x<in->width; x=x+t) {  //ウィンドウをずらす距離を変えとく
						// 改善値が計算済みならSkip
						if(GLOBAL_improved_value_map->data[x][y] != UNCALCULATED) continue;

						diff_sum = break_flag = pnum = 0;

						//ウィンドウの中の差分の合計を取る
						offscrn_count = 0;
						for(xc=-t; xc<=t; xc++) {
							if((x+xc)<0 || (x+xc)>in->width-1) {offscrn_count += 2*t+1;		continue;	}
							for(yc=-t; yc<=t; yc++) {
								if((y+yc)<0 || (y+yc)>in->height-1) {	offscrn_count++;
								}else{
									//diff_sum += abs(nimgV->data[x+xc][y+yc] - cmpr->data[x+xc][y+yc]);
									diff_sum += abs(nimgR->data[x+xc][y+yc] - cmprR->data[x+xc][y+yc]);
									diff_sum += abs(nimgG->data[x+xc][y+yc] - cmprG->data[x+xc][y+yc]);
									diff_sum += abs(nimgB->data[x+xc][y+yc] - cmprB->data[x+xc][y+yc]);
								}
							}
						}
						diff_sum = diff_sum/((2*t+1)*(2*t+1)-offscrn_count)/3;

						//差分の合計平均(画素当たりの差分)が一定以上ならストローク開始位置とする
						if(diff_sum < window_diff_border) {
							GLOBAL_improved_value_map->data[x][y] = SMALL_DIFF;
							continue;
						}
						pnum=1;		//第一点確定
						best_stroke_map[x][y]->p[0].x=x+0.5; best_stroke_map[x][y]->p[0].y=y+0.5;

						//一つ目の描画領域から描画色を取得
						calcu_color_INTEGRATED(best_stroke_map[x][y], in_Lab, cmprC, nimgC, best_stroke_map[x][y]->p[0], t, ColorSet, gauce_filter);


						theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
								gauce_filter, best_stroke_map[x][y]->p[0].x, best_stroke_map[x][y]->p[0].y, t, &break_flag);


						//制御点を方向から計算し代入
						best_stroke_map[x][y]->p[1] = calcu_point(cmpr, best_stroke_map[x][y]->p[0], t, theta);


						//二つ目の描画点周りの色が描画色と一致するか確認する
						if(opt_USE_Lab_ColorDiff){
							sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
						} else{
							sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
						}


						//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
						if(sum < color_diff_border){
							//もう一つの勾配垂直の点を代入
							theta += PI;
							best_stroke_map[x][y]->p[1] = calcu_point(cmpr, best_stroke_map[x][y]->p[0], t, theta);

							//反対方向の第二点の描画点周りの色が描画色と一致するか確認する//点当たりの差異平均
							if(opt_USE_Lab_ColorDiff){
								sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
							} else{
								sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
							}

							//どちらの第二点も不適切なら描画をせず次のループへ
							if( sum < color_diff_border) {
								GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
								continue;
							}
						}

						pnum=2;		//第二点確定


						/*
							POINT2から勾配により次の制御点を探していく
						*/

						while(pnum<max_stroke){
							former_theta=theta;

							//第pnum点周りにおいて、sobelからヒストグラムを作成し最大のものを勾配とする
							theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
													gauce_filter, best_stroke_map[x][y]->p[pnum-1].x, best_stroke_map[x][y]->p[pnum-1].y, t, &break_flag);

							//制御点の為す角が急峻になるようなら逆方向に角度を取る
							if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;}
							best_stroke_map[x][y]->p[pnum] = calcu_point(cmpr, best_stroke_map[x][y]->p[pnum-1], t, theta);

							//pnum+1目の描画点周りの色が描画色と一致するか確認する
							if(opt_USE_Lab_ColorDiff){
								sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[pnum], t, best_stroke_map[x][y]->color, ratio);
							} else{
								sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[pnum], t, best_stroke_map[x][y]->color, ratio);
							}

							/*
								pnum+1目の(次の)制御点周りの色が描画色としきい値以上の差を持つなら
								それまでの制御点を用いて線を描画
							*/
							if( sum < color_diff_border) {break_flag=1; break;}
							else {
								diff_sum += sum;
								pnum++;
							}
						}


						//best_stroke_map[x][y]->color=bright;
						best_stroke_map[x][y]->pnum=pnum;
						//　試しに描いてみて誤差を確認
						if(pnum>=min_stroke){
							// GLOBAL_improved_value_map->data[x][y] = diff_sum;
							GLOBAL_improved_value_map->data[x][y] = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[x][y], t, h, grad_hx, grad_hy, gauce_filter, ratio);
						}else if(pnum<min_stroke){
							GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
						}
					}
					// #ifdef _OPENMP
					// 	FreePPM(test_Canvas);
					// #endif
				}

				// 改善値マップ中の最大値を探索
				best_P = search_max_Point(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height);
				best_x=best_P.x;  best_y=best_P.y;
				// 最適な開始点が前回と一致（ループ）している場合は
				// if(best_P.x==before_P.x && best_P.y==before_P.y){
				// 	GLOBAL_improved_value_map->data[best_x][best_y] = SAME_STROKE;
				// 	printf("LOOP_STROKE:[%d,%d]->",best_x,best_y);
				// 	best_P = search_max_Point(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height);
				// 	best_x=best_P.x;  best_y=best_P.y;
				// 	printf("[%d,%d]\n",best_x,best_y);
				// }
				// before_P.x = best_P.x; before_P.y = best_P.y;
				// diff_stroke_max = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[best_x][best_y], t, h, grad_hx, grad_hy, gauce_filter);
				diff_stroke_max = GLOBAL_improved_value_map->data[best_x][best_y];

				// 直近の50のストロークでの変化がほとんどなければ次の半径ステップに移行
				// int diff_test_stroke_ave[50];
				diff_stroke_max_ave[s_count%50] = diff_stroke_max;
				// diff_test_stroke_ave[s_count%50] = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[best_x][best_y], t, h, grad_hx, grad_hy, gauce_filter);
				if(s_count%50 == 49)
				{
					tmp_num=0;
					for(i=0; i<50; i++) tmp_num+=diff_stroke_max_ave[i];
					p("Stroke_MAX_AVE", tmp_num);
					tmp_num=0;
					// for(i=0; i<50; i++) tmp_num+=diff_test_stroke_ave[i];
					// p("test_STROKE_AVE", tmp_num);
					// if(tmp_num/50 < optimal_improved_value_border) {
					// 	strcat(log_sentence, "\r\n");
					// 	Add_dictionary_to_sentence(log_sentence, "t", t);
					// 	Add_dictionary_to_sentence(log_sentence, "s_count", s_count);
					// 	snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
					// 	strcat(log_sentence, count_name);
					// 	strcat(log_sentence, "\r\n");
					// 	break;
					// }
				}
				if(diff_stroke_max < optimal_improved_value_border) {
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


				//算出したpnum個の制御点を用いてストロークを描画
				Paint_Water_Stroke(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, t, best_stroke_map[best_x][best_y]->color, nimgR->data, nimgG->data, nimgB->data, h, grad_hx, grad_hy, gauce_filter, in->width, in->height, best_stroke_map[best_x][best_y]->ratio);

				// 描画したストロークの周囲のみ改善値マップをリセット
				reset_improved_value_map(GLOBAL_improved_value_map, best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, best_stroke_map, t, max_stroke);


				//一定ストロークごとに途中経過画像を書き出す
				paint_count++;
				nc++;
				// if(nc%1==0)
				if(nc%100==0 || nc<=100)
				{
						strcpy(out_filename, dir_path);
						strcat(out_filename, in_filename);
						// snprintf(count_name, 16, "%02d", t);
						// strcat(out_filename, "_t");
						// strcat(out_filename, count_name);
						snprintf(count_name, 16, "%d", nc);
						strcat(out_filename, "_s");
						strcat(out_filename, count_name);
						strcat(out_filename, ".png");
						out_png = PPM_to_image(nimgC);
						if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
						free_image(out_png);
						printf("%s\n",out_filename);
						printf("%d:",t);
						pd("TIME[s]",my_clock());
				}


				// Greedyアプローチ：最適ストロークを探索するごとに探索位置をずらす
				x_defo += 1;
				if(x_defo>t) {
					x_defo -= t;
					y_defo += 1;
					if(y_defo>t) y_defo -= t;
				}


				/// 制御点を全てとブラシサイズを拡大率に従いスケーリングし、拡大キャンバスに描画
				if(opt_USE_Canvas_Scaling_Method){
					scaling_p = scaling_point(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, opt_canvas_scaling_ratio);
					SINGLE_Paint_Water_Stroke(scaling_p, best_stroke_map[best_x][best_y]->pnum, t_Scaling, best_stroke_map[best_x][best_y]->color, nimgR_Scaling->data, nimgG_Scaling->data, nimgB_Scaling->data, h_Scaling, grad_hx_Scaling, grad_hy_Scaling, gauce_filter_Scaling, nimgC_Scaling->width, nimgC_Scaling->height, ratio);
					free(scaling_p);

					if(nc%100==0 || nc<=100)
					{
						strcpy(out_filename, dir_path);
						strcat(out_filename, in_filename);
						strcat(out_filename, "_SC");
						snprintf(count_name, 16, "%d", nc);
						strcat(out_filename, "_s");
						strcat(out_filename, count_name);
						strcat(out_filename, ".png");
						out_png = PPM_to_image(nimgC_Scaling);
						if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
						free_image(out_png);
						printf("%s\n",out_filename);
						printf("%d:",t);
						pd("TIME[s]",my_clock());
					}
				}
			}

			{
				strcpy(out_filename, dir_path);
				strcat(out_filename, in_filename);
				snprintf(count_name, 16, "%02d", t);
				strcat(out_filename, "__t");
				strcat(out_filename, count_name);
				strcat(out_filename, ".png");
				out_png = PPM_to_image(nimgC);
				if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
				free_image(out_png);
				printf("%s\n",out_filename);
				printf("%d:",t);
				// strcat(log_sentence, "\r\n");
				// Add_dictionary_to_sentence(log_sentence, "t", t);
				// Add_dictionary_to_sentence(log_sentence, "s_count", nc);
				// snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
				// strcat(log_sentence, count_name);
				// strcat(log_sentence, "\r\n");
				pd("TIME[s]",my_clock());
			}

			if(opt_USE_Canvas_Scaling_Method){
				strcpy(out_filename, dir_path);
				strcat(out_filename, in_filename);
				strcat(out_filename, "_SC");
				snprintf(count_name, 16, "%02d", t);
				strcat(out_filename, "__t");
				strcat(out_filename, count_name);
				strcat(out_filename, ".png");
				out_png = PPM_to_image(nimgC_Scaling);
				if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
				free_image(out_png);
				printf("%s\n",out_filename);
				printf("%d:",t);
				pd("TIME[s]",my_clock());
			}

			printf("\n////////////////////\nt%d done.\n////////////////////\n\n\n",t);
			p("Paint_num",paint_count); paint_count=0;
			Free_dally(gauce_filter, 2*t+1);
			x_defo=y_defo=paint_count=0;
		}

	}


	// ラスター順に端からストロークを配置して行く
	if(opt_Stroke_Method==Raster_StrokeOrder)
	{
		printf("Raster_Stroke \n");
		Stroke* tmp_S = create_Stroke(max_stroke);

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

			//ストロークサイズのガウスフィルタを生成
			gauce_filter = create_dally(2*t+1, 2*t+1);
			sum = 0;
			sigma = t/3.0*opt_variance_ratio;
			for(i=0; i<2*t+1; i++){
				for(j=0; j<2*t+1; j++){
					gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
					sum += gauce_filter[i][j];
				}
			}
			for(i=0; i<2*t+1; i++){
				for(j=0; j<2*t+1; j++){
					gauce_filter[i][j] = gauce_filter[i][j] / sum;
				}
			}

			//vecデータ書き込み
			// strcpy(vec_sentence, "t");
			// snprintf(tmp_sentence, 32, "%d", t);
			// strcat(vec_sentence, tmp_sentence);
			// log_print(vec_filename, vec_sentence, "a");


			int x_step=opt_StrokeWindowStep, y_step=opt_StrokeWindowStep;
			// for(y=0; y<in->height; y++) {
			// 	for(x=0; x<in->width; x++) {
			for(y=y_defo; y<in->height; y=y+y_step) {  //ウィンドウの大きさに合わせて
				for(x=x_defo; x<in->width; x=x+x_step) {  //ウィンドウをずらす距離を変えとく

					diff_sum = break_flag = pnum = 0;

					//ウィンドウの中の差分の合計を取る
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
						continue;
					}
					pnum=1;		//第一点確定
					p[0].x=x+0.5; p[0].y=y+0.5;

					//一つ目の描画領域から描画色を取得
					calcu_color_INTEGRATED(tmp_S, in_Lab, cmprC, nimgC, p[0], t, ColorSet, gauce_filter);
					bright = tmp_S->color;
					ratio = tmp_S->ratio;

					theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
							gauce_filter, p[0].x, p[0].y, t, &break_flag);


					//制御点を方向から計算し代入
					p[1] = calcu_point(cmpr, p[0], t, theta);


					//二つ目の描画点周りの色が描画色と一致するか確認する
					if(opt_USE_Lab_ColorDiff){
						sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, ratio);
					} else{
						sum = diffsum_clr_RGB(cmprC, nimgC, p[1], t, bright, ratio);
					}


					//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
					if(sum < color_diff_border){
						//もう一つの勾配垂直の点を代入
						theta += PI;
						p[1] = calcu_point(cmpr, p[0], t, theta);

						//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
						if(opt_USE_Lab_ColorDiff){
							sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, ratio);
						} else{
							sum = diffsum_clr_RGB(cmprC, nimgC, p[1], t, bright, ratio);
						}

						//どちらの第二点も不適切なら描画をせず次のループへ
						if( sum < color_diff_border) {
							continue;
						}
					}

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

						//pnum+1目の描画点周りの色が描画色と一致するか確認する
						if(opt_USE_Lab_ColorDiff){
							sum = diffsum_Lab(in_Lab, nimgC, p[pnum], t, bright, ratio);
						} else{
							sum = diffsum_clr_RGB(cmprC, nimgC, p[pnum], t, bright, ratio);
						}

						/*
							pnum+1目の(次の)制御点周りの色が描画色としきい値以上の差を持つなら
							それまでの制御点を用いて線を描画
						*/
						if( sum < color_diff_border) {break_flag=1; break;}
						else {pnum++;}
					}

					if(pnum>=min_stroke) {
						// double start_PWS = my_clock();
						Paint_Water_Stroke(p, pnum, t, bright, nimgR->data, nimgG->data, nimgB->data, h, grad_hx, grad_hy, gauce_filter, in->width, in->height, ratio);
						// pd("PWS[s]",my_clock()-start_PWS);
						paint_count++;

						if(paint_count%500==0 || paint_count<100)
						{
							strcpy(out_filename, dir_path);
							strcat(out_filename, in_filename);
							snprintf(count_name, 16, "%02d", t);
							strcat(out_filename, "_t");
							strcat(out_filename, count_name);
							snprintf(count_name, 16, "%d", paint_count);
							strcat(out_filename, "_s");
							strcat(out_filename, count_name);
							strcat(out_filename, ".png");
							out_png = PPM_to_image(nimgC);
							if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
							free_image(out_png);
							printf("%s\n",out_filename);
							printf("%d:",t);
							pd("TIME[s]",my_clock());
						}
					}
				}
				x_defo =  (int)(h[0][y]*10000) % t; //適当な乱数らしい値を持ってきて初期位置を散らす
			}

			nc++;


			{
				strcpy(out_filename, dir_path);
				strcat(out_filename, in_filename);
				snprintf(count_name, 16, "%02d", t);
				strcat(out_filename, "__t");
				strcat(out_filename, count_name);
				strcat(out_filename, ".png");
				out_png = PPM_to_image(nimgC);
				if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
				free_image(out_png);
				printf("%s\n",out_filename);
				printf("%d:",t);
				strcat(log_sentence, "\r\n");
				Add_dictionary_to_sentence(log_sentence, "t", t);
				Add_dictionary_to_sentence(log_sentence, "s_count", paint_count);
				Add_dictionary_to_sentence_d(log_sentence, "TIME[s]", my_clock());
				pd("TIME[s]",my_clock());
			}

			printf("\n////////////////////\nt%d done.\n////////////////////\n\n\n",t);
			p("Paint_num",paint_count); paint_count=0;
			Free_dally(gauce_filter, 2*t+1);
			x_defo=y_defo=paint_count=0;
		}

	}


	// 最も誤差の減らせるストロークから順に配置して行く
	if(opt_Stroke_Method==EveryColor_StrokeOrder)
	{
		printf("Every_Color \n");

		for (int c = 0; c < ColorNum; c++){
			RGB PaintColor = ColorSet[c];
		
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


				//ストロークサイズのガウスフィルタを生成
				gauce_filter = create_dally(2*t+1, 2*t+1);
				sigma = t/3.0*opt_variance_ratio;
				sum = 0;
				for(i=0; i<2*t+1; i++){
					for(j=0; j<2*t+1; j++){
						gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
						sum += gauce_filter[i][j];
					}
				}
				for(i=0; i<2*t+1; i++){
					for(j=0; j<2*t+1; j++){
						gauce_filter[i][j] = gauce_filter[i][j] / sum;
					}
				}


				stroke_num=99999;

				// 最適なストロークに関するデータを初期化
				format_ally(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height, UNCALCULATED);

				// #ifdef _OPENMP
				// 	Format_SINGLE_Paint_Water();
				// #endif
				for(s_count=0; s_count<stroke_num; s_count++) {  //ある半径におけるストローク回数
					// 誤差探索の変数を初期化
					diff_stroke_max = UNCALCULATED;

					int x_step=opt_StrokeWindowStep, y_step=opt_StrokeWindowStep;
					// #pragma omp parallel for private(x,xc,yc,test_Canvas,diff_sum,break_flag,pnum,offscrn_count,theta,sum,former_theta) schedule(static, 1)
					for(y=0; y<in->height; y=y+y_step) {
						// #ifdef _OPENMP
						// 	test_Canvas = create_ppm(in->width, in->height, in->bright);
						// #endif
						for(x=0; x<in->width; x=x+x_step) {
							// 改善値が計算済みならSkip
							if(GLOBAL_improved_value_map->data[x][y] != UNCALCULATED) continue;

							diff_sum = break_flag = pnum = 0;

							//ウィンドウの中の差分の合計を取る
							offscrn_count = 0;
							for(xc=-t; xc<=t; xc++) {
								if((x+xc)<0 || (x+xc)>in->width-1) {offscrn_count += 2*t+1;		continue;	}
								for(yc=-t; yc<=t; yc++) {
									if((y+yc)<0 || (y+yc)>in->height-1) {	offscrn_count++;
									}else{
										//diff_sum += abs(nimgV->data[x+xc][y+yc] - cmpr->data[x+xc][y+yc]);
										diff_sum += abs(nimgR->data[x+xc][y+yc] - cmprR->data[x+xc][y+yc]);
										diff_sum += abs(nimgG->data[x+xc][y+yc] - cmprG->data[x+xc][y+yc]);
										diff_sum += abs(nimgB->data[x+xc][y+yc] - cmprB->data[x+xc][y+yc]);
									}
								}
							}
							diff_sum = diff_sum/((2*t+1)*(2*t+1)-offscrn_count)/3;

							//差分の合計平均(画素当たりの差分)が一定以上ならストローク開始位置とする
							if(diff_sum < window_diff_border) {
								GLOBAL_improved_value_map->data[x][y] = SMALL_DIFF;
								continue;
							}
							pnum=1;		//第一点確定
							best_stroke_map[x][y]->p[0].x=x+0.5; best_stroke_map[x][y]->p[0].y=y+0.5;

							// 現在の描画色の最も適切な濃度を計算
							ratio = opt_ratio;
							if(opt_USE_best_ratio){
								double max = -9999999;
								for(double roop_ratio=opt_min_ratio; roop_ratio<=opt_max_ratio; roop_ratio+=opt_ratio_step){
									if(opt_USE_Lab_PaintColorDiff){
										diff_sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[0], t, PaintColor, roop_ratio);
									}else{
										diff_sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[0], t, PaintColor, roop_ratio);
									}
									if(diff_sum>max){
										max = diff_sum;
										ratio = roop_ratio;
									}
								}
							}
							// 描画色と描画濃度を取得
							best_stroke_map[x][y]->color = PaintColor;
							best_stroke_map[x][y]->ratio = ratio;


							theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
									gauce_filter, best_stroke_map[x][y]->p[0].x, best_stroke_map[x][y]->p[0].y, t, &break_flag);


							//制御点を方向から計算し代入
							best_stroke_map[x][y]->p[1] = calcu_point(cmpr, best_stroke_map[x][y]->p[0], t, theta);


							//二つ目の描画点周りの色が描画色と一致するか確認する
							if(opt_USE_Lab_ColorDiff){
								sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
							} else{
								sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
							}


							//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
							if(sum < color_diff_border){
								//もう一つの勾配垂直の点を代入
								theta += PI;
								best_stroke_map[x][y]->p[1] = calcu_point(cmpr, best_stroke_map[x][y]->p[0], t, theta);

								//反対方向の第二点の描画点周りの色が描画色と一致するか確認する//点当たりの差異平均
								if(opt_USE_Lab_ColorDiff){
									sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
								} else{
									sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[1], t, best_stroke_map[x][y]->color, ratio);
								}

								//どちらの第二点も不適切なら描画をせず次のループへ
								if( sum < color_diff_border) {
									GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
									continue;
								}
							}

							pnum=2;		//第二点確定


							/*
								POINT2から勾配により次の制御点を探していく
							*/

							while(pnum<max_stroke){
								former_theta=theta;

								//第pnum点周りにおいて、sobelからヒストグラムを作成し最大のものを勾配とする
								theta =  calcu_histogram(cmpr, sobel_abs, sobel_angle, histogram_partition,
														gauce_filter, best_stroke_map[x][y]->p[pnum-1].x, best_stroke_map[x][y]->p[pnum-1].y, t, &break_flag);

								//制御点の為す角が急峻になるようなら逆方向に角度を取る
								if( (theta < former_theta-PI/2) || (theta > former_theta+PI/2) ) {theta += PI;}
								best_stroke_map[x][y]->p[pnum] = calcu_point(cmpr, best_stroke_map[x][y]->p[pnum-1], t, theta);

								//pnum+1目の描画点周りの色が描画色と一致するか確認する
								if(opt_USE_Lab_ColorDiff){
									sum = diffsum_Lab(in_Lab, nimgC, best_stroke_map[x][y]->p[pnum], t, best_stroke_map[x][y]->color, ratio);
								} else{
									sum = diffsum_clr_RGB(cmprC, nimgC, best_stroke_map[x][y]->p[pnum], t, best_stroke_map[x][y]->color, ratio);
								}

								/*
									pnum+1目の(次の)制御点周りの色が描画色としきい値以上の差を持つなら
									それまでの制御点を用いて線を描画
								*/
								if( sum < color_diff_border) {break_flag=1; break;}
								else {
									diff_sum += sum;
									pnum++;
								}
							}


							//best_stroke_map[x][y]->color=bright;
							best_stroke_map[x][y]->pnum=pnum;
							//　試しに描いてみて誤差を確認
							if(pnum>=min_stroke){
								// GLOBAL_improved_value_map->data[x][y] = diff_sum;
								GLOBAL_improved_value_map->data[x][y] = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[x][y], t, h, grad_hx, grad_hy, gauce_filter, ratio);
							}else if(pnum<min_stroke){
								GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
							}
						}
						// #ifdef _OPENMP
						// 	FreePPM(test_Canvas);
						// #endif
					}

					// 改善値マップ中の最大値を探索
					best_P = search_max_Point(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height);
					best_x=best_P.x;  best_y=best_P.y;
					diff_stroke_max = GLOBAL_improved_value_map->data[best_x][best_y];

					// 直近の50のストロークでの変化を表示
					diff_stroke_max_ave[s_count%50] = diff_stroke_max;
					if(s_count%50 == 49)
					{
						tmp_num=0;
						for(i=0; i<50; i++) tmp_num+=diff_stroke_max_ave[i];
						p("Stroke_MAX_AVE", tmp_num);
						tmp_num=0;
					}
					if(diff_stroke_max < optimal_improved_value_border) {
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


					//算出したpnum個の制御点を用いてストロークを描画
					Paint_Water_Stroke(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, t, best_stroke_map[best_x][best_y]->color, nimgR->data, nimgG->data, nimgB->data, h, grad_hx, grad_hy, gauce_filter, in->width, in->height, best_stroke_map[best_x][best_y]->ratio);

					// 描画したストロークの周囲のみ改善値マップをリセット
					reset_improved_value_map(GLOBAL_improved_value_map, best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, best_stroke_map, t, max_stroke);


					//一定ストロークごとに途中経過画像を書き出す
					paint_count++;
					nc++;
					// if(nc%1==0)
					if(nc%100==0 || nc<=100)
					{
							strcpy(out_filename, dir_path);
							strcat(out_filename, in_filename);
							// snprintf(count_name, 16, "%02d", t);
							// strcat(out_filename, "_t");
							// strcat(out_filename, count_name);
							snprintf(count_name, 16, "%d", nc);
							strcat(out_filename, "_s");
							strcat(out_filename, count_name);
							strcat(out_filename, ".png");
							out_png = PPM_to_image(nimgC);
							if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
							free_image(out_png);
							printf("%s\n",out_filename);
							printf("%d:",t);
							pd("TIME[s]",my_clock());
					}


					// Greedyアプローチ：最適ストロークを探索するごとに探索位置をずらす
					x_defo += 1;
					if(x_defo>t) {
						x_defo -= t;
						y_defo += 1;
						if(y_defo>t) y_defo -= t;
					}
				}

				{
					strcpy(out_filename, dir_path);
					strcat(out_filename, in_filename);
					snprintf(count_name, 16, "%02d", c);
					strcat(out_filename, "__c");
					strcat(out_filename, count_name);
					snprintf(count_name, 16, "%02d", t);
					strcat(out_filename, "_t");
					strcat(out_filename, count_name);
					strcat(out_filename, ".png");
					out_png = PPM_to_image(nimgC);
					if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
					free_image(out_png);
					printf("%s\n",out_filename);
					printf("%d:",t);
					pd("TIME[s]",my_clock());
				}

				printf("\n////////////////////\nt%d done.\n////////////////////\n\n\n",t);
				p("Paint_num",paint_count); paint_count=0;
				Free_dally(gauce_filter, 2*t+1);
				x_defo=y_defo=paint_count=0;
			}
		

			strcat(log_sentence, "\r\n///////////\r\nc[");
			Add_NUM_to_sentence(log_sentence, 2, "%d", c);
			strcat(log_sentence, "]:");
			Add_NUM_to_sentence(log_sentence, 4, "%d", PaintColor.R);
			strcat(log_sentence, ",");
			Add_NUM_to_sentence(log_sentence, 4, "%d", PaintColor.G);
			strcat(log_sentence, ",");
			Add_NUM_to_sentence(log_sentence, 4, "%d", PaintColor.B);
			strcat(log_sentence, "\r\n///////////\r\n");

			printf("\n////////////////////\nColorNUM:%d done.\n////////////////////\n\n\n",c);
		}

	}



	//第一段階描画後の中間画像を出力
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__STEP1");
	strcat(out_filename, ".png");
	out_png = PPM_to_image(nimgC);
	if(write_png_file(out_filename, out_png)){ printf("WRITE JPG ERROR.");}
	free_image(out_png);
	printf("%s\n",out_filename);
	printf("%d",t);
	pd("TIME[s]",my_clock());



	//---------------------------
	//エッジマップを計算し、エッジの複雑な周辺だけに描画を行う
	//---------------------------
	loop_cont = opt2_loop_cont;
	min_stroke = opt2_min_stroke;
	ratio = opt2_ratio;
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
					sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, ratio);
				} else{
					sum = diffsum_clr_RGB(cmprC, nimgC, p[1], t, bright, ratio);
				}

				//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
				if( sum < color_diff_border){
					//もう一つの勾配垂直の点を代入
					theta += PI;
					p[1] = calcu_point(cmpr, p[0], t, theta);

					//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
					if(opt_USE_Lab_ColorDiff){
						sum = diffsum_Lab(in_Lab, nimgC, p[1], t, bright, ratio);
					} else{
						sum = diffsum_clr_RGB(cmprC, nimgC, p[1], t, bright, ratio);
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
						sum = diffsum_Lab(in_Lab, nimgC, p[pnum], t, bright, ratio);
					} else{
						sum = diffsum_clr_RGB(cmprC, nimgC, p[pnum], t, bright, ratio);
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

					// 制御点を全てとブラシサイズを拡大率に従いスケーリングし、拡大キャンバスに描画
					if(opt_USE_Canvas_Scaling_Method){
						scaling_p = scaling_point(p, pnum, opt_canvas_scaling_ratio);
						Paint_Bezier_ex(scaling_p, pnum, nimgR_Scaling, t*opt_canvas_scaling_ratio, bright.R, ratio);
						Paint_Bezier_ex(scaling_p, pnum, nimgG_Scaling, t*opt_canvas_scaling_ratio, bright.G, ratio);
						Paint_Bezier_ex(scaling_p, pnum, nimgB_Scaling, t*opt_canvas_scaling_ratio, bright.B, ratio);
					}
				}

				paint_count++;
				nc++;
				// if(nc%100==0 || nc<=100)
				// {
				// 		strcpy(out_filename, dir_path);
				// 		strcat(out_filename, in_filename);
				// 		// snprintf(count_name, 16, "%02d", t);
				// 		// strcat(out_filename, "_t");
				// 		// strcat(out_filename, count_name);
				// 		snprintf(count_name, 16, "%d", nc);
				// 		strcat(out_filename, "_s");
				// 		strcat(out_filename, count_name);
				// 		strcat(out_filename, ".png");
				// 		out_png = PPM_to_image(nimgC);
				// 		if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
				// 		free_image(out_png);
				// 		printf("%s\n",out_filename);
				// 		printf("%d:",t);
				// 		pd("TIME[s]",my_clock());
				// }
			}
		}


		printf("////////////////////\nt%d done.\n////////////////////\n\n\n",t);
		Free_dally(gauce_filter, 2*t+1);

		{
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			snprintf(count_name, 16, "%02d", t);
			strcat(out_filename, "__st");
			strcat(out_filename, count_name);
			snprintf(count_name, 16, "%02d", lc);
			strcat(out_filename, "_lc");
			strcat(out_filename, count_name);
			strcat(out_filename, ".png");
			out_png = PPM_to_image(nimgC);
			if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
			free_image(out_png);
			printf("%s\n",out_filename);
			printf("%d:",t);
			strcat(log_sentence, "\r\n");
			Add_dictionary_to_sentence(log_sentence, "t", t);
			Add_dictionary_to_sentence(log_sentence, "s_count", paint_count);
			Add_dictionary_to_sentence_d(log_sentence, "TIME[s]", my_clock());
			pd("TIME[s]",my_clock());
		}

		if(opt_USE_Canvas_Scaling_Method){
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			strcat(out_filename, "_SC");
			snprintf(count_name, 16, "%02d", t);
			strcat(out_filename, "__st");
			strcat(out_filename, count_name);
			snprintf(count_name, 16, "%02d", lc);
			strcat(out_filename, "_lc");
			strcat(out_filename, count_name);
			strcat(out_filename, ".png");
			out_png = PPM_to_image(nimgC);
			if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
			free_image(out_png);
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
	Free_dally(h_Scaling, nimgC_Scaling->width);
	Free_dally(grad_hx_Scaling, nimgC_Scaling->width+1);
	Free_dally(grad_hy_Scaling, nimgC_Scaling->width);
	FreePGM(gray);
	FreePGM(inR);
	FreePGM(inG);
	FreePGM(inB);
	FreePGM(nimgV);
	FreePGM(nimgR_Scaling);
	FreePGM(nimgG_Scaling);
	FreePGM(nimgB_Scaling);
	free(nimgC_Scaling);
	free(nimgR);
	free(nimgG);
	free(nimgB);
	FreePPM(test_Canvas);
	if(thick_max){
		FreePGM(canny);
		FreePGM(EdgeMap);
	}
	Free_ally(x_centlabel_2D,in->width);
	Free_ally(in_Lab,in->width);
	free(ColorSet);

    return nimgC;
}