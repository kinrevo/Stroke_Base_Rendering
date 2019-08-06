#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <dirent.h>
#include <png.h>
#include <jpeglib.h>
#include "ImageIO/image.h"
#include "sbr.h"
#include "sbr_opt.h"
#include "water.h"


#define p(s,a) printf("%s:%d\n",s,a)

PPM *c_Illust_brush_Water(PPM *in, char *filename);
PPM *c_Illust_brush_Water_best(PPM *in, char *filename);

PGM* GLOBAL_improved_value_map;


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
	if(opt_USE_Best_Stroke_Method){
		trans_ppm = c_Illust_brush_Water_best(in_ppm, argv[2]);
	}else{
		trans_ppm = c_Illust_brush_Water(in_ppm, argv[2]);
	}

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
		
	
	pd("TOTAL_TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	return 0;
}



//与えられたイメージからブラッシングによる絵画を作る
PPM *c_Illust_brush_Water(PPM *in, char *filename) 
{
	clock_t start = clock();
	int i,j,x,y,xc,yc,t,break_flag,pnum, offscrn_count;
	// int P2a_count, P2b_count, inapp_count, direct_count;
	// P2a_count = P2b_count = inapp_count = direct_count = 0;
	int window_diff_border = opt_window_diff_border; 	//ストローク位置探索のしきい値
	int color_diff_border = opt_color_diff_border;  	//描画色の差異のしきい値
	int max_stroke = opt_max_stroke;
	int min_stroke = opt_min_stroke;
	Point p[max_stroke];
//	Point *scaling_p;
	int stroke_histogram[max_stroke+1];
	for(i=0; i<max_stroke+1; i++){stroke_histogram[i]=0;}
	double ratio=opt_ratio;		//ストロークの濃度
	double theta, former_theta;	
	double **gauce_filter;
	double sigma, diff_sum, sum;
	int histogram_partition=opt_histogram_partition;
	// int histogram_direct[31]={}; 
	int paint_count=0, nc=0, tc=-1;
	int loop_cont=opt_loop_cont, x_defo=0, y_defo=0;
	int lc=0;
//	double canvas_scaling_ratio = 1.0;
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
	Add_dictionary_to_sentence(log_sentence, "width", in->width);
	Add_dictionary_to_sentence(log_sentence, "height", in->height);
	Add_dictionary_to_sentence(log_sentence, "thick_max", thick_max);
	Add_dictionary_to_sentence(log_sentence, "thick_min", thick_min);
	Add_dictionary_to_sentence(log_sentence, "max_stroke", max_stroke);
	Add_dictionary_to_sentence(log_sentence, "min_stroke", min_stroke);
	Add_dictionary_to_sentence(log_sentence, "window_diff_border", window_diff_border);
	Add_dictionary_to_sentence(log_sentence, "color_diff_border", color_diff_border);
	Add_dictionary_to_sentence(log_sentence, "ratio", (int)(ratio*100));
	Add_dictionary_to_sentence(log_sentence, "histogram_partition", histogram_partition);
	Add_dictionary_to_sentence(log_sentence, "loop_cont", loop_cont);
	Add_dictionary_to_sentence(log_sentence, "USE_calcu_color_bi", opt_USE_calcu_color_bi);
	Add_dictionary_to_sentence(log_sentence, "USE_gause_histogram", opt_USE_gause_histogram);
	Add_dictionary_to_sentence(log_sentence, "optimal_improved_value_border", opt_optimal_improved_value_border);
	strcat(log_sentence, "[Water Option]\r\n");
	Add_dictionary_to_sentence(log_sentence, "mhu", (int)(opt_mhu*100));
	Add_dictionary_to_sentence(log_sentence, "kappa", (int)(opt_kappa*100));
	Add_dictionary_to_sentence(log_sentence, "N", opt_N);
	Add_dictionary_to_sentence(log_sentence, "tau", (int)(opt_tau*100));
	Add_dictionary_to_sentence(log_sentence, "xi", (int)(opt_xi*100));
	Add_dictionary_to_sentence(log_sentence, "K", opt_K);
	Add_dictionary_to_sentence(log_sentence, "eta", (int)(opt_eta*100));
	Add_dictionary_to_sentence(log_sentence, "gamma", (int)(opt_gamma*100));
	Add_dictionary_to_sentence(log_sentence, "rho", (int)(opt_rho*100));
	Add_dictionary_to_sentence(log_sentence, "omega", (int)(opt_omega*100));
	Add_dictionary_to_sentence(log_sentence, "SoakTme", opt_SoakTime);
	Add_dictionary_to_sentence(log_sentence, "SoakTimeStep", (int)(opt_SoakTimeStep*100));
	Add_dictionary_to_sentence(log_sentence, "perlin_freq", (int)(opt_perlin_freq*100));
	Add_dictionary_to_sentence(log_sentence, "perlin_depth", opt_perlin_depth);
	if(opt_USE_Backrun){
		Add_dictionary_to_sentence(log_sentence, "alpha", (int)(opt_alpha*100));
		Add_dictionary_to_sentence(log_sentence, "epsilon", (int)(opt_epsilon*100));
		Add_dictionary_to_sentence(log_sentence, "delta", (int)(opt_delta*100));
		Add_dictionary_to_sentence(log_sentence, "sigma", (int)(opt_sigma*100));
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
	// PGM *nimgR_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
	// PGM *nimgG_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
	// PGM *nimgB_Scaling = create_pgm(gray->width*canvas_scaling_ratio, gray->height*canvas_scaling_ratio, gray->bright); 
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
	PPM *nimgC = create_ppm(in->width, in->height, in->bright); //実際に描画するキャンバス
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;
	// PPM *nimgC_Scaling = create_ppm(in->width*canvas_scaling_ratio, in->height*canvas_scaling_ratio, in->bright); //実際に描画するキャンバス（拡縮描画用）
	// nimgC_Scaling->dataR = nimgR_Scaling->data;
	// nimgC_Scaling->dataG = nimgG_Scaling->data;
	// nimgC_Scaling->dataB = nimgB_Scaling->data;
	
	// sobelフィルタを適応した計算結果を予め格納しておく
	double **sobel_abs = create_dally(in->width, in->height);
	double **sobel_angle = create_dally(in->width, in->height);
	sobel_calcu(gray, sobel_abs, sobel_angle);

    //Water実装
    double** h = perlin_img(in->width, in->height, opt_perlin_freq, opt_perlin_depth);
    double** grad_hx = create_dally(in->width+1, in->height); 
    double** grad_hy = create_dally(in->width, in->height+1); 
    calcu_grad_h(h, grad_hx, grad_hy, in->width, in->height);

	
	//太いストロークから順番にストロークを小さくしておおまかに絵の形を取っていく
	for(t=thick_max; t>=thick_min; t--){
		if(t!=10 && t!=6 && t!=3 ) continue;
			
		//ストロークサイズのガウスフィルタを生成
		gauce_filter = create_dally(2*t+1, 2*t+1);
		sigma = t/2.0;	// Water:3sigma -> 2
		for(i=0; i<2*t+1; i++){
	        for(j=0; j<2*t+1; j++){
	        	gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
	        }
	    }
				
		//vecデータ書き込み
		// strcpy(vec_sentence, "t");
		// snprintf(tmp_sentence, 32, "%d", t);
		// strcat(vec_sentence, tmp_sentence);
		// log_print(vec_filename, vec_sentence, "a");


			
		// for(y=0; y<in->height; y++) {  
		// 	for(x=0; x<in->width; x++) {  
		for(y=y_defo; y<in->height; y=y+t) {  //ウィンドウの大きさに合わせて
			for(x=x_defo; x<in->width; x=x+t) {  //ウィンドウをずらす距離を変えとく
			
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
					stroke_histogram[pnum]++;
					continue;
				}
				pnum=1;		//第一点確定
				p[0].x=x+0.5; p[0].y=y+0.5;
				
				//一つ目の描画領域から描画色を取得
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
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
				
				
				//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
				if(sum < color_diff_border){
					//もう一つの勾配垂直の点を代入
					theta += PI;
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //点当たりの差異平均
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					
					//どちらの第二点も不適切なら描画をせず次のループへ
					if( sum < color_diff_border) {
						stroke_histogram[pnum]++;
						continue;
					}
					else{ 
						// reversal_map[x][y]=Reversal_ON;
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
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
					
					/*
						pnum+1目の(次の)制御点周りの色が描画色としきい値以上の差を持つなら
						それまでの制御点を用いて線を描画
					*/
					if( sum < color_diff_border) {break_flag=1; break;}
					else {pnum++;}
				}

				if(pnum>=min_stroke) { 
                    Paint_Water_Stroke(p, pnum, t, bright, nimgR->data, nimgG->data, nimgB->data, h, grad_hx, grad_hy, gauce_filter, in->width, in->height);	
					stroke_histogram[pnum]++;  
					paint_count++;
				}
				if(paint_count%500==0)
       			 // if(t!=tc)
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
					// out_png = PPM_to_image(nimgC_Scaling);
					if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
					free_image(out_png);
					printf("%s\n",out_filename);
					printf("%d:",t);
					pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
				}
			}
		}
			
        nc++;

						
		// if(t!=tc)
		{
			tc=t;	
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			snprintf(count_name, 16, "%02d", t);
			strcat(out_filename, "__t");
			strcat(out_filename, count_name);
			strcat(out_filename, ".png");
			out_png = PPM_to_image(nimgC);
			// out_png = PPM_to_image(nimgC_Scaling);
			if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
			free_image(out_png);
			printf("%s\n",out_filename);
			printf("%d:",t);
			strcat(log_sentence, "\r\n");
			Add_dictionary_to_sentence(log_sentence, "t", t);
			Add_dictionary_to_sentence(log_sentence, "s_count", paint_count);
			snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
			strcat(log_sentence, count_name);
			strcat(log_sentence, "\r\n");
			pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
		}
	
		printf("\n////////////////////\nt%d done.\n////////////////////\n",t);
		p("Paint_num",paint_count); paint_count=0;
		Free_dally(gauce_filter, 2*t+1);
		x_defo=y_defo=paint_count=0;		
	}
	
	//第一段階描画後の中間画像を出力
	snprintf(count_name, 16, "%03d", t);
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__t");
	strcat(out_filename, count_name);
	strcat(out_filename, ".png");
	out_png = PPM_to_image(nimgC);
	// out_png = PPM_to_image(nimgC_Scaling);
	if(write_png_file(out_filename, out_png)){ printf("WRITE JPG ERROR.");}
	free_image(out_png);
	printf("%s\n",out_filename);
	tc=t;
	printf("%d",t);
	pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	


	//---------------------------
	//エッジマップを計算し、エッジの複雑な周辺だけに描画を行う
	//---------------------------
	loop_cont = opt2_loop_cont;
	tc=-1;
	double maxValue=0.30, minValue=0.10;
	PGM *canny;
	PGM *EdgeMap;
	
	thick_max = opt2_thick_max;
	if(thick_max){
		canny = cannyedge_detector(gray, maxValue, minValue, thick_min);
		EdgeMap = calcu_EdgeMap(canny, thick_min, sobel_angle);
		//EdgeMap = expand_Edge(canny, thick_min);
	}
	thick_min = opt2_thick_min;
	strcat(log_sentence, "[EdgeStroke Option]\r\n");
	
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
				sum = 0;
				sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
				sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
				sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
				
				//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
				if( sum < color_diff_border){
					//もう一つの勾配垂直の点を代入
					theta += PI; 
					p[1] = calcu_point(cmpr, p[0], t, theta);
					
					//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
					//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //点当たりの差異平均
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					
					//どちらの第二点も不適切なら描画をせず次のループへ
					if( sum < color_diff_border) {
						// inapp_count++;
						// stroke_histogram[pnum]++;
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
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
					
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
					// stroke_histogram[pnum]++;  
				}

				paint_count++;
				nc++;
			}
		}
		
	
		printf("////////////////////\nt%d done.\n////////////////////\n\n",t);
		Free_dally(gauce_filter, 2*t+1);
		
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
			strcat(out_filename, ".png");
			out_png = PPM_to_image(nimgC);
			// out_png = PPM_to_image(nimgC_Scaling);
			if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
			free_image(out_png);
			printf("%s\n",out_filename);
			printf("%d:",t);
			strcat(log_sentence, "\r\n");
			Add_dictionary_to_sentence(log_sentence, "t", t);
			Add_dictionary_to_sentence(log_sentence, "s_count", paint_count);
			snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
			strcat(log_sentence, count_name);
			strcat(log_sentence, "\r\n");
			pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
		}

	}


	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__EdgeMap");
	strcat(out_filename, ".pgm");
	if(write_pgm(out_filename, EdgeMap)){ printf("WRITE PNG ERROR.");}

	
	double MSE = image_MSE(nimgC, in);
	Add_dictionary_to_sentence(log_sentence, "MSE", MSE);
	
	Add_dictionary_to_sentence(log_sentence, "All_Execution_TIME", (double)(clock()-start)/CLOCKS_PER_SEC);
	printf("%s\n", log_filename);
	if(log_print(log_filename, log_sentence, "w") ){ printf("LOG_PRINTING_FAIL\n"); }
	
	
	Free_dally(sobel_abs, in->width);
	Free_dally(sobel_angle, in->width);
	FreePGM(gray);
	FreePGM(inR);
	FreePGM(inG);
	FreePGM(inB);
	FreePGM(nimgV);
    
    return nimgC;
}




//与えられたイメージからブラッシングによる絵画を作る
PPM *c_Illust_brush_Water_best(PPM *in, char *filename) 
{
	clock_t start = clock();
	int i,j,x,y,xc,yc,t,break_flag,pnum, offscrn_count;
	// int P2a_count, P2b_count, inapp_count, direct_count;
	// P2a_count = P2b_count = inapp_count = direct_count = 0;
	int window_diff_border = opt_window_diff_border; 	//ストローク位置探索のしきい値
	int color_diff_border = opt_color_diff_border;  	//描画色の差異のしきい値
	int max_stroke = opt_max_stroke;
	int min_stroke = opt_min_stroke;
	Point p[max_stroke];
//	Point *scaling_p;
	int stroke_histogram[max_stroke+1];
	for(i=0; i<max_stroke+1; i++){stroke_histogram[i]=0;}
	double ratio=opt_ratio;		//ストロークの濃度
	double theta, former_theta;	
	double **gauce_filter;
	double sigma, diff_sum, sum;
	int histogram_partition=opt_histogram_partition;
	int paint_count=0, nc=0, tc=-1;
	int loop_cont=opt_loop_cont, x_defo=0, y_defo=0;
	// int lc=0;
	// double maxValue, minValue;
//	double canvas_scaling_ratio = 1.0;
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
	Add_dictionary_to_sentence(log_sentence, "width", in->width);
	Add_dictionary_to_sentence(log_sentence, "height", in->height);
	Add_dictionary_to_sentence(log_sentence, "thick_max", thick_max);
	Add_dictionary_to_sentence(log_sentence, "thick_min", thick_min);
	Add_dictionary_to_sentence(log_sentence, "max_stroke", max_stroke);
	Add_dictionary_to_sentence(log_sentence, "min_stroke", min_stroke);
	Add_dictionary_to_sentence(log_sentence, "window_diff_border", window_diff_border);
	Add_dictionary_to_sentence(log_sentence, "color_diff_border", color_diff_border);
	Add_dictionary_to_sentence(log_sentence, "ratio", (int)(ratio*100));
	Add_dictionary_to_sentence(log_sentence, "histogram_partition", histogram_partition);
	Add_dictionary_to_sentence(log_sentence, "loop_cont", loop_cont);
	Add_dictionary_to_sentence(log_sentence, "USE_calcu_color_bi", opt_USE_calcu_color_bi);
	Add_dictionary_to_sentence(log_sentence, "USE_gause_histogram", opt_USE_gause_histogram);
	Add_dictionary_to_sentence(log_sentence, "optimal_improved_value_border", opt_optimal_improved_value_border);
	strcat(log_sentence, "[Water Option]\r\n");
	Add_dictionary_to_sentence(log_sentence, "mhu", (int)(opt_mhu*100));
	Add_dictionary_to_sentence(log_sentence, "kappa", (int)(opt_kappa*100));
	Add_dictionary_to_sentence(log_sentence, "N", opt_N);
	Add_dictionary_to_sentence(log_sentence, "tau", (int)(opt_tau*100));
	Add_dictionary_to_sentence(log_sentence, "xi", (int)(opt_xi*100));
	Add_dictionary_to_sentence(log_sentence, "K", opt_K);
	Add_dictionary_to_sentence(log_sentence, "eta", (int)(opt_eta*100));
	Add_dictionary_to_sentence(log_sentence, "gamma", (int)(opt_gamma*100));
	Add_dictionary_to_sentence(log_sentence, "rho", (int)(opt_rho*100));
	Add_dictionary_to_sentence(log_sentence, "omega", (int)(opt_omega*100));
	Add_dictionary_to_sentence(log_sentence, "SoakTme", opt_SoakTime);
	Add_dictionary_to_sentence(log_sentence, "SoakTimeStep", (int)(opt_SoakTimeStep*100));
	Add_dictionary_to_sentence(log_sentence, "perlin_freq", (int)(opt_perlin_freq*100));
	Add_dictionary_to_sentence(log_sentence, "perlin_depth", opt_perlin_depth);
	if(opt_USE_Backrun){
		Add_dictionary_to_sentence(log_sentence, "alpha", (int)(opt_alpha*100));
		Add_dictionary_to_sentence(log_sentence, "epsilon", (int)(opt_epsilon*100));
		Add_dictionary_to_sentence(log_sentence, "delta", (int)(opt_delta*100));
		Add_dictionary_to_sentence(log_sentence, "sigma", (int)(opt_sigma*100));
	}

		
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
	PPM *nimgC = create_ppm(in->width, in->height, in->bright); //実際に描画するキャンバス
	nimgC->dataR = nimgR->data;
	nimgC->dataG = nimgG->data;
	nimgC->dataB = nimgB->data;
	PPM *test_Canvas = create_ppm(in->width, in->height, in->bright); //誤差予測に使うテストキャンバス
	// PGM *improved_value_map = create_pgm(in->width, in->height, in->bright); //改善値マップ
	GLOBAL_improved_value_map = create_pgm(in->width, in->height, in->bright); //改善値マップ


	// sobelフィルタを適応した計算結果を予め格納しておく
	double **sobel_abs = create_dally(in->width, in->height);
	double **sobel_angle = create_dally(in->width, in->height);
	sobel_calcu(gray, sobel_abs, sobel_angle);
	
	//Greedyアプローチ
	int s_count,stroke_num;
	// int best_pnum;
	int diff_stroke_max=0;
	int tmp_num;
	int diff_stroke_max_ave[50]={};
	int **reversal_map = create_ally(in->width, in->height);
	int best_x=0, best_y=0;
	// int miss_stroke_count=0;
	Stroke*** best_stroke_map = create_Stroke_ally(in->width, in->height, max_stroke);
	Point best_P;
	int optimal_improved_value_border = opt_optimal_improved_value_border;
	// int diff_stroke=0;
	Point before_P={0,0};

    //Water実装
    double** h = perlin_img(in->width, in->height, opt_perlin_freq, opt_perlin_depth);
    double** grad_hx = create_dally(in->width+1, in->height); 
    double** grad_hy = create_dally(in->width, in->height+1); 
    calcu_grad_h(h, grad_hx, grad_hy, in->width, in->height);

	
	//太いストロークから順番にストロークを小さくしておおまかに絵の形を取っていく
	for(t=thick_max; t>=thick_min; t--){
		if(t!=20 && t!=15 && t!=10 && t!=6 && t!=3 ) continue;
			
		//ストロークサイズのガウスフィルタを生成
		gauce_filter = create_dally(2*t+1, 2*t+1);
		sigma = t/2.0;	// Water:3sigma -> 2
		for(i=0; i<2*t+1; i++){
	        for(j=0; j<2*t+1; j++){
	        	gauce_filter[i][j] = gause_func(i-t, j-t, sigma);
	        }
	    }
				

		
		stroke_num=9999;
		// stroke_num = in->width/t * in->height/t / (max_stroke+min_stroke)/2;
		p("stroke_num", stroke_num);
		
		format_ally(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height, UNCALCULATED);
		format_ally(reversal_map, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height, Reversal_OFF);
		
		for(s_count=0; s_count<stroke_num; s_count++) {  //ある半径におけるストローク回数
			// 誤差探索の変数を初期化
			diff_stroke_max = UNCALCULATED;
			
			for(y=0; y<in->height; y++) {  
				for(x=0; x<in->width; x++) {  
			// for(y=y_defo; y<in->height; y=y+t) {  //ウィンドウの大きさに合わせて
			// 	for(x=x_defo; x<in->width; x=x+t) {  //ウィンドウをずらす距離を変えとく
					// 改善値が計算済みならSkip
					if(GLOBAL_improved_value_map->data[x][y] != UNCALCULATED) continue;
				
					diff_sum = break_flag = pnum = 0;
					// diff_stroke=0;	//Greedyアプローチ
					
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
						stroke_histogram[pnum]++;
						continue;
					}
					pnum=1;		//第一点確定
					p[0].x=x+0.5; p[0].y=y+0.5;
					
					//一つ目の描画領域から描画色を取得
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
					sum = 0;
					sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
					sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
					sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
					diff_sum += sum;
					
					
					//二つ目の制御点周りの色が描画色としきい値以上の差を持つなら描画せず反対方向の制御点を見る
					if(sum < color_diff_border){
						//もう一つの勾配垂直の点を代入
						theta += PI;
						p[1] = calcu_point(cmpr, p[0], t, theta);
						
						//反対方向の第二点の描画点周りの色が描画色と一致するか確認する
						//sum = diffsum_clr(cmpr, nimgV, p[1], t, bright);  //点当たりの差異平均
						sum = 0;
						sum += diffsum_clr(cmprR, nimgR, p[1], t, bright.R);
						sum += diffsum_clr(cmprG, nimgG, p[1], t, bright.G);
						sum += diffsum_clr(cmprB, nimgB, p[1], t, bright.B);
						diff_sum += sum;
						
						//どちらの第二点も不適切なら描画をせず次のループへ
						if( sum < color_diff_border) {
							stroke_histogram[pnum]++;
							GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
							continue;
						}
						else{ 
							// reversal_map[x][y]=Reversal_ON;
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
						sum = 0;
						sum += diffsum_clr(cmprR, nimgR, p[pnum], t, bright.R);
						sum += diffsum_clr(cmprG, nimgG, p[pnum], t, bright.G);
						sum += diffsum_clr(cmprB, nimgB, p[pnum], t, bright.B);
						
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
					
					//　試しに描いてみて誤差を確認
					if(pnum>=min_stroke){
						// diff_stroke = test_stroke(test_Canvas, cmprR, cmprG, cmprB, nimgR, nimgG, nimgB, p, pnum, t, bright.R, bright.G, bright.B, ratio);
						GLOBAL_improved_value_map->data[x][y] = diff_sum;
						//現在のストロークが保存している開始点のものより良いなら更新
						for(i=0; i<pnum; i++){
							best_stroke_map[x][y]->p[i] = p[i];
						}
						best_stroke_map[x][y]->color=bright;
						best_stroke_map[x][y]->pnum=pnum;
					}else if(pnum<min_stroke){
						GLOBAL_improved_value_map->data[x][y] = MIN_STROKE;
					}

					// if( (diff_sum>diff_stroke_max) && (pnum>=min_stroke) ){
					// 	best_x=x; best_y=y;
					// 	diff_stroke_max = diff_sum;
					// }
				}
			}
			
			// 改善値マップ中の最大値を探索
			best_P = search_max_Point(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height);
			best_x=best_P.x;  best_y=best_P.y;
			if(best_P.x==before_P.x && best_P.y==before_P.y){
				GLOBAL_improved_value_map->data[best_x][best_y] = MIN_STROKE;
				printf("LOUP_STROKE:[%d,%d]->",best_x,best_y);
				best_P = search_max_Point(GLOBAL_improved_value_map->data, GLOBAL_improved_value_map->width, GLOBAL_improved_value_map->height);
				best_x=best_P.x;  best_y=best_P.y;
				printf("[%d,%d]\n",best_x,best_y);
			}
			before_P.x = best_P.x; before_P.y = best_P.y;
			// diff_stroke_max = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[best_x][best_y], t, h, grad_hx, grad_hy, gauce_filter);
			diff_stroke_max = GLOBAL_improved_value_map->data[best_x][best_y];
			
			// 直近の50のストロークでの変化がほとんどなければ次の半径ステップに移行
			// if(t==1)optimal_improved_value_border = 1;
			int diff_test_stroke_ave[50];
			diff_stroke_max_ave[s_count%50] = diff_stroke_max;
			diff_test_stroke_ave[s_count%50] = test_water_stroke(test_Canvas, in, nimgC, best_stroke_map[best_x][best_y], t, h, grad_hx, grad_hy, gauce_filter);
			if(s_count%50 == 49)
			{
				tmp_num=0;
				for(i=0; i<50; i++) tmp_num+=diff_stroke_max_ave[i];
				p("Stroke_MAX_AVE", tmp_num);
				tmp_num=0;
				for(i=0; i<50; i++) tmp_num+=diff_test_stroke_ave[i];
				p("test_STROKE_AVE", tmp_num);
				// printf("C:[%d,%d,%d]",best_stroke_map[best_x][best_y]->color.R, best_stroke_map[best_x][best_y]->color.G, best_stroke_map[best_x][best_y]->color.B);
				// display_Point_ally(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum);
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
				snprintf(count_name, 16, "%f", (double)(clock()-start)/CLOCKS_PER_SEC);
				strcat(log_sentence, count_name);
				strcat(log_sentence, "\r\n");
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
			Paint_Water_Stroke(best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, t, best_stroke_map[best_x][best_y]->color, nimgR->data, nimgG->data, nimgB->data, h, grad_hx, grad_hy, gauce_filter, in->width, in->height);	
			
			// 描画したストロークの周囲のみ改善値マップをリセット
			// reset_improved_value_map(GLOBAL_improved_value_map, best_stroke_map[best_x][best_y]->p, best_stroke_map[best_x][best_y]->pnum, t, max_stroke);
			
					
			//一定ストロークごとに途中経過画像を書き出す
			paint_count++;
        	nc++;	
			if(nc%100==0 || nc<=100)
			//if(t!=tc)
			{
					strcpy(out_filename, dir_path);
					strcat(out_filename, in_filename);
					snprintf(count_name, 16, "%02d", t);
					strcat(out_filename, "_t");
					strcat(out_filename, count_name);
					snprintf(count_name, 16, "%d", nc);
					strcat(out_filename, "_s");
					strcat(out_filename, count_name);
					strcat(out_filename, ".png");
					out_png = PPM_to_image(nimgC);
					// out_png = PPM_to_image(nimgC_Scaling);
					if(write_png_file(out_filename, out_png)){ printf("WRITE PNG ERROR.");}
					free_image(out_png);
					printf("%s\n",out_filename);
					printf("%d:",t);
					pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
			}
						
				
			// Greedyアプローチ：最適ストロークを探索するごとに探索位置をずらす
			x_defo += 1;
			if(x_defo>t) {
				x_defo -= t;
				y_defo += 1;
				if(y_defo>t) y_defo -= t;
			}
		}
		
		// printf("miss_stroke/t%d:%d\n",t,miss_stroke_count);

		// if(t!=tc)
		{
			// tc=t;	
			strcpy(out_filename, dir_path);
			strcat(out_filename, in_filename);
			snprintf(count_name, 16, "%02d", t);
			strcat(out_filename, "__t");
			strcat(out_filename, count_name);
			strcat(out_filename, ".png");
			out_png = PPM_to_image(nimgC);
			// out_png = PPM_to_image(nimgC_Scaling);
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
			pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
		}

		printf("\n////////////////////\nt%d done.\n////////////////////\n",t);
		p("Paint_num",paint_count); paint_count=0;
		Free_dally(gauce_filter, 2*t+1);
		x_defo=y_defo=paint_count=0;		
	}
	
	//第一段階描画後の中間画像を出力
	snprintf(count_name, 16, "%03d", tc);
	strcpy(out_filename, dir_path);
	strcat(out_filename, in_filename);
	strcat(out_filename, "__t");
	strcat(out_filename, count_name);
	strcat(out_filename, ".png");
	out_png = PPM_to_image(nimgC);
	// out_png = PPM_to_image(nimgC_Scaling);
	if(write_png_file(out_filename, out_png)){ printf("WRITE JPG ERROR.");}
	free_image(out_png);
	printf("%s\n",out_filename);
	tc=t;
	printf("%d",t);
	pd("TIME[s]",(double)(clock()-start)/CLOCKS_PER_SEC);
	

	
	double MSE = image_MSE(nimgC, in);
	Add_dictionary_to_sentence(log_sentence, "MSE", MSE);
	
	Add_dictionary_to_sentence(log_sentence, "All_Execution_TIME", (double)(clock()-start)/CLOCKS_PER_SEC);
	printf("%s\n", log_filename);
	if(log_print(log_filename, log_sentence, "w") ){ printf("LOG_PRINTING_FAIL\n"); }
	
	
	Free_dally(sobel_abs, in->width);
	Free_dally(sobel_angle, in->width);
	FreePGM(gray);
	FreePGM(inR);
	FreePGM(inG);
	FreePGM(inB);
	FreePGM(nimgV);
    
    return nimgC;
}