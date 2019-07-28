#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "water.h"

//　試しに描いてみて誤差を確認(water)
int test_water_stroke(PPM* test_Canvas, PPM* cmpr, PPM* nimgC, Stroke* stroke, int t, double** h, double** grad_hx, double** grad_hy, double** gauce_filter)
{
	int i,j,diffsum=0;
	double left_end,right_end,upper_end,lower_end;//upper,lowerは画像の見かけの上下（値の上下ではない）
	
	//ストローク点を囲む端の座標を特定
	left_end=right_end=stroke->p[0].x;
	upper_end=lower_end=stroke->p[0].y;
	for(i=1; i<stroke->pnum; i++){
		if(stroke->p[i].x < left_end) left_end=stroke->p[i].x;
		if(right_end < stroke->p[i].x) right_end=stroke->p[i].x;
		if(stroke->p[i].y < upper_end) upper_end=stroke->p[i].y;
		if(lower_end < stroke->p[i].y) lower_end=stroke->p[i].y;
	}
	//ストローク半径分、端座標を膨張
	left_end-=2*t; right_end+=2*t; upper_end-=2*t; lower_end+=2*t;
	if(left_end<0)left_end=0; 
	if(cmpr->width <= right_end)right_end=cmpr->width-1; 
	if(upper_end<0)upper_end=0; 
	if(cmpr->height <= lower_end)lower_end=cmpr->height-1; 
	

	//テストキャンバスにおけるストロークによる影響箇所をコピー
	for(i=left_end; i<=right_end; i++) {
		for(j=upper_end; j<=lower_end; j++) {
			//if(i<0 || i>cmprR->width-1 || j<0 || j>cmprR->height-1)continue;
			test_Canvas->dataR[i][j] = nimgC->dataR[i][j];
			test_Canvas->dataG[i][j] = nimgC->dataG[i][j];
			test_Canvas->dataB[i][j] = nimgC->dataB[i][j];
		}
	}
	//試し描き
	Paint_Water_Stroke(stroke->p, stroke->pnum, t, stroke->color, test_Canvas->dataR, test_Canvas->dataG, test_Canvas->dataB, h, grad_hx, grad_hy, gauce_filter, nimgC->width, nimgC->height);
	//試し描きによる変化を確認
	for(i=left_end; i<=right_end; i++) {
        if(i<0 || i>cmpr->width-1)continue;
		for(j=upper_end; j<=lower_end; j++) {
			if(j<0 || j>cmpr->height-1)continue;
			diffsum += abs(cmpr->dataR[i][j]-nimgC->dataR[i][j]) - abs(cmpr->dataR[i][j]-test_Canvas->dataR[i][j]);
			diffsum += abs(cmpr->dataG[i][j]-nimgC->dataG[i][j]) - abs(cmpr->dataG[i][j]-test_Canvas->dataG[i][j]);
			diffsum += abs(cmpr->dataB[i][j]-nimgC->dataB[i][j]) - abs(cmpr->dataB[i][j]-test_Canvas->dataB[i][j]);
		}
	}
	
	return diffsum;
}



void Paint_Water_Stroke(Point StrokeP[], int pnum, int thick, RGB color, int** CanR, int** CanG, int** CanB, 
                            double** h, double** grad_hx, double** grad_hy, double** gauce_filter, int width, int height)
{
    int i,j;
    double** u = create_dally(width+1, height);
    double** v = create_dally(width, height+1);
    format_dally(u, width+1, height, 0);
    format_dally(v, width, height+1, 0);

    int** M = create_ally(width, height);
    double** p = create_dally(width, height);
    format_ally(M, width, height, 0);
    format_dally(p, width, height, 0);

    double** gR = create_dally(width, height);
    double** gG = create_dally(width, height);
    double** gB = create_dally(width, height);
    format_dally(gR, width, height, 0);
    format_dally(gG, width, height, 0);
    format_dally(gB, width, height, 0);

    double** dR = create_dally(width, height);
    double** dG = create_dally(width, height);
    double** dB = create_dally(width, height);
    format_dally(dR, width, height, 0);
    format_dally(dG, width, height, 0);
    format_dally(dB, width, height, 0);


    // キャンバスの色をCMYに変換し堆積顔料を計算
    for(i=0; i<width; i++) {
        for(j=0; j<height; j++) {
            dR[i][j] = 1 - CanR[i][j]/255.0;    //RGB[0,255]->CMY[0,1]
            dG[i][j] = 1 - CanG[i][j]/255.0;
            dB[i][j] = 1 - CanB[i][j]/255.0;
        }
    }

    set_WetStroke(M, p, gR, gG, gB, StrokeP, pnum, thick, color, gauce_filter, width, height);   //ストロークのエリアと水量を計算
    Paint_Water(M, u, v, p, h, grad_hx, grad_hy, gR, gG, gB, dR, dG, dB, width, height);    //水と顔料の移動を計算

    // 堆積顔料をRGBに変換しキャンバスの色を計算
    for(i=0; i<width; i++) {
        for(j=0; j<height; j++) {
            CanR[i][j] = (1 - dR[i][j]) * 255;    //CMY[0,1]->RGB[0,255]
            CanG[i][j] = (1 - dG[i][j]) * 255;
            CanB[i][j] = (1 - dB[i][j]) * 255;
        }
    }

    Free_dally(u, width+1);
    Free_dally(v, width);
    Free_ally(M, width);
    Free_dally(p, width);
    Free_dally(gR, width);
    Free_dally(gG, width);
    Free_dally(gB, width);
    Free_dally(dR, width);
    Free_dally(dG, width);
    Free_dally(dB, width);
}


// ストローク点に従いウェットエリアと水量を計算
void set_WetStroke(int** M, double** p, double** gR, double** gG, double** gB, Point SP[], int pnum, int thick, RGB color, double** gauce_filter, int width, int height)
{
    int i,j,r;
	double t;
	Point temp;  //描画線を通らない制御点

    //二点までしか与えられなければ直線を引く
	if(pnum==2){
	    int partition = abs(SP[0].x-SP[1].x) + abs(SP[0].y-SP[1].y); //線分の分割数
        for(i=0; i <= partition; i++) {
            t = (double)i/partition;
            r = thick*sin((10*t<PI/2 ? 10*t:PI/2));     //初期速度のみをｔに従い減衰
            temp.x = SP[0].x+(SP[1].x-SP[0].x)*t;
            temp.y = SP[0].y+(SP[1].y-SP[0].y)*t;
            Circle_fill_Water(M, p, gR, gG, gB, temp, r, color, gauce_filter, width, height);
        }
		return;
	}
	
	Point p0={2*SP[0].x-SP[1].x, 2*SP[0].y-SP[1].y}, 	//両端の一つ外の制御点を適当に決める
	p_np1={2*SP[pnum-1].x-SP[pnum-2].x, 2*SP[pnum-1].y-SP[pnum-2].y};
	int partition = abs(SP[0].x-SP[1].x)+abs(SP[0].y-SP[1].y); //線分の分割数
	
	Point bp0 = {(SP[1].x-p0.x)/6.0 + SP[0].x, (SP[1].y-p0.y)/6.0 + SP[0].y}
		,bp1 = {(SP[0].x-SP[2].x)/6.0 + SP[1].x, (SP[0].y-SP[2].y)/6.0 + SP[1].y};
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick*sin((10*t<PI/2 ? 10*t:PI/2));
		temp = BezierCurve_P(SP[0], bp0, bp1, SP[1], t);
        Circle_fill_Water(M, p, gR, gG, gB, temp, r, color, gauce_filter, width, height);
	}
	
	for(i=1; i<pnum-2; i++){
		bp0.x = (SP[i+1].x-SP[i-1].x)/6.0 + SP[i].x;	bp0.y = (SP[i+1].y-SP[i-1].y)/6.0 + SP[i].y;
		bp1.x = (SP[i].x-SP[i+2].x)/6.0 + SP[i+1].x;	bp1.y = (SP[i].y-SP[i+2].y)/6.0 + SP[i+1].y;
		partition = abs(SP[i].x-SP[i+1].x)+abs(SP[i].y-SP[i+1].y);
		for(j=0; j <= partition; j++) {
			t = (double)j/partition;
			r = thick;
			temp = BezierCurve_P(SP[i], bp0, bp1, SP[i+1], t);
            Circle_fill_Water(M, p, gR, gG, gB, temp, r, color, gauce_filter, width, height);
		}
	}
	

	bp0.x = (SP[pnum-1].x-SP[pnum-3].x)/6.0 + SP[pnum-2].x;
	bp0.y = (SP[pnum-1].y-SP[pnum-3].y)/6.0 + SP[pnum-2].y;
	bp1.x = (SP[pnum-2].x-p_np1.x)/6.0 + SP[pnum-1].x;
	bp1.y = (SP[pnum-2].y-p_np1.y)/6.0 + SP[pnum-1].y;
	partition = abs(SP[i].x-SP[i+1].x)+abs(SP[i].y-SP[i+1].y);
	for(i=0; i <= partition; i++) {
		t = (double)i/partition;
		r = thick;
		temp = BezierCurve_P(SP[pnum-2], bp0, bp1, SP[pnum-1], t);
        Circle_fill_Water(M, p, gR, gG, gB, temp, r, color, gauce_filter, width, height);
	}
}


//与えられた画像の座標を中心とする円に水を置く
void Circle_fill_Water(int** M, double** p, double** gR, double** gG, double** gB, Point SP, int r, RGB color, double** gauce_filter, int width, int height) {
	int x,y;

	for(x=SP.x-r; x <= SP.x+r; x++) {
		for(y=SP.y-r; y <= SP.y+r; y++) {
			if(x<0 || x>width-1 || y<0 || y>height-1) {}
			else if( (x-SP.x)*(x-SP.x)+(y-SP.y)*(y-SP.y) <= r*r ) {
                M[x][y] = 1;
                p[x][y] += gauce_filter[x-((int)SP.x-r)][y-((int)SP.y-r)]*r;    //筆の入りのときガウスフィルタの半径とｒが合わない
                gR[x][y] = 1-color.R/255.0;     //RGB[0,255]->CMY[0,1]
                gG[x][y] = 1-color.G/255.0;
                gB[x][y] = 1-color.B/255.0;
			}
		}
	}
}


// 水彩筆による顔料の移動をシミュレーション
void Paint_Water(int** M, double** u, double** v, double** p, double** h, double** grad_hx, double** grad_hy, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, int width, int height)
{
    int i,j;
    double t, var_t;
    double max=0;
    double** s = create_dally(width, height);
    format_dally(s, width, height, 0);

    char count_name[8];
    char out_name[32];
    double** dM = create_dally(width, height);
    // int paint_count=0;
    PPM* fig_img = create_ppm(width, height, 255);
    PPM* Canvas_img = create_ppm(width, height, 255);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            u[i][j] = u[i][j] - grad_hx[i][j];
            max = fmax( max, fabs(u[i][j]) );
            v[i][j] = v[i][j] - grad_hy[i][j];
            max = fmax( max, fabs(v[i][j]) );
        }
    }

    var_t = fmin(1/max, opt_SoakTimeStep);  // maxは１以下なのでおかしい・・・おかしくない？
    // pd("var_t", var_t);
    
    // format_ally(M, width, height, 1);  ///////////////////////////////////全体をウェットエリアにしても変わらない

    for ( t = 0; t < opt_SoakTime; t=t+var_t)
    {   
        UpdateVelocities(M, u, v, p, var_t, width, height);	
        if(opt_USE_MoveWater) MoveWater(M, u, v, p, var_t, width, height);
        RelaxDivergence(M, u, v, p, var_t, width, height);
        FlowOutward(M, p, var_t, width, height);
        MovePigment(M, u, v, gR, gG, gB, var_t, width, height);
        TransferPigment(M, h, gR, gG, gB, dR, dG, dB, var_t, width, height);
        if(opt_USE_Backrun) SimulateCapillaryFlow(M, p, h, s, var_t, width, height);

        for (i = 0; i < width; i++) {	
            for (j = 0; j < height; j++) {
                dM[i][j] = M[i][j];
                Canvas_img->dataR[i][j] = (1 - dR[i][j]) * 255;    //CMY[0,1]->RGB[0,255]
                Canvas_img->dataG[i][j] = (1 - dG[i][j]) * 255;
                Canvas_img->dataB[i][j] = (1 - dB[i][j]) * 255;
            }
        }

        // paint_count++;
        // if((int)(t*100)%100==0){
        //     snprintf(count_name, 16, "%02d", (int)t);
        //     strcpy(out_name, "WetArea");
        //     strcat(out_name, count_name);
        //     strcat(out_name, ".ppm");
        //     trans_Vector_img(fig_img, dM, width, height);
        //     write_ppm(out_name, fig_img);
        //     strcpy(out_name, "paperWater");
        //     strcat(out_name, count_name);
        //     strcat(out_name, ".ppm");
        //     trans_Vector_img(fig_img, p, width, height);
        //     write_ppm(out_name, fig_img);
        //     strcpy(out_name, "underWater");
        //     strcat(out_name, count_name);
        //     strcat(out_name, ".ppm");
        //     trans_Vector_img(fig_img, s, width, height);
        //     write_ppm(out_name, fig_img);
        //     strcpy(out_name, "vero_u");
        //     strcat(out_name, count_name);
        //     strcat(out_name, ".ppm");
        //     trans_Vector_img(fig_img, u, width, height);
        //     write_ppm(out_name, fig_img);
        //     strcpy(out_name, "Can");
        //     strcat(out_name, count_name);
        //     strcat(out_name, ".ppm");
        //     write_ppm(out_name, Canvas_img);
        // }
    }
}



// スタッガード格子表現から求まる値を実際の配列の値から計算
double uf(double** u, double x, double y){
    int x_LD = (int)(x*10)%10;  //少数第一位の値を整数で取得
    int y_LD = (int)(y*10)%10;
    // if(x<-0.5 || y<0){
    //    printf("uf_MINUS:%f,%f\n",x,y);
    //    return 0;
    // }
    if( x_LD == 0){
        return ( uf(u, x-0.5, y) + uf(u, x+0.5, y) )/2;
    }else if( y_LD == 5){
        return ( uf(u, x, y-0.5) + uf(u, x, y+0.5) )/2;
    }else {//if( x_LD == 5 || x_LD==-5){
        return u[(int)(x+0.5)][(int)y];
    }
    printf("uf_ERROR:%f\n",x);
    return 0;
}
double vf(double** v, double x, double y){
    int x_LD = (int)(x*10)%10;  //少数第一位の値を整数で取得
    int y_LD = (int)(y*10)%10;
    // if(x<0 || y<-0.5){
    //     printf("vf_MINUS:%f,%f\n",x,y);
    //     return 0;
    // }
    if( y_LD == 0){
        return ( vf(v, x, y-0.5) + vf(v, x, y+0.5) )/2;
    }else if( x_LD == 5){
        return ( vf(v, x-0.5, y) + vf(v, x+0.5, y) )/2;
    }else {//if( y_LD == 5 || y_LD==-5){
        return v[(int)x][(int)(y+0.5)];
    }
    printf("vf_ERROR:%f\n",y);
    return 0;
}



// 一定時間経過後の速度変化を計算
void UpdateVelocities(int** M,  double** u, double** v, double** p, double var_t, int width, int height){
    int i,j;
    double A,B;
    double** new_u = create_dally(width+1, height);
    double** new_v = create_dally(width, height+1);
    double mhu = opt_mhu;
    double kappa = opt_kappa;
    format_dally(new_u, width+1, height, 0);
    format_dally(new_v, width, height+1, 0);

    for (i = 0; i < width-1; i++){    // x:[i-0.5,i+1.5]
        for (j = 1; j < height-1; j++)    // y:[j-1.0,j+1.0]
        {
            // (uv)i+0.5,j-0.5 = uf(u, i+0.5, j)*vf(v, i, j-0.5)
            if(M[i][j]==1){
                A = pow((u[i][j]+u[i+1][j])/2, 2) - pow((u[i+1][j]+u[i+2][j])/2, 2) + (u[i+1][j-1]+u[i+1][j])/2*(v[i][j]+v[i+1][j])/2 - (u[i+1][j]+u[i+1][j+1])/2*(v[i][j+1]+v[i+1][j+1])/2;
                B = u[i+2][j] + u[i][j] + u[i+1][j+1] + u[i+1][j-1] - 4*u[i+1][j];
                new_u[i+1][j] = u[i+1][j] + var_t*(A - mhu*B + p[i][j] - p[i+1][j] - kappa*u[i+1][j]);  //u[i][j]はu[i+0.5][j]
                // double old_A = pow(uf(u,i,j), 2) - pow(uf(u,i+1.0,j), 2) + uf(u, i+0.5, j-0.5)*vf(v, i+0.5, j-0.5) - uf(u, i+0.5, j+0.5)*vf(v, i+0.5, j+0.5);
                // double old_B = uf(u,i+1.5,j) + uf(u, i-0.5, j) + uf(u, i+0.5, j+1.0) + uf(u, i+0.5, j-1.0) - 4*uf(u, i+0.5, j);
                // if(A!=old_A) printf("A:%f,oA:%f\n",A,old_A);
                // if(B!=old_B) printf("B:%f,oB:%f\n",B,old_B);
            }else if(M[i][j]==0){   //ウェットエリア外を速度０に
                new_u[i][j] = 0;
                new_u[i+1][j] = 0;
            }else{
                printf("UV_ERROR:M=%d,x=%d,y=%d\n", M[i][j], i, j);
            }
        }
    }

    for (i = 1; i < width-1; i++){
        for (j = 0; j < height-1; j++)
        {
            if(M[i][j]==1){
                A = pow((v[i][j]+v[i][j+1])/2, 2) - pow((v[i][j+1]+v[i][j+2])/2, 2) + (u[i][j]+u[i][j+1])/2*(v[i-1][j+1]+v[i][j+1])/2 - (u[i+1][j]+u[i+1][j+1])/2*(v[i][j+1]+v[i+1][j+1])/2;
                B = v[i+1][j+1] + v[i-1][j+1] + v[i][j+2] + v[i][j] - 4*v[i][j+1];
                new_v[i][j+1] = v[i][j+1] + var_t*(A - mhu*B + p[i][j] - p[i][j+1] - kappa*v[i][j+1]);  //u[i][j]はu[i+0.5][j]
                // double old_A = pow(vf(v,i,j), 2) - pow(vf(v,i,j+1.0), 2) + uf(u, i-0.5, j+0.5)*vf(v, i-0.5, j+0.5) - uf(u, i+0.5, j+0.5)*vf(v, i+0.5, j+0.5);
                // double old_B = vf(v,i+1.0,j+0.5) + vf(v, i-1.0, j+0.5) + vf(v, i, j+1.5) + vf(v, i, j-0.5) - 4*vf(v, i, j+0.5);
                // if(A!=old_A) printf("A:%f,oA:%f\n",A,old_A);
                // if(B!=old_B) printf("B:%f,oB:%f\n",B,old_B);
            }else if(M[i][j]==0){   //ウェットエリア外を速度０に
                new_v[i][j] = 0;
                new_v[i][j+1] = 0;
            }else{
                printf("UV_ERROR:M=%d,x=%d,y=%d\n", M[i][j], i, j);
            }
        }
    }

    // for (i = 0; i < width-1; i++){  //端付近は計算していないのでwidth-1にすべき？（しないと0が増殖）
    //     for (j = 0; j < height-1; j++) {
    //         if(j!=0) u[i][j] = new_u[i][j];
    //         if(i!=0) v[i][j] = new_v[i][j];
    //     }
    // }
    for (i = 0; i < width-2; i++){    // -2にしないと0が拡がる
        for (j = 1; j < height-1; j++)    // 
            {u[i+1][j] = new_u[i+1][j];}}

    for (i = 1; i < width-1; i++){
        for (j = 0; j < height-2; j++)    // -2にしないと0が拡がる
            {v[i][j+1] = new_v[i][j+1];}}

    Free_dally(new_u, width+1);
    Free_dally(new_v, width);
}


// 速度ベクトルの発散をある許容範囲τ未満になるまで緩和
void RelaxDivergence(int** M, double** u, double** v, double** p, double var_t, int width, int height)
{
    int i,j,t;
    double delta, delta_MAX;
    double delta_SUM;
    double N=opt_N, tau=opt_tau, xi=opt_xi;
    double** new_u = create_dally(width+1, height);
    double** new_v = create_dally(width, height+1);
    format_dally(new_u, width+1, height, 0);
    format_dally(new_v, width, height+1, 0);

    for (t = 0; t < N; t++)
    {
        delta_MAX = delta_SUM = 0;
        copy_dally(u, new_u, width+1, height);
        copy_dally(v, new_v, width, height+1);
        for(i=0; i<width; i++){
            for(j=0; j<height; j++){
                if(M[i][j]==1){
                    delta = xi*(u[i+1][j] - u[i][j] + v[i][j+1] - v[i][j]);
                    p[i][j] =  fmax(0, p[i][j] - delta);     // 計算符号正負不明、fmaxがないと水量が負になる
                    new_u[i+1][j] = new_u[i+1][j] - delta;
                    new_u[i][j] = new_u[i][j] + delta;
                    new_v[i][j+1] = new_v[i][j+1] - delta;
                    new_v[i][j] = new_v[i][j] + delta;
                    delta_MAX = fmax(fabs(delta), delta_MAX);
                    // if(isnan(new_u[i][j]) == 1) 
                        // printf("nan:%d,%d\n", i,j);
                }
            }
        }
        if(delta_MAX<tau) {
            // p("t",t);
            break;
        }
        // pd("deltaMAX",delta_MAX);
        copy_dally(new_u, u, width+1, height);
        copy_dally(new_v, v, width, height+1);
    }

    Free_dally(new_u, width+1);
    Free_dally(new_v, width);
}



void FlowOutward(int** M, double** p, double var_t, int width, int height)
{
    int i,j;
    int K=opt_K;
    double eta=opt_eta;
    
    double** dM = create_dally(width, height);  //関数に渡すためにdouble型に変換
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            dM[i][j] = M[i][j];
        }
    }
    
    double** gauss_M = gaussian_filter_d(dM, K, width, height);
    
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            p[i][j] = fmax(0, p[i][j] - eta*var_t*(1-gauss_M[i][j])*M[i][j]);
        }
    }

    Free_dally(gauss_M, width);
    Free_dally(dM, width);
}


void MovePigment(int** M,  double** u, double** v, double** gR, double** gG, double** gB, double var_t, int width, int height)
{
    int i,j;
    double** new_gR = create_dally(width, height);
    double** new_gG = create_dally(width, height);
    double** new_gB = create_dally(width, height);
    copy_dally(gR, new_gR, width, height);
    copy_dally(gG, new_gG, width, height);
    copy_dally(gB, new_gB, width, height);

    for(i=1; i<width-1; i++){
        for(j=1; j<height-1; j++){
            if(M[i][j]==1){
                // R
                new_gR[i+1][j] = new_gR[i+1][j] + fmax(0, var_t*u[i+1][j]*gR[i][j]);
                new_gR[i-1][j] = new_gR[i-1][j] + fmax(0, var_t*-u[i][j]*gR[i][j]);
                new_gR[i][j+1] = new_gR[i][j+1] + fmax(0, var_t*v[i][j+1]*gR[i][j]);
                new_gR[i][j-1] = new_gR[i][j-1] + fmax(0, var_t*-v[i][j]*gR[i][j]);
                // 負の値になる可能性がありそう
                new_gR[i][j] = new_gR[i][j] - fmax(0, var_t*u[i+1][j]*gR[i][j]) - fmax(0, var_t*-u[i][j]*gR[i][j]) - fmax(0, var_t*v[i][j+1]*gR[i][j]) - fmax(0, var_t*-v[i][j]*gR[i][j]);  

                // G
                new_gG[i+1][j] = new_gG[i+1][j] + fmax(0, var_t*u[i+1][j]*gG[i][j]);
                new_gG[i-1][j] = new_gG[i-1][j] + fmax(0, var_t*-u[i][j]*gG[i][j]);
                new_gG[i][j+1] = new_gG[i][j+1] + fmax(0, var_t*v[i][j+1]*gG[i][j]);
                new_gG[i][j-1] = new_gG[i][j-1] + fmax(0, var_t*-v[i][j]*gG[i][j]);
                // 負の値になる可能性がありそう
                new_gG[i][j] = new_gG[i][j] - fmax(0, var_t*u[i+1][j]*gG[i][j]) - fmax(0, var_t*-u[i][j]*gG[i][j]) - fmax(0, var_t*v[i][j+1]*gG[i][j]) - fmax(0, var_t*-v[i][j]*gG[i][j]);
                
                // B
                new_gB[i+1][j] = new_gB[i+1][j] + fmax(0, var_t*u[i+1][j]*gB[i][j]);
                new_gB[i-1][j] = new_gB[i-1][j] + fmax(0, var_t*-u[i][j]*gB[i][j]);
                new_gB[i][j+1] = new_gB[i][j+1] + fmax(0, var_t*v[i][j+1]*gB[i][j]);
                new_gB[i][j-1] = new_gB[i][j-1] + fmax(0, var_t*-v[i][j]*gB[i][j]);
                // 負の値になる可能性がありそう
                new_gB[i][j] = new_gB[i][j] - fmax(0, var_t*u[i+1][j]*gB[i][j]) - fmax(0, var_t*-u[i][j]*gB[i][j]) - fmax(0, var_t*v[i][j+1]*gB[i][j]) - fmax(0, var_t*-v[i][j]*gB[i][j]);
            }
        }
    }

    copy_dally(new_gR, gR, width, height);
    copy_dally(new_gG, gG, width, height);
    copy_dally(new_gB, gB, width, height);
    Free_dally(new_gR, width);
    Free_dally(new_gG, width);
    Free_dally(new_gB, width);
}


void TransferPigment(int** M, double** h, double** gR, double** gG, double** gB, double** dR, double** dG, double** dB, double var_t, int width, int height)
{
    int i,j;
    double down, up;
    double gamma=opt_gamma; //gammaR=0.5, gammaG=0.5, gammaB=0.5;
    double rho=opt_rho; //rhoR=0.05, rhoG=0.05, rhoB=0.05;
    double omega=opt_omega; //omegaR=1.0, omegaG=1.0, omegaB=1.0;

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            if(M[i][j]==1){
                //R
                down = gR[i][j] * var_t * (1-h[i][j]*gamma) * rho;
                up   = dR[i][j] * var_t * (1+(h[i][j]-1)*gamma) * rho / omega;
                if(dR[i][j]+down>1)
                    down = fmax(0,1-dR[i][j]);
                if(gR[i][j]+up>1)
                    up = fmax(0,1-gR[i][j]);
                dR[i][j] = dR[i][j]+down-up;
                gR[i][j] = gR[i][j]+up-down;

                //G
                down = gG[i][j] * var_t * (1-h[i][j]*gamma) * rho;
                up   = dG[i][j] * var_t * (1+(h[i][j]-1)*gamma) * rho / omega;
                if(dG[i][j]+down>1)
                    down = fmax(0,1-dG[i][j]);
                if(gG[i][j]+up>1)
                    up = fmax(0,1-gG[i][j]);
                dG[i][j] = dG[i][j]+down-up;
                gG[i][j] = gG[i][j]+up-down;

                //B
                down = gB[i][j] * var_t * (1-h[i][j]*gamma) * rho;
                up   = dB[i][j] * var_t * (1+(h[i][j]-1)*gamma) * rho / omega;
                if(dB[i][j]+down>1)
                    down = fmax(0,1-dB[i][j]);
                if(gB[i][j]+up>1)
                    up = fmax(0,1-gB[i][j]);
                dB[i][j] = dB[i][j]+down-up;
                gB[i][j] = gB[i][j]+up-down;
            }
        }
    }
}



// 水の浸透を計算しバックランをシミュレーション
void SimulateCapillaryFlow(int** M, double** p, double** c, double** s, double var_t, int width, int height)
{
    int i,j,k,l;
    double var_s;
    double alpha=opt_alpha, epsilon=opt_epsilon, delta=opt_delta, sigma=opt_sigma;
    double** new_s = create_dally(width, height);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            if(M[i][j]==1){
                s[i][j] = s[i][j] + var_t*fmax(0, fmin(alpha*p[i][j], c[i][j]-s[i][j]) );
            }
        }
    }
    copy_dally(s, new_s, width, height);

    for(i=1; i<width-1; i++){
        for(j=1; j<height-1; j++){
            if(s[i][j] < epsilon) continue;    //水量が０であるほとんどの領域をスキップ
            for(k=i-1; k<=i+1; k++){
                for(l=j-1; l<=j+1; l++){
                    if((k==i-1&&l==j-1) || (k==i-1&&l==j+1) || (k==i&&l==j) || (k==i+1&&l==j-1) || (k==i+1&&l==j+1))continue;   //４近傍以外スキップ
                    if(s[i][j] > s[k][l] && s[k][l]<delta){
                        var_s = fmax(0, fmin(s[i][j]-s[k][l], c[k][l]-s[k][l])/4);
                        new_s[i][j] = new_s[i][j] - var_s;
                        new_s[k][l] = new_s[k][l] + var_s;
                    }
                }
            }
        }
    }
    copy_dally(new_s, s, width, height);

    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            if(s[i][j] > sigma) M[i][j]=1;
        }
    }

    Free_dally(new_s, width);
}


// 自作関数
void MoveWater(int** M,  double** u, double** v, double** p, double var_t, int width, int height)
{
    int i,j;
    double** new_p = create_dally(width, height);
    copy_dally(p, new_p, width, height);

    for(i=1; i<width-1; i++){
        for(j=1; j<height-1; j++){
            if(M[i][j]==1){
                // R
                new_p[i+1][j] = new_p[i+1][j] + fmax(0, var_t*u[i+1][j]*p[i][j]);
                new_p[i-1][j] = new_p[i-1][j] + fmax(0, var_t*-u[i][j]*p[i][j]);
                new_p[i][j+1] = new_p[i][j+1] + fmax(0, var_t*v[i][j+1]*p[i][j]);
                new_p[i][j-1] = new_p[i][j-1] + fmax(0, var_t*-v[i][j]*p[i][j]);
                new_p[i][j] = new_p[i][j] - fmax(0, var_t*u[i+1][j]*p[i][j]) - fmax(0, var_t*-u[i][j]*p[i][j]) - fmax(0, var_t*v[i][j+1]*p[i][j]) - fmax(0, var_t*-v[i][j]*p[i][j]);
            }
        }
    }

    copy_dally(new_p, p, width, height);
    Free_dally(new_p, width);
}



void calcu_grad_h(double** h, double** grad_hx, double** grad_hy, int width, int height)
{
    int i,j;

    for(i=0; i<width+1; i++){
        for(j=0; j<height; j++){
            if(i==0 || i==width) grad_hx[i][j]=0;
            else grad_hx[i][j] = h[i][j] - h[i-1][j];
        }
    }

    for(i=0; i<width; i++){
        for(j=0; j<height+1; j++){
            if(j==0 || j==height) grad_hy[i][j]=0;
            else grad_hy[i][j] = h[i][j] - h[i][j-1];            
        }
    }
}



// ベクトルの正負範囲外を赤青緑の二次元で表現(画像出力はしない)
void trans_Vector_img(PPM* img, double** u, int width, int height)
{
    int i,j;
	format_ally(img->dataR, img->width, img->height, 0);
	format_ally(img->dataG, img->width, img->height, 0);
	format_ally(img->dataB, img->width, img->height, 0);
    for(i=0; i<width; i++){
        for(j=0; j<height; j++){
            if(fabs(u[i][j])>1 ) img->dataG[i][j] = 255;
            if(u[i][j]>=0) img->dataR[i][j] = u[i][j]*255;
            else if(u[i][j]<0) img->dataB[i][j] = -u[i][j]*255;
            //else printf("UnpredictionValue_ERROR:%f\n", u[i][j]);
            // if(img->dataR[i][j]<0 || img->dataB[i][j]<0) printf("%f\n%f\n",u[i][j], fabs(u[i][j]));
        }
    }
}


/////////////////////////////
//perlin関数
/////////////////////////////
static int SEED = 0;
static int hash[] = {208,34,231,213,32,248,233,56,161,78,24,140,71,48,140,254,245,255,247,247,40,
                     185,248,251,245,28,124,204,204,76,36,1,107,28,234,163,202,224,245,128,167,204,
                     9,92,217,54,239,174,173,102,193,189,190,121,100,108,167,44,43,77,180,204,8,81,
                     70,223,11,38,24,254,210,210,177,32,81,195,243,125,8,169,112,32,97,53,195,13,
                     203,9,47,104,125,117,114,124,165,203,181,235,193,206,70,180,174,0,167,181,41,
                     164,30,116,127,198,245,146,87,224,149,206,57,4,192,210,65,210,129,240,178,105,
                     228,108,245,148,140,40,35,195,38,58,65,207,215,253,65,85,208,76,62,3,237,55,89,
                     232,50,217,64,244,157,199,121,252,90,17,212,203,149,152,140,187,234,177,73,174,
                     193,100,192,143,97,53,145,135,19,103,13,90,135,151,199,91,239,247,33,39,145,
                     101,120,99,3,186,86,99,41,237,203,111,79,220,135,158,42,30,154,120,67,87,167,
                     135,176,183,191,253,115,184,21,233,58,129,233,142,39,128,211,118,137,139,255,
                     114,20,218,113,154,27,127,246,250,1,8,198,250,209,92,222,173,21,88,102,219};

int noise2(int x, int y){
    int tmp = hash[(y + SEED) % 256];
    return hash[(tmp + x) % 256];
}
float lin_inter(float x, float y, float s){
    return x + s * (y-x);
}
float smooth_inter(float x, float y, float s){
    return lin_inter(x, y, s * s * (3-2*s));
}
float noise2d(float x, float y)
{
    int x_int = x;
    int y_int = y;
    float x_frac = x - x_int;
    float y_frac = y - y_int;
    int s = noise2(x_int, y_int);
    int t = noise2(x_int+1, y_int);
    int u = noise2(x_int, y_int+1);
    int v = noise2(x_int+1, y_int+1);
    float low = smooth_inter(s, t, x_frac);
    float high = smooth_inter(u, v, x_frac);
    return smooth_inter(low, high, y_frac);
}
float perlin2d(float x, float y, float freq, int depth)
{
    float xa = x*freq;
    float ya = y*freq;
    float amp = 1.0;
    float fin = 0;
    float div = 0.0;

    int i;
    for(i=0; i<depth; i++)
    {
        div += 256 * amp;
        fin += noise2d(xa, ya) * amp;
        amp /= 2;
        xa *= 2;
        ya *= 2;
    }

    return fin/div;
}
// 本体
double** perlin_img(int width, int height, double freq, int depth)
{
    int x,y;
    double** perlin = create_dally(width, height);

    for(x=0; x<width; x++){
        for(y=0; y<height; y++){
            perlin[x][y] = perlin2d(x, y, freq, depth);//ここでノイズを調整
        }
    }

    return perlin;
}

