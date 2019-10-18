#ifndef __SBR_OPT_H_INCLUDED__
#define __SBR_OPT_H_INCLUDED__

////////////////////////////////////////
//	Stroke Base Renderingのパラメータ設定
////////////////////////////////////////

////////////////	第一段階のパラメータ	/////////////////

// ストローク半径の最大小値。この大きさの範囲でストロークを行う
#define opt_thick_max 10
#define opt_thick_min 5
// ストローク半径の種類の数（指定しないとき０）と指定の半径
#define opt_num_thick 3
#define opt_thick_assignment {10,7,5}

// ストロークの最大小長
#define opt_max_stroke 4
#define opt_min_stroke 2

// 最適ストローク手法を用いるかどうか
#define opt_USE_Best_Stroke_Method 1
    // 最適ストローク手法においてのストローク半径切り替えに用いるしきい値
    #define opt_optimal_improved_value_border 100

// ストローク開始位置のウィンドウステップ幅（x,y共通）(Raster:1.2*t,Best:1)
#define opt_StrokeWindowStep 1.2*t

// ストロークの描画の濃さ。この値が大きいほど下地の色を無視して塗りつぶす
#define opt_ratio 0.2

// キャンバスを描画途中画像から始めるかどうかとそのアドレス
#define opt_USE_input_progress_image 0
    #define opt_progress_image_address ".png"

// キャンバスのスケーリングを行うか
#define opt_USE_Canvas_Scaling_Method 1
    // キャンバスの拡大率
    #define opt_canvas_scaling_ratio 4.0

// ストロークの半径ごとの繰り返し数。多いと丁寧に塗りつぶす
#define opt_loop_cont 1

// ストローク開始位置探索における閾値。入力画像（目的画像）とキャンバスの間にこの値以上の誤差がある位置を開始位置とする
#define	opt_window_diff_border 1

// 描画改善値の計算にLabユークリッド距離を用いるか（用いなければRGBマンハッタン距離）
#define opt_USE_Lab_ColorDiff 1
    // ストローク停止位置決定における閾値。ある色で描画した際この値以下の改善値しかなければストロークを停止
    #define opt_color_diff_border 0

// 描画色の計算法においてバイラテラル距離を用いるか
#define opt_USE_calcu_color_bi 0
// 描画色の計算法においてKmeanによるカラーセットを用いるか
#define opt_USE_calcu_Kmean_ColorSet 1
    #define opt_Kmean_ClusterNum 12

// ストローク方向を決定する際の、方向の種類の数。小さいと大雑把な方向に線を引く
#define opt_histogram_partition 31
// ストローク方向の決定においてガウス重みによるヒストグラムを用いるか
#define opt_USE_gause_histogram 0



////////////////	第二段階のパラメータ	/////////////////

// ストローク半径の最大小値。この大きさの範囲でストロークを行う（最大0のとき第二段階は用いられない）
#define opt2_thick_max 3
#define opt2_thick_min 1

// ストロークの最小長
#define opt2_min_stroke 3

// ストロークの描画の濃さ。この値が大きいほど下地の色を無視して塗りつぶす
#define opt2_ratio 0.2

// ストロークの半径ごとの繰り返し数。多いと丁寧に塗りつぶす
#define opt2_loop_cont 1

#endif