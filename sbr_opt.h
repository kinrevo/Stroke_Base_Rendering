#ifndef __SBR_OPT_H_INCLUDED__
#define __SBR_OPT_H_INCLUDED__

////////////////////
//	Stroke Base Renderingのパラメータ設定
////////////////////

// ストローク開始位置探索における閾値。入力画像（目的画像）とキャンバスの間にこの値以上の誤差がある位置を開始位置とする
#define	opt_window_diff_border 1

// ストローク停止位置決定における閾値。入力画像（目的画像）と描画色との間にこの値以上の誤差があればストロークを停止
#define opt_color_diff_border 0

// ストロークの最大長
#define opt_max_stroke 10

// ストロークの最小長
#define opt_min_stroke 4

// ストロークの描画の濃さ。この値が大きいほど下地の色を無視して塗りつぶす
#define opt_ratio 0.6

// ストローク方向を決定する際の、方向の種類の数。小さいと大雑把な方向に線を引く
#define opt_histogram_partition 31

// ストロークの半径ごとの繰り返し数。多いと丁寧に塗りつぶす
#define opt_loop_cont 1

// ストローク半径の最大値。この大きさからストロークを始める
#define opt_thick_max 12

// ストローク半径の最小値。この大きさでストロークを終える
#define opt_thick_min 4


#endif