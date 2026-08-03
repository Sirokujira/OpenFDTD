/*
ofd_hdf5.h

HDF5 出力 (time_series_data.h5) の共有 API。

従来は sol/solve.c, mpi/solve.c, cuda/solve.cu, cuda_mpi/solve.cu の
4 箇所に同じ書き込みコードが重複していた。ここに集約して 1 箇所にする。

=== 画面表示側 (OpenFDTD-X) 向けのファイル構成 ===

/geometry/
    Xn {Nx+1}, Yn {Ny+1}, Zn {Nz+1}   節点座標 [m] (f8)
    Xc {Nx},   Yc {Ny},   Zc {Nz}     セル中心座標 [m] (f8)
    Gline {NGline,2,3}                形状ワイヤフレームの線分 [m] (f8)
/timeseries/                          時間領域アニメーション用 (瞬時値)
    itime  {nsnap}                    時間ステップ番号 (i4, 追記)
    time   {nsnap}                    E の時刻 [s] = (itime+1)*Dt (f8, 追記)
    time_H {nsnap}                    H の時刻 [s] = time - Dt/2  (f8, 追記)
                                      (leapfrog なので E と H は半ステップずれる)
    E     {nsnap, Nx+1, Ny+1, Nz+1, 3}  Ex,Ey,Ez の瞬時値 [V/m] (f4, 追記)
    H     {nsnap, Nx+1, Ny+1, Nz+1, 3}  Hx,Hy,Hz の瞬時値 [A/m] (f4, 追記)
/freqdomain/                          周波数領域の分布図用 (最終 DFT 結果)
    freq  {NFreq2}                    周波数 [Hz] (f8)
    E     {NFreq2, Nx+1, Ny+1, Nz+1, 3, 2}  複素振幅 (re,im) [V/m] (f4)
    H     {NFreq2, Nx+1, Ny+1, Nz+1, 3, 2}  複素振幅 (re,im) [A/m] (f4)
/loss/
    P_loss {NFreq2, Nx+1, Ny+1, Nz+1} 損失電力密度 (f4)
        p = (1/2)(σe|E|^2 + σm|H|^2)。ただし E/H は近傍界 DFT の値で、
        sol/setupDft.c が入射スペクトルで正規化しているため
        **入射振幅 1 あたりの相対量**であり W/m^3 の絶対値ではない。
        絶対値が要るときは入射電力で換算すること。
/convergence/
    iter {Niter}, E {Niter}, H {Niter}   収束履歴 (平均電磁界)
/metadata/                            スカラーと解析条件 (従来どおり)

配列はいずれも (i,j,k) の自然な 3 次元形状で書く。ソルバー内部の
1 次元添字 (n = i*Ni + j*Nj + k*Nk + N0) は冗長領域 (PML/Mur のマージン) を
含むため表示側が扱いにくい。ここで物理領域だけに詰め替える。

電磁界成分は Yee 格子上の生の値。節点 (i,j,k) に格納される成分の定義位置は
    Ex : (i,j,k)-(i+1,j,k) の稜線中点     Hx : (i,j,k) を含む x 面の中心
のように成分ごとに半セルずれる。表示側で必要なら補間すること
(周波数領域については sol/nearfield_c.c の NodeE_c/NodeH_c が節点補間の実装)。

=== 実装ごとの対応状況 ===

4 実装 (CPU / MPI / CUDA / CUDA+MPI) すべてで全グループを出力する。
MPI 版の /timeseries は mpi/comm.c の comm_snapshot() が時間ループ内で
全域を rank 0 に集めてから書く (dipole で CPU 版とビット一致を確認済み)。
/freqdomain は comm_near3d() の集約後に書く。

出力は hdf5 キーで制御できる (sol/input_data.c):
    hdf5 = <output> [interval]
      output   : 0 = 出力しない / 1 = 出力する (既定 1)
      interval : 瞬時値スナップショットの間隔 [ステップ]
                 (既定 0 = solver の nout に従う。nout の倍数に切り上げる)
*/

#ifndef _OFD_HDF5_H_
#define _OFD_HDF5_H_

/* 実装は C (sol/outputHdf5.c) だが、cuda/solve.cu と cuda_mpi/solve.cu は
   nvcc が C++ としてコンパイルする。ガードが無いと呼び出し側だけ C++ 名前修飾に
   なり "undefined reference to hdf5_open(int)" でリンクに失敗する
   (ofd.h も同じ理由で同じガードを持つ)。 */
#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

#define OFD_H5_FILE "time_series_data.h5"

/* ファイルを作成し /geometry を書く。失敗時は 0 を返し、
   以降の hdf5_* 呼び出しは何もしない。
   with_timeseries : 0 なら /timeseries を作らない (MPI 版。各ランクが
   部分領域しか持たないため、瞬時値スナップショットを出せない) */
int  hdf5_open(int with_timeseries);

/* このステップでスナップショットを出すか (hdf5 キーの output / interval)。
   MPI では全ランクが同じ答を返すので、集約の通信をこの判定で括ってよい。 */
int  hdf5_snapshot_enabled(int itime);

/* 瞬時値スナップショットを /timeseries に追記する (時間ループ内)。
   ex..hz はソルバー内部の 1 次元配列 (NA(i,j,k) で添字付け)。
   全域を 1 プロセスで持っている CPU / CUDA 版はこれを使う。 */
void hdf5_write_snapshot(int itime, double t,
	const real_t *ex, const real_t *ey, const real_t *ez,
	const real_t *hx, const real_t *hy, const real_t *hz);

/* 逐次ストリーミング集約 (MPI 版)。
   rank 0 が全域配列を抱えずに済むよう、begin -> put を繰り返す -> commit の
   形で少しずつ出力バッファへ書き込む。保持するのは出力と同じ
   (Nx+1)(Ny+1)(Nz+1)*3 の 2 本だけ。
     field : 0 = E, 1 = H
     comp  : 0 = x, 1 = y, 2 = z
     src   : src[(ni*i) + (nj*j) + (nk*k) + n0] で引ける配列
             (ソルバーのローカル配列でも、受信した箱でも同じ形で渡せる) */
int  hdf5_snapshot_begin(void);
void hdf5_snapshot_put(int field, int comp,
	int i0, int i1, int j0, int j1, int k0, int k1,
	const real_t *src, int64_t ni, int64_t nj, int64_t nk, int64_t n0);
void hdf5_snapshot_commit(int itime, double t);

/* 周波数領域の最終結果を /freqdomain に書く (時間ループ後)。
   MPI では comm_near3d() で全域に集約した g_c* を渡すこと。 */
void hdf5_write_freqdomain(
	const float *cex_r, const float *cex_i,
	const float *cey_r, const float *cey_i,
	const float *cez_r, const float *cez_i,
	const float *chx_r, const float *chx_i,
	const float *chy_r, const float *chy_i,
	const float *chz_r, const float *chz_i);

/* 損失電力密度を /loss に書く。ploss は {NFreq2, NN} の配列。
   値は入射振幅 1 あたりの相対量 (上のファイル構成の注記を参照)。
   NULL 可 (熱解析なしの場合)。 */
void hdf5_write_loss(const double *ploss);

/* 収束履歴を /convergence に書く。 */
void hdf5_write_convergence(int niter, const double *eiter, const double *hiter);

/* /metadata のスカラー類を書いてファイルを閉じる。 */
void hdf5_close(void);

#ifdef __cplusplus
}
#endif

#endif  /* _OFD_HDF5_H_ */
