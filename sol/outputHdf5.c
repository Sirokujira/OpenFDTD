/*
outputHdf5.c

HDF5 出力 (time_series_data.h5) の共有実装。
ファイル構成と実装ごとの対応状況は include/ofd_hdf5.h を参照。

従来の実装からの変更点 (いずれも表示用途で問題になっていた点):

1. 格子点 1 点ごとに H5Sselect_hyperslab + H5Dwrite を呼んでいた。
   dipole (NN=37026) でも 220 万回の H5Dwrite になり、これが出力サイズと
   実行時間の主因だった。データセット単位の一括書き込みにする。
2. double で無圧縮だった。表示用途に倍精度は不要なので float32 + gzip に
   する (チャンクはスナップショット 1 枚単位)。
3. 時間ステップごとのグループ (/data000050 等) に書いていたのは
   cEx_r 等の「DFT の累積途中の値」で、瞬時値ではなかった。
   途中経過の部分和には物理的な意味がなく、アニメーション表示に使えない。
   瞬時値 Ex..Hz を /timeseries に追記する形に改め、DFT の結果は
   収束後に 1 度だけ /freqdomain に書く。
4. 冗長領域 (PML/Mur のマージン) を含む 1 次元配列をそのまま出していた。
   物理領域だけを (i,j,k) の 3 次元形状に詰め替える。
*/

#include "ofd.h"
#include "ofd_hdf5.h"
#include "hdf5.h"

/* gzip 圧縮レベル (6 以上は時間の割に縮まない) */
#define H5_GZIP_LEVEL (4)

static hid_t h5file = -1;   /* < 0 なら出力しない */
static int   nsnap = 0;     /* 追記したスナップショット数 */

/* 物理領域の節点数 */
static int64_t node_count(void)
{
	return (int64_t)(Nx + 1) * (Ny + 1) * (Nz + 1);
}

/* 複素 3 成分を (i,j,k,component,reim) に詰め替える (周波数 ifreq 分) */
static void pack_complex(float *dst, int64_t off,
	const float *xr, const float *xi,
	const float *yr, const float *yi,
	const float *zr, const float *zi)
{
	int i, j, k;
	int64_t m = 0;

	for (i = 0; i <= Nx; i++) {
	for (j = 0; j <= Ny; j++) {
	for (k = 0; k <= Nz; k++) {
		const int64_t n = off + NA(i, j, k);
		dst[m++] = xr[n];  dst[m++] = xi[n];
		dst[m++] = yr[n];  dst[m++] = yi[n];
		dst[m++] = zr[n];  dst[m++] = zi[n];
	}
	}
	}
}

/* 固定長データセットを一括で書く (圧縮あり) */
static void write_dataset(hid_t loc, const char *name, int rank,
	const hsize_t *dims, hid_t type, const void *data)
{
	hid_t space, dset, plist;
	int r;
	hsize_t nelem = 1;

	if (h5file < 0) return;

	for (r = 0; r < rank; r++) nelem *= dims[r];
	if (nelem <= 0) return;

	space = H5Screate_simple(rank, dims, NULL);
	plist = H5Pcreate(H5P_DATASET_CREATE);
	/* 要素数が十分あるときだけ圧縮する (小さい配列はチャンク化の方が損) */
	if (nelem >= 1024) {
		H5Pset_chunk(plist, rank, dims);
		H5Pset_deflate(plist, H5_GZIP_LEVEL);
	}
	dset = H5Dcreate(loc, name, type, space, H5P_DEFAULT, plist, H5P_DEFAULT);
	if (dset >= 0) {
		H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		H5Dclose(dset);
	}
	H5Pclose(plist);
	H5Sclose(space);
}

static void write_scalar_int(hid_t loc, const char *name, int64_t v)
{
	hid_t space, dset;
	if (h5file < 0) return;
	space = H5Screate(H5S_SCALAR);
	dset = H5Dcreate(loc, name, H5T_NATIVE_INT64, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset >= 0) {
		H5Dwrite(dset, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
		H5Dclose(dset);
	}
	H5Sclose(space);
}

static void write_scalar_double(hid_t loc, const char *name, double v)
{
	hid_t space, dset;
	if (h5file < 0) return;
	space = H5Screate(H5S_SCALAR);
	dset = H5Dcreate(loc, name, H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset >= 0) {
		H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &v);
		H5Dclose(dset);
	}
	H5Sclose(space);
}

/* 追記可能 (H5S_UNLIMITED) なデータセットを作る */
static hid_t create_extendible(hid_t loc, const char *name, int rank,
	const hsize_t *chunk, hid_t type)
{
	hsize_t dims[8], maxdims[8];
	hid_t space, plist, dset;
	int r;
	hsize_t nelem = 1;

	for (r = 0; r < rank; r++) {
		dims[r] = (r == 0) ? 0 : chunk[r];
		maxdims[r] = (r == 0) ? H5S_UNLIMITED : chunk[r];
		nelem *= chunk[r];
	}

	space = H5Screate_simple(rank, dims, maxdims);
	plist = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist, rank, chunk);
	if (nelem >= 1024) {
		H5Pset_deflate(plist, H5_GZIP_LEVEL);
	}
	dset = H5Dcreate(loc, name, type, space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Pclose(plist);
	H5Sclose(space);

	return dset;
}

/* 追記データセットの末尾に 1 枚 (先頭次元 1 個分) 足す */
static void append_slab(hid_t dset, int rank, const hsize_t *chunk,
	int index, hid_t type, const void *data)
{
	hsize_t dims[8], offset[8];
	hid_t space, memspace;
	int r;

	if (dset < 0) return;

	for (r = 0; r < rank; r++) {
		dims[r] = (r == 0) ? (hsize_t)(index + 1) : chunk[r];
		offset[r] = (r == 0) ? (hsize_t)index : 0;
	}

	H5Dset_extent(dset, dims);
	space = H5Dget_space(dset);
	H5Sselect_hyperslab(space, H5S_SELECT_SET, offset, NULL, chunk, NULL);
	memspace = H5Screate_simple(rank, chunk, NULL);
	H5Dwrite(dset, type, memspace, space, H5P_DEFAULT, data);
	H5Sclose(memspace);
	H5Sclose(space);
}

/* /timeseries の追記データセット (open 時に作り close 時に閉じる) */
static hid_t ts_group = -1;
static hid_t ts_itime = -1, ts_time = -1, ts_time_h = -1, ts_e = -1, ts_h = -1;
static hsize_t ts_field_chunk[5];
static hsize_t ts_scalar_chunk[1];

int hdf5_open(int with_timeseries)
{
	hid_t gid;

	/* hdf5 キーで出力そのものを止められる (既定は出力する)。
	   h5file < 0 のままにすることで、以降の hdf5_* はすべて何もしない。 */
	if (!Hdf5Output) {
		h5file = -1;
		return 0;
	}

	h5file = H5Fcreate(OFD_H5_FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (h5file < 0) {
		fprintf(stderr, "*** cannot create %s\n", OFD_H5_FILE);
		return 0;
	}
	nsnap = 0;

	/* /geometry : 表示に必要な格子座標と形状 */
	gid = H5Gcreate(h5file, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (gid >= 0) {
		hsize_t d[3];
		d[0] = (hsize_t)(Nx + 1); write_dataset(gid, "Xn", 1, d, H5T_NATIVE_DOUBLE, Xn);
		d[0] = (hsize_t)(Ny + 1); write_dataset(gid, "Yn", 1, d, H5T_NATIVE_DOUBLE, Yn);
		d[0] = (hsize_t)(Nz + 1); write_dataset(gid, "Zn", 1, d, H5T_NATIVE_DOUBLE, Zn);
		d[0] = (hsize_t)Nx;       write_dataset(gid, "Xc", 1, d, H5T_NATIVE_DOUBLE, Xc);
		d[0] = (hsize_t)Ny;       write_dataset(gid, "Yc", 1, d, H5T_NATIVE_DOUBLE, Yc);
		d[0] = (hsize_t)Nz;       write_dataset(gid, "Zc", 1, d, H5T_NATIVE_DOUBLE, Zc);
		if ((NGline > 0) && (Gline != NULL)) {
			d[0] = (hsize_t)NGline; d[1] = 2; d[2] = 3;
			write_dataset(gid, "Gline", 3, d, H5T_NATIVE_DOUBLE, Gline);
		}
		/* 線分ごとの材料 id。表示側が材料別に色を変えるのに使う。
		   id_t はビルド構成で幅が変わる (_ID16 なら unsigned short、
		   既定は unsigned char) ので、そのまま書くとファイルの型が
		   ビルドに依存してしまう。必ず int32 に直してから書く。 */
		if ((NGline > 0) && (MGline != NULL)) {
			int32_t *buf = (int32_t *)malloc((size_t)NGline * sizeof(int32_t));
			if (buf != NULL) {
				for (int64_t n = 0; n < NGline; n++) {
					buf[n] = (int32_t)MGline[n];
				}
				d[0] = (hsize_t)NGline;
				write_dataset(gid, "MGline", 1, d, H5T_NATIVE_INT32, buf);
				free(buf);
			}
		}
		H5Gclose(gid);
	}

	/* /timeseries : 瞬時値の追記先を用意する */
	ts_group = with_timeseries
		? H5Gcreate(h5file, "/timeseries", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)
		: -1;
	if (ts_group >= 0) {
		ts_field_chunk[0] = 1;
		ts_field_chunk[1] = (hsize_t)(Nx + 1);
		ts_field_chunk[2] = (hsize_t)(Ny + 1);
		ts_field_chunk[3] = (hsize_t)(Nz + 1);
		ts_field_chunk[4] = 3;
		ts_scalar_chunk[0] = 1;

		ts_itime = create_extendible(ts_group, "itime", 1, ts_scalar_chunk, H5T_NATIVE_INT);
		ts_time  = create_extendible(ts_group, "time",  1, ts_scalar_chunk, H5T_NATIVE_DOUBLE);
		ts_time_h = create_extendible(ts_group, "time_H", 1, ts_scalar_chunk, H5T_NATIVE_DOUBLE);
		ts_e     = create_extendible(ts_group, "E", 5, ts_field_chunk, H5T_NATIVE_FLOAT);
		ts_h     = create_extendible(ts_group, "H", 5, ts_field_chunk, H5T_NATIVE_FLOAT);
	}

	return 1;
}

/* このステップでスナップショットを出すか。
   Hdf5Output / Hdf5Interval はどちらも comm_broadcast() で全ランクに配られる
   ので、MPI でも全ランクが同じ答を返す。集約の通信をこの判定で括ることで、
   間引いているときに無駄な転送が起きない。 */
int hdf5_snapshot_enabled(int itime)
{
	if (!Hdf5Output) return 0;
	if ((Hdf5Interval > 0) && (itime % Hdf5Interval != 0)) return 0;
	return 1;
}

/* === 逐次ストリーミング集約 ===
   MPI では rank 0 が全域配列 (6 成分ぶん) を抱えると、大きなモデルで
   rank 0 のノードだけメモリを食い潰す。begin -> put を繰り返す -> commit の
   形にして、ランクから届いた分をそのまま出力用バッファへ書き込む。
   保持するのは出力そのものと同じ (Nx+1)(Ny+1)(Nz+1)*3 の 2 本だけになる。 */

static float *snap_e = NULL;   /* 出力用パックバッファ (E) */
static float *snap_h = NULL;   /* 同 (H) */

int hdf5_snapshot_begin(void)
{
	const size_t n = (size_t)node_count() * 3;

	if ((h5file < 0) || (ts_group < 0)) return 0;

	if (snap_e == NULL) {
		snap_e = (float *)malloc(n * sizeof(float));
		snap_h = (float *)malloc(n * sizeof(float));
		if ((snap_e == NULL) || (snap_h == NULL)) {
			free(snap_e); free(snap_h);
			snap_e = snap_h = NULL;
			fprintf(stderr, "*** snapshot buffer malloc error\n");
			return 0;
		}
	}

	/* 更新対象外の成分 (Ex の i=Nx 等) は書かれないので 0 のままにする */
	memset(snap_e, 0, n * sizeof(float));
	memset(snap_h, 0, n * sizeof(float));

	return 1;
}

/* 1 成分の直方体を出力バッファへ書き込む。
   src は src[(ni*i) + (nj*j) + (nk*k) + n0] で引ける配列
   (ソルバーのローカル配列でも、受信した箱でも同じ形で渡せる)。 */
void hdf5_snapshot_put(int field, int comp,
	int i0, int i1, int j0, int j1, int k0, int k1,
	const real_t *src, int64_t ni, int64_t nj, int64_t nk, int64_t n0)
{
	float *dst = (field == 0) ? snap_e : snap_h;
	const int64_t sj = (int64_t)(Ny + 1);
	const int64_t sk = (int64_t)(Nz + 1);
	int i, j, k;

	if ((dst == NULL) || (src == NULL)) return;

	/* 物理領域だけを書く (冗長領域は出力に含めない) */
	if (i0 < 0) i0 = 0;
	if (j0 < 0) j0 = 0;
	if (k0 < 0) k0 = 0;
	if (i1 > Nx) i1 = Nx;
	if (j1 > Ny) j1 = Ny;
	if (k1 > Nz) k1 = Nz;

	for (i = i0; i <= i1; i++) {
	for (j = j0; j <= j1; j++) {
	for (k = k0; k <= k1; k++) {
		const int64_t d = (((int64_t)i * sj + j) * sk + k) * 3 + comp;
		dst[d] = (float)src[(ni * i) + (nj * j) + (nk * k) + n0];
	}
	}
	}
}

/* 集めたスナップショットを 1 枚書き出す */
void hdf5_snapshot_commit(int itime, double t)
{
	if ((h5file < 0) || (ts_group < 0)) return;
	if ((snap_e == NULL) || (snap_h == NULL)) return;

	append_slab(ts_e, 5, ts_field_chunk, nsnap, H5T_NATIVE_FLOAT, snap_e);
	append_slab(ts_h, 5, ts_field_chunk, nsnap, H5T_NATIVE_FLOAT, snap_h);

	append_slab(ts_itime, 1, ts_scalar_chunk, nsnap, H5T_NATIVE_INT, &itime);
	append_slab(ts_time,  1, ts_scalar_chunk, nsnap, H5T_NATIVE_DOUBLE, &t);

	/* FDTD の leapfrog では H は E より半ステップ前の時刻の値。
	   アニメーションで E と H を重ねるときに効くので明示する。 */
	{
		const double th = t - (0.5 * Dt);
		append_slab(ts_time_h, 1, ts_scalar_chunk, nsnap, H5T_NATIVE_DOUBLE, &th);
	}

	nsnap++;
}

/* 単一プロセス版 (CPU / CUDA) : 全域を自分で持っているのでそのまま詰める。
   成分ごとの範囲は Yee 格子で実際に更新される範囲。 */
void hdf5_write_snapshot(int itime, double t,
	const real_t *ex, const real_t *ey, const real_t *ez,
	const real_t *hx, const real_t *hy, const real_t *hz)
{
	if (!hdf5_snapshot_enabled(itime)) return;
	if ((ex == NULL) || (hx == NULL)) return;
	if (!hdf5_snapshot_begin()) return;

	hdf5_snapshot_put(0, 0, 0, Nx, 0, Ny, 0, Nz, ex, Ni, Nj, Nk, N0);
	hdf5_snapshot_put(0, 1, 0, Nx, 0, Ny, 0, Nz, ey, Ni, Nj, Nk, N0);
	hdf5_snapshot_put(0, 2, 0, Nx, 0, Ny, 0, Nz, ez, Ni, Nj, Nk, N0);
	hdf5_snapshot_put(1, 0, 0, Nx, 0, Ny, 0, Nz, hx, Ni, Nj, Nk, N0);
	hdf5_snapshot_put(1, 1, 0, Nx, 0, Ny, 0, Nz, hy, Ni, Nj, Nk, N0);
	hdf5_snapshot_put(1, 2, 0, Nx, 0, Ny, 0, Nz, hz, Ni, Nj, Nk, N0);

	hdf5_snapshot_commit(itime, t);
}

void hdf5_write_freqdomain(
	const float *cex_r, const float *cex_i,
	const float *cey_r, const float *cey_i,
	const float *cez_r, const float *cez_i,
	const float *chx_r, const float *chx_i,
	const float *chy_r, const float *chy_i,
	const float *chz_r, const float *chz_i)
{
	hid_t gid;
	float *buf;
	size_t per_freq;
	int ifreq;
	hsize_t dims[6];

	if (h5file < 0) return;
	if ((NFreq2 <= 0) || (cex_r == NULL) || (chx_r == NULL)) return;

	gid = H5Gcreate(h5file, "/freqdomain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (gid < 0) return;

	dims[0] = (hsize_t)NFreq2;
	write_dataset(gid, "freq", 1, dims, H5T_NATIVE_DOUBLE, Freq2);

	/* {NFreq2, Nx+1, Ny+1, Nz+1, 3, 2} を周波数ごとに組み立てて一括で書く */
	per_freq = (size_t)node_count() * 3 * 2;
	buf = (float *)malloc((size_t)NFreq2 * per_freq * sizeof(float));
	if (buf == NULL) {
		H5Gclose(gid);
		return;
	}

	dims[0] = (hsize_t)NFreq2;
	dims[1] = (hsize_t)(Nx + 1);
	dims[2] = (hsize_t)(Ny + 1);
	dims[3] = (hsize_t)(Nz + 1);
	dims[4] = 3;
	dims[5] = 2;

	for (ifreq = 0; ifreq < NFreq2; ifreq++) {
		pack_complex(buf + (size_t)ifreq * per_freq, (int64_t)ifreq * NN,
			cex_r, cex_i, cey_r, cey_i, cez_r, cez_i);
	}
	write_dataset(gid, "E", 6, dims, H5T_NATIVE_FLOAT, buf);

	for (ifreq = 0; ifreq < NFreq2; ifreq++) {
		pack_complex(buf + (size_t)ifreq * per_freq, (int64_t)ifreq * NN,
			chx_r, chx_i, chy_r, chy_i, chz_r, chz_i);
	}
	write_dataset(gid, "H", 6, dims, H5T_NATIVE_FLOAT, buf);

	free(buf);
	H5Gclose(gid);

	/* /post : 後処理 (ofd_post -h5) 用に、ofd.out と同じ値をソルバー内部の
	   配置 (余白込みの NN) のまま書く。/freqdomain は物理領域 0..N だけだが、
	   ofd_post の節点補間 NodeH_c は境界の節点で余白 (添字 -1) の H を読むので、
	   /freqdomain からは境界の値が復元できない */
	{
		static const char *name[12] = {
			"cEx_r", "cEx_i", "cEy_r", "cEy_i", "cEz_r", "cEz_i",
			"cHx_r", "cHx_i", "cHy_r", "cHy_i", "cHz_r", "cHz_i"};
		const float *src[12] = {cex_r, cex_i, cey_r, cey_i, cez_r, cez_i,
		                        chx_r, chx_i, chy_r, chy_i, chz_r, chz_i};
		hid_t pg = H5Gcreate(h5file, OFD_H5_POST, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (pg >= 0) {
			hsize_t d2[2], ch[2];
			d2[0] = (hsize_t)NFreq2;
			d2[1] = (hsize_t)NN;
			/* チャンクは 1 周波数ぶん (HDF5 のチャンクは 4 GiB 未満の制限があるので
			   大きな問題では分ける) */
			ch[0] = 1;
			ch[1] = (d2[1] < ((hsize_t)1 << 28)) ? d2[1] : ((hsize_t)1 << 28);
			for (int b = 0; b < 12; b++) {
				if (src[b] == NULL) continue;
				hid_t space = H5Screate_simple(2, d2, NULL);
				hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
				H5Pset_chunk(plist, 2, ch);
				H5Pset_deflate(plist, H5_GZIP_LEVEL);
				hid_t dset = H5Dcreate(pg, name[b], H5T_NATIVE_FLOAT, space, H5P_DEFAULT, plist, H5P_DEFAULT);
				if (dset >= 0) {
					H5Dwrite(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, src[b]);
					H5Dclose(dset);
				}
				H5Pclose(plist);
				H5Sclose(space);
			}
			H5Gclose(pg);
		}
	}
}

/* 給電の構造体 Feed と入力電力 Pin[0] / Pin[1] を /post に追記する。
   Pin は outputChars() の calcPin() で求まる (ファイルを閉じた後) ので、
   閉じたファイルを開き直して書く (sol/outputZin.c の input_impedance と同じ方式)。
   hdf5 = 0 で出力していなければ何もしない */
void hdf5_write_post_feed(void)
{
	if (!Hdf5Output || (NFeed <= 0) || (Feed == NULL)) return;

	/* ファイルが無ければ (書き出しに失敗した等) 黙って何もしない
	   (H5Fopen は無いファイルでエラースタックを表示する) */
	FILE *probe = fopen(OFD_H5_FILE, "rb");
	if (probe == NULL) return;
	fclose(probe);

	hid_t file = H5Fopen(OFD_H5_FILE, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file < 0) return;
	hid_t pg = (H5Lexists(file, OFD_H5_POST, H5P_DEFAULT) > 0)
		? H5Gopen(file, OFD_H5_POST, H5P_DEFAULT)
		: H5Gcreate(file, OFD_H5_POST, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (pg < 0) {
		H5Fclose(file);
		return;
	}

	/* feed_t の dir は char なので、型の幅を固定するため int に直した複合型で書く */
	typedef struct {
		int dir, i, j, k;
		double dx, dy, dz, volt, delay, z0;
	} h5feed_t;
	h5feed_t *f = (h5feed_t *)malloc((size_t)NFeed * sizeof(h5feed_t));
	if (f != NULL) {
		for (int n = 0; n < NFeed; n++) {
			f[n].dir = Feed[n].dir;
			f[n].i = Feed[n].i; f[n].j = Feed[n].j; f[n].k = Feed[n].k;
			f[n].dx = Feed[n].dx; f[n].dy = Feed[n].dy; f[n].dz = Feed[n].dz;
			f[n].volt = Feed[n].volt; f[n].delay = Feed[n].delay; f[n].z0 = Feed[n].z0;
		}
		hid_t t = H5Tcreate(H5T_COMPOUND, sizeof(h5feed_t));
		H5Tinsert(t, "dir",   HOFFSET(h5feed_t, dir),   H5T_NATIVE_INT);
		H5Tinsert(t, "i",     HOFFSET(h5feed_t, i),     H5T_NATIVE_INT);
		H5Tinsert(t, "j",     HOFFSET(h5feed_t, j),     H5T_NATIVE_INT);
		H5Tinsert(t, "k",     HOFFSET(h5feed_t, k),     H5T_NATIVE_INT);
		H5Tinsert(t, "dx",    HOFFSET(h5feed_t, dx),    H5T_NATIVE_DOUBLE);
		H5Tinsert(t, "dy",    HOFFSET(h5feed_t, dy),    H5T_NATIVE_DOUBLE);
		H5Tinsert(t, "dz",    HOFFSET(h5feed_t, dz),    H5T_NATIVE_DOUBLE);
		H5Tinsert(t, "volt",  HOFFSET(h5feed_t, volt),  H5T_NATIVE_DOUBLE);
		H5Tinsert(t, "delay", HOFFSET(h5feed_t, delay), H5T_NATIVE_DOUBLE);
		H5Tinsert(t, "z0",    HOFFSET(h5feed_t, z0),    H5T_NATIVE_DOUBLE);
		hsize_t d = (hsize_t)NFeed;
		hid_t space = H5Screate_simple(1, &d, NULL);
		if (H5Lexists(pg, "Feed", H5P_DEFAULT) <= 0) {
			hid_t dset = H5Dcreate(pg, "Feed", t, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if (dset >= 0) {
				H5Dwrite(dset, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, f);
				H5Dclose(dset);
			}
		}
		H5Sclose(space);
		H5Tclose(t);
		free(f);
	}

	/* Pin[0] / Pin[1] : {NFeed * NFreq2} (calcPin が求めていれば) */
	if ((NFreq2 > 0) && (Pin[0] != NULL) && (Pin[1] != NULL)) {
		hsize_t d = (hsize_t)NFeed * (hsize_t)NFreq2;
		const char *nm[2] = {"Pin0", "Pin1"};
		for (int m = 0; m < 2; m++) {
			if (H5Lexists(pg, nm[m], H5P_DEFAULT) > 0) continue;
			hid_t space = H5Screate_simple(1, &d, NULL);
			hid_t dset = H5Dcreate(pg, nm[m], H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if (dset >= 0) {
				H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Pin[m]);
				H5Dclose(dset);
			}
			H5Sclose(space);
		}
	}

	H5Gclose(pg);
	H5Fclose(file);
}

void hdf5_write_loss(const double *ploss)
{
	hid_t gid;
	float *buf;
	int ifreq, i, j, k;
	int64_t m = 0;
	hsize_t dims[4];

	if ((h5file < 0) || (ploss == NULL) || (NFreq2 <= 0)) return;

	gid = H5Gcreate(h5file, "/loss", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (gid < 0) return;

	buf = (float *)malloc((size_t)NFreq2 * node_count() * sizeof(float));
	if (buf == NULL) {
		H5Gclose(gid);
		return;
	}

	/* ploss は {NFreq2, NN} (NN はソルバー内部の 1 次元添字) */
	for (ifreq = 0; ifreq < NFreq2; ifreq++) {
		const int64_t off = (int64_t)ifreq * NN;
		for (i = 0; i <= Nx; i++) {
		for (j = 0; j <= Ny; j++) {
		for (k = 0; k <= Nz; k++) {
			buf[m++] = (float)ploss[off + NA(i, j, k)];
		}
		}
		}
	}

	dims[0] = (hsize_t)NFreq2;
	dims[1] = (hsize_t)(Nx + 1);
	dims[2] = (hsize_t)(Ny + 1);
	dims[3] = (hsize_t)(Nz + 1);
	write_dataset(gid, "P_loss", 4, dims, H5T_NATIVE_FLOAT, buf);

	free(buf);
	H5Gclose(gid);
}

void hdf5_write_convergence(int niter, const double *eiter, const double *hiter)
{
	hid_t gid;
	int *iter;
	int n;
	hsize_t dims[1];

	if ((h5file < 0) || (niter <= 0)) return;

	gid = H5Gcreate(h5file, "/convergence", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (gid < 0) return;

	iter = (int *)malloc((size_t)niter * sizeof(int));
	if (iter != NULL) {
		for (n = 0; n < niter; n++) {
			iter[n] = n * Solver.nout;
		}
		dims[0] = (hsize_t)niter;
		write_dataset(gid, "iter", 1, dims, H5T_NATIVE_INT, iter);
		free(iter);
	}

	dims[0] = (hsize_t)niter;
	if (eiter != NULL) write_dataset(gid, "E", 1, dims, H5T_NATIVE_DOUBLE, eiter);
	if (hiter != NULL) write_dataset(gid, "H", 1, dims, H5T_NATIVE_DOUBLE, hiter);

	H5Gclose(gid);
}

void hdf5_close(void)
{
	hid_t gid;

	if (h5file < 0) return;

	/* スナップショットの詰め替えバッファを解放する。大きさは最初の
	   hdf5_snapshot_begin() の格子で決まるので、残しておくと次のファイルを
	   別の格子で書いたときにあふれる (ソルバーは 1 プロセスで 1 回しか書かないが、
	   ユニットテストは格子を変えて何度も書く) */
	free(snap_e);
	free(snap_h);
	snap_e = snap_h = NULL;

	/* /timeseries の追記データセットを閉じる */
	if (ts_itime >= 0) { H5Dclose(ts_itime); ts_itime = -1; }
	if (ts_time  >= 0) { H5Dclose(ts_time);  ts_time  = -1; }
	if (ts_time_h >= 0) { H5Dclose(ts_time_h); ts_time_h = -1; }
	if (ts_e     >= 0) { H5Dclose(ts_e);     ts_e     = -1; }
	if (ts_h     >= 0) { H5Dclose(ts_h);     ts_h     = -1; }
	if (ts_group >= 0) { H5Gclose(ts_group); ts_group = -1; }

	/* /metadata : スカラーと解析条件 */
	gid = H5Gcreate(h5file, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (gid >= 0) {
		hsize_t dims[1];

		write_scalar_int(gid, "Nx", Nx);
		write_scalar_int(gid, "Ny", Ny);
		write_scalar_int(gid, "Nz", Nz);
		write_scalar_int(gid, "NN", NN);
		write_scalar_int(gid, "Ni", Ni);
		write_scalar_int(gid, "Nj", Nj);
		write_scalar_int(gid, "Nk", Nk);
		write_scalar_int(gid, "N0", N0);
		write_scalar_int(gid, "Ntime", Ntime);
		write_scalar_int(gid, "Nsnapshot", nsnap);
		write_scalar_int(gid, "Solver_maxiter", Solver.maxiter);
		write_scalar_int(gid, "Solver_nout", Solver.nout);
		write_scalar_int(gid, "NFreq1", NFreq1);
		write_scalar_int(gid, "NFreq2", NFreq2);
		write_scalar_int(gid, "NFeed", NFeed);
		write_scalar_int(gid, "NPoint", NPoint);
		write_scalar_int(gid, "IPlanewave", IPlanewave);
		write_scalar_int(gid, "Niter", Niter);
		write_scalar_int(gid, "NGline", NGline);
		write_scalar_double(gid, "Dt", Dt);

		/* 平面波の波源パラメータ (post/readhdf5.c が復元する)。
		   planewave_t をそのまま書くと int の pol が混ざって型が壊れるので、
		   double の並びと pol を分けて書く。 */
		{
			double pw[15];
			int n = 0, m;
			pw[n++] = Planewave.theta;
			pw[n++] = Planewave.phi;
			for (m = 0; m < 3; m++) pw[n++] = Planewave.ei[m];
			for (m = 0; m < 3; m++) pw[n++] = Planewave.hi[m];
			for (m = 0; m < 3; m++) pw[n++] = Planewave.ri[m];
			for (m = 0; m < 3; m++) pw[n++] = Planewave.r0[m];
			pw[n++] = Planewave.ai;
			dims[0] = 15;
			write_dataset(gid, "Planewave", 1, dims, H5T_NATIVE_DOUBLE, pw);
			write_scalar_int(gid, "Planewave_pol", Planewave.pol);
		}

		if (NFreq1 > 0) {
			dims[0] = (hsize_t)NFreq1;
			write_dataset(gid, "Freq1", 1, dims, H5T_NATIVE_DOUBLE, Freq1);
		}
		if (NFreq2 > 0) {
			dims[0] = (hsize_t)NFreq2;
			write_dataset(gid, "Freq2", 1, dims, H5T_NATIVE_DOUBLE, Freq2);
		}
		if ((NPoint > 0) && (VPoint != NULL)) {
			dims[0] = (hsize_t)((size_t)NPoint * (Solver.maxiter + 1));
			write_dataset(gid, "VPoint", 1, dims, H5T_NATIVE_DOUBLE, VPoint);
		}

		/* 遠方界の等価面 (ofd_post / 表示側が参照する) */
		write_scalar_int(gid, "NSurface", NSurface);
		if ((NSurface > 0) && (Surface != NULL)) {
			hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
			H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "x",  HOFFSET(surface_t, x),  H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "y",  HOFFSET(surface_t, y),  H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "z",  HOFFSET(surface_t, z),  H5T_NATIVE_DOUBLE);
			H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);
			dims[0] = (hsize_t)NSurface;
			write_dataset(gid, "Surface", 1, dims, memtype, Surface);
			H5Tclose(memtype);
		}
		/* 配列の確保長は Solver.maxiter+1 (Ntime は収束で早期終了した実長)。
		   読み手が確保長で読めるよう、確保長そのままの大きさで書く。 */
		if ((NFeed > 0) && (Solver.maxiter >= 0)) {
			dims[0] = (hsize_t)((size_t)NFeed * (Solver.maxiter + 1));
			if (VFeed != NULL) write_dataset(gid, "VFeed", 1, dims, H5T_NATIVE_DOUBLE, VFeed);
			if (IFeed != NULL) write_dataset(gid, "IFeed", 1, dims, H5T_NATIVE_DOUBLE, IFeed);
		}

		/* Title は固定長文字列 */
		{
			hid_t space, type, dset;
			space = H5Screate(H5S_SCALAR);
			type = H5Tcopy(H5T_C_S1);
			H5Tset_size(type, sizeof(Title));
			dset = H5Dcreate(gid, "Title", type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
			if (dset >= 0) {
				H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
				H5Dclose(dset);
			}
			H5Tclose(type);
			H5Sclose(space);
		}

		H5Gclose(gid);
	}

	H5Fclose(h5file);
	h5file = -1;
}
