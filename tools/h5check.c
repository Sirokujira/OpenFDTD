/*
h5check.c

HDF5 ファイルに指定したデータセットが存在するかを調べる小さな道具。

    ofd_h5check <ファイル> < パス一覧   (標準入力から 1 行 1 パス)
    ofd_h5check <ファイル> <パス> [<パス> ...]
    ofd_h5check -d <ファイル> < パス     (データセットの中身を 1 行 1 値で出す)

各パスについて "OK   <パス>" か "MISS <パス>" を出力し、
1 つでも欠けていれば終了コード 1 を返す。

-d は数値データセットを double に変換して "%.10e" で 1 行ずつ出力する
(多次元でも格納順に平坦化)。data/sample/pml_reflection_check.sh が
観測点の時系列 /metadata/VPoint を読むのに使う。h5dump はツール類が
入っている環境にしか無いので、ここでも Windows / macOS を考えて自前で読む。

## パスは標準入力から渡すこと (Windows で必須)

Git Bash (MSYS) は、**先頭が '/' の引数をネイティブ Windows 実行ファイルへ
渡すとき、勝手に Windows パスへ書き換える**。HDF5 のデータセットパスは
"/metadata/Niter" のように '/' で始まるので、コマンドライン引数で渡すと

    MISS C:/Program Files/Git/metadata/Niter

のように MSYS のルートを前置された形で届き、当然すべて MISS になる。
実際に CI の Windows ジョブでこれを踏んだ。標準入力の中身は変換されない
ので、一覧は stdin から渡す。引数指定も残してあるが、これは Linux/macOS
での手動確認用と考えること (MSYS_NO_PATHCONV=1 のような環境変数に
頼る手もあるが、シェルの版によって効かないことがある)。

## なぜ h5ls ではなくこれが要るか

data/sample/hdf5_layout_check.sh は HDF5 の構成を h5ls で調べていたが、
h5ls は hdf5-tools が入っている環境にしか無い。CI の Windows ジョブは
vcpkg の静的ビルド (hdf5[core,zlib]:x64-windows-static-md) を使うので
ツール類が PATH に来ず、**Windows と macOS では HDF5 の中身が
まったく検証されていなかった**。表示側 (OpenFDTD-X) が読むデータが
特定の OS だけで欠けていても、ソルバーは normal end で終わるため
気付けない状態だった。

この道具はソルバーと同じ HDF5 ライブラリをリンクするだけなので、
HDF5 でビルドできる環境なら必ず動く。

## 実装の注意

データセットの探索に H5Ovisit / H5Literate を使わないこと。
これらは HDF5 1.10 と 1.12 以降で関数名とシグネチャが変わっており
(H5Ovisit / H5Ovisit2 / H5Ovisit3)、Linux (1.10) と vcpkg の新しい版で
ビルドが割れる。ここでは受け取ったパスを 1 つずつ H5Lexists で
辿るだけにしてある。この API はバージョン間で安定している。
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5.h"

/* "/a/b/c" が存在するか。途中のグループも 1 段ずつ確認する
   (H5Lexists は中間のリンクが無いとエラーを出すため) */
static int path_exists(hid_t file, const char *path)
{
	char buf[1024];
	size_t i;
	size_t n = 0;

	if ((path == NULL) || (path[0] != '/')) return 0;
	if (strlen(path) >= sizeof(buf)) return 0;

	/* 先頭の '/' から順に、区切りごとに部分パスを作って確認する */
	for (i = 1; ; i++) {
		if ((path[i] == '/') || (path[i] == '\0')) {
			memcpy(buf, path, i);
			buf[i] = '\0';
			if (H5Lexists(file, buf, H5P_DEFAULT) <= 0) return 0;
			n++;
			if (path[i] == '\0') break;
		}
	}

	return (n > 0);
}

/* 1 件報告して、欠けていれば 1 を返す */
static int report(hid_t file, const char *path)
{
	if (path_exists(file, path)) {
		printf("  OK   %s\n", path);
		return 0;
	}
	printf("  MISS %s\n", path);
	return 1;
}

/* 数値データセットを 1 行 1 値で出力する。成功で 0 */
static int dump(hid_t file, const char *path)
{
	if (!path_exists(file, path)) {
		fprintf(stderr, "*** no dataset : %s\n", path);
		return 1;
	}

	hid_t dset = H5Dopen2(file, path, H5P_DEFAULT);
	if (dset < 0) {
		fprintf(stderr, "*** not a dataset : %s\n", path);
		return 1;
	}

	hid_t space = H5Dget_space(dset);
	const hssize_t n = H5Sget_simple_extent_npoints(space);
	int status = 1;
	if (n >= 0) {
		double *buf = (double *)malloc(((size_t)n + 1) * sizeof(double));
		if ((buf != NULL) && (H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf) >= 0)) {
			for (hssize_t i = 0; i < n; i++) {
				printf("%.10e\n", buf[i]);
			}
			status = 0;
		}
		else {
			fprintf(stderr, "*** cannot read as double : %s\n", path);
		}
		free(buf);
	}
	H5Sclose(space);
	H5Dclose(dset);

	return status;
}

/* 標準入力から 1 行読み、行末の改行と CR を落とす。空行・EOF なら 0 */
static int read_line(char *line, int size)
{
	if (fgets(line, size, stdin) == NULL) return 0;
	size_t n = strlen(line);
	while ((n > 0) && ((line[n - 1] == '\n') || (line[n - 1] == '\r'))) {
		line[--n] = '\0';
	}
	return (n > 0);
}


int main(int argc, char **argv)
{
	hid_t file;
	int i;
	int status = 0;

	if (argc < 2) {
		fprintf(stderr, "Usage: ofd_h5check <file> < path-list\n");
		fprintf(stderr, "       ofd_h5check <file> <dataset-path> [<dataset-path> ...]\n");
		fprintf(stderr, "       ofd_h5check -d <file> < dataset-path\n");
		return 2;
	}

	/* 見つからない項目は自分で報告するので、HDF5 のエラー表示は止める */
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	/* -d : データセットの中身を出力する (パスは標準入力から。理由は上記) */
	if (!strcmp(argv[1], "-d")) {
		char line[1024];
		if (argc < 3) {
			fprintf(stderr, "Usage: ofd_h5check -d <file> < dataset-path\n");
			return 2;
		}
		file = H5Fopen(argv[2], H5F_ACC_RDONLY, H5P_DEFAULT);
		if (file < 0) {
			fprintf(stderr, "*** cannot open : %s\n", argv[2]);
			return 1;
		}
		if (read_line(line, (int)sizeof(line))) {
			status = dump(file, line);
		}
		else {
			fprintf(stderr, "*** no dataset path on stdin\n");
			status = 2;
		}
		H5Fclose(file);
		return status;
	}

	file = H5Fopen(argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file < 0) {
		fprintf(stderr, "*** cannot open : %s\n", argv[1]);
		return 1;
	}

	if (argc >= 3) {
		/* 引数で渡された場合 (Linux/macOS での手動確認用) */
		for (i = 2; i < argc; i++) {
			status |= report(file, argv[i]);
		}
	}
	else {
		/* 標準入力から 1 行 1 パス。Windows ではこちらを使うこと
		   (上の「パスは標準入力から渡すこと」を参照) */
		char line[1024];
		while (!feof(stdin)) {
			if (!read_line(line, (int)sizeof(line))) continue;
			status |= report(file, line);
		}
	}

	H5Fclose(file);

	return status;
}
