#include "ofd.h"
#include "complex.h"
#include "ofd_prototype.h"
#include "ev.h"

#include "hdf5.h"

#define FILE_NAME "time_series_data.h5"
#define GROUP_NAME "/metadata"
#define DATA_GROUP_PREFIX "/data"
//#define TIME_STEP 200 // Solver.nout
//#define MAX_STEP 1000 // Solver.maxiter

#define SPARA_DATASET_NAME "s_parameters"
#define ZIN_DATASET_NAME "input_impedance"
typedef struct {
    double magnitude_dB;
    double phase_deg;
} spara_data_t;

typedef struct {
    double frequency;
    double rin;
    double xin;
    double gin;
    double bin;
    double ref;
    double vswr;
} input_impedance_data_t;

int readhdf5(void) {
    // HDF5ファイルのオープン
    hid_t file_id = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    herr_t status;
    if (file_id < 0) {
        fprintf(stderr, "*** %s : cannot open\n", FILE_NAME);
        return 1;
    }

    // /metadata グループのオープン
    hid_t metadata_group_id = H5Gopen(file_id, "/metadata", H5P_DEFAULT);

    // Titleの読み込み
    // 書き手 (sol/outputHdf5.c) は sizeof(Title) の固定長文字列で書く。
    // メモリ側も同じ長さの文字列型で読む (以前は H5T_NATIVE_CHAR で読んでいて、
    // 文字列型から整数型への変換は無いので H5Dread が失敗し Title が空のままだった)
    hid_t dataset_id = H5Dopen(metadata_group_id, "Title", H5P_DEFAULT);
    {
        hid_t strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(strtype, sizeof(Title));
        status = H5Dread(dataset_id, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
        H5Tclose(strtype);
        Title[sizeof(Title) - 1] = '\0';
    }
    H5Dclose(dataset_id);

    // 各種整数型メタデータの読み込み
    struct {
        const char *name;
        void *value;
        hid_t type;
    } metadata[] = {
        // Ni/Nj/Nk/N0/NN/NGline は int64_t なので INT64 で読むこと
        // (INT で読むと 4 バイトしか埋まらず上位が壊れる)
        {"Nx", &Nx, H5T_NATIVE_INT},
        {"Ny", &Ny, H5T_NATIVE_INT},
        {"Nz", &Nz, H5T_NATIVE_INT},
        {"Ni", &Ni, H5T_NATIVE_INT64},
        {"Nj", &Nj, H5T_NATIVE_INT64},
        {"Nk", &Nk, H5T_NATIVE_INT64},
        {"N0", &N0, H5T_NATIVE_INT64},
        {"NN", &NN, H5T_NATIVE_INT64},
        {"NFreq1", &NFreq1, H5T_NATIVE_INT},
        {"NFreq2", &NFreq2, H5T_NATIVE_INT},
        {"NFeed", &NFeed, H5T_NATIVE_INT},
        {"NPoint", &NPoint, H5T_NATIVE_INT},
        {"Niter", &Niter, H5T_NATIVE_INT},
        {"Ntime", &Ntime, H5T_NATIVE_INT},
        {"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
        {"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
        {"NGline", &NGline, H5T_NATIVE_INT64},
        {"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
    };

    for (int i = 0; i < (int)(sizeof(metadata) / sizeof(metadata[0])); i++) {
        dataset_id = H5Dopen(metadata_group_id, metadata[i].name, H5P_DEFAULT);
        status = H5Dread(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
    }

    // Dtの読み込み
    dataset_id = H5Dopen(metadata_group_id, "Dt", H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);

    // Planewave の読み込み
    // planewave_t をそのまま読むと int の pol が混ざって型が壊れるので、
    // double の並び (15 個) と pol を分けて持たせている (sol/outputHdf5.c)
    {
        double pw[15];
        int n = 0, m;
        dataset_id = H5Dopen(metadata_group_id, "Planewave", H5P_DEFAULT);
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pw);
        H5Dclose(dataset_id);
        Planewave.theta = pw[n++];
        Planewave.phi   = pw[n++];
        for (m = 0; m < 3; m++) Planewave.ei[m] = pw[n++];
        for (m = 0; m < 3; m++) Planewave.hi[m] = pw[n++];
        for (m = 0; m < 3; m++) Planewave.ri[m] = pw[n++];
        for (m = 0; m < 3; m++) Planewave.r0[m] = pw[n++];
        Planewave.ai = pw[n++];

        dataset_id = H5Dopen(metadata_group_id, "Planewave_pol", H5P_DEFAULT);
        status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave.pol);
        H5Dclose(dataset_id);
    }

    // 配列データの読み込み
    //
    // 表示側が扱いやすいようグループを用途別に分けてある (include/ofd_hdf5.h)。
    // 格子座標と形状は /geometry、収束履歴は /convergence、
    // 解析条件は /metadata。ここは書き手 (sol/outputHdf5.c) と対になっている
    // ので、片方を変えたら必ず両方直すこと。
    hid_t geometry_group_id = H5Gopen(file_id, "/geometry", H5P_DEFAULT);
    hid_t convergence_group_id = H5Gopen(file_id, "/convergence", H5P_DEFAULT);

    struct {
        hid_t group;
        const char *name;
        double **data;
        size_t size;
    } arrays[] = {
        {geometry_group_id,    "Xn",     &Xn,     Nx + 1},
        {geometry_group_id,    "Yn",     &Yn,     Ny + 1},
        {geometry_group_id,    "Zn",     &Zn,     Nz + 1},
        {geometry_group_id,    "Xc",     &Xc,     Nx},
        {geometry_group_id,    "Yc",     &Yc,     Ny},
        {geometry_group_id,    "Zc",     &Zc,     Nz},
        {convergence_group_id, "E",      &Eiter,  Niter},
        {convergence_group_id, "H",      &Hiter,  Niter},
        {metadata_group_id,    "VFeed",  &VFeed,  NFeed * (Solver.maxiter + 1)},
        {metadata_group_id,    "IFeed",  &IFeed,  NFeed * (Solver.maxiter + 1)},
        {metadata_group_id,    "VPoint", &VPoint, NPoint * (Solver.maxiter + 1)},
        {metadata_group_id,    "Freq1",  &Freq1,  NFreq1},
        {metadata_group_id,    "Freq2",  &Freq2,  NFreq2}
    };

    for (int i = 0; i < (int)(sizeof(arrays) / sizeof(arrays[0])); i++) {
        if ((arrays[i].group < 0) || (arrays[i].size <= 0)) continue;
        *arrays[i].data = (double *)malloc(sizeof(double) * arrays[i].size);
        dataset_id = H5Dopen(arrays[i].group, arrays[i].name, H5P_DEFAULT);
        if (dataset_id < 0) {
            fprintf(stderr, "*** %s : dataset not found\n", arrays[i].name);
            continue;
        }
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *arrays[i].data);
        H5Dclose(dataset_id);
    }

    // Gline は double[2][3] の配列なので上の表 (double**) に入れられない
    if ((geometry_group_id >= 0) && (NGline > 0)) {
        Gline = (double (*)[2][3])malloc(sizeof(double) * NGline * 2 * 3);
        dataset_id = H5Dopen(geometry_group_id, "Gline", H5P_DEFAULT);
        if (dataset_id >= 0) {
            status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Gline);
            H5Dclose(dataset_id);
        }

        // MGline (線分ごとの材料 id) はファイル上は i4。
        // ソルバー内部の id_t はビルド構成で幅が変わるので、
        // i4 で受けてから id_t へ詰め替える (直接 H5Dread しない)。
        MGline = (id_t *)malloc(sizeof(id_t) * NGline);
        dataset_id = H5Dopen(geometry_group_id, "MGline", H5P_DEFAULT);
        if ((dataset_id >= 0) && (MGline != NULL)) {
            int32_t *buf = (int32_t *)malloc(sizeof(int32_t) * NGline);
            if (buf != NULL) {
                status = H5Dread(dataset_id, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
                for (int64_t n = 0; n < NGline; n++) {
                    MGline[n] = (id_t)buf[n];
                }
                free(buf);
            }
            H5Dclose(dataset_id);
        }
        else if (MGline != NULL) {
            // 古いファイルには MGline が無い。全部真空扱いにして落ちないようにする
            memset(MGline, 0, sizeof(id_t) * NGline);
        }
    }

    if (geometry_group_id >= 0) H5Gclose(geometry_group_id);
    if (convergence_group_id >= 0) H5Gclose(convergence_group_id);

    //サイズ確保
    Feed   =         (feed_t *)malloc(sizeof(feed_t)      * NFeed);
    Zin    =    (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
    Ref    =         (double *)malloc(sizeof(double)      * NFeed * NFreq1);
    Pin[0] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
    Pin[1] =         (double *)malloc(sizeof(double)      * NFeed * NFreq2);
    Spara  =    (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
    // Gline / MGline はここで確保しない。
    // 上の /geometry の読み込みで確保済みで、ここで malloc し直すと
    // **読んだ線分データを捨てて未初期化のバッファに差し替えてしまう**
    // (かつ元の領域が漏れる)。形状が無い (NGline == 0) 場合は
    // どちらも NULL のままで、描画側は NGline で判定する。

	//plot2d?
	size_t size;
    size = NN * NFreq2 * sizeof(float);
    cEx_r = (float *)malloc(size);
    cEx_i = (float *)malloc(size);
    cEy_r = (float *)malloc(size);
    cEy_i = (float *)malloc(size);
    cEz_r = (float *)malloc(size);
    cEz_i = (float *)malloc(size);
    cHx_r = (float *)malloc(size);
    cHx_i = (float *)malloc(size);
    cHy_r = (float *)malloc(size);
    cHy_i = (float *)malloc(size);
    cHz_r = (float *)malloc(size);
    cHz_i = (float *)malloc(size);

    // 遠方界の等価面
    //
    // ソルバーは /metadata を書く時点 (hdf5_close、時間ループの直後) ではまだ
    // 遠方界を設定していない (alloc_farfield は outputChars の中) ので、通常の
    // ファイルの NSurface は 0 で Surface も無い。面の幾何と場は ofd_post が
    // 近傍界 DFT から setup_farfield() で作り直すので、ここでは ofd.out の読み手
    // (readout) と同じ大きさ (alloc_farfield の式) で確保だけする。
    // NSurface > 0 と Surface が書かれているファイルではそれを読む。
    NSurface = 0;
    if (H5Lexists(metadata_group_id, "NSurface", H5P_DEFAULT) > 0) {
        int64_t ns = 0;
        dataset_id = H5Dopen(metadata_group_id, "NSurface", H5P_DEFAULT);
        status = H5Dread(dataset_id, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ns);
        H5Dclose(dataset_id);
        NSurface = ns;
    }
    if ((NSurface > 0) && (H5Lexists(metadata_group_id, "Surface", H5P_DEFAULT) > 0)) {
        // surface_t構造体に対応する複合データ型を再定義
        hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
        H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
        H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);
        if (Surface == NULL) {
            Surface = (surface_t *)malloc(NSurface * sizeof(surface_t));
        }
        dataset_id = H5Dopen(metadata_group_id, "Surface", H5P_DEFAULT);
        status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
        if (status < 0) {
            fprintf(stderr, "Error reading dataset: Surface\n");
        }
        H5Dclose(dataset_id);
        H5Tclose(memtype);

        size = NFreq2 * sizeof(d_complex_t *);
        SurfaceEx = (d_complex_t **)malloc(size);
        SurfaceEy = (d_complex_t **)malloc(size);
        SurfaceEz = (d_complex_t **)malloc(size);
        SurfaceHx = (d_complex_t **)malloc(size);
        SurfaceHy = (d_complex_t **)malloc(size);
        SurfaceHz = (d_complex_t **)malloc(size);
        for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
            size = NSurface * sizeof(d_complex_t);
            SurfaceEx[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
            SurfaceEy[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
            SurfaceEz[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
            SurfaceHx[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
            SurfaceHy[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
            SurfaceHz[ifreq] = (d_complex_t *)calloc(NSurface, sizeof(d_complex_t));
        }
    }
    else if ((Nx > 0) && (Ny > 0) && (Nz > 0)) {
        alloc_farfield();
    }

//
    // coupling (/metadata/coupling) と散乱断面積 (/metadata/cross_section) は
    // 読まない。ofd_post は coupling の保存値を使わず給電・観測点の時系列から
    // その場で計算し直し (post/plot2dFreq.c が sol/outputCoupling.c の coupling()
    // を呼ぶ)、散乱断面積はそもそも使わない (ソルバーの ofd.log と表示側向け)。
    // 以前はここで coupling を読んで何にも使わずに捨てていた

    // s-parameters データの読み込み
    spara_data_t *spara_data =
        (spara_data_t *)malloc((size_t)NFreq1 * NPoint * sizeof(spara_data_t));
    if (spara_data == NULL) {
        fprintf(stderr, "Error: malloc failed (spara_data)\n");
        return 1;
    }

    // データタイプを作成
    hid_t datatype_id = H5Tcreate(H5T_COMPOUND, sizeof(spara_data_t));
    H5Tinsert(datatype_id, "magnitude_dB", HOFFSET(spara_data_t, magnitude_dB), H5T_NATIVE_DOUBLE);
    H5Tinsert(datatype_id, "phase_deg", HOFFSET(spara_data_t, phase_deg), H5T_NATIVE_DOUBLE);

    // 観測点が無い (または frequency1 が無い) と書き手は作らない
    const int has_spara = (NPoint > 0) && (NFreq1 > 0) && (H5Lexists(metadata_group_id, SPARA_DATASET_NAME, H5P_DEFAULT) > 0);
    dataset_id = has_spara ? H5Dopen(metadata_group_id, SPARA_DATASET_NAME, H5P_DEFAULT) : -1;
    status = has_spara ? H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, spara_data) : 0;
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", SPARA_DATASET_NAME);
    }
    // spara_data から Spara へ復元 (dB/deg → 複素数)
    if (Spara == NULL) {
        Spara = (d_complex_t *)malloc(sizeof(d_complex_t) * NPoint * NFreq1);
    }
    for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
        for (int ipoint = 0; ipoint < NPoint; ipoint++) {
            const int id = (ipoint * NFreq1) + ifreq;
            const double mag = pow(10, spara_data[(size_t)ifreq * NPoint + ipoint].magnitude_dB / 20);
            const double ph  = spara_data[(size_t)ifreq * NPoint + ipoint].phase_deg * (PI / 180.0);
            Spara[id] = d_complex(mag * cos(ph), mag * sin(ph));
        }
    }
    if (dataset_id >= 0) H5Dclose(dataset_id);
    H5Tclose(datatype_id);

    free(spara_data);


    // input_impedance データの読み込み
    fprintf(stdout, "input_impedance_data (start)\n");
    input_impedance_data_t *input_impedance_data =
        (input_impedance_data_t *)malloc((size_t)NFeed * NFreq1 * sizeof(input_impedance_data_t));
    if (input_impedance_data == NULL) {
        fprintf(stderr, "Error: malloc failed (input_impedance_data)\n");
        return 1;
    }
    
    hid_t memtype_id;

    // メモリタイプを定義
    memtype_id = H5Tcreate(H5T_COMPOUND, sizeof(input_impedance_data_t));
    H5Tinsert(memtype_id, "frequency", HOFFSET(input_impedance_data_t, frequency), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "rin", HOFFSET(input_impedance_data_t, rin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "xin", HOFFSET(input_impedance_data_t, xin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "gin", HOFFSET(input_impedance_data_t, gin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "bin", HOFFSET(input_impedance_data_t, bin), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "ref", HOFFSET(input_impedance_data_t, ref), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype_id, "vswr", HOFFSET(input_impedance_data_t, vswr), H5T_NATIVE_DOUBLE);

    // 給電が無い (または frequency1 が無い) と書き手は作らない
    const int has_zin = (NFeed > 0) && (NFreq1 > 0) && (H5Lexists(metadata_group_id, ZIN_DATASET_NAME, H5P_DEFAULT) > 0);
    dataset_id = has_zin ? H5Dopen(metadata_group_id, ZIN_DATASET_NAME, H5P_DEFAULT) : -1;
    status = has_zin ? H5Dread(dataset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, input_impedance_data) : 0;
    if (status < 0) {
        fprintf(stderr, "Error reading dataset: %s\n", ZIN_DATASET_NAME);
    }
    //input_impedance_data から Zin への加工
    //Zin = (d_complex_t *)malloc(sizeof(d_complex_t) * NFeed * NFreq1);
    //Ref = (double *)malloc(sizeof(double)      * NFeed * NFreq1);
    for (int ifeed = 0; ifeed < NFeed; ifeed++) {
        for (int ifreq = 0; ifreq < NFreq1; ifreq++) {
            const int id = (ifeed * NFreq1) + ifreq;
            Zin[id].r = input_impedance_data[(size_t)ifeed * NFreq1 + ifreq].rin;
            Zin[id].i = input_impedance_data[(size_t)ifeed * NFreq1 + ifreq].xin;
            Ref[id] = input_impedance_data[(size_t)ifeed * NFreq1 + ifreq].ref;
        }
    }
    fprintf(stdout, "input_impedance_data (end)\n");
    if (dataset_id >= 0) H5Dclose(dataset_id);
    H5Tclose(memtype_id);

    free(input_impedance_data);

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    // /post : ofd.out と同じ値 (sol/outputHdf5.c)。給電・入力電力・余白込みの
    // 近傍界 DFT。無ければ (古いファイル) 後処理には使えないので 1 を返す
    int ok = (H5Lexists(file_id, "/post", H5P_DEFAULT) > 0);
    hid_t post_id = ok ? H5Gopen(file_id, "/post", H5P_DEFAULT) : -1;
    if (post_id >= 0) {
        // 近傍界 DFT {NFreq2, NN}
        float *dst[12] = {cEx_r, cEx_i, cEy_r, cEy_i, cEz_r, cEz_i,
                          cHx_r, cHx_i, cHy_r, cHy_i, cHz_r, cHz_i};
        static const char *nm[12] = {
            "cEx_r", "cEx_i", "cEy_r", "cEy_i", "cEz_r", "cEz_i",
            "cHx_r", "cHx_i", "cHy_r", "cHy_i", "cHz_r", "cHz_i"};
        if (NFreq2 > 0) {
            for (int b = 0; b < 12; b++) {
                dataset_id = H5Dopen(post_id, nm[b], H5P_DEFAULT);
                if ((dataset_id < 0) || (dst[b] == NULL)) { ok = 0; continue; }
                status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst[b]);
                if (status < 0) ok = 0;
                H5Dclose(dataset_id);
            }
        }

        // 給電 (dir は int で書いてある)
        if (NFeed > 0) {
            typedef struct {
                int dir, i, j, k;
                double dx, dy, dz, volt, delay, z0;
            } h5feed_t;
            h5feed_t *f = (h5feed_t *)malloc((size_t)NFeed * sizeof(h5feed_t));
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
            dataset_id = H5Dopen(post_id, "Feed", H5P_DEFAULT);
            if ((dataset_id >= 0) && (f != NULL) && (Feed != NULL)
             && (H5Dread(dataset_id, t, H5S_ALL, H5S_ALL, H5P_DEFAULT, f) >= 0)) {
                for (int n = 0; n < NFeed; n++) {
                    memset(&Feed[n], 0, sizeof(feed_t));
                    Feed[n].dir = (char)f[n].dir;
                    Feed[n].i = f[n].i; Feed[n].j = f[n].j; Feed[n].k = f[n].k;
                    Feed[n].dx = f[n].dx; Feed[n].dy = f[n].dy; Feed[n].dz = f[n].dz;
                    Feed[n].volt = f[n].volt; Feed[n].delay = f[n].delay; Feed[n].z0 = f[n].z0;
                }
            }
            else {
                ok = 0;
            }
            if (dataset_id >= 0) H5Dclose(dataset_id);
            H5Tclose(t);
            free(f);

            // 入力電力
            if (NFreq2 > 0) {
                double *pin[2] = {Pin[0], Pin[1]};
                const char *pn[2] = {"Pin0", "Pin1"};
                for (int m = 0; m < 2; m++) {
                    dataset_id = H5Dopen(post_id, pn[m], H5P_DEFAULT);
                    if ((dataset_id < 0) || (pin[m] == NULL)
                     || (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pin[m]) < 0)) {
                        ok = 0;
                    }
                    if (dataset_id >= 0) H5Dclose(dataset_id);
                }
            }
        }
        H5Gclose(post_id);
    }
    if (!ok) {
        fprintf(stderr, "*** %s : /post (near field / feed / input power) not found or incomplete\n", FILE_NAME);
    }

    // HDF5ファイルのクローズ
    H5Fclose(file_id);

    return ok ? 0 : 1;
}
