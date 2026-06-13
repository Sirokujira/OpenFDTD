// I18n.cpp
#include "I18n.h"

using namespace ofd;

I18n &I18n::instance() { static I18n inst; return inst; }

I18n::I18n() { loadTables(); }

void I18n::setLanguage(const QString &lang) { m_lang = lang; }

QString I18n::tr(const QString &key) {
    auto &I = instance();
    const QString ja = I.m_ja.value(key, key);
    const QString en = I.m_en.value(key, key);
    if (I.m_lang == "en")   return en;
    if (I.m_lang == "both") return (ja == en) ? ja : (ja + " / " + en);
    return ja;
}

// All strings declared once here. Add a key → it shows up everywhere.
void I18n::loadTables() {
    auto add = [this](const char *key, const char *ja, const char *en) {
        m_ja[key] = QString::fromUtf8(ja);
        m_en[key] = QString::fromUtf8(en);
    };

    // App
    add("app_name", "OpenFDTD-X", "OpenFDTD-X");
    add("app_subtitle", "マルチドメインFDTDシミュレーター", "Multi-Domain FDTD Simulator");
    add("untitled", "無題", "Untitled");

    // Menus & toolbar
    add("m_file", "ファイル(&F)", "&File");
    add("m_edit", "編集(&E)", "&Edit");
    add("m_view", "表示(&V)", "&View");
    add("m_tools", "ツール(&T)", "&Tools");
    add("m_run", "実行(&R)", "&Run");
    add("m_post", "ポスト処理(&P)", "&Post");
    add("m_help", "ヘルプ(&H)", "&Help");

    add("tb_new", "新規", "New");
    add("tb_open", "開く", "Open");
    add("tb_save", "保存", "Save");
    add("tb_saveas", "名前を付けて保存", "Save As");
    add("tb_calc", "計算", "Compute");
    add("tb_post", "ポスト処理", "Post");
    add("tb_plot2d", "図形表示2D", "Plot 2D");
    add("tb_plot3d", "図形表示3D", "Plot 3D");
    add("tb_stop", "停止", "Stop");
    add("tb_cloud", "クラウド送信", "Cloud");
    add("tb_cloud_optical_only", "クラウド送信 (光のみ)", "Cloud (optical only)");
    add("tb_export", "エクスポート", "Export");
    add("m_exit", "終了", "Exit");
    add("m_lang", "言語", "Language");
    add("m_lang_restart", "言語設定を保存しました。再起動後に反映されます。",
        "Language saved. It takes effect after restart.");

    // Domains
    add("d_em", "電磁波", "Electromagnetic");
    add("d_optical", "光", "Optical");
    add("d_acoustic", "室内音響", "Room Acoustic");
    add("d_underwater", "水中音響", "Underwater");

    // Left tabs
    add("t_general", "全般", "General");
    add("t_mesh", "メッシュ", "Mesh");
    add("t_material", "物性値・集中定数", "Materials && Lumped");
    add("t_geometry", "物体形状", "Geometry");
    add("t_source", "波源・観測点", "Source && Probe");
    add("t_post1", "ポスト処理(1)", "Post-Proc (1)");
    add("t_post2", "ポスト処理(2)", "Post-Proc (2)");
    add("t_optical", "光解析", "Optical");
    add("t_acoustic", "音響解析", "Acoustic");
    add("t_underwater", "水中音響", "Underwater");
    add("t_tidy3d", "tidy3d連携", "tidy3d");

    // General tab
    add("g_title", "タイトル", "Title");
    add("g_solver", "計算条件", "Solver");
    add("g_max_iter", "最大反復回数", "Max iterations");
    add("g_nout", "出力間隔", "Output interval");
    add("g_conv_tol", "収束判定値", "Tolerance");
    add("g_abc", "吸収境界条件", "Absorbing BC");
    add("g_mur1", "Mur 1次", "Mur 1st");
    add("g_pml", "PML", "PML");
    add("g_pml_layers", "PML層数", "PML layers");
    add("g_pml_order", "次数 m", "Order m");
    add("g_pml_r0", "反射係数 R0", "Reflection R0");
    add("g_periodic", "周期境界条件 (PBC)", "Periodic BC");
    add("g_freq1", "給電点・観測点周波数 [Hz]", "Feed/Probe frequency [Hz]");
    add("g_freq2", "遠方界・近傍界周波数 [Hz]", "Far/Near field frequency [Hz]");
    add("g_freq_min", "最小", "min");
    add("g_freq_max", "最大", "max");
    add("g_freq_div", "分割数", "divisions");
    add("g_advanced", "詳細設定", "Advanced");
    add("g_timestep", "タイムステップ Δt [s] (0=自動)", "Timestep Δt [s] (0=auto)");
    add("g_pulsewidth", "パルス幅 Tw [s] (0=自動)", "Pulse width Tw [s] (0=auto)");
    add("g_rfeed", "給電点抵抗 rfeed [Ω]", "Feed resistance rfeed [Ω]");
    add("g_plot3dgeom", "形状確認モード (plot3dgeom)", "Geometry check (plot3dgeom)");

    // Mesh tab
    add("me_axis_x", "Xメッシュ", "X mesh");
    add("me_axis_y", "Yメッシュ", "Y mesh");
    add("me_axis_z", "Zメッシュ", "Z mesh");
    add("me_coord", "座標 [m]", "Coord [m]");
    add("me_div", "分割数", "Divisions");
    add("me_add", "行を追加", "Add row");
    add("me_del", "行を削除", "Delete row");
    add("me_cells", "セル数", "cells");
    add("me_dmin", "最小間隔", "min spacing");
    add("me_total", "総セル数", "Total cells");

    // Material tab
    add("ma_section", "物性値 (材質番号は2から)", "Materials (IDs start at 2)");
    add("ma_builtin", "0=真空, 1=PEC (組込み)", "0=vacuum, 1=PEC (built-in)");
    add("ma_id", "番号", "ID");
    add("ma_type", "種別", "Type");
    add("ma_normal", "通常", "Normal");
    add("ma_dispersive", "分散性", "Dispersive");
    add("ma_name", "名前", "Name");
    add("ma_add", "材質を追加", "Add material");
    add("ma_del", "材質を削除", "Delete material");
    add("ma_lumped", "集中定数 (load)", "Lumped elements (load)");
    add("ma_kind", "種類", "Kind");
    add("ma_value", "値 [Ω/H/F]", "Value [Ω/H/F]");
    add("ma_dir", "方向", "Dir");
    add("ma_add_load", "素子を追加", "Add element");
    add("ma_del_load", "素子を削除", "Delete element");

    // Geometry tab
    add("ge_section", "物体形状 (ユニット番号順, 後優先)", "Geometry units (later wins)");
    add("ge_mat", "材質", "Material");
    add("ge_shape", "形状", "Shape");
    add("ge_add", "ユニットを追加", "Add unit");
    add("ge_del", "ユニットを削除", "Delete unit");
    add("ge_unit", "ユニット", "Unit");
    add("ge_params", "形状パラメータ", "Shape parameters");
    add("ge_shape_1", "直方体", "Brick");
    add("ge_shape_2", "楕円体/球", "Ellipsoid/Sphere");
    add("ge_shape_11", "円柱 (X軸)", "Cylinder (X)");
    add("ge_shape_12", "円柱 (Y軸)", "Cylinder (Y)");
    add("ge_shape_13", "円柱 (Z軸)", "Cylinder (Z)");
    add("ge_shape_31", "三角柱 (X軸)", "Triangle pillar (X)");
    add("ge_shape_32", "三角柱 (Y軸)", "Triangle pillar (Y)");
    add("ge_shape_33", "三角柱 (Z軸)", "Triangle pillar (Z)");
    add("ge_shape_41", "角錐台 (X軸)", "Pyramid (X)");
    add("ge_shape_42", "角錐台 (Y軸)", "Pyramid (Y)");
    add("ge_shape_43", "角錐台 (Z軸)", "Pyramid (Z)");
    add("ge_shape_51", "円錐台 (X軸)", "Cone (X)");
    add("ge_shape_52", "円錐台 (Y軸)", "Cone (Y)");
    add("ge_shape_53", "円錐台 (Z軸)", "Cone (Z)");
    add("ge_import", "3Dモデル取込 (STL)", "3D model import (STL)");
    add("ge_import_btn", "STLを取り込む…", "Import STL…");
    add("ge_import_hint",
        "STL→Yee格子ボクセル化は libigl 連携でビルドすると有効になります (-DUSE_LIBIGL=ON)。",
        "Build with -DUSE_LIBIGL=ON to enable STL→Yee voxelization.");

    // Source tab
    add("so_feeds", "給電点 (feed)", "Feeds");
    add("so_volt", "電圧 [V]", "Volt [V]");
    add("so_delay", "遅延 [deg]", "Delay [deg]");
    add("so_z0", "内部抵抗 Z0 [Ω]", "Z0 [Ω]");
    add("so_add_feed", "給電点を追加", "Add feed");
    add("so_del_feed", "給電点を削除", "Delete feed");
    add("so_planewave", "平面波入射 (planewave)", "Plane wave");
    add("so_pw_enable", "平面波を使用", "Enable plane wave");
    add("so_theta", "θ [deg]", "θ [deg]");
    add("so_phi", "φ [deg]", "φ [deg]");
    add("so_pol", "偏波", "Polarization");
    add("so_pol_v", "1: 垂直 (V)", "1: Vertical");
    add("so_pol_h", "2: 水平 (H)", "2: Horizontal");
    add("so_points", "観測点 (point)", "Observation points");
    add("so_prop", "伝搬方向 (#1のみ)", "Propagation (#1 only)");
    add("so_add_point", "観測点を追加", "Add point");
    add("so_del_point", "観測点を削除", "Delete point");
    add("so_exclusive", "給電点と平面波は同時に指定できません (本家仕様)",
        "Feed and plane wave are mutually exclusive");

    // Post1 tab
    add("p1_freq_section", "周波数特性プロット", "Frequency-domain plots");
    add("p1_iter", "反復回数特性 (plotiter)", "Iteration history (plotiter)");
    add("p1_feed", "給電点波形 (plotfeed)", "Feed waveform (plotfeed)");
    add("p1_point", "観測点波形 (plotpoint)", "Point waveform (plotpoint)");
    add("p1_smith", "スミスチャート (plotsmith)", "Smith chart (plotsmith)");
    add("p1_zin", "入力インピーダンス (plotzin)", "Input impedance (plotzin)");
    add("p1_yin", "入力アドミタンス (plotyin)", "Input admittance (plotyin)");
    add("p1_ref", "反射係数 (plotref)", "Reflection (plotref)");
    add("p1_spara", "Sパラメータ (plotspara)", "S-parameters (plotspara)");
    add("p1_coupling", "結合係数 (plotcoupling)", "Coupling (plotcoupling)");
    add("p1_matching", "整合損を含む (matchingloss)", "Include matching loss");
    add("p1_freqdiv", "周波数目盛分割 (freqdiv)", "Frequency divisions (freqdiv)");
    add("p1_user_scale", "スケール指定", "User scale");
    add("p1_min", "最小", "min");
    add("p1_max", "最大", "max");
    add("p1_div", "分割", "div");

    // Post2 tab
    add("p2_far0d", "遠方界周波数特性 (plotfar0d)", "Far field vs freq (plotfar0d)");
    add("p2_far1d", "遠方界指向性 (plotfar1d)", "Far-field pattern (plotfar1d)");
    add("p2_far2d", "遠方界全方向 (plotfar2d)", "Far-field 3D (plotfar2d)");
    add("p2_near1d", "近傍界1D (plotnear1d)", "Near field 1D (plotnear1d)");
    add("p2_near2d", "近傍界2D (plotnear2d)", "Near field 2D (plotnear2d)");
    add("p2_dir", "方向/面", "Dir/Plane");
    add("p2_division", "分割数", "Divisions");
    add("p2_angle", "角度 [deg]", "Angle [deg]");
    add("p2_style", "スタイル (far1dstyle)", "Style (far1dstyle)");
    add("p2_db", "dB表示", "dB scale");
    add("p2_norm", "正規化 (far1dnorm)", "Normalize (far1dnorm)");
    add("p2_component", "成分", "Components");
    add("p2_obj", "物体表示 (far2dobj)", "Object overlay (far2dobj)");
    add("p2_cmp", "成分", "Component");
    add("p2_pos", "位置", "Position");
    add("p2_add", "追加", "Add");
    add("p2_del", "削除", "Delete");
    add("p2_noinc", "入射波を除く (noinc)", "Exclude incident (noinc)");
    add("p2_contour", "等高線 (contour)", "Contour");
    add("p2_theta_div", "θ分割", "θ div");
    add("p2_phi_div", "φ分割", "φ div");

    // Optical tab
    add("opt_solver", "光ソルバー選択", "Optical solver");
    add("opt_solver_fdtd", "FDTD — 広帯域・任意形状・分散材料",
        "FDTD — broadband, arbitrary shapes, dispersive media");
    add("opt_solver_rcwa", "RCWA — 周期格子・薄膜・回折次数 (OpenRCWA)",
        "RCWA — periodic gratings, thin films (OpenRCWA)");
    add("opt_solver_bpm", "BPM — 導波路の高速伝搬解析 (OpenBPM)",
        "BPM — fast waveguide propagation (OpenBPM)");
    add("opt_solver_fmm", "FMM — 多層周期構造のS行列カスケード",
        "FMM — layered periodic structures, S-matrix cascade");
    add("opt_mode", "光解析モード", "Optical mode");
    add("opt_bpf", "バンドパスフィルタ (BPF)", "Band-pass filter");
    add("opt_wg", "導波路モード解析", "Waveguide mode");
    add("opt_ring", "リング共振器", "Ring resonator");
    add("opt_mzi", "MZI干渉計", "MZI interferometer");
    add("opt_meta", "メタサーフェス", "Metasurface");
    add("opt_phc", "フォトニック結晶", "Photonic crystal");
    add("opt_nf2ff", "近傍界→遠方界変換", "NF→FF transform");
    add("opt_spara", "Sパラメータ抽出", "S-parameter extraction");
    add("opt_wavelength", "波長範囲", "Wavelength range");
    add("opt_lambda_min", "λmin [nm]", "λmin [nm]");
    add("opt_lambda_max", "λmax [nm]", "λmax [nm]");
    add("opt_lambda_div", "分割数", "divisions");
    add("opt_rcwa_section", "RCWAパラメータ", "RCWA parameters");
    add("opt_rcwa_orders", "Fourier次数 (Nx × Ny)", "Fourier orders (Nx × Ny)");
    add("opt_rcwa_period", "格子周期 Λ [nm]", "Grating period Λ [nm]");
    add("opt_rcwa_layers", "層分割数", "Layer slices");
    add("opt_bpm_section", "BPMパラメータ", "BPM parameters");
    add("opt_bpm_algo", "アルゴリズム", "Algorithm");
    add("opt_bpm_dz", "Δz ステップ [nm]", "Δz step [nm]");
    add("opt_bpm_n0", "参照屈折率 n0", "Reference index n0");
    add("opt_bpm_input", "入射モード", "Input mode");
    add("opt_fmm_section", "FMMパラメータ", "FMM parameters");
    add("opt_fmm_harmonics", "ハーモニクス数", "Harmonics");
    add("opt_fmm_li", "Liの因数化規則", "Li's factorization rules");
    add("opt_bpf_section", "BPF目標仕様", "BPF target spec");
    add("opt_band", "通過帯域 [nm]", "Passband [nm]");
    add("opt_q", "Q値", "Q factor");
    add("opt_ring_section", "リング共振器", "Ring resonator");
    add("opt_radius", "半径 [μm]", "Radius [μm]");
    add("opt_gap", "ギャップ [nm]", "Gap [nm]");
    add("opt_kernel_hint",
        "RCWA は OpenRCWA (orcwa)、BPM は OpenBPM (obpm) のカーネルを実行します。",
        "RCWA runs the OpenRCWA kernel (orcwa); BPM runs OpenBPM (obpm).");

    // Acoustic tab
    add("ac_metrics", "室内音響指標", "Room-acoustic metrics");
    add("ac_rt60", "残響時間 RT60", "Reverb time RT60");
    add("ac_c80", "クラリティ C80", "Clarity C80");
    add("ac_d50", "明瞭度 D50", "Definition D50");
    add("ac_sti", "STI (音声明瞭度)", "STI");
    add("ac_edt", "初期減衰時間 EDT", "Early decay time EDT");
    add("ac_irf", "インパルス応答 (IRF)", "Impulse response");
    add("ac_aurora", "可聴化 (オーラリゼーション)", "Auralization");
    add("ac_sample_rate", "サンプリング周波数 [Hz]", "Sample rate [Hz]");
    add("ac_source", "音源", "Source");
    add("ac_directivity", "指向性", "Directivity");
    add("ac_omni", "無指向", "Omni");
    add("ac_cardioid", "カーディオイド", "Cardioid");
    add("ac_speaker", "スピーカー", "Speaker");
    add("ac_spl", "音圧レベル [dB SPL]", "SPL [dB]");
    add("ac_mics", "受音点 / マイクアレイ", "Mic array");
    add("ac_mic_count", "マイク数", "Mic count");
    add("ac_mapping_hint",
        "音響カーネルは .ofd の εr↔(ρ,c) マッピングで実行されます (.ofdx 参照)。",
        "Acoustic kernels map εr↔(ρ,c) on the same .ofd schema (see .ofdx).");

    // Underwater tab
    add("uw_env", "海洋環境", "Ocean environment");
    add("uw_temp", "水温 [°C]", "Temperature [°C]");
    add("uw_salinity", "塩分 [psu]", "Salinity [psu]");
    add("uw_ssp", "音速プロファイル (SSP)", "Sound speed profile");
    add("uw_depth", "深度 [m]", "Depth [m]");
    add("uw_speed", "音速 [m/s]", "c [m/s]");
    add("uw_sofar", "SOFARチャネル", "SOFAR channel");
    add("uw_bottom", "海底特性", "Seabed");
    add("uw_bottom_type", "底質", "Bottom type");
    add("uw_bottom_c", "音速 [m/s]", "c [m/s]");
    add("uw_bottom_rho", "密度 [kg/m³]", "ρ [kg/m³]");
    add("uw_sonar", "ソナー", "Sonar");
    add("uw_freq", "周波数 [kHz]", "Frequency [kHz]");
    add("uw_sl", "音源レベル SL [dB]", "Source level [dB]");
    add("uw_range", "最大距離 [km]", "Max range [km]");

    // tidy3d tab
    add("t3_section", "tidy3d クラウド計算 (光FDTD専用)", "tidy3d Cloud (optical FDTD only)");
    add("t3_hint",
        "tidy3d は Flexcompute 社の光FDTD専用クラウドです。光ドメインのプロジェクトを "
        "Python スクリプトへ変換して送信します。",
        "tidy3d is Flexcompute's photonic-FDTD cloud. The optical project is "
        "exported as a Python script for submission.");
    add("t3_apikey", "APIキー", "API key");
    add("t3_project", "プロジェクト名", "Project name");
    add("t3_resolution", "解像度", "Resolution");
    add("t3_auto_pml", "PML自動設定", "Auto PML");
    add("t3_export", "Pythonスクリプト生成…", "Generate Python script…");
    add("t3_mapping", "自動変換マッピング", "Conversion mapping");

    // Run / status
    add("run_engine", "実行エンジン", "Engine");
    add("run_cpu", "CPU+OpenMP", "CPU+OpenMP");
    add("run_cpu_mpi", "CPU+MPI", "CPU+MPI");
    add("run_gpu", "GPU (CUDA)", "GPU (CUDA)");
    add("run_gpu_mpi", "GPU+MPI", "GPU+MPI");
    add("run_threads", "スレッド数", "Threads");
    add("run_solver_only", "計算のみ (-solver相当)", "Solver only");
    add("run_post_only", "ポストのみ", "Post only");
    add("run_both", "一括 (計算+ポスト)", "Solver + post");

    add("sb_ready", "準備完了", "Ready");
    add("sb_running", "実行中…", "Running…");
    add("sb_failed", "失敗", "Failed");
    add("sb_done", "完了", "Done");

    // Right dock
    add("rd_project", "プロジェクト", "Project");
    add("rd_log", "ログ", "Log");
    add("rd_tree_mesh", "メッシュ", "Mesh");
    add("rd_tree_materials", "物性値", "Materials");
    add("rd_tree_geometry", "物体形状", "Geometry");
    add("rd_tree_sources", "波源", "Sources");
    add("rd_tree_points", "観測点", "Points");
    add("rd_tree_loads", "集中定数", "Lumped");

    // Plot panel / viewport
    add("vp_wire", "ワイヤ", "Wire");
    add("vp_solid", "ソリッド", "Solid");
    add("vp_fit", "フィット", "Fit");
    add("pp_waveform", "波源波形プレビュー", "Source waveform preview");
    add("pp_convergence", "収束特性 (実行ログ)", "Convergence (run log)");
    add("pp_export_csv", "CSV出力", "Export CSV");
    add("pp_export_h5", "HDF5出力", "Export HDF5");
    add("pp_export_s2p", "Touchstone出力 (.s2p)", "Export Touchstone (.s2p)");

    // ev viewer
    add("ev_backend", "表示バックエンド", "Viewer backend");
    add("ev_html", "HTML出力 (-html) をブラウザで開く", "Open -html output in browser");
    add("ev_process", "外部ビューワー (ev2d/ev3d) を起動", "Launch external viewer");
    add("ev_native", "ネイティブ描画 (未実装)", "Native rendering (TODO)");
    add("ev_open2d", "ev2d を開く", "Open ev2d");
    add("ev_open3d", "ev3d を開く", "Open ev3d");
    add("ev_nofile", "出力ファイルが見つかりません。先にポスト処理を実行してください。",
        "No output file found. Run post-processing first.");
}
