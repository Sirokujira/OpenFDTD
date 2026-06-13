// OfdIO.h — read/write the OpenFDTD .ofd text format (本家完全互換).
//
// 文法 (sol/input_data.c, post/post_data.c が読む形式):
//   OpenFDTD 4 2                      ← ヘッダ (プログラム名 major minor)
//   title = <text>
//   xmesh = x0 d1 x1 d2 x2 ...        ← 座標と分割数の交互列
//   material = 1 εr σ μr σm           ← 通常媒質 (2=分散性: ε∞ a b c)
//   geometry = <mat> <shape> <g...>   ← 形状コード 1,2,11-13,31-33,41-43,51-53
//   feed = <dir> x y z volt delay z0
//   planewave = theta phi pol
//   point = <dir> x y z [prop]
//   load = <dir> x y z {R|L|C} value
//   abc = 0 | 1 <layers> <m> <r0>
//   pbc = <x> <y> <z>
//   frequency1 = f0 f1 div            ← [Hz]
//   frequency2 = f0 f1 div
//   solver = maxiter nout converg
//   timestep / pulsewidth / rfeed / plot3dgeom
//   plotiter / plotzin / plotfar1d / plotnear2d ... ← ポスト処理キー
//   end
//
// 未知のキーは Project::extraLines() に保存し、保存時にそのまま書き戻す
// (手編集ファイルが GUI ラウンドトリップで壊れない)。
#pragma once
#include <QString>
#include "../core/Project.h"

namespace ofd {

class OfdIO {
public:
    static bool load(const QString &path, Project &project, QString *err = nullptr);
    static bool save(const QString &path, const Project &project, QString *err = nullptr);

    // Serialize to a QString (for preview / tests).
    static QString serialize(const Project &project);

    // Parse from a QString (for tests).
    static bool parse(const QString &text, Project &project, QString *err = nullptr);
};

// .ofdx — JSON sidecar with the extension-domain settings the legacy kernel
// ignores (optical / acoustic / underwater / tidy3d). Backward compatible:
// a .ofd without sidecar is a plain EM project.
class OfdxIO {
public:
    static bool load(const QString &path, Project &project, QString *err = nullptr);
    static bool save(const QString &path, const Project &project, QString *err = nullptr);
};

} // namespace ofd
