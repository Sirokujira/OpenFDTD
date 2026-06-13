// Tidy3dExporter.h — generate a tidy3d Python script from the project.
//
// tidy3d は Flexcompute 社の光FDTD専用クラウド。光ドメインのプロジェクトを
// td.Simulation(...) の Python スクリプトへ変換し、ユーザーが
// `python script.py` でジョブ送信できる形にする (APIキーは tidy3d CLI の
// 標準設定 ~/.tidy3d を使用し、スクリプトには埋め込まない)。
//
// 変換マッピング:
//   geometry (直方体)     → td.Box
//   geometry (楕円体)     → td.Sphere (等方近似)
//   geometry (円柱)       → td.Cylinder
//   material (type 1)     → td.Medium(permittivity, conductivity)
//   material (type 2)     → コメントで明示 (PoleResidue 要手動調整)
//   feed                  → td.PointDipole
//   planewave             → td.PlaneWave
//   frequency2 / 波長範囲 → td.GaussianPulse + freqs
//   音響系パラメータ      → 非対応 (光のみ)
#pragma once
#include <QString>

namespace ofd {

class Project;

class Tidy3dExporter {
public:
    static QString generatePython(const Project &project);
    static bool exportTo(const QString &path, const Project &project,
                         QString *err = nullptr);
};

} // namespace ofd
