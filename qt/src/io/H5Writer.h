// H5Writer.h — HDF5 export of the project description + convergence data.
//
// このフォークのソルバーは ofd.out と並行して HDF5 時系列を直接出力できる
// (sol/solve.c)。GUI 側の H5Writer は補助的に「プロジェクト定義+収束履歴」を
// .h5 へ書き出す。-DUSE_HDF5=ON のときのみ有効。
#pragma once
#include <QString>
#include <QVector>

namespace ofd {

class Project;

class H5Writer {
public:
    static bool available();
    static bool write(const QString &path, const Project &project,
                      const QVector<int> &steps,
                      const QVector<double> &eAvg,
                      const QVector<double> &hAvg,
                      QString *err = nullptr);
};

} // namespace ofd
