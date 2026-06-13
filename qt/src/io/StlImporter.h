// StlImporter.h — load STL meshes for the Geometry tab's 3D import.
//
// バイナリ/ASCII STL を自前でパースし、診断情報 (三角形数・バウンディング
// ボックス・面積) を返す。Yee格子へのボクセル化 (winding number 判定・
// 共形メッシュ) は libigl 連携 (-DUSE_LIBIGL=ON) で拡張する設計
// (docs/libigl-integration.md 参照)。
#pragma once
#include <QString>
#include <QVector>

namespace ofd {

struct ImportedMesh {
    QString          name;
    QString          sourcePath;
    QVector<float>   vertices;     // 9 floats per triangle (x1 y1 z1 x2 ...)
    int              numTriangles = 0;
    double           bbox[6] = {0,0,0,0,0,0};   // xmin ymin zmin xmax ymax zmax
    double           surfaceArea = 0.0;
};

class StlImporter {
public:
    static bool load(const QString &path, ImportedMesh &mesh, QString *err = nullptr);
    static QStringList supportedExtensions();   // for the file dialog filter

private:
    static bool loadBinary(const QByteArray &data, ImportedMesh &mesh, QString *err);
    static bool loadAscii(const QByteArray &data, ImportedMesh &mesh, QString *err);
};

} // namespace ofd
