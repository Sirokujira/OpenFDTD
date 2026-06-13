// StlImporter.cpp
#include "StlImporter.h"

#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <cmath>
#include <cstring>
#include <limits>

using namespace ofd;

QStringList StlImporter::supportedExtensions()
{
    return { "stl" };
}

static void accumulate(ImportedMesh &mesh, const float v[9])
{
    for (int i = 0; i < 9; ++i) mesh.vertices.push_back(v[i]);
    for (int k = 0; k < 3; ++k) {
        for (int c = 0; c < 3; ++c) {
            const double val = v[3*k + c];
            if (mesh.numTriangles == 0 && k == 0) {
                mesh.bbox[c] = mesh.bbox[3 + c] = val;
            } else {
                mesh.bbox[c]     = std::min(mesh.bbox[c], val);
                mesh.bbox[3 + c] = std::max(mesh.bbox[3 + c], val);
            }
        }
    }
    // area += |AB × AC| / 2
    const double ab[3] = { v[3]-v[0], v[4]-v[1], v[5]-v[2] };
    const double ac[3] = { v[6]-v[0], v[7]-v[1], v[8]-v[2] };
    const double cx = ab[1]*ac[2] - ab[2]*ac[1];
    const double cy = ab[2]*ac[0] - ab[0]*ac[2];
    const double cz = ab[0]*ac[1] - ab[1]*ac[0];
    mesh.surfaceArea += 0.5 * std::sqrt(cx*cx + cy*cy + cz*cz);
    ++mesh.numTriangles;
}

bool StlImporter::load(const QString &path, ImportedMesh &mesh, QString *err)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        if (err) *err = f.errorString();
        return false;
    }
    const QByteArray data = f.readAll();
    if (data.size() < 84) {
        if (err) *err = "file too small for STL";
        return false;
    }

    mesh = ImportedMesh{};
    mesh.sourcePath = path;
    mesh.name = QFileInfo(path).completeBaseName();

    // binary STL: 80-byte header + uint32 count + 50 bytes per triangle.
    // "solid" prefix alone is not reliable — check the size equation.
    quint32 count = 0;
    std::memcpy(&count, data.constData() + 80, 4);
    const bool sizeMatches =
        (qint64(data.size()) == 84 + qint64(count) * 50);

    if (sizeMatches)
        return loadBinary(data, mesh, err);
    return loadAscii(data, mesh, err);
}

bool StlImporter::loadBinary(const QByteArray &data, ImportedMesh &mesh, QString *err)
{
    quint32 count = 0;
    std::memcpy(&count, data.constData() + 80, 4);
    const char *p = data.constData() + 84;
    mesh.vertices.reserve(int(count) * 9);
    for (quint32 i = 0; i < count; ++i) {
        float v[12];                       // normal + 3 vertices
        std::memcpy(v, p, 48);
        accumulate(mesh, v + 3);
        p += 50;                           // 48 + 2 attribute bytes
    }
    if (mesh.numTriangles == 0 && err) *err = "no triangles";
    return mesh.numTriangles > 0;
}

bool StlImporter::loadAscii(const QByteArray &data, ImportedMesh &mesh, QString *err)
{
    QTextStream in(data);
    float v[9];
    int idx = 0;
    while (!in.atEnd()) {
        QString word;
        in >> word;
        if (word == "vertex") {
            in >> v[idx] >> v[idx+1] >> v[idx+2];
            idx += 3;
            if (idx == 9) {
                accumulate(mesh, v);
                idx = 0;
            }
        }
    }
    if (mesh.numTriangles == 0) {
        if (err) *err = "no triangles found (not an STL file?)";
        return false;
    }
    return true;
}
