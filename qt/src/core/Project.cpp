// Project.cpp
#include "Project.h"
#include "../io/OfdIO.h"

#include <QFileInfo>
#include <cmath>

using namespace ofd;

Project::Project(QObject *parent) : QObject(parent)
{
    clear();
}

void Project::clear()
{
    m_general = GeneralOpts{};
    for (int a = 0; a < 3; ++a) {
        m_mesh[a].nodes = { -0.05, 0.05 };
        m_mesh[a].divs  = { 20 };
    }
    m_materials.clear();
    m_loads.clear();
    m_geometries.clear();
    m_feeds.clear();
    m_planewave = PlaneWave{};
    m_probes.clear();
    m_post = PostOpts{};
    m_optical = OpticalOpts{};
    m_acoustic = AcousticOpts{};
    m_underwater = UnderwaterOpts{};
    m_underwater.ssp = { {0, 1525}, {100, 1510}, {500, 1490}, {1000, 1485},
                         {1500, 1488}, {3000, 1510}, {5000, 1540} };
    m_tidy3d = Tidy3dOpts{};
    m_extraLines.clear();
    m_filePath.clear();
}

void Project::setActiveDomain(Domain d)
{
    if (m_domain == d) return;
    m_domain = d;
    emit domainChanged(d);
    emit changed();
}

qint64 Project::totalCells() const
{
    qint64 n = 1;
    for (int a = 0; a < 3; ++a) n *= qMax(1, m_mesh[a].totalCells());
    return n;
}

double Project::estimatedMemoryMB() const
{
    // E/H 6成分 (double) + 媒質ID等の補助配列 ≈ 60 byte/cell
    return totalCells() * 60.0 / (1024.0 * 1024.0);
}

double Project::courantDt() const
{
    const double c0 = 2.99792458e8;
    double s = 0;
    for (int a = 0; a < 3; ++a) {
        const double d = m_mesh[a].minSpacing();
        if (d <= 0 || d >= 1e308) return 0;
        s += 1.0 / (d * d);
    }
    return (s > 0) ? 1.0 / (c0 * std::sqrt(s)) : 0;
}

bool Project::load(const QString &path, QString *err)
{
    clear();
    if (!OfdIO::load(path, *this, err)) return false;

    const QString ofdx = QFileInfo(path).path() + "/" +
                         QFileInfo(path).completeBaseName() + ".ofdx";
    if (QFileInfo::exists(ofdx))
        OfdxIO::load(ofdx, *this, nullptr);   // sidecar is optional

    m_filePath = path;
    emit loaded();
    emit domainChanged(m_domain);
    emit changed();
    return true;
}

bool Project::save(const QString &path, QString *err)
{
    if (!OfdIO::save(path, *this, err)) return false;

    const QString ofdx = QFileInfo(path).path() + "/" +
                         QFileInfo(path).completeBaseName() + ".ofdx";
    if (!OfdxIO::save(ofdx, *this, err)) return false;

    m_filePath = path;
    return true;
}
