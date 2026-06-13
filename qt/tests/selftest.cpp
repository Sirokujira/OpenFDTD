// selftest.cpp — .ofd ラウンドトリップ検証.
//
// data/sample/*.ofd を全件ロード → シリアライズ → 再パースし、
// 構造 (メッシュ/材質/形状/波源/周波数/ポスト設定) が一致することを確認する。
// 失敗が1件でもあれば非0で終了 (CI 用)。
//
//   ./ofdx_selftest [sample_dir]    (default: ../../data/sample)
#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <cmath>
#include <cstdio>

#include "core/Project.h"
#include "io/OfdIO.h"

using namespace ofd;

static int g_checks = 0;
static int g_failures = 0;
static QString g_file;

static void check(bool cond, const char *what)
{
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::fprintf(stderr, "FAIL %s: %s\n", qPrintable(g_file), what);
    }
}

static bool nearlyEq(double a, double b)
{
    const double m = std::max(std::fabs(a), std::fabs(b));
    return std::fabs(a - b) <= 1e-9 * std::max(m, 1.0);
}

static void compareProjects(const Project &a, const Project &b)
{
    check(a.general().title == b.general().title, "title");
    check(a.general().maxiter == b.general().maxiter, "solver maxiter");
    check(a.general().nout == b.general().nout, "solver nout");
    check(nearlyEq(a.general().converg, b.general().converg), "solver converg");
    check(a.general().abc == b.general().abc, "abc");
    check(a.general().pbcX == b.general().pbcX &&
          a.general().pbcY == b.general().pbcY &&
          a.general().pbcZ == b.general().pbcZ, "pbc");
    check(nearlyEq(a.general().f1min, b.general().f1min) &&
          nearlyEq(a.general().f1max, b.general().f1max) &&
          a.general().f1div == b.general().f1div, "frequency1");
    check(nearlyEq(a.general().f2min, b.general().f2min) &&
          nearlyEq(a.general().f2max, b.general().f2max) &&
          a.general().f2div == b.general().f2div, "frequency2");

    for (int ax = 0; ax < 3; ++ax) {
        const MeshAxis &ma = a.mesh(ax), &mb = b.mesh(ax);
        check(ma.nodes.size() == mb.nodes.size(), "mesh node count");
        check(ma.divs == mb.divs, "mesh divisions");
        for (int i = 0; i < qMin(ma.nodes.size(), mb.nodes.size()); ++i)
            check(nearlyEq(ma.nodes[i], mb.nodes[i]), "mesh node value");
    }

    check(a.materials().size() == b.materials().size(), "material count");
    for (int i = 0; i < qMin(a.materials().size(), b.materials().size()); ++i) {
        const Material &x = a.materials()[i], &y = b.materials()[i];
        check(x.type == y.type, "material type");
        check(nearlyEq(x.epsr, y.epsr) && nearlyEq(x.esgm, y.esgm) &&
              nearlyEq(x.amur, y.amur) && nearlyEq(x.msgm, y.msgm), "material values");
        check(nearlyEq(x.einf, y.einf) && nearlyEq(x.ae, y.ae) &&
              nearlyEq(x.be, y.be) && nearlyEq(x.ce, y.ce), "dispersive values");
    }

    check(a.geometries().size() == b.geometries().size(), "geometry count");
    for (int i = 0; i < qMin(a.geometries().size(), b.geometries().size()); ++i) {
        const Geometry &x = a.geometries()[i], &y = b.geometries()[i];
        check(x.materialId == y.materialId, "geometry material");
        check(x.shape == y.shape, "geometry shape");
        for (int k = 0; k < Geometry::paramCount(x.shape); ++k)
            check(nearlyEq(x.g[k], y.g[k]), "geometry coords");
    }

    check(a.feeds().size() == b.feeds().size(), "feed count");
    for (int i = 0; i < qMin(a.feeds().size(), b.feeds().size()); ++i) {
        const Feed &x = a.feeds()[i], &y = b.feeds()[i];
        check(x.dir == y.dir, "feed dir");
        check(nearlyEq(x.x, y.x) && nearlyEq(x.y, y.y) && nearlyEq(x.z, y.z),
              "feed position");
        check(nearlyEq(x.volt, y.volt) && nearlyEq(x.delay, y.delay) &&
              nearlyEq(x.z0, y.z0), "feed params");
    }
    check(a.planewave().enabled == b.planewave().enabled, "planewave");
    if (a.planewave().enabled && b.planewave().enabled) {
        check(nearlyEq(a.planewave().theta, b.planewave().theta) &&
              nearlyEq(a.planewave().phi, b.planewave().phi) &&
              a.planewave().pol == b.planewave().pol, "planewave params");
    }
    check(a.probes().size() == b.probes().size(), "point count");
    check(a.loads().size() == b.loads().size(), "load count");

    const PostOpts &pa = a.post(), &pb = b.post();
    check(pa.plotiter == pb.plotiter, "plotiter");
    check(pa.plotsmith == pb.plotsmith, "plotsmith");
    check(pa.zin.enabled == pb.zin.enabled, "plotzin");
    check(pa.ref.enabled == pb.ref.enabled, "plotref");
    check(pa.far1d.size() == pb.far1d.size(), "plotfar1d count");
    check(pa.far2d == pb.far2d, "plotfar2d");
    check(pa.near1d.size() == pb.near1d.size(), "plotnear1d count");
    check(pa.near2d.size() == pb.near2d.size(), "plotnear2d count");
    check(pa.far1dDb == pb.far1dDb, "far1ddb");
    check(a.extraLines() == b.extraLines(), "extra lines round-trip");
}

int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    QString dir = argc > 1 ? argv[1]
        : QFileInfo(QString::fromLocal8Bit(argv[0])).path() + "/../../data/sample";
    if (!QDir(dir).exists()) dir = "data/sample";

    const QStringList files =
        QDir(dir).entryList({ "*.ofd" }, QDir::Files, QDir::Name);
    if (files.isEmpty()) {
        std::fprintf(stderr, "no .ofd samples found under %s\n", qPrintable(dir));
        return 2;
    }

    int loaded = 0;
    for (const QString &name : files) {
        g_file = name;
        const QString path = QDir(dir).filePath(name);

        Project p1;
        QString err;
        if (!OfdIO::load(path, p1, &err)) {
            ++g_failures;
            std::fprintf(stderr, "FAIL %s: load: %s\n",
                         qPrintable(name), qPrintable(err));
            continue;
        }
        ++loaded;

        const QString text = OfdIO::serialize(p1);
        Project p2;
        if (!OfdIO::parse(text, p2, &err)) {
            ++g_failures;
            std::fprintf(stderr, "FAIL %s: reparse: %s\n",
                         qPrintable(name), qPrintable(err));
            continue;
        }
        compareProjects(p1, p2);
    }

    std::printf("%d files loaded, %d checks, %d failures\n",
                loaded, g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
