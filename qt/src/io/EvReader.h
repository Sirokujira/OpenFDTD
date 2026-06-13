// EvReader.h — minimal parser for .ev2 / .ev3 (native rendering backend).
//
// The format is undocumented — `ev2d_*` / `ev3d_*` are part of an "ev"
// drawing library shipped with OpenFDTD. The library writes either:
//   - HTML (type=0, text)
//   - binary .ev2/.ev3 (type=1, binary=1)
//   - text .ev2/.ev3   (type=1, binary=0)
//
// What we know from the public API (chap5/chap6 reference):
//   ev3d_init() / ev3d_newPage() / drawing… / ev3d_output()
// → page-oriented display list, OpenGL-style coordinate system.
//
// This header sketches the data structures we'll fill once we reverse the
// exact binary layout. The viewer code is structured so that adding new
// drawing primitives is a matter of adding a case to the switch in render().
//
// To complete this: read `Sirokujira/OpenFDTD` GitHub repo, files:
//   include/ev2d_*.h, include/ev3d_*.h, post/ev2d_*.c, post/ev3d_*.c
#pragma once
#include <QString>
#include <QVector>
#include <QColor>
#include <QVariant>

class QPainter;

namespace ofd {

enum class EvCmd {
    // Page lifecycle
    BeginPage, EndPage,
    // State
    SetColor, SetWidth, SetFont,
    // 2D primitives
    Line2D, Polyline2D, Polygon2D, Text2D, Rect2D, Arc2D,
    // 3D primitives
    Line3D, Triangle3D, Polygon3D, Text3D,
    // Color-mapped data
    Contour2D, Vector2D, Surface3D,
    // Camera (3D)
    SetView, SetProjection
};

struct EvCommand {
    EvCmd            kind;
    QColor           color = Qt::black;
    double           width = 1.0;
    QVector<double>  v;       // coords (length depends on kind)
    QString          text;
    QVariantMap      props;   // extra (font size, label, etc.)
};

struct EvPage {
    QString          title;
    QVector<EvCommand> commands;
};

struct EvDocument {
    int              version = 0;
    bool             is3D = false;
    QVector<EvPage>  pages;
};

class EvReader {
public:
    // Auto-detects binary vs text from the first byte (binary starts with
    // a magic number, text starts with an ASCII keyword).
    static bool load(const QString &path, EvDocument &doc, QString *err = nullptr);

private:
    static bool loadBinary(const QString &path, EvDocument &doc, QString *err);
    static bool loadText  (const QString &path, EvDocument &doc, QString *err);
};

// EvRenderer — replay an EvDocument into a QPainter (2D) or OpenGL (3D).
class EvRenderer {
public:
    // 2D: draw page n into a QPainter at the given pixel rect.
    static void render2D(QPainter &p, const QRectF &rect, const EvPage &page);

    // 3D: replay as immediate-mode OpenGL commands. Called from a
    // QOpenGLWidget::paintGL(). Camera params are read from SetView cmds.
    static void render3D(const EvPage &page);
};

} // namespace ofd
