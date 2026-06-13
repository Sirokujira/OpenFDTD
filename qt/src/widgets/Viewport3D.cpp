// Viewport3D.cpp
#include "Viewport3D.h"
#include "../core/Project.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>
#include <cmath>

using namespace ofd;

Viewport3D::Viewport3D(Project *project, QWidget *parent)
    : QWidget(parent), m_project(project)
{
    setObjectName("Viewport3D");
    setMinimumSize(320, 240);
    setMouseTracking(false);
    setAutoFillBackground(false);
    connect(project, &Project::changed, this, qOverload<>(&QWidget::update));
    connect(project, &Project::loaded, this, qOverload<>(&QWidget::update));
}

void Viewport3D::fitView()
{
    m_zoom = 1.0;
    m_panPx = QPointF();
    update();
}

QPointF Viewport3D::projectPoint(double x, double y, double z) const
{
    // center + rotate (azimuth around Z, then elevation around screen-X)
    const double az = m_azimuthDeg  * M_PI / 180.0;
    const double el = m_elevationDeg * M_PI / 180.0;
    const double dx = (x - m_cx) * m_scale;
    const double dy = (y - m_cy) * m_scale;
    const double dz = (z - m_cz) * m_scale;

    const double x1 =  dx * std::cos(az) + dy * std::sin(az);
    const double y1 = -dx * std::sin(az) + dy * std::cos(az);
    const double y2 =  y1 * std::cos(el) - dz * std::sin(el);
    // screen: x right, y down
    return QPointF(width()  / 2.0 + m_panPx.x() + x1,
                   height() / 2.0 + m_panPx.y() + y2);
}

void Viewport3D::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.fillRect(rect(), QColor("#1d2430"));

    // scene extents from the mesh
    double lo[3], hi[3];
    bool any = false;
    for (int a = 0; a < 3; ++a) {
        const MeshAxis &ax = m_project->mesh(a);
        lo[a] = ax.min(); hi[a] = ax.max();
        if (hi[a] > lo[a]) any = true;
    }
    if (!any) { lo[0]=lo[1]=lo[2]=-0.5; hi[0]=hi[1]=hi[2]=0.5; }

    m_cx = (lo[0] + hi[0]) / 2;
    m_cy = (lo[1] + hi[1]) / 2;
    m_cz = (lo[2] + hi[2]) / 2;
    const double ext = std::max({ hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2], 1e-12 });
    m_scale = 0.55 * std::min(width(), height()) / ext * m_zoom;

    const QColor accent(accentColor(m_domain));

    // axis triad (bottom-left)
    {
        const double L = ext * 0.18;
        const QPointF o  = projectPoint(lo[0], lo[1], lo[2]);
        const QPointF px = projectPoint(lo[0] + L, lo[1], lo[2]);
        const QPointF py = projectPoint(lo[0], lo[1] + L, lo[2]);
        const QPointF pz = projectPoint(lo[0], lo[1], lo[2] + L);
        p.setPen(QPen(QColor("#e05555"), 2)); p.drawLine(o, px); p.drawText(px, "X");
        p.setPen(QPen(QColor("#4fb24f"), 2)); p.drawLine(o, py); p.drawText(py, "Y");
        p.setPen(QPen(QColor("#5b8fd9"), 2)); p.drawLine(o, pz); p.drawText(pz, "Z");
    }

    // mesh region box
    auto drawBox = [&](const double a[3], const double b[3], const QPen &pen) {
        const QPointF v[8] = {
            projectPoint(a[0],a[1],a[2]), projectPoint(b[0],a[1],a[2]),
            projectPoint(b[0],b[1],a[2]), projectPoint(a[0],b[1],a[2]),
            projectPoint(a[0],a[1],b[2]), projectPoint(b[0],a[1],b[2]),
            projectPoint(b[0],b[1],b[2]), projectPoint(a[0],b[1],b[2]),
        };
        static const int e[12][2] = { {0,1},{1,2},{2,3},{3,0},
                                      {4,5},{5,6},{6,7},{7,4},
                                      {0,4},{1,5},{2,6},{3,7} };
        p.setPen(pen);
        for (auto &ed : e) p.drawLine(v[ed[0]], v[ed[1]]);
    };
    drawBox(lo, hi, QPen(QColor(255,255,255,70), 1, Qt::DashLine));

    // mesh grid ticks on the bottom face (z = lo[2])
    {
        p.setPen(QPen(QColor(255,255,255,28), 1));
        const MeshAxis &mx = m_project->mesh(0);
        const MeshAxis &my = m_project->mesh(1);
        for (int i = 0; i < mx.divs.size(); ++i) {
            const double x0 = mx.nodes[i], x1 = mx.nodes[i+1];
            for (int k = 0; k <= mx.divs[i]; ++k) {
                const double x = x0 + (x1 - x0) * k / mx.divs[i];
                p.drawLine(projectPoint(x, lo[1], lo[2]),
                           projectPoint(x, hi[1], lo[2]));
            }
        }
        for (int i = 0; i < my.divs.size(); ++i) {
            const double y0 = my.nodes[i], y1 = my.nodes[i+1];
            for (int k = 0; k <= my.divs[i]; ++k) {
                const double y = y0 + (y1 - y0) * k / my.divs[i];
                p.drawLine(projectPoint(lo[0], y, lo[2]),
                           projectPoint(hi[0], y, lo[2]));
            }
        }
    }

    // geometry units
    int unit = 0;
    for (const Geometry &g : m_project->geometries()) {
        ++unit;
        QColor col = accent;
        col.setAlpha(m_solid ? 110 : 230);
        const QPen pen(col.lighter(120), 1.4);

        // all 6-parameter shapes are drawn from their bounding box; the
        // ellipsoid/cylinder shapes additionally show an inscribed outline
        double a[3] = { g.g[0], g.g[2], g.g[4] };
        double b[3] = { g.g[1], g.g[3], g.g[5] };
        if (Geometry::paramCount(g.shape) == 8) {
            // 8-param shapes: use min/max of the coordinate list as a hull
            a[0] = std::min({g.g[0], g.g[1]}); b[0] = std::max({g.g[0], g.g[1]});
            a[1] = std::min({g.g[2], g.g[3]}); b[1] = std::max({g.g[2], g.g[3]});
            a[2] = std::min({g.g[4], g.g[5], g.g[6], g.g[7]});
            b[2] = std::max({g.g[4], g.g[5], g.g[6], g.g[7]});
        }

        if (m_solid) {
            // shade the top face
            QPainterPath path;
            path.moveTo(projectPoint(a[0], a[1], b[2]));
            path.lineTo(projectPoint(b[0], a[1], b[2]));
            path.lineTo(projectPoint(b[0], b[1], b[2]));
            path.lineTo(projectPoint(a[0], b[1], b[2]));
            path.closeSubpath();
            p.fillPath(path, col);
        }
        drawBox(a, b, pen);

        if (g.shape == 2 || (g.shape >= 11 && g.shape <= 13)) {
            // inscribed ellipse outline on the mid plane
            p.setPen(pen);
            const int N = 36;
            QPolygonF poly;
            for (int k = 0; k <= N; ++k) {
                const double t = 2 * M_PI * k / N;
                double x = (a[0]+b[0])/2, y = (a[1]+b[1])/2, z = (a[2]+b[2])/2;
                const double rx = (b[0]-a[0])/2, ry = (b[1]-a[1])/2,
                             rz = (b[2]-a[2])/2;
                switch (g.shape) {
                    case 11: y += ry*std::cos(t); z += rz*std::sin(t); break;
                    case 12: x += rx*std::cos(t); z += rz*std::sin(t); break;
                    default: x += rx*std::cos(t); y += ry*std::sin(t); break;
                }
                poly << projectPoint(x, y, z);
            }
            p.drawPolyline(poly);
        }

        p.setPen(QColor(255,255,255,140));
        p.drawText(projectPoint(b[0], b[1], b[2]) + QPointF(3,-3),
                   g.name.isEmpty() ? QStringLiteral("#%1").arg(unit) : g.name);
    }

    // feeds (red diamonds) and probes (green circles)
    p.setPen(Qt::NoPen);
    for (const Feed &f : m_project->feeds()) {
        const QPointF c = projectPoint(f.x, f.y, f.z);
        QPolygonF d; d << c+QPointF(0,-5) << c+QPointF(5,0)
                       << c+QPointF(0,5)  << c+QPointF(-5,0);
        p.setBrush(QColor("#ff5252"));
        p.drawPolygon(d);
    }
    for (const Probe &pr : m_project->probes()) {
        const QPointF c = projectPoint(pr.x, pr.y, pr.z);
        p.setBrush(QColor("#69d069"));
        p.drawEllipse(c, 4, 4);
    }
    if (m_project->planewave().enabled) {
        // incident direction arrow from outside the box
        const double th = m_project->planewave().theta * M_PI / 180.0;
        const double ph = m_project->planewave().phi   * M_PI / 180.0;
        const double R = ext * 0.75;
        const QPointF from = projectPoint(m_cx + R*std::sin(th)*std::cos(ph),
                                          m_cy + R*std::sin(th)*std::sin(ph),
                                          m_cz + R*std::cos(th));
        const QPointF to = projectPoint(m_cx, m_cy, m_cz);
        p.setPen(QPen(QColor("#ffd24d"), 2));
        p.drawLine(from, to);
        p.setBrush(QColor("#ffd24d"));
        p.drawEllipse(to, 3, 3);
    }

    // overlay text
    p.setPen(QColor(255,255,255,150));
    p.drawText(8, height() - 10,
               QStringLiteral("%1   cells: %L2   az %3°  el %4°")
               .arg(domainKey(m_domain))
               .arg(m_project->totalCells())
               .arg(int(m_azimuthDeg)).arg(int(m_elevationDeg)));
}

void Viewport3D::mousePressEvent(QMouseEvent *e)
{
    m_lastPos = e->position();
    m_dragButton = e->button();
}

void Viewport3D::mouseMoveEvent(QMouseEvent *e)
{
    const QPointF d = e->position() - m_lastPos;
    m_lastPos = e->position();
    if (m_dragButton == Qt::LeftButton) {
        m_azimuthDeg  += d.x() * 0.5;
        m_elevationDeg = qBound(-89.0, m_elevationDeg + d.y() * 0.5, 89.0);
        update();
    } else if (m_dragButton == Qt::MiddleButton) {
        m_panPx += d;
        update();
    }
}

void Viewport3D::wheelEvent(QWheelEvent *e)
{
    const double f = std::pow(1.0015, e->angleDelta().y());
    m_zoom = qBound(0.05, m_zoom * f, 50.0);
    update();
}

void Viewport3D::mouseDoubleClickEvent(QMouseEvent *)
{
    fitView();
}
