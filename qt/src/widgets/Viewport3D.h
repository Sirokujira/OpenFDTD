// Viewport3D.h — central 3D view of the project (mesh region, geometry
// units, feeds, probes). QPainter-based orthographic wireframe so that no
// OpenGL context is required (works in headless / remote sessions too).
//
// Mouse: left-drag = orbit, middle-drag = pan, wheel = zoom, double = fit.
#pragma once
#include <QWidget>
#include <QPointF>
#include "../core/Domain.h"

namespace ofd {

class Project;

class Viewport3D : public QWidget {
    Q_OBJECT
public:
    explicit Viewport3D(Project *project, QWidget *parent = nullptr);

    void setDomain(Domain d) { m_domain = d; update(); }
    void setSolidMode(bool solid) { m_solid = solid; update(); }
    bool solidMode() const { return m_solid; }

public slots:
    void fitView();

protected:
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void wheelEvent(QWheelEvent *) override;
    void mouseDoubleClickEvent(QMouseEvent *) override;

private:
    QPointF projectPoint(double x, double y, double z) const;

    Project *m_project;
    Domain   m_domain = Domain::EM;
    bool     m_solid = false;

    double   m_azimuthDeg = -60;
    double   m_elevationDeg = 25;
    double   m_zoom = 1.0;
    QPointF  m_panPx;
    QPointF  m_lastPos;
    Qt::MouseButton m_dragButton = Qt::NoButton;

    // cached scene transform (set in paintEvent)
    mutable double m_cx = 0, m_cy = 0, m_cz = 0, m_scale = 1.0;
};

} // namespace ofd
