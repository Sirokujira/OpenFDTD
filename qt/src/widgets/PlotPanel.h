// PlotPanel.h — 2D result/preview plot (QPainter, no QtCharts dependency).
//
// Two data sources, both honest:
//   - source waveform preview: gaussian pulse computed from the project's
//     Δt/Tw settings (what the kernel will inject)
//   - convergence history: "<step> <Eavg> <Havg>" lines parsed live from the
//     running kernel's stdout (Runner::logLine → addConvergencePoint)
#pragma once
#include <QWidget>
#include <QVector>
#include <QPointF>
#include "../core/Domain.h"

namespace ofd {

class Project;

class PlotPanel : public QWidget {
    Q_OBJECT
public:
    explicit PlotPanel(Project *project, QWidget *parent = nullptr);

    void setDomain(Domain d) { m_domain = d; update(); }

public slots:
    void showWaveform()    { m_mode = Waveform;    update(); }
    void showConvergence() { m_mode = Convergence; update(); }
    void clearConvergence();
    void addConvergencePoint(int step, double e, double h);
    bool exportCsv(const QString &path) const;

protected:
    void paintEvent(QPaintEvent *) override;

private:
    enum Mode { Waveform, Convergence };

    Project *m_project;
    Domain   m_domain = Domain::EM;
    Mode     m_mode = Waveform;

    QVector<int>    m_steps;
    QVector<double> m_eAvg, m_hAvg;
};

} // namespace ofd
