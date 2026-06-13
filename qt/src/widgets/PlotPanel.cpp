// PlotPanel.cpp
#include "PlotPanel.h"
#include "../core/Project.h"
#include "../I18n.h"

#include <QFile>
#include <QPainter>
#include <QPainterPath>
#include <QTextStream>
#include <cmath>

using namespace ofd;

PlotPanel::PlotPanel(Project *project, QWidget *parent)
    : QWidget(parent), m_project(project)
{
    setObjectName("PlotPanel");
    setMinimumSize(320, 200);
    connect(project, &Project::changed, this, qOverload<>(&QWidget::update));
}

void PlotPanel::clearConvergence()
{
    m_steps.clear(); m_eAvg.clear(); m_hAvg.clear();
    update();
}

void PlotPanel::addConvergencePoint(int step, double e, double h)
{
    m_steps.push_back(step);
    m_eAvg.push_back(e);
    m_hAvg.push_back(h);
    if (m_mode == Convergence) update();
}

bool PlotPanel::exportCsv(const QString &path) const
{
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    QTextStream out(&f);
    out << "step,Eavg,Havg\n";
    for (int i = 0; i < m_steps.size(); ++i)
        out << m_steps[i] << ',' << m_eAvg[i] << ',' << m_hAvg[i] << '\n';
    return true;
}

void PlotPanel::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    p.fillRect(rect(), palette().base());

    const QRectF plot(56, 28, width() - 76, height() - 64);
    const QColor accent(accentColor(m_domain));

    p.setPen(QPen(palette().mid().color(), 1));
    p.drawRect(plot);

    // grid
    p.setPen(QPen(palette().midlight().color(), 1, Qt::DotLine));
    for (int i = 1; i < 10; ++i) {
        const double x = plot.left() + plot.width() * i / 10.0;
        p.drawLine(QPointF(x, plot.top()), QPointF(x, plot.bottom()));
    }
    for (int i = 1; i < 5; ++i) {
        const double y = plot.top() + plot.height() * i / 5.0;
        p.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y));
    }

    p.setPen(palette().text().color());

    if (m_mode == Waveform) {
        p.drawText(QPointF(plot.left(), 18), I18n::tr("pp_waveform"));

        // gaussian pulse exactly as the kernel computes it:
        //   v(t) = exp(-((t - 4σ)/σ)²·π) 系の正規化パルス。
        //   Tw (pulsewidth) 未指定時はカーネル既定 (周波数帯域から自動)。
        const GeneralOpts &g = m_project->general();
        double dt = g.dt > 0 ? g.dt : m_project->courantDt();
        if (dt <= 0) dt = 1e-12;
        const double tw = g.tw > 0 ? g.tw
                          : (g.f1max > 0 ? 1.27 / g.f1max : 100 * dt);
        const int N = qMax(64, qMin(2048, int(4 * tw / dt)));

        QPainterPath path;
        for (int i = 0; i <= N; ++i) {
            const double t = 4.0 * tw * i / N;
            const double arg = (t - 2.0 * tw) / (tw / 2.0);
            const double v = std::exp(-arg * arg);
            const double x = plot.left() + plot.width() * i / double(N);
            const double y = plot.bottom() - plot.height() * v * 0.92;
            if (i == 0) path.moveTo(x, y); else path.lineTo(x, y);
        }
        p.setPen(QPen(accent, 2));
        p.drawPath(path);

        p.setPen(palette().text().color());
        p.drawText(QPointF(plot.left(), plot.bottom() + 16),
                   QStringLiteral("t: 0 … %1 s   (Δt=%2 s, Tw=%3 s)")
                       .arg(QString::number(4 * tw, 'g', 3),
                            QString::number(dt, 'g', 3),
                            QString::number(tw, 'g', 3)));
    } else {
        p.drawText(QPointF(plot.left(), 18), I18n::tr("pp_convergence"));
        if (m_steps.isEmpty()) {
            p.drawText(plot, Qt::AlignCenter,
                       QStringLiteral("no data — run the solver"));
            return;
        }

        double vmax = 1e-300;
        for (double v : m_eAvg) vmax = std::max(vmax, v);
        for (double v : m_hAvg) vmax = std::max(vmax, v);
        const int smax = qMax(1, m_steps.last());

        auto drawSeries = [&](const QVector<double> &v, const QColor &c) {
            QPainterPath path;
            for (int i = 0; i < v.size(); ++i) {
                const double x = plot.left() + plot.width() * m_steps[i] / double(smax);
                const double y = plot.bottom() - plot.height() * (v[i] / vmax) * 0.92;
                if (i == 0) path.moveTo(x, y); else path.lineTo(x, y);
            }
            p.setPen(QPen(c, 2));
            p.drawPath(path);
        };
        drawSeries(m_eAvg, accent);
        drawSeries(m_hAvg, QColor("#888888"));

        p.setPen(accent);
        p.drawText(QPointF(plot.right() - 110, plot.top() + 16), "⟨E⟩");
        p.setPen(QColor("#888888"));
        p.drawText(QPointF(plot.right() - 70, plot.top() + 16), "⟨H⟩");
        p.setPen(palette().text().color());
        p.drawText(QPointF(plot.left(), plot.bottom() + 16),
                   QStringLiteral("step: 0 … %1").arg(smax));
    }
}
