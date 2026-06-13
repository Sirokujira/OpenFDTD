// SourceTab.cpp
#include "SourceTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QCheckBox>
#include <QComboBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

SourceTab::SourceTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    m_warning = new QLabel(I18n::tr("so_exclusive"), body);
    m_warning->setStyleSheet("color: #c25400; font-weight: 600;");
    m_warning->setVisible(false);
    v->addWidget(m_warning);

    // feeds
    auto *sf = new SectionBox(I18n::tr("so_feeds"), body);
    m_feeds = new QTableWidget(0, 7, sf);
    m_feeds->setHorizontalHeaderLabels({
        I18n::tr("ma_dir"), "X [m]", "Y [m]", "Z [m]",
        I18n::tr("so_volt"), I18n::tr("so_delay"), I18n::tr("so_z0") });
    m_feeds->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_feeds->verticalHeader()->setDefaultSectionSize(24);
    m_feeds->setMinimumHeight(110);
    sf->vbox()->addWidget(m_feeds);

    auto *frow = new QHBoxLayout();
    auto *fadd = new QPushButton(I18n::tr("so_add_feed"), sf);
    auto *fdel = new QPushButton(I18n::tr("so_del_feed"), sf);
    frow->addWidget(fadd);
    frow->addWidget(fdel);
    frow->addStretch(1);
    sf->vbox()->addLayout(frow);
    v->addWidget(sf);

    // plane wave
    auto *sp = new SectionBox(I18n::tr("so_planewave"), body);
    m_pwEnable = new QCheckBox(I18n::tr("so_pw_enable"), sp);
    m_pwTheta = new QLineEdit(sp); m_pwTheta->setMaximumWidth(90);
    m_pwPhi   = new QLineEdit(sp); m_pwPhi->setMaximumWidth(90);
    m_pwPol = new QComboBox(sp);
    m_pwPol->addItem(I18n::tr("so_pol_v"));
    m_pwPol->addItem(I18n::tr("so_pol_h"));
    sp->form()->addRow(m_pwEnable);
    sp->form()->addRow(I18n::tr("so_theta"), m_pwTheta);
    sp->form()->addRow(I18n::tr("so_phi"), m_pwPhi);
    sp->form()->addRow(I18n::tr("so_pol"), m_pwPol);
    v->addWidget(sp);

    // observation points
    auto *so = new SectionBox(I18n::tr("so_points"), body);
    m_points = new QTableWidget(0, 5, so);
    m_points->setHorizontalHeaderLabels({
        I18n::tr("ma_dir"), "X [m]", "Y [m]", "Z [m]", I18n::tr("so_prop") });
    m_points->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_points->verticalHeader()->setDefaultSectionSize(24);
    m_points->setMinimumHeight(110);
    so->vbox()->addWidget(m_points);

    auto *prow = new QHBoxLayout();
    auto *padd = new QPushButton(I18n::tr("so_add_point"), so);
    auto *pdel = new QPushButton(I18n::tr("so_del_point"), so);
    prow->addWidget(padd);
    prow->addWidget(pdel);
    prow->addStretch(1);
    so->vbox()->addLayout(prow);
    v->addWidget(so);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    connect(fadd, &QPushButton::clicked, this, [this] {
        m_p->feeds().push_back(Feed{});
        refresh();
        m_p->touch();
    });
    connect(fdel, &QPushButton::clicked, this, [this] {
        const int r = m_feeds->currentRow();
        if (r >= 0 && r < m_p->feeds().size()) {
            m_p->feeds().removeAt(r);
            refresh();
            m_p->touch();
        }
    });
    connect(m_feeds, &QTableWidget::cellChanged, this, [this] {
        if (m_updating) return;
        applyFeeds();
        m_p->touch();
    });

    auto applyPw = [this] {
        if (m_updating) return;
        PlaneWave &pw = m_p->planewave();
        pw.enabled = m_pwEnable->isChecked();
        pw.theta = m_pwTheta->text().toDouble();
        pw.phi   = m_pwPhi->text().toDouble();
        pw.pol   = m_pwPol->currentIndex() + 1;
        updateExclusiveWarning();
        m_p->touch();
    };
    connect(m_pwEnable, &QCheckBox::toggled, this, applyPw);
    connect(m_pwTheta, &QLineEdit::editingFinished, this, applyPw);
    connect(m_pwPhi, &QLineEdit::editingFinished, this, applyPw);
    connect(m_pwPol, &QComboBox::currentIndexChanged, this, applyPw);

    connect(padd, &QPushButton::clicked, this, [this] {
        Probe pr;
        if (m_p->probes().isEmpty()) pr.propagation = "+X";
        m_p->probes().push_back(pr);
        refresh();
        m_p->touch();
    });
    connect(pdel, &QPushButton::clicked, this, [this] {
        const int r = m_points->currentRow();
        if (r >= 0 && r < m_p->probes().size()) {
            m_p->probes().removeAt(r);
            refresh();
            m_p->touch();
        }
    });
    connect(m_points, &QTableWidget::cellChanged, this, [this] {
        if (m_updating) return;
        applyPoints();
        m_p->touch();
    });

    connect(project, &Project::loaded, this, &SourceTab::refresh);
    refresh();
}

void SourceTab::applyFeeds()
{
    auto &feeds = m_p->feeds();
    for (int r = 0; r < m_feeds->rowCount() && r < feeds.size(); ++r) {
        Feed &f = feeds[r];
        auto cell = [this, r](int c) {
            auto *it = m_feeds->item(r, c);
            return it ? it->text() : QString();
        };
        if (auto *cb = qobject_cast<QComboBox *>(m_feeds->cellWidget(r, 0)))
            f.dir = "XYZ"[cb->currentIndex()];
        f.x = cell(1).toDouble();
        f.y = cell(2).toDouble();
        f.z = cell(3).toDouble();
        f.volt  = cell(4).toDouble();
        f.delay = cell(5).toDouble();
        f.z0    = cell(6).toDouble();
    }
    updateExclusiveWarning();
}

void SourceTab::applyPoints()
{
    auto &probes = m_p->probes();
    for (int r = 0; r < m_points->rowCount() && r < probes.size(); ++r) {
        Probe &pr = probes[r];
        auto cell = [this, r](int c) {
            auto *it = m_points->item(r, c);
            return it ? it->text() : QString();
        };
        if (auto *cb = qobject_cast<QComboBox *>(m_points->cellWidget(r, 0)))
            pr.dir = "XYZ"[cb->currentIndex()];
        pr.x = cell(1).toDouble();
        pr.y = cell(2).toDouble();
        pr.z = cell(3).toDouble();
        if (r == 0) pr.propagation = cell(4);
    }
}

void SourceTab::updateExclusiveWarning()
{
    m_warning->setVisible(!m_p->feeds().isEmpty()
                          && m_p->planewave().enabled);
}

void SourceTab::refresh()
{
    m_updating = true;

    const auto &feeds = m_p->feeds();
    m_feeds->setRowCount(feeds.size());
    for (int r = 0; r < feeds.size(); ++r) {
        const Feed &f = feeds[r];
        auto *dir = new QComboBox(m_feeds);
        dir->addItems({ "X", "Y", "Z" });
        dir->setCurrentIndex(f.dir == 'X' ? 0 : f.dir == 'Y' ? 1 : 2);
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyFeeds();
            m_p->touch();
        });
        m_feeds->setCellWidget(r, 0, dir);
        const double vals[6] = { f.x, f.y, f.z, f.volt, f.delay, f.z0 };
        for (int c = 0; c < 6; ++c)
            m_feeds->setItem(r, 1 + c, new QTableWidgetItem(
                QString::number(vals[c], 'g', 8)));
    }

    const PlaneWave &pw = m_p->planewave();
    m_pwEnable->setChecked(pw.enabled);
    m_pwTheta->setText(QString::number(pw.theta, 'g', 8));
    m_pwPhi->setText(QString::number(pw.phi, 'g', 8));
    m_pwPol->setCurrentIndex(pw.pol == 2 ? 1 : 0);

    const auto &probes = m_p->probes();
    m_points->setRowCount(probes.size());
    for (int r = 0; r < probes.size(); ++r) {
        const Probe &pr = probes[r];
        auto *dir = new QComboBox(m_points);
        dir->addItems({ "X", "Y", "Z" });
        dir->setCurrentIndex(pr.dir == 'X' ? 0 : pr.dir == 'Y' ? 1 : 2);
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyPoints();
            m_p->touch();
        });
        m_points->setCellWidget(r, 0, dir);
        const double vals[3] = { pr.x, pr.y, pr.z };
        for (int c = 0; c < 3; ++c)
            m_points->setItem(r, 1 + c, new QTableWidgetItem(
                QString::number(vals[c], 'g', 8)));
        auto *prop = new QTableWidgetItem(r == 0 ? pr.propagation : QString());
        if (r != 0) prop->setFlags(prop->flags() & ~Qt::ItemIsEditable);
        m_points->setItem(r, 4, prop);
    }

    updateExclusiveWarning();
    m_updating = false;
}
