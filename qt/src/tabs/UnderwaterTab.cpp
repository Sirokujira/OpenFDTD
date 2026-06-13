// UnderwaterTab.cpp
#include "UnderwaterTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

UnderwaterTab::UnderwaterTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    // environment
    auto *se = new SectionBox(I18n::tr("uw_env"), body);
    m_temp = new QDoubleSpinBox(se);
    m_temp->setRange(-2, 40);
    m_salinity = new QDoubleSpinBox(se);
    m_salinity->setRange(0, 45);
    m_sofar = new QCheckBox(se);
    se->form()->addRow(I18n::tr("uw_temp"), m_temp);
    se->form()->addRow(I18n::tr("uw_salinity"), m_salinity);
    se->form()->addRow(I18n::tr("uw_sofar"), m_sofar);
    v->addWidget(se);

    // SSP table
    auto *sp = new SectionBox(I18n::tr("uw_ssp"), body);
    m_ssp = new QTableWidget(0, 2, sp);
    m_ssp->setHorizontalHeaderLabels({ I18n::tr("uw_depth"), I18n::tr("uw_speed") });
    m_ssp->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_ssp->verticalHeader()->setDefaultSectionSize(22);
    m_ssp->setMinimumHeight(140);
    sp->vbox()->addWidget(m_ssp);
    auto *row = new QHBoxLayout();
    auto *add = new QPushButton(I18n::tr("p2_add"), sp);
    auto *del = new QPushButton(I18n::tr("p2_del"), sp);
    row->addWidget(add); row->addWidget(del); row->addStretch(1);
    sp->vbox()->addLayout(row);
    v->addWidget(sp);

    // seabed
    auto *sb = new SectionBox(I18n::tr("uw_bottom"), body);
    m_bottomType = new QComboBox(sb);
    m_bottomType->addItems({ "sand", "mud", "gravel", "rock" });
    m_bottomC = new QDoubleSpinBox(sb);
    m_bottomC->setRange(1000, 6000);
    m_bottomRho = new QDoubleSpinBox(sb);
    m_bottomRho->setRange(1000, 4000);
    sb->form()->addRow(I18n::tr("uw_bottom_type"), m_bottomType);
    sb->form()->addRow(I18n::tr("uw_bottom_c"), m_bottomC);
    sb->form()->addRow(I18n::tr("uw_bottom_rho"), m_bottomRho);
    v->addWidget(sb);

    // sonar
    auto *ss = new SectionBox(I18n::tr("uw_sonar"), body);
    m_sonarFreq = new QDoubleSpinBox(ss);
    m_sonarFreq->setRange(0.01, 1000);
    m_sonarSL = new QDoubleSpinBox(ss);
    m_sonarSL->setRange(0, 300);
    m_rangeMax = new QDoubleSpinBox(ss);
    m_rangeMax->setRange(0.1, 10000);
    ss->form()->addRow(I18n::tr("uw_freq"), m_sonarFreq);
    ss->form()->addRow(I18n::tr("uw_sl"), m_sonarSL);
    ss->form()->addRow(I18n::tr("uw_range"), m_rangeMax);
    v->addWidget(ss);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    auto applyCb = [this] { apply(); };
    for (auto *s : { m_temp, m_salinity, m_bottomC, m_bottomRho,
                     m_sonarFreq, m_sonarSL, m_rangeMax })
        connect(s, &QDoubleSpinBox::valueChanged, this, applyCb);
    connect(m_sofar, &QCheckBox::toggled, this, applyCb);
    connect(m_bottomType, &QComboBox::currentIndexChanged, this, applyCb);

    connect(add, &QPushButton::clicked, this, [this] {
        auto &ssp = m_p->underwater().ssp;
        const double depth = ssp.isEmpty() ? 0 : ssp.last().depth_m + 500;
        ssp.push_back({ depth, 1500 });
        refresh();
        m_p->touch();
    });
    connect(del, &QPushButton::clicked, this, [this] {
        auto &ssp = m_p->underwater().ssp;
        const int r = m_ssp->currentRow();
        if (r >= 0 && r < ssp.size()) {
            ssp.removeAt(r);
            refresh();
            m_p->touch();
        }
    });
    connect(m_ssp, &QTableWidget::cellChanged, this, [this] {
        if (!m_updating) { applySsp(); m_p->touch(); }
    });

    connect(project, &Project::loaded, this, &UnderwaterTab::refresh);
    refresh();
}

void UnderwaterTab::apply()
{
    if (m_updating) return;
    UnderwaterOpts &u = m_p->underwater();
    u.waterTemp_C = m_temp->value();
    u.salinity_psu = m_salinity->value();
    u.sofar = m_sofar->isChecked();
    u.bottomType = m_bottomType->currentText();
    u.bottomC_mps = m_bottomC->value();
    u.bottomRho_kgm3 = m_bottomRho->value();
    u.sonarFreq_kHz = m_sonarFreq->value();
    u.sonarSL_dB = m_sonarSL->value();
    u.rangeMax_km = m_rangeMax->value();
    m_p->touch();
}

void UnderwaterTab::applySsp()
{
    auto &ssp = m_p->underwater().ssp;
    for (int r = 0; r < m_ssp->rowCount() && r < ssp.size(); ++r) {
        if (auto *it = m_ssp->item(r, 0)) ssp[r].depth_m = it->text().toDouble();
        if (auto *it = m_ssp->item(r, 1)) ssp[r].c_mps = it->text().toDouble();
    }
}

void UnderwaterTab::refresh()
{
    m_updating = true;
    const UnderwaterOpts &u = m_p->underwater();
    m_temp->setValue(u.waterTemp_C);
    m_salinity->setValue(u.salinity_psu);
    m_sofar->setChecked(u.sofar);
    m_bottomType->setCurrentText(u.bottomType);
    m_bottomC->setValue(u.bottomC_mps);
    m_bottomRho->setValue(u.bottomRho_kgm3);
    m_sonarFreq->setValue(u.sonarFreq_kHz);
    m_sonarSL->setValue(u.sonarSL_dB);
    m_rangeMax->setValue(u.rangeMax_km);

    m_ssp->setRowCount(u.ssp.size());
    for (int r = 0; r < u.ssp.size(); ++r) {
        m_ssp->setItem(r, 0, new QTableWidgetItem(
            QString::number(u.ssp[r].depth_m, 'g', 8)));
        m_ssp->setItem(r, 1, new QTableWidgetItem(
            QString::number(u.ssp[r].c_mps, 'g', 8)));
    }
    m_updating = false;
}
