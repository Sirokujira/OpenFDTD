// MaterialTab.cpp
#include "MaterialTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QComboBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

MaterialTab::MaterialTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    // materials
    auto *sm = new SectionBox(I18n::tr("ma_section"), body);
    sm->vbox()->addWidget(new QLabel(I18n::tr("ma_builtin"), sm));

    m_mats = new QTableWidget(0, 7, sm);
    m_mats->setHorizontalHeaderLabels({
        I18n::tr("ma_type"),
        "εr / ε∞", "σ / a", "μr / b", "σm / c",
        I18n::tr("ma_name"), I18n::tr("ma_id") });
    m_mats->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_mats->horizontalHeader()->setStretchLastSection(true);
    m_mats->verticalHeader()->setDefaultSectionSize(24);
    m_mats->setMinimumHeight(160);
    sm->vbox()->addWidget(m_mats);

    auto *mrow = new QHBoxLayout();
    auto *madd = new QPushButton(I18n::tr("ma_add"), sm);
    auto *mdel = new QPushButton(I18n::tr("ma_del"), sm);
    mrow->addWidget(madd);
    mrow->addWidget(mdel);
    mrow->addStretch(1);
    sm->vbox()->addLayout(mrow);
    v->addWidget(sm);

    // lumped elements
    auto *sl = new SectionBox(I18n::tr("ma_lumped"), body);
    m_loads = new QTableWidget(0, 6, sl);
    m_loads->setHorizontalHeaderLabels({
        I18n::tr("ma_dir"), "X [m]", "Y [m]", "Z [m]",
        I18n::tr("ma_kind"), I18n::tr("ma_value") });
    m_loads->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_loads->verticalHeader()->setDefaultSectionSize(24);
    m_loads->setMinimumHeight(110);
    sl->vbox()->addWidget(m_loads);

    auto *lrow = new QHBoxLayout();
    auto *ladd = new QPushButton(I18n::tr("ma_add_load"), sl);
    auto *ldel = new QPushButton(I18n::tr("ma_del_load"), sl);
    lrow->addWidget(ladd);
    lrow->addWidget(ldel);
    lrow->addStretch(1);
    sl->vbox()->addLayout(lrow);
    v->addWidget(sl);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    connect(madd, &QPushButton::clicked, this, [this] {
        Material m;
        m.epsr = 2.0;
        m_p->materials().push_back(m);
        refresh();
        m_p->touch();
    });
    connect(mdel, &QPushButton::clicked, this, [this] {
        const int row = m_mats->currentRow();
        auto &mats = m_p->materials();
        if (row >= 0 && row < mats.size()) {
            mats.removeAt(row);
            refresh();
            m_p->touch();
        }
    });
    connect(m_mats, &QTableWidget::cellChanged, this, [this] {
        if (m_updating) return;
        applyMaterials();
        m_p->touch();
    });

    connect(ladd, &QPushButton::clicked, this, [this] {
        m_p->loads().push_back(Load{});
        refresh();
        m_p->touch();
    });
    connect(ldel, &QPushButton::clicked, this, [this] {
        const int row = m_loads->currentRow();
        auto &loads = m_p->loads();
        if (row >= 0 && row < loads.size()) {
            loads.removeAt(row);
            refresh();
            m_p->touch();
        }
    });
    connect(m_loads, &QTableWidget::cellChanged, this, [this] {
        if (m_updating) return;
        applyLoads();
        m_p->touch();
    });

    connect(project, &Project::loaded, this, &MaterialTab::refresh);
    refresh();
}

void MaterialTab::applyMaterials()
{
    auto &mats = m_p->materials();
    for (int r = 0; r < m_mats->rowCount() && r < mats.size(); ++r) {
        Material &m = mats[r];
        if (auto *cb = qobject_cast<QComboBox *>(m_mats->cellWidget(r, 0)))
            m.type = cb->currentIndex() + 1;
        auto cell = [this, r](int c) {
            auto *it = m_mats->item(r, c);
            return it ? it->text() : QString();
        };
        if (m.type == 2) {
            m.einf = cell(1).toDouble();
            m.ae   = cell(2).toDouble();
            m.be   = cell(3).toDouble();
            m.ce   = cell(4).toDouble();
        } else {
            m.epsr = cell(1).toDouble();
            m.esgm = cell(2).toDouble();
            m.amur = cell(3).toDouble();
            m.msgm = cell(4).toDouble();
        }
        m.name = cell(5);
    }
}

void MaterialTab::applyLoads()
{
    auto &loads = m_p->loads();
    for (int r = 0; r < m_loads->rowCount() && r < loads.size(); ++r) {
        Load &l = loads[r];
        auto cell = [this, r](int c) {
            auto *it = m_loads->item(r, c);
            return it ? it->text() : QString();
        };
        if (auto *cb = qobject_cast<QComboBox *>(m_loads->cellWidget(r, 0)))
            l.dir = "XYZ"[cb->currentIndex()];
        l.x = cell(1).toDouble();
        l.y = cell(2).toDouble();
        l.z = cell(3).toDouble();
        if (auto *cb = qobject_cast<QComboBox *>(m_loads->cellWidget(r, 4)))
            l.kind = "RLC"[cb->currentIndex()];
        l.value = cell(5).toDouble();
    }
}

void MaterialTab::refresh()
{
    m_updating = true;

    const auto &mats = m_p->materials();
    m_mats->setRowCount(mats.size());
    for (int r = 0; r < mats.size(); ++r) {
        const Material &m = mats[r];
        auto *type = new QComboBox(m_mats);
        type->addItem(I18n::tr("ma_normal"));
        type->addItem(I18n::tr("ma_dispersive"));
        type->setCurrentIndex(m.type == 2 ? 1 : 0);
        connect(type, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyMaterials();
            refresh();          // re-label value columns
            m_p->touch();
        });
        m_mats->setCellWidget(r, 0, type);

        const double vals[4] = {
            m.type == 2 ? m.einf : m.epsr,
            m.type == 2 ? m.ae   : m.esgm,
            m.type == 2 ? m.be   : m.amur,
            m.type == 2 ? m.ce   : m.msgm };
        for (int c = 0; c < 4; ++c)
            m_mats->setItem(r, 1 + c, new QTableWidgetItem(
                QString::number(vals[c], 'g', 8)));
        m_mats->setItem(r, 5, new QTableWidgetItem(m.name));
        auto *id = new QTableWidgetItem(QString::number(r + 2));
        id->setFlags(id->flags() & ~Qt::ItemIsEditable);
        m_mats->setItem(r, 6, id);
    }

    const auto &loads = m_p->loads();
    m_loads->setRowCount(loads.size());
    for (int r = 0; r < loads.size(); ++r) {
        const Load &l = loads[r];
        auto *dir = new QComboBox(m_loads);
        dir->addItems({ "X", "Y", "Z" });
        dir->setCurrentIndex(l.dir == 'X' ? 0 : l.dir == 'Y' ? 1 : 2);
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyLoads();
            m_p->touch();
        });
        m_loads->setCellWidget(r, 0, dir);
        m_loads->setItem(r, 1, new QTableWidgetItem(QString::number(l.x, 'g', 8)));
        m_loads->setItem(r, 2, new QTableWidgetItem(QString::number(l.y, 'g', 8)));
        m_loads->setItem(r, 3, new QTableWidgetItem(QString::number(l.z, 'g', 8)));
        auto *kind = new QComboBox(m_loads);
        kind->addItems({ "R", "L", "C" });
        kind->setCurrentIndex(l.kind == 'R' ? 0 : l.kind == 'L' ? 1 : 2);
        connect(kind, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyLoads();
            m_p->touch();
        });
        m_loads->setCellWidget(r, 4, kind);
        m_loads->setItem(r, 5, new QTableWidgetItem(QString::number(l.value, 'g', 8)));
    }

    m_updating = false;
}
