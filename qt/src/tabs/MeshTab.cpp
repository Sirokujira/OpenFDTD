// MeshTab.cpp
#include "MeshTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

MeshTab::MeshTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    m_total = new QLabel(body);
    v->addWidget(m_total);

    static const char *secKey[3] = { "me_axis_x", "me_axis_y", "me_axis_z" };
    for (int a = 0; a < 3; ++a) {
        auto *s = new SectionBox(I18n::tr(secKey[a]), body);

        m_table[a] = new QTableWidget(0, 2, s);
        m_table[a]->setHorizontalHeaderLabels(
            { I18n::tr("me_coord"), I18n::tr("me_div") });
        m_table[a]->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
        m_table[a]->verticalHeader()->setDefaultSectionSize(22);
        m_table[a]->setMinimumHeight(110);
        s->vbox()->addWidget(m_table[a]);

        auto *btnRow = new QHBoxLayout();
        auto *add = new QPushButton(I18n::tr("me_add"), s);
        auto *del = new QPushButton(I18n::tr("me_del"), s);
        m_info[a] = new QLabel(s);
        btnRow->addWidget(add);
        btnRow->addWidget(del);
        btnRow->addStretch(1);
        btnRow->addWidget(m_info[a]);
        s->vbox()->addLayout(btnRow);
        v->addWidget(s);

        connect(add, &QPushButton::clicked, this, [this, a] {
            MeshAxis &ax = m_p->mesh(a);
            const double last = ax.nodes.isEmpty() ? 0.0 : ax.nodes.last();
            const double step = ax.nodes.size() >= 2
                ? ax.nodes.last() - ax.nodes[ax.nodes.size()-2] : 0.05;
            ax.nodes.push_back(last + (step > 0 ? step : 0.05));
            ax.divs.push_back(10);
            refresh();
            m_p->touch();
        });
        connect(del, &QPushButton::clicked, this, [this, a] {
            MeshAxis &ax = m_p->mesh(a);
            if (ax.nodes.size() <= 2) return;
            const int row = m_table[a]->currentRow();
            const int i = (row >= 0 && row < ax.nodes.size())
                          ? row : ax.nodes.size() - 1;
            ax.nodes.removeAt(i);
            ax.divs.removeAt(qMin(i, ax.divs.size() - 1));
            refresh();
            m_p->touch();
        });
        connect(m_table[a], &QTableWidget::cellChanged, this, [this, a] {
            if (m_updating) return;
            applyAxis(a);
            refreshAxisInfo(a);
            m_p->touch();
        });
    }

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    connect(project, &Project::loaded, this, &MeshTab::refresh);
    refresh();
}

void MeshTab::applyAxis(int a)
{
    MeshAxis &ax = m_p->mesh(a);
    const int rows = m_table[a]->rowCount();
    ax.nodes.resize(rows);
    ax.divs.resize(qMax(0, rows - 1));
    for (int r = 0; r < rows; ++r) {
        if (auto *it = m_table[a]->item(r, 0))
            ax.nodes[r] = it->text().toDouble();
        if (r < rows - 1) {
            if (auto *it = m_table[a]->item(r, 1))
                ax.divs[r] = qMax(1, it->text().toInt());
        }
    }
}

void MeshTab::refreshAxisInfo(int a)
{
    const MeshAxis &ax = m_p->mesh(a);
    m_info[a]->setText(QStringLiteral("%1: %2   %3: %4")
        .arg(I18n::tr("me_cells")).arg(ax.totalCells())
        .arg(I18n::tr("me_dmin"),
             ax.isValid() ? QString::number(ax.minSpacing(), 'g', 4) : "—"));
    m_total->setText(QStringLiteral("%1: %L2")
        .arg(I18n::tr("me_total")).arg(m_p->totalCells()));
}

void MeshTab::refresh()
{
    m_updating = true;
    for (int a = 0; a < 3; ++a) {
        const MeshAxis &ax = m_p->mesh(a);
        m_table[a]->setRowCount(ax.nodes.size());
        for (int r = 0; r < ax.nodes.size(); ++r) {
            m_table[a]->setItem(r, 0, new QTableWidgetItem(
                QString::number(ax.nodes[r], 'g', 10)));
            auto *divItem = new QTableWidgetItem(
                r < ax.divs.size() ? QString::number(ax.divs[r]) : QString("—"));
            if (r >= ax.divs.size())
                divItem->setFlags(divItem->flags() & ~Qt::ItemIsEditable);
            m_table[a]->setItem(r, 1, divItem);
        }
        refreshAxisInfo(a);
    }
    m_updating = false;
}
