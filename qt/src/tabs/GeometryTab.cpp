// GeometryTab.cpp
#include "GeometryTab.h"
#include "../core/Project.h"
#include "../io/StlImporter.h"
#include "../widgets/SectionBox.h"
#include "../widgets/UnitNav.h"
#include "../I18n.h"

#include <QComboBox>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPushButton>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

static const int kShapeCodes[] = { 1, 2, 11, 12, 13, 31, 32, 33,
                                   41, 42, 43, 51, 52, 53 };

static int shapeIndex(int code)
{
    for (size_t i = 0; i < sizeof(kShapeCodes)/sizeof(int); ++i)
        if (kShapeCodes[i] == code) return int(i);
    return 0;
}

GeometryTab::GeometryTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    auto *s = new SectionBox(I18n::tr("ge_section"), body);

    auto *navRow = new QHBoxLayout();
    navRow->addWidget(new QLabel(I18n::tr("ge_unit"), s));
    m_nav = new UnitNav(s);
    navRow->addWidget(m_nav);
    navRow->addStretch(1);
    s->vbox()->addLayout(navRow);

    m_table = new QTableWidget(0, 11, s);
    QStringList headers { I18n::tr("ge_mat"), I18n::tr("ge_shape") };
    for (int i = 1; i <= 8; ++i) headers << QStringLiteral("g%1").arg(i);
    headers << I18n::tr("ma_name");
    m_table->setHorizontalHeaderLabels(headers);
    m_table->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_table->horizontalHeader()->setStretchLastSection(true);
    m_table->verticalHeader()->setDefaultSectionSize(24);
    m_table->setMinimumHeight(200);
    s->vbox()->addWidget(m_table);

    auto *row = new QHBoxLayout();
    auto *add = new QPushButton(I18n::tr("ge_add"), s);
    auto *del = new QPushButton(I18n::tr("ge_del"), s);
    row->addWidget(add);
    row->addWidget(del);
    row->addStretch(1);
    s->vbox()->addLayout(row);
    v->addWidget(s);

    // STL import
    auto *si = new SectionBox(I18n::tr("ge_import"), body);
    auto *importBtn = new QPushButton(I18n::tr("ge_import_btn"), si);
    m_importInfo = new QLabel(I18n::tr("ge_import_hint"), si);
    m_importInfo->setWordWrap(true);
    si->vbox()->addWidget(importBtn);
    si->vbox()->addWidget(m_importInfo);
    v->addWidget(si);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    connect(add, &QPushButton::clicked, this, [this] {
        Geometry g;
        // default to the mesh region so the new unit is visible
        for (int a = 0; a < 3; ++a) {
            g.g[2*a]   = m_p->mesh(a).min();
            g.g[2*a+1] = m_p->mesh(a).max();
        }
        m_p->geometries().push_back(g);
        refresh();
        m_p->touch();
    });
    connect(del, &QPushButton::clicked, this, [this] {
        const int r = m_table->currentRow();
        auto &gs = m_p->geometries();
        if (r >= 0 && r < gs.size()) {
            gs.removeAt(r);
            refresh();
            m_p->touch();
        }
    });
    connect(m_table, &QTableWidget::cellChanged, this, [this] {
        if (m_updating) return;
        applyTable();
        m_p->touch();
    });
    connect(m_table, &QTableWidget::currentCellChanged, this,
            [this](int r, int, int, int) { m_nav->setCurrent(r); });
    connect(m_nav, &UnitNav::currentChanged, this, [this](int i) {
        m_table->selectRow(i);
    });
    connect(importBtn, &QPushButton::clicked, this, &GeometryTab::importStl);

    connect(project, &Project::loaded, this, &GeometryTab::refresh);
    refresh();
}

void GeometryTab::importStl()
{
    const QString path = QFileDialog::getOpenFileName(
        this, I18n::tr("ge_import_btn"), {}, "STL (*.stl);;All files (*)");
    if (path.isEmpty()) return;

    ImportedMesh mesh;
    QString err;
    if (!StlImporter::load(path, mesh, &err)) {
        m_importInfo->setText("error: " + err);
        return;
    }

    // Add the bounding box as a brick unit so the model footprint shows up
    // in the project; true voxelization needs the libigl build (see hint).
    Geometry g;
    g.shape = 1;
    g.g[0] = mesh.bbox[0]; g.g[1] = mesh.bbox[3];
    g.g[2] = mesh.bbox[1]; g.g[3] = mesh.bbox[4];
    g.g[4] = mesh.bbox[2]; g.g[5] = mesh.bbox[5];
    g.name = mesh.name + " (STL bbox)";
    m_p->geometries().push_back(g);
    refresh();
    m_p->touch();

    m_importInfo->setText(QStringLiteral(
        "%1 — %2 triangles, area %3 m², bbox [%4, %5]×[%6, %7]×[%8, %9]\n%10")
        .arg(mesh.name).arg(mesh.numTriangles)
        .arg(QString::number(mesh.surfaceArea, 'g', 4))
        .arg(QString::number(mesh.bbox[0], 'g', 4), QString::number(mesh.bbox[3], 'g', 4),
             QString::number(mesh.bbox[1], 'g', 4), QString::number(mesh.bbox[4], 'g', 4),
             QString::number(mesh.bbox[2], 'g', 4), QString::number(mesh.bbox[5], 'g', 4))
        .arg(I18n::tr("ge_import_hint")));
}

void GeometryTab::applyTable()
{
    auto &gs = m_p->geometries();
    for (int r = 0; r < m_table->rowCount() && r < gs.size(); ++r) {
        Geometry &g = gs[r];
        auto cell = [this, r](int c) {
            auto *it = m_table->item(r, c);
            return it ? it->text() : QString();
        };
        g.materialId = cell(0).toInt();
        if (auto *cb = qobject_cast<QComboBox *>(m_table->cellWidget(r, 1)))
            g.shape = kShapeCodes[cb->currentIndex()];
        for (int i = 0; i < 8; ++i)
            g.g[i] = cell(2 + i).toDouble();
        g.name = cell(10);
    }
}

void GeometryTab::refresh()
{
    m_updating = true;
    const auto &gs = m_p->geometries();
    m_table->setRowCount(gs.size());
    for (int r = 0; r < gs.size(); ++r) {
        const Geometry &g = gs[r];
        m_table->setItem(r, 0, new QTableWidgetItem(QString::number(g.materialId)));

        auto *shape = new QComboBox(m_table);
        for (int code : kShapeCodes)
            shape->addItem(I18n::tr("ge_shape_" + QString::number(code)));
        shape->setCurrentIndex(shapeIndex(g.shape));
        connect(shape, &QComboBox::currentIndexChanged, this, [this] {
            if (m_updating) return;
            applyTable();
            m_p->touch();
        });
        m_table->setCellWidget(r, 1, shape);

        const int np = Geometry::paramCount(g.shape);
        for (int i = 0; i < 8; ++i) {
            auto *it = new QTableWidgetItem(
                i < np ? QString::number(g.g[i], 'g', 8) : QString());
            if (i >= np) it->setFlags(it->flags() & ~Qt::ItemIsEditable);
            m_table->setItem(r, 2 + i, it);
        }
        m_table->setItem(r, 10, new QTableWidgetItem(g.name));
    }
    m_nav->setRange(gs.size());
    m_updating = false;
}
