// RightDock.cpp
#include "RightDock.h"
#include "I18n.h"
#include "core/Project.h"
#include "widgets/LogConsole.h"

#include <QSplitter>
#include <QTreeWidget>
#include <QVBoxLayout>

using namespace ofd;

RightDock::RightDock(Project *project, QWidget *parent)
    : QWidget(parent), m_project(project)
{
    auto *v = new QVBoxLayout(this);
    v->setContentsMargins(0, 0, 0, 0);

    auto *split = new QSplitter(Qt::Vertical, this);

    m_tree = new QTreeWidget(split);
    m_tree->setHeaderLabels({ I18n::tr("rd_project"), "" });
    m_tree->setColumnWidth(0, 160);
    m_tree->setRootIsDecorated(true);

    m_log = new LogConsole(split);

    split->addWidget(m_tree);
    split->addWidget(m_log);
    split->setStretchFactor(0, 3);
    split->setStretchFactor(1, 2);
    v->addWidget(split);

    connect(project, &Project::changed, this, &RightDock::rebuildTree);
    connect(project, &Project::loaded,  this, &RightDock::rebuildTree);
    rebuildTree();
}

void RightDock::appendLog(const QString &line)
{
    m_log->appendLine(line);
}

void RightDock::rebuildTree()
{
    m_tree->clear();

    auto *root = new QTreeWidgetItem(m_tree,
        { m_project->general().title.isEmpty() ? I18n::tr("untitled")
                                               : m_project->general().title });
    root->setExpanded(true);

    auto *mesh = new QTreeWidgetItem(root, { I18n::tr("rd_tree_mesh"),
        QStringLiteral("%L1 cells").arg(m_project->totalCells()) });
    static const char *axisName[3] = { "X", "Y", "Z" };
    for (int a = 0; a < 3; ++a) {
        const MeshAxis &ax = m_project->mesh(a);
        new QTreeWidgetItem(mesh, { axisName[a],
            QStringLiteral("%1 cells [%2, %3]")
                .arg(ax.totalCells())
                .arg(QString::number(ax.min(), 'g', 4),
                     QString::number(ax.max(), 'g', 4)) });
    }

    auto *mats = new QTreeWidgetItem(root, { I18n::tr("rd_tree_materials"),
        QString::number(m_project->materials().size()) });
    int id = 2;
    for (const Material &m : m_project->materials()) {
        const QString desc = (m.type == 2)
            ? QStringLiteral("disp ε∞=%1").arg(m.einf)
            : QStringLiteral("εr=%1 σ=%2").arg(m.epsr).arg(m.esgm);
        new QTreeWidgetItem(mats, {
            QStringLiteral("#%1 %2").arg(id++).arg(m.name), desc });
    }

    auto *geom = new QTreeWidgetItem(root, { I18n::tr("rd_tree_geometry"),
        QString::number(m_project->geometries().size()) });
    int unit = 1;
    for (const Geometry &g : m_project->geometries())
        new QTreeWidgetItem(geom, {
            QStringLiteral("#%1 %2").arg(unit++).arg(
                g.name.isEmpty() ? I18n::tr("ge_shape_" + QString::number(g.shape))
                                 : g.name),
            QStringLiteral("mat %1").arg(g.materialId) });

    auto *srcs = new QTreeWidgetItem(root, { I18n::tr("rd_tree_sources"),
        QString::number(m_project->feeds().size()
                        + (m_project->planewave().enabled ? 1 : 0)) });
    for (const Feed &f : m_project->feeds())
        new QTreeWidgetItem(srcs, { QStringLiteral("feed %1").arg(f.dir),
            QStringLiteral("(%1, %2, %3)").arg(f.x).arg(f.y).arg(f.z) });
    if (m_project->planewave().enabled)
        new QTreeWidgetItem(srcs, { "planewave",
            QStringLiteral("θ=%1 φ=%2").arg(m_project->planewave().theta)
                                       .arg(m_project->planewave().phi) });

    auto *pts = new QTreeWidgetItem(root, { I18n::tr("rd_tree_points"),
        QString::number(m_project->probes().size()) });
    for (const Probe &pr : m_project->probes())
        new QTreeWidgetItem(pts, { QStringLiteral("point %1").arg(pr.dir),
            QStringLiteral("(%1, %2, %3)").arg(pr.x).arg(pr.y).arg(pr.z) });

    if (!m_project->loads().isEmpty()) {
        auto *lds = new QTreeWidgetItem(root, { I18n::tr("rd_tree_loads"),
            QString::number(m_project->loads().size()) });
        for (const Load &l : m_project->loads())
            new QTreeWidgetItem(lds, { QStringLiteral("%1").arg(l.kind),
                QStringLiteral("%1 @ (%2, %3, %4)")
                    .arg(l.value).arg(l.x).arg(l.y).arg(l.z) });
    }
}
