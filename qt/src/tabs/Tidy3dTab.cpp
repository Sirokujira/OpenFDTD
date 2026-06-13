// Tidy3dTab.cpp
#include "Tidy3dTab.h"
#include "../core/Project.h"
#include "../io/Tidy3dExporter.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QSettings>
#include <QTableWidget>
#include <QHeaderView>
#include <QVBoxLayout>

using namespace ofd;

Tidy3dTab::Tidy3dTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    auto *s = new SectionBox(I18n::tr("t3_section"), body);
    auto *hint = new QLabel(I18n::tr("t3_hint"), s);
    hint->setWordWrap(true);
    s->vbox()->addWidget(hint);

    m_apiKey = new QLineEdit(s);
    m_apiKey->setEchoMode(QLineEdit::Password);
    m_project = new QLineEdit(s);
    m_resolution = new QComboBox(s);
    m_resolution->addItems({ "coarse", "medium", "fine" });
    m_autoPml = new QCheckBox(s);
    s->form()->addRow(I18n::tr("t3_apikey"), m_apiKey);
    s->form()->addRow(I18n::tr("t3_project"), m_project);
    s->form()->addRow(I18n::tr("t3_resolution"), m_resolution);
    s->form()->addRow(I18n::tr("t3_auto_pml"), m_autoPml);

    auto *exportBtn = new QPushButton(I18n::tr("t3_export"), s);
    s->vbox()->addWidget(exportBtn);
    m_status = new QLabel(s);
    m_status->setWordWrap(true);
    s->vbox()->addWidget(m_status);
    v->addWidget(s);

    // conversion mapping (informational, mirrors the mock's table)
    auto *sm = new SectionBox(I18n::tr("t3_mapping"), body);
    auto *table = new QTableWidget(7, 2, sm);
    table->setHorizontalHeaderLabels({ "OpenFDTD-X", "tidy3d API" });
    table->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    table->verticalHeader()->setVisible(false);
    table->setEditTriggers(QAbstractItemView::NoEditTriggers);
    const char *rows[7][2] = {
        { "geometry (直方体)",      "td.Box" },
        { "geometry (楕円体/円柱)", "td.Sphere / td.Cylinder" },
        { "material (type 1)",      "td.Medium" },
        { "material (type 2 分散)", "td.PoleResidue (要手動調整)" },
        { "feed / planewave",       "td.PointDipole / td.PlaneWave" },
        { "音響パラメータ",          "非対応 (光のみ)" },
        { "集中定数 (R/L/C)",       "部分対応 (td.LumpedElement)" },
    };
    for (int r = 0; r < 7; ++r) {
        table->setItem(r, 0, new QTableWidgetItem(QString::fromUtf8(rows[r][0])));
        table->setItem(r, 1, new QTableWidgetItem(QString::fromUtf8(rows[r][1])));
    }
    table->resizeRowsToContents();
    table->setMinimumHeight(220);
    sm->vbox()->addWidget(table);
    v->addWidget(sm);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    // API key lives in QSettings, not in the project files
    connect(m_apiKey, &QLineEdit::editingFinished, this, [this] {
        QSettings().setValue("tidy3d/apiKey", m_apiKey->text());
    });
    auto applyCb = [this] { apply(); };
    connect(m_project, &QLineEdit::editingFinished, this, applyCb);
    connect(m_resolution, &QComboBox::currentIndexChanged, this, applyCb);
    connect(m_autoPml, &QCheckBox::toggled, this, applyCb);
    connect(exportBtn, &QPushButton::clicked, this, &Tidy3dTab::exportScript);

    connect(project, &Project::loaded, this, &Tidy3dTab::refresh);
    refresh();
}

void Tidy3dTab::apply()
{
    if (m_updating) return;
    Tidy3dOpts &t = m_p->tidy3d();
    t.projectName = m_project->text();
    t.resolution = m_resolution->currentText();
    t.autoPml = m_autoPml->isChecked();
    m_p->touch();
}

void Tidy3dTab::exportScript()
{
    const QString path = QFileDialog::getSaveFileName(
        this, I18n::tr("t3_export"),
        m_p->tidy3d().projectName + ".py", "Python (*.py)");
    if (path.isEmpty()) return;
    QString err;
    if (Tidy3dExporter::exportTo(path, *m_p, &err))
        m_status->setText("OK: " + path);
    else
        m_status->setText("error: " + err);
}

void Tidy3dTab::refresh()
{
    m_updating = true;
    m_apiKey->setText(QSettings().value("tidy3d/apiKey").toString());
    m_project->setText(m_p->tidy3d().projectName);
    m_resolution->setCurrentText(m_p->tidy3d().resolution);
    m_autoPml->setChecked(m_p->tidy3d().autoPml);
    m_updating = false;
}
