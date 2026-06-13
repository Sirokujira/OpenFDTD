// Post2Tab.cpp
#include "Post2Tab.h"
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
#include <QSpinBox>
#include <QTableWidget>
#include <QVBoxLayout>

using namespace ofd;

Post2Tab::Post2Tab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    auto applyCb = [this] { apply(); };

    // far0d
    auto *s0 = new SectionBox(I18n::tr("p2_far0d"), body);
    m_far0d = new QCheckBox(I18n::tr("p2_far0d"), s0);
    m_far0dTheta = new QLineEdit(s0); m_far0dTheta->setMaximumWidth(80);
    m_far0dPhi   = new QLineEdit(s0); m_far0dPhi->setMaximumWidth(80);
    auto *r0 = new QHBoxLayout();
    r0->addWidget(m_far0d, 1);
    r0->addWidget(new QLabel("θ", s0));
    r0->addWidget(m_far0dTheta);
    r0->addWidget(new QLabel("φ", s0));
    r0->addWidget(m_far0dPhi);
    s0->vbox()->addLayout(r0);
    v->addWidget(s0);

    // far1d
    auto *s1 = new SectionBox(I18n::tr("p2_far1d"), body);
    m_far1d = new QTableWidget(0, 3, s1);
    m_far1d->setHorizontalHeaderLabels({
        I18n::tr("p2_dir"), I18n::tr("p2_division"), I18n::tr("p2_angle") });
    m_far1d->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_far1d->verticalHeader()->setDefaultSectionSize(24);
    m_far1d->setMinimumHeight(90);
    s1->vbox()->addWidget(m_far1d);
    auto *r1 = new QHBoxLayout();
    auto *add1 = new QPushButton(I18n::tr("p2_add"), s1);
    auto *del1 = new QPushButton(I18n::tr("p2_del"), s1);
    r1->addWidget(add1); r1->addWidget(del1); r1->addStretch(1);
    s1->vbox()->addLayout(r1);

    auto *r1b = new QHBoxLayout();
    r1b->addWidget(new QLabel(I18n::tr("p2_style"), s1));
    m_far1dStyle = new QComboBox(s1);
    m_far1dStyle->addItems({ "0", "1", "2" });
    r1b->addWidget(m_far1dStyle);
    m_far1dDb = new QCheckBox(I18n::tr("p2_db"), s1);
    r1b->addWidget(m_far1dDb);
    m_far1dNorm = new QCheckBox(I18n::tr("p2_norm"), s1);
    r1b->addWidget(m_far1dNorm);
    r1b->addStretch(1);
    s1->vbox()->addLayout(r1b);

    auto *r1c = new QHBoxLayout();
    r1c->addWidget(new QLabel(I18n::tr("p2_component"), s1));
    m_far1dCompE     = new QCheckBox("E", s1);
    m_far1dCompTheta = new QCheckBox("Eθ", s1);
    m_far1dCompPhi   = new QCheckBox("Eφ", s1);
    r1c->addWidget(m_far1dCompE);
    r1c->addWidget(m_far1dCompTheta);
    r1c->addWidget(m_far1dCompPhi);
    r1c->addStretch(1);
    s1->vbox()->addLayout(r1c);
    v->addWidget(s1);

    // far2d
    auto *s2 = new SectionBox(I18n::tr("p2_far2d"), body);
    m_far2d = new QCheckBox(I18n::tr("p2_far2d"), s2);
    m_far2dTheta = new QSpinBox(s2); m_far2dTheta->setRange(1, 3600);
    m_far2dPhi   = new QSpinBox(s2); m_far2dPhi->setRange(1, 3600);
    m_far2dDb = new QCheckBox(I18n::tr("p2_db"), s2);
    m_far2dObj = new QLineEdit(s2); m_far2dObj->setMaximumWidth(80);
    auto *r2 = new QHBoxLayout();
    r2->addWidget(m_far2d, 1);
    r2->addWidget(new QLabel(I18n::tr("p2_theta_div"), s2));
    r2->addWidget(m_far2dTheta);
    r2->addWidget(new QLabel(I18n::tr("p2_phi_div"), s2));
    r2->addWidget(m_far2dPhi);
    r2->addWidget(m_far2dDb);
    r2->addWidget(new QLabel(I18n::tr("p2_obj"), s2));
    r2->addWidget(m_far2dObj);
    s2->vbox()->addLayout(r2);
    v->addWidget(s2);

    // near1d
    auto *s3 = new SectionBox(I18n::tr("p2_near1d"), body);
    m_near1d = new QTableWidget(0, 4, s3);
    m_near1d->setHorizontalHeaderLabels({
        I18n::tr("p2_cmp"), I18n::tr("p2_dir"),
        I18n::tr("p2_pos") + " 1", I18n::tr("p2_pos") + " 2" });
    m_near1d->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_near1d->verticalHeader()->setDefaultSectionSize(24);
    m_near1d->setMinimumHeight(90);
    s3->vbox()->addWidget(m_near1d);
    auto *r3 = new QHBoxLayout();
    auto *add3 = new QPushButton(I18n::tr("p2_add"), s3);
    auto *del3 = new QPushButton(I18n::tr("p2_del"), s3);
    m_near1dDb = new QCheckBox(I18n::tr("p2_db"), s3);
    m_near1dNoinc = new QCheckBox(I18n::tr("p2_noinc"), s3);
    r3->addWidget(add3); r3->addWidget(del3);
    r3->addWidget(m_near1dDb); r3->addWidget(m_near1dNoinc);
    r3->addStretch(1);
    s3->vbox()->addLayout(r3);
    v->addWidget(s3);

    // near2d
    auto *s4 = new SectionBox(I18n::tr("p2_near2d"), body);
    m_near2d = new QTableWidget(0, 3, s4);
    m_near2d->setHorizontalHeaderLabels({
        I18n::tr("p2_cmp"), I18n::tr("p2_dir"), I18n::tr("p2_pos") });
    m_near2d->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    m_near2d->verticalHeader()->setDefaultSectionSize(24);
    m_near2d->setMinimumHeight(90);
    s4->vbox()->addWidget(m_near2d);
    auto *r4 = new QHBoxLayout();
    auto *add4 = new QPushButton(I18n::tr("p2_add"), s4);
    auto *del4 = new QPushButton(I18n::tr("p2_del"), s4);
    r4->addWidget(add4); r4->addWidget(del4); r4->addStretch(1);
    s4->vbox()->addLayout(r4);
    auto *r4b = new QHBoxLayout();
    r4b->addWidget(new QLabel("dim", s4));
    m_near2dDim0 = new QSpinBox(s4); m_near2dDim0->setRange(0, 1);
    m_near2dDim1 = new QSpinBox(s4); m_near2dDim1->setRange(0, 1);
    m_near2dDb = new QCheckBox(I18n::tr("p2_db"), s4);
    m_near2dContour = new QCheckBox(I18n::tr("p2_contour"), s4);
    m_near2dNoinc = new QCheckBox(I18n::tr("p2_noinc"), s4);
    r4b->addWidget(m_near2dDim0);
    r4b->addWidget(m_near2dDim1);
    r4b->addWidget(m_near2dDb);
    r4b->addWidget(m_near2dContour);
    r4b->addWidget(m_near2dNoinc);
    r4b->addStretch(1);
    s4->vbox()->addLayout(r4b);
    v->addWidget(s4);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    // wiring
    for (auto *c : { m_far0d, m_far1dDb, m_far1dNorm, m_far1dCompE,
                     m_far1dCompTheta, m_far1dCompPhi, m_far2d, m_far2dDb,
                     m_near1dDb, m_near1dNoinc, m_near2dDb, m_near2dContour,
                     m_near2dNoinc })
        connect(c, &QCheckBox::toggled, this, applyCb);
    for (auto *e : { m_far0dTheta, m_far0dPhi, m_far2dObj })
        connect(e, &QLineEdit::editingFinished, this, applyCb);
    for (auto *sp : { m_far2dTheta, m_far2dPhi, m_near2dDim0, m_near2dDim1 })
        connect(sp, &QSpinBox::valueChanged, this, applyCb);
    connect(m_far1dStyle, &QComboBox::currentIndexChanged, this, applyCb);

    connect(add1, &QPushButton::clicked, this, [this] {
        m_p->post().far1d.push_back(Far1d{});
        refresh(); m_p->touch();
    });
    connect(del1, &QPushButton::clicked, this, [this] {
        auto &v = m_p->post().far1d;
        const int r = m_far1d->currentRow();
        if (r >= 0 && r < v.size()) { v.removeAt(r); refresh(); m_p->touch(); }
    });
    connect(m_far1d, &QTableWidget::cellChanged, this, [this] {
        if (!m_updating) { applyFar1dTable(); m_p->touch(); }
    });

    connect(add3, &QPushButton::clicked, this, [this] {
        m_p->post().near1d.push_back(Near1d{});
        refresh(); m_p->touch();
    });
    connect(del3, &QPushButton::clicked, this, [this] {
        auto &v = m_p->post().near1d;
        const int r = m_near1d->currentRow();
        if (r >= 0 && r < v.size()) { v.removeAt(r); refresh(); m_p->touch(); }
    });
    connect(m_near1d, &QTableWidget::cellChanged, this, [this] {
        if (!m_updating) { applyNear1dTable(); m_p->touch(); }
    });

    connect(add4, &QPushButton::clicked, this, [this] {
        m_p->post().near2d.push_back(Near2d{});
        refresh(); m_p->touch();
    });
    connect(del4, &QPushButton::clicked, this, [this] {
        auto &v = m_p->post().near2d;
        const int r = m_near2d->currentRow();
        if (r >= 0 && r < v.size()) { v.removeAt(r); refresh(); m_p->touch(); }
    });
    connect(m_near2d, &QTableWidget::cellChanged, this, [this] {
        if (!m_updating) { applyNear2dTable(); m_p->touch(); }
    });

    connect(project, &Project::loaded, this, &Post2Tab::refresh);
    refresh();
}

void Post2Tab::apply()
{
    if (m_updating) return;
    PostOpts &po = m_p->post();
    po.far0d = m_far0d->isChecked();
    po.far0dTheta = m_far0dTheta->text().toDouble();
    po.far0dPhi   = m_far0dPhi->text().toDouble();
    po.far1dStyle = m_far1dStyle->currentIndex();
    po.far1dDb = m_far1dDb->isChecked();
    po.far1dNorm = m_far1dNorm->isChecked();
    po.far1dComp[0] = m_far1dCompE->isChecked();
    po.far1dComp[1] = m_far1dCompTheta->isChecked();
    po.far1dComp[2] = m_far1dCompPhi->isChecked();
    po.far2d = m_far2d->isChecked();
    po.far2dDivTheta = m_far2dTheta->value();
    po.far2dDivPhi   = m_far2dPhi->value();
    po.far2dDb = m_far2dDb->isChecked();
    po.far2dObj = m_far2dObj->text().toDouble();
    po.near1dDb = m_near1dDb->isChecked();
    po.near1dNoinc = m_near1dNoinc->isChecked();
    po.near2dDim[0] = m_near2dDim0->value();
    po.near2dDim[1] = m_near2dDim1->value();
    po.near2dDb = m_near2dDb->isChecked();
    po.near2dContour = m_near2dContour->isChecked();
    po.near2dNoinc = m_near2dNoinc->isChecked();
    m_p->touch();
}

void Post2Tab::applyFar1dTable()
{
    auto &v = m_p->post().far1d;
    for (int r = 0; r < m_far1d->rowCount() && r < v.size(); ++r) {
        if (auto *cb = qobject_cast<QComboBox *>(m_far1d->cellWidget(r, 0)))
            v[r].dir = "XYZVH"[cb->currentIndex()];
        if (auto *it = m_far1d->item(r, 1)) v[r].div = it->text().toInt();
        if (auto *it = m_far1d->item(r, 2)) v[r].angle = it->text().toDouble();
    }
}

void Post2Tab::applyNear1dTable()
{
    auto &v = m_p->post().near1d;
    for (int r = 0; r < m_near1d->rowCount() && r < v.size(); ++r) {
        if (auto *it = m_near1d->item(r, 0)) v[r].cmp = it->text();
        if (auto *cb = qobject_cast<QComboBox *>(m_near1d->cellWidget(r, 1)))
            v[r].dir = "XYZ"[cb->currentIndex()];
        if (auto *it = m_near1d->item(r, 2)) v[r].pos1 = it->text().toDouble();
        if (auto *it = m_near1d->item(r, 3)) v[r].pos2 = it->text().toDouble();
    }
}

void Post2Tab::applyNear2dTable()
{
    auto &v = m_p->post().near2d;
    for (int r = 0; r < m_near2d->rowCount() && r < v.size(); ++r) {
        if (auto *it = m_near2d->item(r, 0)) v[r].cmp = it->text();
        if (auto *cb = qobject_cast<QComboBox *>(m_near2d->cellWidget(r, 1)))
            v[r].dir = "XYZ"[cb->currentIndex()];
        if (auto *it = m_near2d->item(r, 2)) v[r].pos = it->text().toDouble();
    }
}

void Post2Tab::refresh()
{
    m_updating = true;
    const PostOpts &po = m_p->post();

    m_far0d->setChecked(po.far0d);
    m_far0dTheta->setText(QString::number(po.far0dTheta, 'g', 6));
    m_far0dPhi->setText(QString::number(po.far0dPhi, 'g', 6));

    m_far1d->setRowCount(po.far1d.size());
    for (int r = 0; r < po.far1d.size(); ++r) {
        const Far1d &f = po.far1d[r];
        auto *dir = new QComboBox(m_far1d);
        dir->addItems({ "X", "Y", "Z", "V", "H" });
        const int di = QString("XYZVH").indexOf(f.dir);
        dir->setCurrentIndex(qMax(0, di));
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (!m_updating) { applyFar1dTable(); m_p->touch(); }
        });
        m_far1d->setCellWidget(r, 0, dir);
        m_far1d->setItem(r, 1, new QTableWidgetItem(QString::number(f.div)));
        m_far1d->setItem(r, 2, new QTableWidgetItem(
            QString::number(f.angle, 'g', 6)));
    }
    m_far1dStyle->setCurrentIndex(qBound(0, po.far1dStyle, 2));
    m_far1dDb->setChecked(po.far1dDb);
    m_far1dNorm->setChecked(po.far1dNorm);
    m_far1dCompE->setChecked(po.far1dComp[0]);
    m_far1dCompTheta->setChecked(po.far1dComp[1]);
    m_far1dCompPhi->setChecked(po.far1dComp[2]);

    m_far2d->setChecked(po.far2d);
    m_far2dTheta->setValue(po.far2dDivTheta);
    m_far2dPhi->setValue(po.far2dDivPhi);
    m_far2dDb->setChecked(po.far2dDb);
    m_far2dObj->setText(QString::number(po.far2dObj, 'g', 6));

    m_near1d->setRowCount(po.near1d.size());
    for (int r = 0; r < po.near1d.size(); ++r) {
        const Near1d &n = po.near1d[r];
        m_near1d->setItem(r, 0, new QTableWidgetItem(n.cmp));
        auto *dir = new QComboBox(m_near1d);
        dir->addItems({ "X", "Y", "Z" });
        dir->setCurrentIndex(n.dir == 'X' ? 0 : n.dir == 'Y' ? 1 : 2);
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (!m_updating) { applyNear1dTable(); m_p->touch(); }
        });
        m_near1d->setCellWidget(r, 1, dir);
        m_near1d->setItem(r, 2, new QTableWidgetItem(QString::number(n.pos1, 'g', 6)));
        m_near1d->setItem(r, 3, new QTableWidgetItem(QString::number(n.pos2, 'g', 6)));
    }
    m_near1dDb->setChecked(po.near1dDb);
    m_near1dNoinc->setChecked(po.near1dNoinc);

    m_near2d->setRowCount(po.near2d.size());
    for (int r = 0; r < po.near2d.size(); ++r) {
        const Near2d &n = po.near2d[r];
        m_near2d->setItem(r, 0, new QTableWidgetItem(n.cmp));
        auto *dir = new QComboBox(m_near2d);
        dir->addItems({ "X", "Y", "Z" });
        dir->setCurrentIndex(n.dir == 'X' ? 0 : n.dir == 'Y' ? 1 : 2);
        connect(dir, &QComboBox::currentIndexChanged, this, [this] {
            if (!m_updating) { applyNear2dTable(); m_p->touch(); }
        });
        m_near2d->setCellWidget(r, 1, dir);
        m_near2d->setItem(r, 2, new QTableWidgetItem(QString::number(n.pos, 'g', 6)));
    }
    m_near2dDim0->setValue(po.near2dDim[0]);
    m_near2dDim1->setValue(po.near2dDim[1]);
    m_near2dDb->setChecked(po.near2dDb);
    m_near2dContour->setChecked(po.near2dContour);
    m_near2dNoinc->setChecked(po.near2dNoinc);

    m_updating = false;
}
