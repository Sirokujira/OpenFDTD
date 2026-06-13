// OpticalTab.cpp
#include "OpticalTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QStackedWidget>
#include <QVBoxLayout>

using namespace ofd;

OpticalTab::OpticalTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    // solver method selection (FDTD / RCWA / BPM / FMM)
    auto *ss = new SectionBox(I18n::tr("opt_solver"), body);
    m_solver = new QComboBox(ss);
    m_solver->addItem(I18n::tr("opt_solver_fdtd"));
    m_solver->addItem(I18n::tr("opt_solver_rcwa"));
    m_solver->addItem(I18n::tr("opt_solver_bpm"));
    m_solver->addItem(I18n::tr("opt_solver_fmm"));
    ss->vbox()->addWidget(m_solver);
    auto *hint = new QLabel(I18n::tr("opt_kernel_hint"), ss);
    hint->setWordWrap(true);
    ss->vbox()->addWidget(hint);

    // per-method parameter pages
    m_solverStack = new QStackedWidget(ss);

    // [0] FDTD — uses the General/Mesh tabs; nothing extra here
    {
        auto *page = new QLabel(I18n::tr("opt_solver_fdtd"), m_solverStack);
        page->setWordWrap(true);
        page->setAlignment(Qt::AlignTop | Qt::AlignLeft);
        m_solverStack->addWidget(page);
    }
    // [1] RCWA
    {
        auto *page = new SectionBox(I18n::tr("opt_rcwa_section"), m_solverStack);
        m_rcwaNx = new QSpinBox(page); m_rcwaNx->setRange(1, 101);
        m_rcwaNy = new QSpinBox(page); m_rcwaNy->setRange(1, 101);
        m_rcwaPx = new QLineEdit(page); m_rcwaPx->setMaximumWidth(100);
        m_rcwaPy = new QLineEdit(page); m_rcwaPy->setMaximumWidth(100);
        m_rcwaLayers = new QSpinBox(page); m_rcwaLayers->setRange(1, 1000);
        page->form()->addRow(I18n::tr("opt_rcwa_orders") + " Nx", m_rcwaNx);
        page->form()->addRow("Ny", m_rcwaNy);
        page->form()->addRow(I18n::tr("opt_rcwa_period") + " Λx", m_rcwaPx);
        page->form()->addRow("Λy", m_rcwaPy);
        page->form()->addRow(I18n::tr("opt_rcwa_layers"), m_rcwaLayers);
        m_solverStack->addWidget(page);
    }
    // [2] BPM
    {
        auto *page = new SectionBox(I18n::tr("opt_bpm_section"), m_solverStack);
        m_bpmAlgo = new QComboBox(page);
        m_bpmAlgo->addItems({ "FFT-BPM", "FDM-BPM", "Wide-Angle Padé(1,1)" });
        m_bpmDz = new QLineEdit(page); m_bpmDz->setMaximumWidth(100);
        m_bpmN0 = new QLineEdit(page); m_bpmN0->setMaximumWidth(100);
        m_bpmInput = new QComboBox(page);
        m_bpmInput->addItems({ "TE₀", "TE₁", "TM₀", "Gaussian" });
        page->form()->addRow(I18n::tr("opt_bpm_algo"), m_bpmAlgo);
        page->form()->addRow(I18n::tr("opt_bpm_dz"), m_bpmDz);
        page->form()->addRow(I18n::tr("opt_bpm_n0"), m_bpmN0);
        page->form()->addRow(I18n::tr("opt_bpm_input"), m_bpmInput);
        m_solverStack->addWidget(page);
    }
    // [3] FMM
    {
        auto *page = new SectionBox(I18n::tr("opt_fmm_section"), m_solverStack);
        m_fmmHarmonics = new QSpinBox(page);
        m_fmmHarmonics->setRange(1, 201);
        m_fmmLi = new QCheckBox(page);
        page->form()->addRow(I18n::tr("opt_fmm_harmonics"), m_fmmHarmonics);
        page->form()->addRow(I18n::tr("opt_fmm_li"), m_fmmLi);
        m_solverStack->addWidget(page);
    }
    ss->vbox()->addWidget(m_solverStack);
    v->addWidget(ss);

    // wavelength range
    auto *sw = new SectionBox(I18n::tr("opt_wavelength"), body);
    m_lambdaMin = new QLineEdit(sw); m_lambdaMin->setMaximumWidth(100);
    m_lambdaMax = new QLineEdit(sw); m_lambdaMax->setMaximumWidth(100);
    m_lambdaDiv = new QSpinBox(sw);  m_lambdaDiv->setRange(2, 100000);
    sw->form()->addRow(I18n::tr("opt_lambda_min"), m_lambdaMin);
    sw->form()->addRow(I18n::tr("opt_lambda_max"), m_lambdaMax);
    sw->form()->addRow(I18n::tr("opt_lambda_div"), m_lambdaDiv);
    v->addWidget(sw);

    // optical mode
    auto *sm = new SectionBox(I18n::tr("opt_mode"), body);
    m_mode = new QComboBox(sm);
    m_mode->addItem(I18n::tr("opt_bpf"));
    m_mode->addItem(I18n::tr("opt_wg"));
    m_mode->addItem(I18n::tr("opt_ring"));
    m_mode->addItem(I18n::tr("opt_mzi"));
    m_mode->addItem(I18n::tr("opt_meta"));
    m_mode->addItem(I18n::tr("opt_phc"));
    m_mode->addItem(I18n::tr("opt_nf2ff"));
    m_mode->addItem(I18n::tr("opt_spara"));
    sm->vbox()->addWidget(m_mode);
    v->addWidget(sm);

    // BPF spec
    auto *sb = new SectionBox(I18n::tr("opt_bpf_section"), body);
    m_bpfMin = new QLineEdit(sb); m_bpfMin->setMaximumWidth(100);
    m_bpfMax = new QLineEdit(sb); m_bpfMax->setMaximumWidth(100);
    m_bpfQ = new QLineEdit(sb); m_bpfQ->setMaximumWidth(100);
    sb->form()->addRow(I18n::tr("opt_band") + " min", m_bpfMin);
    sb->form()->addRow(I18n::tr("opt_band") + " max", m_bpfMax);
    sb->form()->addRow(I18n::tr("opt_q"), m_bpfQ);
    v->addWidget(sb);

    // Ring spec
    auto *sr = new SectionBox(I18n::tr("opt_ring_section"), body);
    m_ringR = new QLineEdit(sr); m_ringR->setMaximumWidth(100);
    m_ringGap = new QLineEdit(sr); m_ringGap->setMaximumWidth(100);
    sr->form()->addRow(I18n::tr("opt_radius"), m_ringR);
    sr->form()->addRow(I18n::tr("opt_gap"), m_ringGap);
    v->addWidget(sr);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    auto applyCb = [this] { apply(); };
    connect(m_solver, &QComboBox::currentIndexChanged, this, [this](int i) {
        m_solverStack->setCurrentIndex(i);
        apply();
    });
    connect(m_mode, &QComboBox::currentIndexChanged, this, applyCb);
    connect(m_bpmAlgo, &QComboBox::currentIndexChanged, this, applyCb);
    connect(m_bpmInput, &QComboBox::currentIndexChanged, this, applyCb);
    for (auto *e : { m_lambdaMin, m_lambdaMax, m_rcwaPx, m_rcwaPy, m_bpmDz,
                     m_bpmN0, m_bpfMin, m_bpfMax, m_bpfQ, m_ringR, m_ringGap })
        connect(e, &QLineEdit::editingFinished, this, applyCb);
    for (auto *s : { m_lambdaDiv, m_rcwaNx, m_rcwaNy, m_rcwaLayers,
                     m_fmmHarmonics })
        connect(s, &QSpinBox::valueChanged, this, applyCb);
    connect(m_fmmLi, &QCheckBox::toggled, this, applyCb);

    connect(project, &Project::loaded, this, &OpticalTab::refresh);
    refresh();
}

void OpticalTab::apply()
{
    if (m_updating) return;
    OpticalOpts &o = m_p->optical();
    o.solver = OpticalSolver(m_solver->currentIndex());
    o.mode   = OpticalMode(m_mode->currentIndex());
    o.lambdaMin = m_lambdaMin->text().toDouble();
    o.lambdaMax = m_lambdaMax->text().toDouble();
    o.lambdaDiv = m_lambdaDiv->value();
    o.rcwaNx = m_rcwaNx->value();
    o.rcwaNy = m_rcwaNy->value();
    o.rcwaPeriodX = m_rcwaPx->text().toDouble();
    o.rcwaPeriodY = m_rcwaPy->text().toDouble();
    o.rcwaLayers = m_rcwaLayers->value();
    o.bpmAlgorithm = m_bpmAlgo->currentIndex();
    o.bpmDz = m_bpmDz->text().toDouble();
    o.bpmRefIndex = m_bpmN0->text().toDouble();
    o.bpmInputMode = m_bpmInput->currentIndex();
    o.fmmHarmonics = m_fmmHarmonics->value();
    o.fmmLiRules = m_fmmLi->isChecked();
    o.bpfBandMin = m_bpfMin->text().toDouble();
    o.bpfBandMax = m_bpfMax->text().toDouble();
    o.bpfQ = m_bpfQ->text().toDouble();
    o.ringRadius_um = m_ringR->text().toDouble();
    o.ringGap_nm = m_ringGap->text().toDouble();
    m_p->touch();
}

void OpticalTab::refresh()
{
    m_updating = true;
    const OpticalOpts &o = m_p->optical();
    m_solver->setCurrentIndex(int(o.solver));
    m_solverStack->setCurrentIndex(int(o.solver));
    m_mode->setCurrentIndex(int(o.mode));
    m_lambdaMin->setText(QString::number(o.lambdaMin, 'g', 8));
    m_lambdaMax->setText(QString::number(o.lambdaMax, 'g', 8));
    m_lambdaDiv->setValue(o.lambdaDiv);
    m_rcwaNx->setValue(o.rcwaNx);
    m_rcwaNy->setValue(o.rcwaNy);
    m_rcwaPx->setText(QString::number(o.rcwaPeriodX, 'g', 8));
    m_rcwaPy->setText(QString::number(o.rcwaPeriodY, 'g', 8));
    m_rcwaLayers->setValue(o.rcwaLayers);
    m_bpmAlgo->setCurrentIndex(o.bpmAlgorithm);
    m_bpmDz->setText(QString::number(o.bpmDz, 'g', 8));
    m_bpmN0->setText(QString::number(o.bpmRefIndex, 'g', 8));
    m_bpmInput->setCurrentIndex(o.bpmInputMode);
    m_fmmHarmonics->setValue(o.fmmHarmonics);
    m_fmmLi->setChecked(o.fmmLiRules);
    m_bpfMin->setText(QString::number(o.bpfBandMin, 'g', 8));
    m_bpfMax->setText(QString::number(o.bpfBandMax, 'g', 8));
    m_bpfQ->setText(QString::number(o.bpfQ, 'g', 8));
    m_ringR->setText(QString::number(o.ringRadius_um, 'g', 8));
    m_ringGap->setText(QString::number(o.ringGap_nm, 'g', 8));
    m_updating = false;
}
