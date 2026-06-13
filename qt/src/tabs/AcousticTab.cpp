// AcousticTab.cpp
#include "AcousticTab.h"
#include "../core/Project.h"
#include "../widgets/SectionBox.h"
#include "../I18n.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QSpinBox>
#include <QVBoxLayout>

using namespace ofd;

AcousticTab::AcousticTab(Project *project, QWidget *parent)
    : QScrollArea(parent), m_p(project)
{
    auto *body = new QWidget(this);
    auto *v = new QVBoxLayout(body);
    v->setContentsMargins(8, 8, 8, 8);
    v->setSpacing(8);

    auto *hint = new QLabel(I18n::tr("ac_mapping_hint"), body);
    hint->setWordWrap(true);
    v->addWidget(hint);

    auto *sm = new SectionBox(I18n::tr("ac_metrics"), body);
    m_rt60 = new QCheckBox(I18n::tr("ac_rt60"), sm);
    m_c80  = new QCheckBox(I18n::tr("ac_c80"), sm);
    m_d50  = new QCheckBox(I18n::tr("ac_d50"), sm);
    m_sti  = new QCheckBox(I18n::tr("ac_sti"), sm);
    m_edt  = new QCheckBox(I18n::tr("ac_edt"), sm);
    m_irf  = new QCheckBox(I18n::tr("ac_irf"), sm);
    m_aural = new QCheckBox(I18n::tr("ac_aurora"), sm);
    for (auto *c : { m_rt60, m_c80, m_d50, m_sti, m_edt, m_irf, m_aural })
        sm->vbox()->addWidget(c);
    m_sampleRate = new QComboBox(sm);
    m_sampleRate->addItems({ "44100", "48000", "96000" });
    sm->form()->addRow(I18n::tr("ac_sample_rate"), m_sampleRate);
    v->addWidget(sm);

    auto *ss = new SectionBox(I18n::tr("ac_source"), body);
    m_directivity = new QComboBox(ss);
    m_directivity->addItem(I18n::tr("ac_omni"));      // omni
    m_directivity->addItem(I18n::tr("ac_cardioid"));  // cardioid
    m_directivity->addItem(I18n::tr("ac_speaker"));   // speaker
    m_spl = new QDoubleSpinBox(ss);
    m_spl->setRange(0, 200);
    m_spl->setSuffix(" dB");
    ss->form()->addRow(I18n::tr("ac_directivity"), m_directivity);
    ss->form()->addRow(I18n::tr("ac_spl"), m_spl);
    v->addWidget(ss);

    auto *sr = new SectionBox(I18n::tr("ac_mics"), body);
    m_micCount = new QSpinBox(sr);
    m_micCount->setRange(1, 256);
    sr->form()->addRow(I18n::tr("ac_mic_count"), m_micCount);
    v->addWidget(sr);

    v->addStretch(1);
    setWidget(body);
    setWidgetResizable(true);
    setFrameShape(QFrame::NoFrame);

    auto applyCb = [this] { apply(); };
    for (auto *c : { m_rt60, m_c80, m_d50, m_sti, m_edt, m_irf, m_aural })
        connect(c, &QCheckBox::toggled, this, applyCb);
    connect(m_sampleRate, &QComboBox::currentIndexChanged, this, applyCb);
    connect(m_directivity, &QComboBox::currentIndexChanged, this, applyCb);
    connect(m_spl, &QDoubleSpinBox::valueChanged, this, applyCb);
    connect(m_micCount, &QSpinBox::valueChanged, this, applyCb);

    connect(project, &Project::loaded, this, &AcousticTab::refresh);
    refresh();
}

void AcousticTab::apply()
{
    if (m_updating) return;
    AcousticOpts &a = m_p->acoustic();
    a.rt60 = m_rt60->isChecked();
    a.c80  = m_c80->isChecked();
    a.d50  = m_d50->isChecked();
    a.sti  = m_sti->isChecked();
    a.edt  = m_edt->isChecked();
    a.impulseResponse = m_irf->isChecked();
    a.auralization = m_aural->isChecked();
    a.sampleRate = m_sampleRate->currentText().toInt();
    static const char *dirs[] = { "omni", "cardioid", "speaker" };
    a.srcDirectivity = dirs[qBound(0, m_directivity->currentIndex(), 2)];
    a.srcSPL_dB = m_spl->value();
    a.micCount = m_micCount->value();
    m_p->touch();
}

void AcousticTab::refresh()
{
    m_updating = true;
    const AcousticOpts &a = m_p->acoustic();
    m_rt60->setChecked(a.rt60);
    m_c80->setChecked(a.c80);
    m_d50->setChecked(a.d50);
    m_sti->setChecked(a.sti);
    m_edt->setChecked(a.edt);
    m_irf->setChecked(a.impulseResponse);
    m_aural->setChecked(a.auralization);
    m_sampleRate->setCurrentText(QString::number(a.sampleRate));
    const int di = (a.srcDirectivity == "cardioid") ? 1
                 : (a.srcDirectivity == "speaker")  ? 2 : 0;
    m_directivity->setCurrentIndex(di);
    m_spl->setValue(a.srcSPL_dB);
    m_micCount->setValue(a.micCount);
    m_updating = false;
}
