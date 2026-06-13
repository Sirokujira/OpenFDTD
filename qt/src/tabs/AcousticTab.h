// AcousticTab.h — 室内音響タブ (RT60/C80/D50/STI, マイクアレイ, 可聴化).
// Settings persist in the .ofdx sidecar.
#pragma once
#include <QScrollArea>

class QCheckBox;
class QComboBox;
class QSpinBox;
class QDoubleSpinBox;

namespace ofd {

class Project;

class AcousticTab : public QScrollArea {
    Q_OBJECT
public:
    explicit AcousticTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void apply();

    Project   *m_p;
    bool       m_updating = false;

    QCheckBox *m_rt60, *m_c80, *m_d50, *m_sti, *m_edt, *m_irf, *m_aural;
    QComboBox *m_sampleRate;
    QComboBox *m_directivity;
    QDoubleSpinBox *m_spl;
    QSpinBox  *m_micCount;
};

} // namespace ofd
