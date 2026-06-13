// OpticalTab.h — 光解析タブ: FDTD/RCWA/BPM/FMM ソルバー切替 + 光モード設定.
// Settings persist in the .ofdx sidecar; RCWA/BPM run the OpenRCWA (orcwa) /
// OpenBPM (obpm) sister kernels through Runner.
#pragma once
#include <QScrollArea>

class QComboBox;
class QLineEdit;
class QSpinBox;
class QCheckBox;
class QStackedWidget;

namespace ofd {

class Project;

class OpticalTab : public QScrollArea {
    Q_OBJECT
public:
    explicit OpticalTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void apply();

    Project   *m_p;
    bool       m_updating = false;

    QComboBox *m_solver;
    QStackedWidget *m_solverStack;
    QComboBox *m_mode;
    QLineEdit *m_lambdaMin, *m_lambdaMax;
    QSpinBox  *m_lambdaDiv;

    // RCWA
    QSpinBox  *m_rcwaNx, *m_rcwaNy, *m_rcwaLayers;
    QLineEdit *m_rcwaPx, *m_rcwaPy;
    // BPM
    QComboBox *m_bpmAlgo, *m_bpmInput;
    QLineEdit *m_bpmDz, *m_bpmN0;
    // FMM
    QSpinBox  *m_fmmHarmonics;
    QCheckBox *m_fmmLi;
    // BPF / Ring
    QLineEdit *m_bpfMin, *m_bpfMax, *m_bpfQ;
    QLineEdit *m_ringR, *m_ringGap;
};

} // namespace ofd
