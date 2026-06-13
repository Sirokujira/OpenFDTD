// UnderwaterTab.h — 水中音響タブ (SSP, SOFAR, 海底底質, ソナー).
// Settings persist in the .ofdx sidecar.
#pragma once
#include <QScrollArea>

class QDoubleSpinBox;
class QCheckBox;
class QComboBox;
class QTableWidget;

namespace ofd {

class Project;

class UnderwaterTab : public QScrollArea {
    Q_OBJECT
public:
    explicit UnderwaterTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void apply();
    void applySsp();

    Project        *m_p;
    bool            m_updating = false;

    QDoubleSpinBox *m_temp, *m_salinity;
    QTableWidget   *m_ssp;
    QCheckBox      *m_sofar;
    QComboBox      *m_bottomType;
    QDoubleSpinBox *m_bottomC, *m_bottomRho;
    QDoubleSpinBox *m_sonarFreq, *m_sonarSL, *m_rangeMax;
};

} // namespace ofd
