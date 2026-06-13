// GeneralTab.h — solver / ABC / PBC / frequency settings (全般タブ).
// Maps 1:1 to the .ofd keys: title, solver, abc, pbc, frequency1/2,
// timestep, pulsewidth, rfeed, plot3dgeom.
#pragma once
#include <QScrollArea>

class QLineEdit;
class QSpinBox;
class QComboBox;
class QCheckBox;
class QDoubleSpinBox;

namespace ofd {

class Project;

class GeneralTab : public QScrollArea {
    Q_OBJECT
public:
    explicit GeneralTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();   // model → widgets

private:
    Project *m_p;
    bool m_updating = false;

    QLineEdit *m_title;
    QSpinBox  *m_maxiter;
    QSpinBox  *m_nout;
    QLineEdit *m_converg;
    QComboBox *m_abc;
    QSpinBox  *m_pmlL;
    QDoubleSpinBox *m_pmlM;
    QLineEdit *m_pmlR0;
    QCheckBox *m_pbc[3];
    QLineEdit *m_f1min, *m_f1max; QSpinBox *m_f1div;
    QLineEdit *m_f2min, *m_f2max; QSpinBox *m_f2div;
    QLineEdit *m_dt, *m_tw, *m_rfeed;
    QCheckBox *m_plot3dgeom;
};

} // namespace ofd
