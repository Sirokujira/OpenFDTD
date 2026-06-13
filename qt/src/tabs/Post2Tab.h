// Post2Tab.h — ポスト処理制御(2): 遠方界・近傍界プロット類.
// Maps 1:1 to plotfar0d/plotfar1d/plotfar2d/plotnear1d/plotnear2d
// and their style keys (far1dstyle, far1dcomponent, near2ddim, ...).
#pragma once
#include <QScrollArea>

class QCheckBox;
class QSpinBox;
class QLineEdit;
class QComboBox;
class QTableWidget;

namespace ofd {

class Project;

class Post2Tab : public QScrollArea {
    Q_OBJECT
public:
    explicit Post2Tab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void apply();
    void applyFar1dTable();
    void applyNear1dTable();
    void applyNear2dTable();

    Project   *m_p;
    bool       m_updating = false;

    // far0d
    QCheckBox *m_far0d;
    QLineEdit *m_far0dTheta, *m_far0dPhi;

    // far1d
    QTableWidget *m_far1d;
    QComboBox *m_far1dStyle;
    QCheckBox *m_far1dDb, *m_far1dNorm;
    QCheckBox *m_far1dCompE, *m_far1dCompTheta, *m_far1dCompPhi;

    // far2d
    QCheckBox *m_far2d;
    QSpinBox  *m_far2dTheta, *m_far2dPhi;
    QCheckBox *m_far2dDb;
    QLineEdit *m_far2dObj;

    // near1d / near2d
    QTableWidget *m_near1d;
    QCheckBox *m_near1dDb, *m_near1dNoinc;
    QTableWidget *m_near2d;
    QSpinBox  *m_near2dDim0, *m_near2dDim1;
    QCheckBox *m_near2dDb, *m_near2dContour, *m_near2dNoinc;
};

} // namespace ofd
