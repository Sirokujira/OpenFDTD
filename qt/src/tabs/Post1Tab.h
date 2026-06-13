// Post1Tab.h — ポスト処理制御(1): 周波数特性プロット類.
// Maps 1:1 to plotiter/plotfeed/plotpoint/plotsmith/plotzin/plotyin/
// plotref/plotspara/plotcoupling/matchingloss/freqdiv.
#pragma once
#include <QScrollArea>
#include <QVector>

class QCheckBox;
class QSpinBox;
class QLineEdit;

namespace ofd {

class Project;
struct FreqPlot;

class Post1Tab : public QScrollArea {
    Q_OBJECT
public:
    explicit Post1Tab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    struct FreqRow {
        QCheckBox *enabled;
        QCheckBox *userScale;
        QLineEdit *min, *max;
        QSpinBox  *div;
        FreqPlot  *target;
    };
    void addFreqRow(QWidget *parent, class SectionBox *s,
                    const QString &label, FreqPlot *target);
    void apply();

    Project   *m_p;
    bool       m_updating = false;
    QCheckBox *m_iter, *m_feed, *m_point, *m_smith, *m_matching;
    QSpinBox  *m_freqdiv;
    QVector<FreqRow> m_rows;
};

} // namespace ofd
