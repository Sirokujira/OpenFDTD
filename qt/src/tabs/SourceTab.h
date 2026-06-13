// SourceTab.h — feeds, plane wave, observation points (波源・観測点タブ).
// Maps 1:1 to the "feed =", "planewave =", "point =" lines.
#pragma once
#include <QScrollArea>

class QTableWidget;
class QCheckBox;
class QLineEdit;
class QComboBox;
class QLabel;

namespace ofd {

class Project;

class SourceTab : public QScrollArea {
    Q_OBJECT
public:
    explicit SourceTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void applyFeeds();
    void applyPoints();
    void updateExclusiveWarning();

    Project      *m_p;
    bool          m_updating = false;
    QTableWidget *m_feeds;
    QCheckBox    *m_pwEnable;
    QLineEdit    *m_pwTheta, *m_pwPhi;
    QComboBox    *m_pwPol;
    QTableWidget *m_points;
    QLabel       *m_warning;
};

} // namespace ofd
