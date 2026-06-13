// GeometryTab.h — geometry unit editor (物体形状タブ) + STL import.
// Maps 1:1 to the "geometry =" lines; unit order = ユニット番号 (later wins).
#pragma once
#include <QScrollArea>

class QTableWidget;
class QLabel;

namespace ofd {

class Project;
class UnitNav;

class GeometryTab : public QScrollArea {
    Q_OBJECT
public:
    explicit GeometryTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();
    void importStl();

private:
    void applyTable();

    Project      *m_p;
    bool          m_updating = false;
    QTableWidget *m_table;
    UnitNav      *m_nav;
    QLabel       *m_importInfo;
};

} // namespace ofd
