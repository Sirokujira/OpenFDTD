// MeshTab.h — non-uniform mesh editor (メッシュタブ).
// Three coord/division tables map 1:1 to the xmesh/ymesh/zmesh lines.
#pragma once
#include <QScrollArea>

class QTableWidget;
class QLabel;

namespace ofd {

class Project;

class MeshTab : public QScrollArea {
    Q_OBJECT
public:
    explicit MeshTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void applyAxis(int axis);
    void refreshAxisInfo(int axis);

    Project      *m_p;
    bool          m_updating = false;
    QTableWidget *m_table[3];
    QLabel       *m_info[3];
    QLabel       *m_total;
};

} // namespace ofd
