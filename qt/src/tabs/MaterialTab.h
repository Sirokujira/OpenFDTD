// MaterialTab.h — materials & lumped elements (物性値・集中定数タブ).
// Maps 1:1 to the "material =" and "load =" lines (+ rfeed on GeneralTab).
#pragma once
#include <QScrollArea>

class QTableWidget;

namespace ofd {

class Project;

class MaterialTab : public QScrollArea {
    Q_OBJECT
public:
    explicit MaterialTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();

private:
    void applyMaterials();
    void applyLoads();

    Project      *m_p;
    bool          m_updating = false;
    QTableWidget *m_mats;
    QTableWidget *m_loads;
};

} // namespace ofd
