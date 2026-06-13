// RightDock.h — project tree + run log (right side of the main window).
#pragma once
#include <QWidget>

class QTreeWidget;

namespace ofd {

class Project;
class LogConsole;

class RightDock : public QWidget {
    Q_OBJECT
public:
    explicit RightDock(Project *project, QWidget *parent = nullptr);

    void appendLog(const QString &line);

private slots:
    void rebuildTree();

private:
    Project     *m_project;
    QTreeWidget *m_tree;
    LogConsole  *m_log;
};

} // namespace ofd
