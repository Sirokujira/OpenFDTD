// EvViewer.h — handles the kernel's ev2/ev3 figure output via 3 strategies
// (see docs/ev-format.md):
//   (A) Html    : run post with -html → ev2d.htm / ev3d.htm → default browser
//   (B) Process : launch the standalone ev2d / ev3d viewer executables
//   (C) Native  : parse .ev2/.ev3 → Qt rendering (skeleton, see io/EvReader.h)
#pragma once
#include <QWidget>
#include "../core/Domain.h"

class QComboBox;
class QLabel;

namespace ofd {

enum class EvBackend { Html, Process, Native };

class EvViewer : public QWidget {
    Q_OBJECT
public:
    explicit EvViewer(QWidget *parent = nullptr);

    EvBackend backend() const;
    void setWorkdir(const QString &dir) { m_workdir = dir; }

public slots:
    void open2D();
    void open3D();

private:
    void open(bool threeD);

    QComboBox *m_backendBox;
    QLabel    *m_status;
    QString    m_workdir;
};

} // namespace ofd
