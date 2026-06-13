// LogConsole.h — monospace log pane fed by Runner::logLine.
#pragma once
#include <QPlainTextEdit>

namespace ofd {

class LogConsole : public QPlainTextEdit {
    Q_OBJECT
public:
    explicit LogConsole(QWidget *parent = nullptr);

public slots:
    void appendLine(const QString &line);
};

} // namespace ofd
