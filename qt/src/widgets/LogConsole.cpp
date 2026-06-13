// LogConsole.cpp
#include "LogConsole.h"

#include <QFontDatabase>
#include <QScrollBar>

using namespace ofd;

LogConsole::LogConsole(QWidget *parent)
    : QPlainTextEdit(parent)
{
    setObjectName("LogConsole");
    setReadOnly(true);
    setMaximumBlockCount(5000);
    setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    setLineWrapMode(QPlainTextEdit::NoWrap);
}

void LogConsole::appendLine(const QString &line)
{
    appendPlainText(line);
    verticalScrollBar()->setValue(verticalScrollBar()->maximum());
}
