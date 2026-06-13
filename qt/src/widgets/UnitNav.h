// UnitNav.h — "ユニット番号" navigator (◀ n / total ▶), as in the 本家 GUI.
#pragma once
#include <QWidget>

class QToolButton;
class QLabel;

namespace ofd {

class UnitNav : public QWidget {
    Q_OBJECT
public:
    explicit UnitNav(QWidget *parent = nullptr);

    void setRange(int count);          // 1..count, empty when count == 0
    void setCurrent(int index);        // 0-based
    int  current() const { return m_index; }

signals:
    void currentChanged(int index);    // 0-based

private:
    void refresh();

    QToolButton *m_prev;
    QToolButton *m_next;
    QLabel      *m_label;
    int          m_index = -1;
    int          m_count = 0;
};

} // namespace ofd
