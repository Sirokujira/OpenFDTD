// DomainBar.h — horizontal tab switcher for the active physics domain.
// Equivalent of the React <DomainTabs> in the HTML mock.
#pragma once
#include <QWidget>
#include <QVector>
#include "core/Domain.h"

class QButtonGroup;
class QToolButton;

namespace ofd {

class DomainBar : public QWidget {
    Q_OBJECT
public:
    explicit DomainBar(QWidget *parent = nullptr);

    void setActiveDomain(Domain d);

signals:
    void domainSelected(ofd::Domain d);

private:
    struct Entry { Domain d; QToolButton *btn; };
    QVector<Entry> m_buttons;
    QButtonGroup  *m_group;
};

} // namespace ofd
