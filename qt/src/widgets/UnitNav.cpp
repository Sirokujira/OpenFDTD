// UnitNav.cpp
#include "UnitNav.h"

#include <QHBoxLayout>
#include <QToolButton>
#include <QLabel>

using namespace ofd;

UnitNav::UnitNav(QWidget *parent)
    : QWidget(parent)
{
    auto *h = new QHBoxLayout(this);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(4);

    m_prev = new QToolButton(this);
    m_prev->setArrowType(Qt::LeftArrow);
    m_next = new QToolButton(this);
    m_next->setArrowType(Qt::RightArrow);
    m_label = new QLabel("- / -", this);
    m_label->setMinimumWidth(56);
    m_label->setAlignment(Qt::AlignCenter);

    h->addWidget(m_prev);
    h->addWidget(m_label);
    h->addWidget(m_next);

    connect(m_prev, &QToolButton::clicked, this, [this] {
        if (m_index > 0) { setCurrent(m_index - 1); emit currentChanged(m_index); }
    });
    connect(m_next, &QToolButton::clicked, this, [this] {
        if (m_index + 1 < m_count) { setCurrent(m_index + 1); emit currentChanged(m_index); }
    });
    refresh();
}

void UnitNav::setRange(int count)
{
    m_count = count;
    if (m_index >= count) m_index = count - 1;
    if (m_index < 0 && count > 0) m_index = 0;
    refresh();
}

void UnitNav::setCurrent(int index)
{
    m_index = qBound(-1, index, m_count - 1);
    refresh();
}

void UnitNav::refresh()
{
    m_label->setText(m_count > 0
        ? QStringLiteral("%1 / %2").arg(m_index + 1).arg(m_count)
        : QStringLiteral("- / -"));
    m_prev->setEnabled(m_index > 0);
    m_next->setEnabled(m_index + 1 < m_count);
}
