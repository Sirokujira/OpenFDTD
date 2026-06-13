// DomainBar.cpp
#include "DomainBar.h"
#include "I18n.h"

#include <QHBoxLayout>
#include <QToolButton>
#include <QButtonGroup>

using namespace ofd;

DomainBar::DomainBar(QWidget *parent)
    : QWidget(parent)
    , m_group(new QButtonGroup(this))
{
    setObjectName("DomainBar");
    setFixedHeight(30);
    setAutoFillBackground(true);

    auto *h = new QHBoxLayout(this);
    h->setContentsMargins(8, 0, 8, 0);
    h->setSpacing(0);

    const struct { Domain d; const char *labelKey; const char *glyph; } items[] = {
        { Domain::EM,         "d_em",         "⚡" },
        { Domain::Optical,    "d_optical",    "✦" },
        { Domain::Acoustic,   "d_acoustic",   "♪" },
        { Domain::Underwater, "d_underwater", "≋" },
    };

    int id = 0;
    for (const auto &it : items) {
        auto *btn = new QToolButton(this);
        btn->setText(QStringLiteral("%1 %2")
                     .arg(QString::fromUtf8(it.glyph), I18n::tr(it.labelKey)));
        btn->setCheckable(true);
        btn->setAutoRaise(true);
        btn->setCursor(Qt::PointingHandCursor);
        btn->setStyleSheet(QStringLiteral(
            "QToolButton { padding: 5px 14px; border: none;"
            "  border-top: 2px solid transparent; }"
            "QToolButton:hover { background: palette(midlight); }"
            "QToolButton:checked { background: palette(base);"
            "  border-top-color: %1; font-weight: 600; }")
            .arg(accentColor(it.d)));
        h->addWidget(btn);
        m_group->addButton(btn, id++);
        m_buttons.push_back({ it.d, btn });
        connect(btn, &QToolButton::clicked, this, [this, d = it.d] {
            emit domainSelected(d);
        });
    }
    h->addStretch(1);

    if (!m_buttons.isEmpty()) m_buttons.first().btn->setChecked(true);
}

void DomainBar::setActiveDomain(Domain d) {
    for (auto &e : m_buttons) e.btn->setChecked(e.d == d);
}
