// SectionBox.cpp
#include "SectionBox.h"

using namespace ofd;

SectionBox::SectionBox(const QString &title, QWidget *parent)
    : QGroupBox(title, parent)
{
    setObjectName("SectionBox");

    m_vbox = new QVBoxLayout(this);
    m_vbox->setContentsMargins(10, 12, 10, 10);
    m_vbox->setSpacing(4);
}

QFormLayout *SectionBox::form()
{
    if (!m_form) {
        m_form = new QFormLayout();
        m_form->setRowWrapPolicy(QFormLayout::DontWrapRows);
        m_form->setFieldGrowthPolicy(QFormLayout::ExpandingFieldsGrow);
        m_form->setFormAlignment(Qt::AlignLeft | Qt::AlignTop);
        m_form->setLabelAlignment(Qt::AlignLeft);
        m_form->setHorizontalSpacing(8);
        m_form->setVerticalSpacing(4);
        m_vbox->addLayout(m_form);
    }
    return m_form;
}
