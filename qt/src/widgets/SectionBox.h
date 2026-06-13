// SectionBox.h — equivalent of the React <Section title>…</Section>.
//
// A QGroupBox with the OpenFDTD-style flat border and floated title.
// Use it as a generic container inside any tab to mirror the mock layout.
#pragma once
#include <QGroupBox>
#include <QFormLayout>
#include <QVBoxLayout>

namespace ofd {

class SectionBox : public QGroupBox {
    Q_OBJECT
public:
    explicit SectionBox(const QString &title, QWidget *parent = nullptr);

    // Shortcut: form layout for label/widget rows (the typical case).
    QFormLayout *form();
    // Raw vertical layout for non-form content (tables, plots…).
    QVBoxLayout *vbox() const { return m_vbox; }

private:
    QFormLayout *m_form = nullptr;
    QVBoxLayout *m_vbox;
};

} // namespace ofd
