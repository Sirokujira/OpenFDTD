// EvViewer.cpp
#include "EvViewer.h"
#include "../I18n.h"

#include <QComboBox>
#include <QDesktopServices>
#include <QDir>
#include <QFileInfo>
#include <QHBoxLayout>
#include <QLabel>
#include <QProcess>
#include <QPushButton>
#include <QUrl>

using namespace ofd;

EvViewer::EvViewer(QWidget *parent)
    : QWidget(parent)
{
    auto *h = new QHBoxLayout(this);
    h->setContentsMargins(0, 0, 0, 0);
    h->setSpacing(6);

    h->addWidget(new QLabel(I18n::tr("ev_backend") + ":", this));

    m_backendBox = new QComboBox(this);
    m_backendBox->addItem(I18n::tr("ev_html"));      // Html
    m_backendBox->addItem(I18n::tr("ev_process"));   // Process
    m_backendBox->addItem(I18n::tr("ev_native"));    // Native (TODO)
    h->addWidget(m_backendBox, 1);

    auto *b2 = new QPushButton(I18n::tr("ev_open2d"), this);
    auto *b3 = new QPushButton(I18n::tr("ev_open3d"), this);
    h->addWidget(b2);
    h->addWidget(b3);

    m_status = new QLabel(this);
    h->addWidget(m_status, 1);

    connect(b2, &QPushButton::clicked, this, &EvViewer::open2D);
    connect(b3, &QPushButton::clicked, this, &EvViewer::open3D);
}

EvBackend EvViewer::backend() const
{
    return EvBackend(m_backendBox->currentIndex());
}

void EvViewer::open2D() { open(false); }
void EvViewer::open3D() { open(true); }

void EvViewer::open(bool threeD)
{
    const QDir dir(m_workdir.isEmpty() ? QDir::currentPath() : m_workdir);
    m_status->clear();

    switch (backend()) {
    case EvBackend::Html: {
        const QString htm = dir.filePath(threeD ? "ev3d.htm" : "ev2d.htm");
        if (!QFileInfo::exists(htm)) {
            m_status->setText(I18n::tr("ev_nofile"));
            return;
        }
        QDesktopServices::openUrl(QUrl::fromLocalFile(htm));
        break;
    }
    case EvBackend::Process: {
        const QString ev = dir.filePath(threeD ? "ev.ev3" : "ev.ev2");
        if (!QFileInfo::exists(ev)) {
            m_status->setText(I18n::tr("ev_nofile"));
            return;
        }
        const QString exe = threeD ? "ev3d" : "ev2d";
        if (!QProcess::startDetached(exe, { ev }, dir.path()))
            m_status->setText("error: cannot launch " + exe);
        break;
    }
    case EvBackend::Native:
        // EvReader (io/EvReader.h) はパーサ骨格のみ — docs/ev-format.md 参照
        m_status->setText(I18n::tr("ev_native"));
        break;
    }
}
