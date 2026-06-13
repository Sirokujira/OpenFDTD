// Tidy3dTab.h — tidy3d クラウド連携タブ (光ドメイン専用).
// 設計判断: tidy3d は物理ドメインではなく光FDTDのクラウドバックエンド —
// このタブは光ドメイン選択時のみ表示される (MainWindow::onDomainChanged)。
#pragma once
#include <QScrollArea>

class QLineEdit;
class QComboBox;
class QCheckBox;
class QLabel;

namespace ofd {

class Project;

class Tidy3dTab : public QScrollArea {
    Q_OBJECT
public:
    explicit Tidy3dTab(Project *project, QWidget *parent = nullptr);

private slots:
    void refresh();
    void exportScript();

private:
    void apply();

    Project   *m_p;
    bool       m_updating = false;

    QLineEdit *m_apiKey;
    QLineEdit *m_project;
    QComboBox *m_resolution;
    QCheckBox *m_autoPml;
    QLabel    *m_status;
};

} // namespace ofd
