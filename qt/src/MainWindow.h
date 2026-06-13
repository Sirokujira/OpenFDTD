// MainWindow.h — main application shell.
//
// Mirrors the HTML mock's structure:
//   ┌──────────────────────────────────────────────────┐
//   │ menubar                                          │
//   │ toolbar [New][Open][Save] [Run][Post][2D][3D]    │
//   │         [Cloud]  …  engine / threads             │
//   │ domain tabs [電磁 | 光 | 室内音響 | 水中]          │
//   ├─────────┬───────────────────────────────┬────────┤
//   │ left    │ center                        │ right  │
//   │ QTabW.  │ Viewport3D / PlotPanel        │ tree + │
//   │ 7+域タブ │ + EvViewer bar                │ log    │
//   └─────────┴───────────────────────────────┴────────┘
//   statusbar: state | cells | mem | Δt | step | progress
#pragma once
#include <QMainWindow>
#include "core/Domain.h"
#include "kernel/Runner.h"

class QTabWidget;
class QStackedWidget;
class QLabel;
class QProgressBar;
class QComboBox;
class QSpinBox;
class QAction;

namespace ofd {

class Project;
class DomainBar;
class RightDock;
class Viewport3D;
class PlotPanel;
class EvViewer;
class GeneralTab;
class MeshTab;
class MaterialTab;
class GeometryTab;
class SourceTab;
class Post1Tab;
class Post2Tab;
class OpticalTab;
class AcousticTab;
class UnderwaterTab;
class Tidy3dTab;

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

public slots:
    void newProject();
    void openProject(const QString &path = {});
    void saveProject();
    void saveProjectAs();
    void runSimulation();
    void runPostProcess();
    void show2DPlot();
    void show3DPlot();
    void exportHdf5();
    void exportTouchstone();
    void exportTidy3d();
    void setDomain(ofd::Domain d);

private slots:
    void onDomainChanged(ofd::Domain d);
    void onProjectChanged();
    void onRunnerProgress(int step, int total);
    void onRunnerLog(const QString &line);
    void onRunnerFinished(bool ok);

private:
    void buildMenu();
    void buildToolbar();
    void buildCentral();
    void buildDocks();
    void buildStatusBar();
    RunConfig currentRunConfig() const;
    void updateWindowTitle();

    Project *m_project = nullptr;
    Runner  *m_runner  = nullptr;

    DomainBar      *m_domainBar = nullptr;
    QTabWidget     *m_leftTabs = nullptr;
    RightDock      *m_rightDock = nullptr;
    Viewport3D     *m_viewport = nullptr;
    PlotPanel      *m_plotPanel = nullptr;
    EvViewer       *m_evViewer = nullptr;
    QStackedWidget *m_centerStack = nullptr;

    GeneralTab    *m_tabGeneral = nullptr;
    MeshTab       *m_tabMesh = nullptr;
    MaterialTab   *m_tabMaterial = nullptr;
    GeometryTab   *m_tabGeometry = nullptr;
    SourceTab     *m_tabSource = nullptr;
    Post1Tab      *m_tabPost1 = nullptr;
    Post2Tab      *m_tabPost2 = nullptr;
    OpticalTab    *m_tabOptical = nullptr;
    AcousticTab   *m_tabAcoustic = nullptr;
    UnderwaterTab *m_tabUnderwater = nullptr;
    Tidy3dTab     *m_tabTidy3d = nullptr;

    QComboBox *m_engineBox = nullptr;
    QComboBox *m_modeBox = nullptr;
    QSpinBox  *m_threadsBox = nullptr;
    QAction   *m_cloudAction = nullptr;

    QLabel       *m_sbState = nullptr;
    QLabel       *m_sbCells = nullptr;
    QLabel       *m_sbMem = nullptr;
    QLabel       *m_sbDt = nullptr;
    QLabel       *m_sbStep = nullptr;
    QProgressBar *m_sbProgress = nullptr;
};

} // namespace ofd
