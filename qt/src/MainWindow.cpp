// MainWindow.cpp
#include "MainWindow.h"
#include "DomainBar.h"
#include "RightDock.h"
#include "I18n.h"

#include "core/Project.h"
#include "io/H5Writer.h"
#include "io/Tidy3dExporter.h"
#include "io/Touchstone.h"

#include "widgets/EvViewer.h"
#include "widgets/PlotPanel.h"
#include "widgets/Viewport3D.h"

#include "tabs/GeneralTab.h"
#include "tabs/MeshTab.h"
#include "tabs/MaterialTab.h"
#include "tabs/GeometryTab.h"
#include "tabs/SourceTab.h"
#include "tabs/Post1Tab.h"
#include "tabs/Post2Tab.h"
#include "tabs/OpticalTab.h"
#include "tabs/AcousticTab.h"
#include "tabs/UnderwaterTab.h"
#include "tabs/Tidy3dTab.h"

#include <QApplication>
#include <QComboBox>
#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QProgressBar>
#include <QRegularExpression>
#include <QSettings>
#include <QSpinBox>
#include <QSplitter>
#include <QStackedWidget>
#include <QStatusBar>
#include <QStyle>
#include <QTabWidget>
#include <QToolBar>
#include <QToolButton>
#include <QVBoxLayout>

using namespace ofd;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , m_project(new Project(this))
    , m_runner(new Runner(this))
{
    setObjectName("OpenFDTD_MainWindow");
    resize(1500, 940);
    setMinimumSize(1100, 700);

    buildMenu();
    buildToolbar();
    buildCentral();
    buildDocks();
    buildStatusBar();

    connect(m_project, &Project::domainChanged, this, &MainWindow::onDomainChanged);
    connect(m_project, &Project::changed, this, &MainWindow::onProjectChanged);

    connect(m_runner, &Runner::progress, this, &MainWindow::onRunnerProgress);
    connect(m_runner, &Runner::logLine, this, &MainWindow::onRunnerLog);
    connect(m_runner, &Runner::finished, this, &MainWindow::onRunnerFinished);

    onDomainChanged(m_project->activeDomain());
    onProjectChanged();
    updateWindowTitle();
}

MainWindow::~MainWindow() = default;

// ── Menus ───────────────────────────────────────────────────────────────────
void MainWindow::buildMenu()
{
    auto *mb = menuBar();
    auto *mFile = mb->addMenu(I18n::tr("m_file"));
    auto *mView = mb->addMenu(I18n::tr("m_view"));
    auto *mRun  = mb->addMenu(I18n::tr("m_run"));
    auto *mPost = mb->addMenu(I18n::tr("m_post"));
    auto *mTools= mb->addMenu(I18n::tr("m_tools"));
    auto *mHelp = mb->addMenu(I18n::tr("m_help"));

    mFile->addAction(I18n::tr("tb_new"), QKeySequence::New,
                     this, &MainWindow::newProject);
    mFile->addAction(I18n::tr("tb_open"), QKeySequence::Open,
                     this, [this] { openProject(); });
    mFile->addAction(I18n::tr("tb_save"), QKeySequence::Save,
                     this, &MainWindow::saveProject);
    mFile->addAction(I18n::tr("tb_saveas"), QKeySequence::SaveAs,
                     this, &MainWindow::saveProjectAs);
    mFile->addSeparator();
    mFile->addAction(I18n::tr("m_exit"), this, [] { qApp->quit(); });

    auto *langMenu = mView->addMenu(I18n::tr("m_lang"));
    for (const auto &[code, label] : std::initializer_list<
             std::pair<const char *, const char *>>{
             { "ja", "日本語" }, { "en", "English" }, { "both", "日英 / JA+EN" } }) {
        langMenu->addAction(QString::fromUtf8(label), this, [this, code = QString(code)] {
            QSettings().setValue("ui/language", code);
            QMessageBox::information(this, "OpenFDTD-X", I18n::tr("m_lang_restart"));
        });
    }

    mRun->addAction(I18n::tr("tb_calc"), QKeySequence(Qt::Key_F5),
                    this, &MainWindow::runSimulation);
    mRun->addAction(I18n::tr("tb_stop"), QKeySequence(Qt::Key_Escape),
                    this, [this] { m_runner->stop(); });

    mPost->addAction(I18n::tr("tb_post"), this, &MainWindow::runPostProcess);
    mPost->addAction(I18n::tr("tb_plot2d"), this, &MainWindow::show2DPlot);
    mPost->addAction(I18n::tr("tb_plot3d"), this, &MainWindow::show3DPlot);
    mPost->addSeparator();
    mPost->addAction(I18n::tr("pp_export_csv"), this, [this] {
        const QString p = QFileDialog::getSaveFileName(
            this, I18n::tr("pp_export_csv"), "convergence.csv", "CSV (*.csv)");
        if (!p.isEmpty()) m_plotPanel->exportCsv(p);
    });
    mPost->addAction(I18n::tr("pp_export_h5"), this, &MainWindow::exportHdf5);
    mPost->addAction(I18n::tr("pp_export_s2p"), this, &MainWindow::exportTouchstone);

    mTools->addAction(I18n::tr("tb_cloud"), this, &MainWindow::exportTidy3d);

    mHelp->addAction("About OpenFDTD-X…", this, [this] {
        QMessageBox::about(this, "OpenFDTD-X",
            QStringLiteral("<b>OpenFDTD-X</b><br>%1<br><br>"
                           "OpenFDTD / OpenRCWA / OpenBPM GUI front-end")
                .arg(I18n::tr("app_subtitle")));
    });
}

// ── Toolbar ─────────────────────────────────────────────────────────────────
void MainWindow::buildToolbar()
{
    auto *tb = addToolBar("Main");
    tb->setObjectName("MainToolBar");
    tb->setMovable(false);
    tb->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
    tb->setIconSize({ 16, 16 });
    auto icon = [this](QStyle::StandardPixmap sp) { return style()->standardIcon(sp); };

    tb->addAction(icon(QStyle::SP_FileIcon), I18n::tr("tb_new"),
                  this, &MainWindow::newProject);
    tb->addAction(icon(QStyle::SP_DialogOpenButton), I18n::tr("tb_open"),
                  this, [this] { openProject(); });
    tb->addAction(icon(QStyle::SP_DialogSaveButton), I18n::tr("tb_save"),
                  this, &MainWindow::saveProject);
    tb->addSeparator();

    auto *runAct = tb->addAction(icon(QStyle::SP_MediaPlay), I18n::tr("tb_calc"),
                                 this, &MainWindow::runSimulation);
    if (auto *btn = qobject_cast<QToolButton *>(tb->widgetForAction(runAct)))
        btn->setObjectName("primaryAction");
    tb->addAction(icon(QStyle::SP_MediaSeekForward), I18n::tr("tb_post"),
                  this, &MainWindow::runPostProcess);
    tb->addAction(icon(QStyle::SP_FileDialogContentsView), I18n::tr("tb_plot2d"),
                  this, &MainWindow::show2DPlot);
    tb->addAction(icon(QStyle::SP_FileDialogListView), I18n::tr("tb_plot3d"),
                  this, &MainWindow::show3DPlot);
    tb->addSeparator();

    m_cloudAction = tb->addAction(icon(QStyle::SP_ArrowUp), I18n::tr("tb_cloud"),
                                  this, &MainWindow::exportTidy3d);

    auto *spacer = new QWidget(tb);
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    tb->addWidget(spacer);

    tb->addWidget(new QLabel(I18n::tr("run_engine") + ": ", tb));
    m_engineBox = new QComboBox(tb);
    m_engineBox->addItem(I18n::tr("run_cpu"));      // Engine::CPU
    m_engineBox->addItem(I18n::tr("run_cpu_mpi"));  // Engine::CPU_MPI
    m_engineBox->addItem(I18n::tr("run_gpu"));      // Engine::GPU
    m_engineBox->addItem(I18n::tr("run_gpu_mpi"));  // Engine::GPU_MPI
    tb->addWidget(m_engineBox);

    m_modeBox = new QComboBox(tb);
    m_modeBox->addItem(I18n::tr("run_both"));        // RunMode::Both
    m_modeBox->addItem(I18n::tr("run_solver_only")); // RunMode::Solver
    m_modeBox->addItem(I18n::tr("run_post_only"));   // RunMode::Post
    tb->addWidget(m_modeBox);

    tb->addWidget(new QLabel(" " + I18n::tr("run_threads") + ": ", tb));
    m_threadsBox = new QSpinBox(tb);
    m_threadsBox->setRange(1, 256);
    m_threadsBox->setValue(QSettings().value("run/threads", 4).toInt());
    tb->addWidget(m_threadsBox);
}

// ── Central widget ──────────────────────────────────────────────────────────
void MainWindow::buildCentral()
{
    auto *central = new QWidget(this);
    auto *v = new QVBoxLayout(central);
    v->setContentsMargins(0, 0, 0, 0);
    v->setSpacing(0);

    m_domainBar = new DomainBar(central);
    v->addWidget(m_domainBar);
    connect(m_domainBar, &DomainBar::domainSelected,
            m_project, &Project::setActiveDomain);

    auto *split = new QSplitter(Qt::Horizontal, central);
    split->setChildrenCollapsible(false);

    // left: tabs
    m_leftTabs = new QTabWidget(split);
    m_leftTabs->setDocumentMode(true);
    m_leftTabs->setUsesScrollButtons(true);

    m_tabGeneral   = new GeneralTab(m_project);
    m_tabMesh      = new MeshTab(m_project);
    m_tabMaterial  = new MaterialTab(m_project);
    m_tabGeometry  = new GeometryTab(m_project);
    m_tabSource    = new SourceTab(m_project);
    m_tabPost1     = new Post1Tab(m_project);
    m_tabPost2     = new Post2Tab(m_project);
    m_tabOptical   = new OpticalTab(m_project);
    m_tabAcoustic  = new AcousticTab(m_project);
    m_tabUnderwater= new UnderwaterTab(m_project);
    m_tabTidy3d    = new Tidy3dTab(m_project);

    m_leftTabs->addTab(m_tabGeneral, I18n::tr("t_general"));
    m_leftTabs->addTab(m_tabMesh, I18n::tr("t_mesh"));
    m_leftTabs->addTab(m_tabMaterial, I18n::tr("t_material"));
    m_leftTabs->addTab(m_tabGeometry, I18n::tr("t_geometry"));
    m_leftTabs->addTab(m_tabSource, I18n::tr("t_source"));
    m_leftTabs->addTab(m_tabPost1, I18n::tr("t_post1"));
    m_leftTabs->addTab(m_tabPost2, I18n::tr("t_post2"));
    // domain tabs are added/removed in onDomainChanged()

    // center: viewport / plot stack + ev viewer bar
    auto *centerWrap = new QWidget(split);
    auto *cv = new QVBoxLayout(centerWrap);
    cv->setContentsMargins(0, 0, 0, 0);
    cv->setSpacing(2);

    m_centerStack = new QStackedWidget(centerWrap);
    m_viewport = new Viewport3D(m_project, m_centerStack);
    m_plotPanel = new PlotPanel(m_project, m_centerStack);
    m_centerStack->addWidget(m_viewport);
    m_centerStack->addWidget(m_plotPanel);
    cv->addWidget(m_centerStack, 1);

    m_evViewer = new EvViewer(centerWrap);
    cv->addWidget(m_evViewer);

    split->addWidget(m_leftTabs);
    split->addWidget(centerWrap);
    split->setStretchFactor(0, 0);
    split->setStretchFactor(1, 1);
    split->setSizes({ 480, 720 });

    v->addWidget(split, 1);
    setCentralWidget(central);
}

// ── Docks ───────────────────────────────────────────────────────────────────
void MainWindow::buildDocks()
{
    m_rightDock = new RightDock(m_project, this);
    auto *dock = new QDockWidget(I18n::tr("rd_project"), this);
    dock->setObjectName("rightDock");
    dock->setAllowedAreas(Qt::RightDockWidgetArea | Qt::LeftDockWidgetArea);
    dock->setFeatures(QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable);
    dock->setWidget(m_rightDock);
    dock->setMinimumWidth(260);
    addDockWidget(Qt::RightDockWidgetArea, dock);
}

// ── Statusbar ───────────────────────────────────────────────────────────────
void MainWindow::buildStatusBar()
{
    auto *sb = statusBar();
    m_sbState = new QLabel("● " + I18n::tr("sb_ready"));
    m_sbCells = new QLabel;
    m_sbMem = new QLabel;
    m_sbDt = new QLabel;
    m_sbStep = new QLabel;
    m_sbProgress = new QProgressBar;
    m_sbProgress->setRange(0, 100);
    m_sbProgress->setFixedWidth(140);
    m_sbProgress->setVisible(false);

    sb->addWidget(m_sbState);
    sb->addPermanentWidget(m_sbCells);
    sb->addPermanentWidget(m_sbMem);
    sb->addPermanentWidget(m_sbDt);
    sb->addPermanentWidget(m_sbStep);
    sb->addPermanentWidget(m_sbProgress);
    sb->addPermanentWidget(new QLabel("Qt " QT_VERSION_STR));
}

// ── Domain switching ────────────────────────────────────────────────────────
void MainWindow::onDomainChanged(Domain d)
{
    auto removeTab = [this](QWidget *w) {
        const int idx = m_leftTabs->indexOf(w);
        if (idx >= 0) m_leftTabs->removeTab(idx);
    };
    removeTab(m_tabOptical);
    removeTab(m_tabAcoustic);
    removeTab(m_tabUnderwater);
    removeTab(m_tabTidy3d);

    switch (d) {
    case Domain::Optical:
        m_leftTabs->addTab(m_tabOptical, I18n::tr("t_optical"));
        // tidy3d は光FDTD専用のクラウドバックエンド (独立ドメインではない)
        m_leftTabs->addTab(m_tabTidy3d, I18n::tr("t_tidy3d"));
        break;
    case Domain::Acoustic:
        m_leftTabs->addTab(m_tabAcoustic, I18n::tr("t_acoustic"));
        break;
    case Domain::Underwater:
        m_leftTabs->addTab(m_tabUnderwater, I18n::tr("t_underwater"));
        break;
    case Domain::EM:
        break;
    }

    // cloud submission is optical-only
    m_cloudAction->setEnabled(d == Domain::Optical);
    m_cloudAction->setText(d == Domain::Optical
        ? I18n::tr("tb_cloud") : I18n::tr("tb_cloud_optical_only"));

    m_domainBar->setActiveDomain(d);
    m_viewport->setDomain(d);
    m_plotPanel->setDomain(d);

    setStyleSheet(QStringLiteral(
        "QToolButton#primaryAction { background: %1; color: white;"
        " font-weight: 600; border-radius: 3px; padding: 3px 10px; }")
        .arg(accentColor(d)));
}

void MainWindow::onProjectChanged()
{
    m_sbCells->setText(QStringLiteral("cells: %L1").arg(m_project->totalCells()));
    m_sbMem->setText(QStringLiteral("mem: %1 MB")
        .arg(m_project->estimatedMemoryMB(), 0, 'f', 1));
    const double dt = m_project->courantDt();
    m_sbDt->setText(dt > 0
        ? QStringLiteral("Δt(CFL): %1 s").arg(QString::number(dt, 'g', 3))
        : QStringLiteral("Δt: -"));
}

// ── File actions ────────────────────────────────────────────────────────────
void MainWindow::updateWindowTitle()
{
    const QString file = m_project->filePath().isEmpty()
        ? I18n::tr("untitled")
        : QFileInfo(m_project->filePath()).fileName();
    setWindowTitle(QStringLiteral("OpenFDTD-X — %1").arg(file));
}

void MainWindow::newProject()
{
    m_project->clear();
    emit m_project->loaded();
    emit m_project->changed();
    updateWindowTitle();
}

void MainWindow::openProject(const QString &path)
{
    QString p = path;
    if (p.isEmpty()) {
        p = QFileDialog::getOpenFileName(this, I18n::tr("tb_open"), {},
            "OpenFDTD (*.ofd);;All files (*)");
        if (p.isEmpty()) return;
    }
    QString err;
    if (!m_project->load(p, &err)) {
        QMessageBox::warning(this, I18n::tr("tb_open"), err);
        return;
    }
    m_evViewer->setWorkdir(QFileInfo(p).path());
    updateWindowTitle();
}

void MainWindow::saveProject()
{
    if (m_project->filePath().isEmpty()) {
        saveProjectAs();
        return;
    }
    QString err;
    if (!m_project->save(m_project->filePath(), &err))
        QMessageBox::warning(this, I18n::tr("tb_save"), err);
}

void MainWindow::saveProjectAs()
{
    const QString p = QFileDialog::getSaveFileName(this, I18n::tr("tb_saveas"),
        m_project->general().title.isEmpty() ? "project.ofd"
            : m_project->general().title + ".ofd",
        "OpenFDTD (*.ofd)");
    if (p.isEmpty()) return;
    QString err;
    if (!m_project->save(p, &err)) {
        QMessageBox::warning(this, I18n::tr("tb_save"), err);
        return;
    }
    m_evViewer->setWorkdir(QFileInfo(p).path());
    updateWindowTitle();
}

// ── Run ─────────────────────────────────────────────────────────────────────
RunConfig MainWindow::currentRunConfig() const
{
    RunConfig cfg;
    cfg.engine = Engine(m_engineBox->currentIndex());
    cfg.mode = (m_modeBox->currentIndex() == 1) ? RunMode::Solver
             : (m_modeBox->currentIndex() == 2) ? RunMode::Post
                                                : RunMode::Both;
    cfg.threads = m_threadsBox->value();
    QSettings().setValue("run/threads", cfg.threads);

    // 光ドメイン: RCWA/BPM は姉妹カーネル (orcwa / obpm) を使う
    if (m_project->activeDomain() == Domain::Optical) {
        switch (m_project->optical().solver) {
            case OpticalSolver::RCWA: cfg.kernel = Kernel::RCWA; break;
            case OpticalSolver::BPM:  cfg.kernel = Kernel::BPM;  break;
            default:                  cfg.kernel = Kernel::FDTD; break;
        }
    }
    return cfg;
}

void MainWindow::runSimulation()
{
    if (m_runner->isRunning()) {
        m_runner->stop();
        return;
    }
    m_plotPanel->clearConvergence();
    m_sbProgress->setVisible(true);
    m_sbProgress->setValue(0);
    m_sbState->setText("● " + I18n::tr("sb_running"));
    m_runner->start(m_project, currentRunConfig());
    m_evViewer->setWorkdir(m_runner->workingDir());
}

void MainWindow::runPostProcess()
{
    if (m_runner->isRunning()) return;
    RunConfig cfg = currentRunConfig();
    cfg.mode = RunMode::Post;
    m_sbState->setText("● " + I18n::tr("sb_running"));
    m_runner->start(m_project, cfg);
    m_evViewer->setWorkdir(m_runner->workingDir());
}

void MainWindow::setDomain(Domain d)
{
    m_project->setActiveDomain(d);
}

void MainWindow::show2DPlot()
{
    m_centerStack->setCurrentWidget(m_plotPanel);
}

void MainWindow::show3DPlot()
{
    m_centerStack->setCurrentWidget(m_viewport);
}

void MainWindow::exportHdf5()
{
    const QString p = QFileDialog::getSaveFileName(this, I18n::tr("pp_export_h5"),
        (m_project->general().title.isEmpty() ? "project"
            : m_project->general().title) + ".h5",
        "HDF5 (*.h5)");
    if (p.isEmpty()) return;
    QString err;
    QVector<int> steps; QVector<double> e, h;   // PlotPanel keeps the live copy
    if (!H5Writer::write(p, *m_project, steps, e, h, &err))
        QMessageBox::warning(this, I18n::tr("pp_export_h5"), err);
}

void MainWindow::exportTouchstone()
{
    // S-parameters come from the post processor's frequency output; the
    // Touchstone writer is wired here so post-run data can be exported once
    // a parser for the kernel's spara output lands.
    QMessageBox::information(this, I18n::tr("pp_export_s2p"),
        "Run the post processor with plotspara enabled, then export from "
        "the produced frequency data (io/Touchstone).");
}

void MainWindow::exportTidy3d()
{
    if (m_project->activeDomain() != Domain::Optical) return;
    const QString p = QFileDialog::getSaveFileName(this, I18n::tr("t3_export"),
        m_project->tidy3d().projectName + ".py", "Python (*.py)");
    if (p.isEmpty()) return;
    QString err;
    if (!Tidy3dExporter::exportTo(p, *m_project, &err))
        QMessageBox::warning(this, I18n::tr("t3_export"), err);
}

// ── Runner feedback ─────────────────────────────────────────────────────────
void MainWindow::onRunnerProgress(int step, int total)
{
    if (total > 0)
        m_sbProgress->setValue(int(100.0 * step / total));
    m_sbStep->setText(QStringLiteral("step: %1 / %2").arg(step).arg(total));
}

void MainWindow::onRunnerLog(const QString &line)
{
    m_rightDock->appendLog(line);
    // feed convergence points to the plot ("%7d %f %f")
    static const QRegularExpression stepRe(
        "^\\s*(\\d+)\\s+([-+0-9.eE]+)\\s+([-+0-9.eE]+)\\s*$");
    const auto m = stepRe.match(line);
    if (m.hasMatch())
        m_plotPanel->addConvergencePoint(m.captured(1).toInt(),
                                         m.captured(2).toDouble(),
                                         m.captured(3).toDouble());
}

void MainWindow::onRunnerFinished(bool ok)
{
    m_sbProgress->setVisible(false);
    m_sbState->setText("● " + (ok ? I18n::tr("sb_done") : I18n::tr("sb_failed")));
}
