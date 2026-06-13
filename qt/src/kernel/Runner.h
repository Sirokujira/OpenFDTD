// Runner.h — drives the OpenFDTD-family CLI kernels as subprocesses.
//
// 本家のパイプライン:
//   solver:  ofd [-n <thread>] [-out <outfile>] <datafile>   → ofd.out + ofd.log
//   post:    ofd_post [-n <thread>] [-html] <datafile>       → ev.ev2/ev.ev3 or *.htm
//   MPI:     mpiexec -n <process> ofd_mpi [-p x y z] [-n <thread>] <datafile>
//   CUDA:    ofd_cuda / ofd_cuda_mpi
//
// 姉妹ソルバー (光ドメインの RCWA / BPM):
//   OpenRCWA: orcwa / orcwa_mpi / orcwa_cuda + orcwa_post
//   OpenBPM:  obpm  / obpm_mpi  / obpm_cuda  + obpm_post
//
// バイナリ探索順: cfg.binaryDir → $OPENFDTD_HOME (orcwa は $OPENRCWA_HOME,
// obpm は $OPENBPM_HOME) → アプリ実行ディレクトリ → PATH。
#pragma once
#include <QObject>
#include <QProcess>
#include <QString>

namespace ofd {

class Project;

enum class Engine  { CPU, CPU_MPI, GPU, GPU_MPI };
enum class Kernel  { FDTD, RCWA, BPM };          // solver family
enum class RunMode { Solver, Post, Both };

struct RunConfig {
    Engine  engine    = Engine::CPU;
    Kernel  kernel    = Kernel::FDTD;
    RunMode mode      = RunMode::Both;
    int     threads   = 4;       // OpenMP threads
    int     processes = 2;       // MPI processes
    bool    evHtml    = false;   // post: -html → ev2d.htm / ev3d.htm
    QString workingDir;          // where .ofd lives — outputs land here too
    QString binaryDir;           // explicit kernel location (optional)
};

class Runner : public QObject {
    Q_OBJECT
public:
    explicit Runner(QObject *parent = nullptr);
    ~Runner() override;

    bool isRunning() const;
    QString workingDir() const { return m_cfg.workingDir; }

    // Serialize project → .ofd/.ofdx in the working dir, then launch.
    void start(Project *project, const RunConfig &cfg = {});
    void stop();

    static QString solverBinary(const RunConfig &cfg);
    static QString postBinary(const RunConfig &cfg);
    static QString resolveBinary(const RunConfig &cfg, const QString &name);

signals:
    void started();
    void logLine(const QString &line);     // forwarded stdout/stderr
    void progress(int step, int total);    // parsed from "%7d %f %f" lines
    void finished(bool ok);

private slots:
    void onReadyRead();
    void onFinished(int exitCode, QProcess::ExitStatus status);

private:
    void launch(bool solverPhase);

    QProcess  *m_proc = nullptr;
    RunConfig  m_cfg;
    QString    m_ofdPath;
    bool       m_postPending = false;   // Both mode: post runs after solver
    int        m_totalSteps = 1000;
};

} // namespace ofd
