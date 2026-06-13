// Runner.cpp
#include "Runner.h"
#include "../core/Project.h"
#include "../io/OfdIO.h"

#include <QCoreApplication>
#include <QDir>
#include <QFileInfo>
#include <QRegularExpression>
#include <QStandardPaths>

using namespace ofd;

Runner::Runner(QObject *parent) : QObject(parent) {}
Runner::~Runner() { if (m_proc) m_proc->kill(); }

bool Runner::isRunning() const {
    return m_proc && m_proc->state() != QProcess::NotRunning;
}

static QString kernelPrefix(Kernel k) {
    switch (k) {
        case Kernel::FDTD: return "ofd";
        case Kernel::RCWA: return "orcwa";
        case Kernel::BPM:  return "obpm";
    }
    return "ofd";
}

QString Runner::solverBinary(const RunConfig &cfg) {
    QString base = kernelPrefix(cfg.kernel);
    switch (cfg.engine) {
        case Engine::CPU:     break;
        case Engine::CPU_MPI: base += "_mpi";      break;
        case Engine::GPU:     base += "_cuda";     break;
        case Engine::GPU_MPI: base += "_cuda_mpi"; break;
    }
    return resolveBinary(cfg, base);
}

QString Runner::postBinary(const RunConfig &cfg) {
    return resolveBinary(cfg, kernelPrefix(cfg.kernel) + "_post");
}

QString Runner::resolveBinary(const RunConfig &cfg, const QString &name) {
    QString base = name;
#ifdef Q_OS_WIN
    base += ".exe";
#endif
    const char *homeVar = (cfg.kernel == Kernel::RCWA) ? "OPENRCWA_HOME"
                        : (cfg.kernel == Kernel::BPM)  ? "OPENBPM_HOME"
                                                       : "OPENFDTD_HOME";
    const QString dirs[] = {
        cfg.binaryDir,
        qEnvironmentVariable(homeVar),
        QCoreApplication::applicationDirPath() + "/kernel",
        QCoreApplication::applicationDirPath(),
    };
    for (const QString &d : dirs) {
        if (d.isEmpty()) continue;
        const QString full = QDir(d).absoluteFilePath(base);
        if (QFileInfo::exists(full)) return full;
    }
    return base;   // let PATH resolve it
}

void Runner::start(Project *project, const RunConfig &cfg)
{
    if (isRunning() || !project) return;
    m_cfg = cfg;

    if (m_cfg.workingDir.isEmpty()) {
        m_cfg.workingDir = project->filePath().isEmpty()
            ? QStandardPaths::writableLocation(QStandardPaths::TempLocation)
              + "/openfdtd-x"
            : QFileInfo(project->filePath()).path();
    }
    QDir().mkpath(m_cfg.workingDir);

    const QString baseName = project->filePath().isEmpty()
        ? QStringLiteral("project")
        : QFileInfo(project->filePath()).completeBaseName();
    m_ofdPath = QDir(m_cfg.workingDir).filePath(baseName + ".ofd");

    QString err;
    if (!project->save(m_ofdPath, &err)) {
        emit logLine("error: cannot write .ofd: " + err);
        emit finished(false);
        return;
    }

    m_totalSteps = qMax(1, project->general().maxiter);
    m_postPending = (m_cfg.mode == RunMode::Both);
    launch(m_cfg.mode != RunMode::Post);
}

void Runner::launch(bool solverPhase)
{
    QString program;
    QStringList args;

    if (solverPhase) {
        const QString bin = solverBinary(m_cfg);
        if (m_cfg.engine == Engine::CPU_MPI || m_cfg.engine == Engine::GPU_MPI) {
            program = "mpiexec";
            args << "-n" << QString::number(m_cfg.processes) << bin;
        } else {
            program = bin;
        }
        args << "-n" << QString::number(m_cfg.threads);
        args << m_ofdPath;
    } else {
        program = postBinary(m_cfg);
        args << "-n" << QString::number(m_cfg.threads);
        if (m_cfg.evHtml) args << "-html";
        args << m_ofdPath;
    }

    m_proc = new QProcess(this);
    m_proc->setWorkingDirectory(m_cfg.workingDir);
    m_proc->setProcessChannelMode(QProcess::MergedChannels);

    auto env = QProcessEnvironment::systemEnvironment();
    env.insert("OMP_NUM_THREADS", QString::number(m_cfg.threads));
    m_proc->setProcessEnvironment(env);

    connect(m_proc, &QProcess::readyRead, this, &Runner::onReadyRead);
    connect(m_proc, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
            this, &Runner::onFinished);
    connect(m_proc, &QProcess::errorOccurred, this, [this](QProcess::ProcessError) {
        emit logLine("error: " + m_proc->errorString()
                     + " (" + m_proc->program() + ")");
    });

    emit logLine(QStringLiteral("$ cd %1").arg(m_cfg.workingDir));
    emit logLine(QStringLiteral("$ %1 %2").arg(program, args.join(' ')));
    m_proc->start(program, args);
    emit started();
}

void Runner::stop()
{
    m_postPending = false;
    if (!isRunning()) return;
    m_proc->terminate();
    if (!m_proc->waitForFinished(2000))
        m_proc->kill();
}

void Runner::onReadyRead()
{
    if (!m_proc) return;
    // solver iteration lines: "%7d %.6f %.6f" (sol/solve.c)
    static const QRegularExpression stepRe(
        "^\\s*(\\d+)\\s+([-+0-9.eE]+)\\s+([-+0-9.eE]+)\\s*$");
    while (m_proc->canReadLine()) {
        const QString line = QString::fromUtf8(m_proc->readLine())
                                 .remove('\r').trimmed();
        if (line.isEmpty()) continue;
        emit logLine(line);
        const auto m = stepRe.match(line);
        if (m.hasMatch())
            emit progress(m.captured(1).toInt(), m_totalSteps);
    }
}

void Runner::onFinished(int exitCode, QProcess::ExitStatus status)
{
    const bool ok = (status == QProcess::NormalExit && exitCode == 0);
    m_proc->deleteLater();
    m_proc = nullptr;

    if (ok && m_postPending) {
        m_postPending = false;
        emit logLine("=== solver done, running post ===");
        launch(false);
        return;
    }
    m_postPending = false;
    emit logLine(ok ? "=== normal end ==="
                    : QStringLiteral("=== failed (exit %1) ===").arg(exitCode));
    emit finished(ok);
}
