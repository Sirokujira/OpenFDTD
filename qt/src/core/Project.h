// Project.h — top-level data model.
//
// One Project holds everything the GUI edits. Tabs are *views* of this model;
// when the user clicks Compute we serialize Project to a .ofd file (compatible
// with the original OpenFDTD kernel) plus a .ofdx JSON sidecar carrying the
// extension-domain settings, then hand both to Runner.
#pragma once

#include <QObject>
#include <QString>
#include <QStringList>
#include <QVector>
#include "Domain.h"
#include "Material.h"
#include "Geometry.h"
#include "Source.h"
#include "MeshAxis.h"
#include "PostOpts.h"

namespace ofd {

// ── 全般タブ (solver / abc / pbc / frequency) ───────────────────────────────
struct GeneralOpts {
    QString title;
    int     maxiter = 3000;
    int     nout    = 50;
    double  converg = 1e-3;

    int     abc = 0;          // 0 = Mur-1st, 1 = PML
    int     pmlL = 5;         // PML layers
    double  pmlM = 2.0;       // PML order m
    double  pmlR0 = 1e-5;     // PML reflection coefficient

    bool    pbcX = false, pbcY = false, pbcZ = false;

    double  f1min = 2e9, f1max = 3e9;  int f1div = 10;  // frequency1 [Hz]
    double  f2min = 3e9, f2max = 3e9;  int f2div = 0;   // frequency2 [Hz]
    bool    hasF1 = true, hasF2 = true;  // false = key absent in loaded file

    double  dt = 0;           // timestep   (0 = auto)
    double  tw = 0;           // pulsewidth (0 = auto)
    double  rfeed = 0;        // feed resistance correction
    int     plot3dgeom = 0;
};

// ── 光ドメイン拡張 (.ofdx) ──────────────────────────────────────────────────
enum class OpticalSolver { FDTD, RCWA, BPM, FMM };
enum class OpticalMode { BPF, Waveguide, Ring, MZI, Metasurface, PhC, NF2FF, SParam };

struct OpticalOpts {
    OpticalSolver solver = OpticalSolver::FDTD;
    OpticalMode   mode   = OpticalMode::BPF;
    double  lambdaMin = 1500.0;   // nm
    double  lambdaMax = 1600.0;
    int     lambdaDiv = 201;

    // RCWA
    int     rcwaNx = 11, rcwaNy = 11;   // Fourier orders
    double  rcwaPeriodX = 600.0, rcwaPeriodY = 600.0;  // nm
    int     rcwaLayers = 8;

    // BPM
    int     bpmAlgorithm = 0;     // 0=FFT, 1=FDM, 2=Wide-Angle Padé
    double  bpmDz = 50.0;         // nm
    double  bpmRefIndex = 1.45;
    int     bpmInputMode = 0;     // 0=TE0, 1=TE1, 2=TM0, 3=Gaussian

    // FMM
    int     fmmHarmonics = 15;
    bool    fmmLiRules = true;

    // BPF
    double  bpfBandMin = 1540.0, bpfBandMax = 1560.0;  // nm
    double  bpfQ = 10000.0;

    // Ring
    double  ringRadius_um = 5.0;
    double  ringGap_nm = 200.0;
};

// ── 室内音響ドメイン拡張 (.ofdx) ────────────────────────────────────────────
struct AcousticOpts {
    bool    rt60 = true, c80 = true, d50 = false, sti = false, edt = false;
    bool    impulseResponse = true;
    bool    auralization = false;
    int     sampleRate = 48000;     // WAV出力サンプリング周波数
    QString srcDirectivity = "omni"; // omni / cardioid / speaker
    double  srcSPL_dB = 94.0;
    int     micCount = 1;
};

// ── 水中音響ドメイン拡張 (.ofdx) ────────────────────────────────────────────
struct SSPPoint { double depth_m; double c_mps; };

struct UnderwaterOpts {
    double  waterTemp_C = 15.0;
    double  salinity_psu = 34.5;
    QVector<SSPPoint> ssp;
    bool    sofar = false;
    QString bottomType = "sand";
    double  bottomC_mps = 1650.0;
    double  bottomRho_kgm3 = 1900.0;
    double  sonarFreq_kHz = 3.5;
    double  sonarSL_dB = 220.0;
    double  rangeMax_km = 50.0;
};

// ── tidy3d クラウドバックエンド (光ドメイン専用, .ofdx) ─────────────────────
struct Tidy3dOpts {
    QString projectName = "openfdtd-x";
    QString resolution  = "medium";   // coarse / medium / fine
    bool    autoPml     = true;
};

// ────────────────────────────────────────────────────────────────────────────
class Project : public QObject {
    Q_OBJECT
public:
    explicit Project(QObject *parent = nullptr);

    Domain activeDomain() const { return m_domain; }

    GeneralOpts        &general()     { return m_general; }
    MeshAxis           &mesh(int axis) { return m_mesh[axis]; }   // 0=x 1=y 2=z
    QVector<Material>  &materials()   { return m_materials; }
    QVector<Load>      &loads()       { return m_loads; }
    QVector<Geometry>  &geometries()  { return m_geometries; }
    QVector<Feed>      &feeds()       { return m_feeds; }
    PlaneWave          &planewave()   { return m_planewave; }
    QVector<Probe>     &probes()      { return m_probes; }
    PostOpts           &post()        { return m_post; }
    OpticalOpts        &optical()     { return m_optical; }
    AcousticOpts       &acoustic()    { return m_acoustic; }
    UnderwaterOpts     &underwater()  { return m_underwater; }
    Tidy3dOpts         &tidy3d()      { return m_tidy3d; }

    const GeneralOpts       &general()    const { return m_general; }
    const MeshAxis          &mesh(int axis) const { return m_mesh[axis]; }
    const QVector<Material> &materials()  const { return m_materials; }
    const QVector<Load>     &loads()      const { return m_loads; }
    const QVector<Geometry> &geometries() const { return m_geometries; }
    const QVector<Feed>     &feeds()      const { return m_feeds; }
    const PlaneWave         &planewave()  const { return m_planewave; }
    const QVector<Probe>    &probes()     const { return m_probes; }
    const PostOpts          &post()       const { return m_post; }
    const OpticalOpts       &optical()    const { return m_optical; }
    const AcousticOpts      &acoustic()   const { return m_acoustic; }
    const UnderwaterOpts    &underwater() const { return m_underwater; }
    const Tidy3dOpts        &tidy3d()     const { return m_tidy3d; }

    // Lines from a loaded .ofd that the GUI doesn't model — preserved verbatim
    // on save so a hand-edited file survives a GUI round trip.
    QStringList &extraLines() { return m_extraLines; }
    const QStringList &extraLines() const { return m_extraLines; }

    QString filePath() const { return m_filePath; }
    void    setFilePath(const QString &p) { m_filePath = p; }

    qint64 totalCells() const;
    double estimatedMemoryMB() const;
    double courantDt() const;    // CFL timestep estimate from min spacing

    void clear();                // reset to the default new-project state

    // Persistence (implemented in io/OfdIO.cpp; these are thin wrappers)
    bool load(const QString &path, QString *err = nullptr);  // .ofd + .ofdx
    bool save(const QString &path, QString *err = nullptr);

public slots:
    void setActiveDomain(ofd::Domain d);
    void touch() { emit changed(); }   // call after editing through references

signals:
    void domainChanged(ofd::Domain);
    void changed();    // structural change → viewport / tree / statusbar refresh
    void loaded();     // file (re)loaded → all tabs re-read their widgets

private:
    Domain  m_domain = Domain::EM;
    QString m_filePath;

    GeneralOpts        m_general;
    MeshAxis           m_mesh[3];
    QVector<Material>  m_materials;
    QVector<Load>      m_loads;
    QVector<Geometry>  m_geometries;
    QVector<Feed>      m_feeds;
    PlaneWave          m_planewave;
    QVector<Probe>     m_probes;
    PostOpts           m_post;
    OpticalOpts        m_optical;
    AcousticOpts       m_acoustic;
    UnderwaterOpts     m_underwater;
    Tidy3dOpts         m_tidy3d;
    QStringList        m_extraLines;
};

} // namespace ofd
