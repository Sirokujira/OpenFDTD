// OfdIO.cpp
#include "OfdIO.h"

#include <QFile>
#include <QTextStream>
#include <QRegularExpression>
#include <QStringList>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>

using namespace ofd;

// ── helpers ─────────────────────────────────────────────────────────────────
static QString num(double v)            { return QString::number(v, 'g', 10); }
static QString joinNums(const double *v, int n) {
    QString s;
    for (int i = 0; i < n; ++i) { s += ' '; s += num(v[i]); }
    return s;
}

// Split off a trailing " # name" comment (GUI metadata the kernel tokenizer
// simply ignores as extra tokens). Returns the name, shortens the line.
static QString takeTrailingName(QString &line) {
    const int hash = line.indexOf('#');
    if (hash < 0) return {};
    const QString name = line.mid(hash + 1).trimmed();
    line = line.left(hash).trimmed();
    return name;
}

// ── Serializer ──────────────────────────────────────────────────────────────
QString OfdIO::serialize(const Project &p)
{
    QString text;
    QTextStream out(&text);

    out << "OpenFDTD 4 2\n";
    const GeneralOpts &g = p.general();
    if (!g.title.isEmpty())
        out << "title = " << g.title << "\n";

    // mesh: x0 d1 x1 d2 x2 ...
    static const char *meshKey[3] = { "xmesh", "ymesh", "zmesh" };
    for (int a = 0; a < 3; ++a) {
        const MeshAxis &ax = p.mesh(a);
        if (ax.nodes.isEmpty()) continue;
        out << meshKey[a] << " = " << num(ax.nodes[0]);
        for (int i = 0; i < ax.divs.size(); ++i)
            out << " " << ax.divs[i] << " " << num(ax.nodes[i+1]);
        out << "\n";
    }

    // materials (ID 0/1 are built-in; user materials start at 2)
    for (const Material &m : p.materials()) {
        out << "material = " << m.type;
        if (m.type == 2)
            out << " " << num(m.einf) << " " << num(m.ae)
                << " " << num(m.be)   << " " << num(m.ce);
        else
            out << " " << num(m.epsr) << " " << num(m.esgm)
                << " " << num(m.amur) << " " << num(m.msgm);
        if (!m.name.isEmpty()) out << " # " << m.name;
        out << "\n";
    }

    // geometries
    for (const Geometry &ge : p.geometries()) {
        out << "geometry = " << ge.materialId << " " << ge.shape
            << joinNums(ge.g, Geometry::paramCount(ge.shape));
        if (!ge.name.isEmpty()) out << " # " << ge.name;
        out << "\n";
    }

    // sources
    for (const Feed &f : p.feeds())
        out << "feed = " << f.dir << " " << num(f.x) << " " << num(f.y)
            << " " << num(f.z) << " " << num(f.volt) << " " << num(f.delay)
            << " " << num(f.z0) << "\n";
    if (p.planewave().enabled)
        out << "planewave = " << num(p.planewave().theta) << " "
            << num(p.planewave().phi) << " " << p.planewave().pol << "\n";
    for (int i = 0; i < p.probes().size(); ++i) {
        const Probe &pr = p.probes()[i];
        out << "point = " << pr.dir << " " << num(pr.x) << " " << num(pr.y)
            << " " << num(pr.z);
        if (i == 0 && !pr.propagation.isEmpty())
            out << " " << pr.propagation;
        out << "\n";
    }
    for (const Load &l : p.loads())
        out << "load = " << l.dir << " " << num(l.x) << " " << num(l.y)
            << " " << num(l.z) << " " << l.kind << " " << num(l.value) << "\n";
    if (g.rfeed != 0)
        out << "rfeed = " << num(g.rfeed) << "\n";

    // boundary conditions
    if (g.abc == 1)
        out << "abc = 1 " << g.pmlL << " " << num(g.pmlM)
            << " " << num(g.pmlR0) << "\n";
    if (g.pbcX || g.pbcY || g.pbcZ)
        out << "pbc = " << int(g.pbcX) << " " << int(g.pbcY)
            << " " << int(g.pbcZ) << "\n";

    // frequencies / solver
    if (g.hasF1)
        out << "frequency1 = " << num(g.f1min) << " " << num(g.f1max)
            << " " << g.f1div << "\n";
    if (g.hasF2)
        out << "frequency2 = " << num(g.f2min) << " " << num(g.f2max)
            << " " << g.f2div << "\n";
    out << "solver = " << g.maxiter << " " << g.nout
        << " " << num(g.converg) << "\n";
    if (g.dt > 0) out << "timestep = "   << num(g.dt) << "\n";
    if (g.tw > 0) out << "pulsewidth = " << num(g.tw) << "\n";
    if (g.plot3dgeom) out << "plot3dgeom = " << g.plot3dgeom << "\n";

    // ── post-processing keys ────────────────────────────────────────────
    const PostOpts &po = p.post();
    auto freqPlot = [&out](const char *key, const FreqPlot &fp) {
        if (!fp.enabled) return;
        out << key << " = ";
        if (fp.userScale)
            out << "2 " << num(fp.min) << " " << num(fp.max) << " " << fp.div;
        else
            out << "1";
        out << "\n";
    };

    if (po.matchingloss) out << "matchingloss = 1\n";
    if (po.plotiter)     out << "plotiter = 1\n";
    if (po.plotfeed)     out << "plotfeed = 1\n";
    if (po.plotpoint)    out << "plotpoint = 1\n";
    if (po.plotsmith)    out << "plotsmith = 1\n";
    freqPlot("plotzin",      po.zin);
    freqPlot("plotyin",      po.yin);
    freqPlot("plotref",      po.ref);
    freqPlot("plotspara",    po.spara);
    freqPlot("plotcoupling", po.coupling);
    if (po.freqdiv != 10) out << "freqdiv = " << po.freqdiv << "\n";

    if (po.far0d) {
        out << "plotfar0d = " << num(po.far0dTheta) << " " << num(po.far0dPhi);
        if (po.far0dUserScale)
            out << " 2 " << num(po.far0dMin) << " " << num(po.far0dMax)
                << " " << po.far0dDiv;
        else
            out << " 1";
        out << "\n";
    }

    for (const Far1d &f : po.far1d) {
        out << "plotfar1d = " << f.dir << " " << f.div;
        if (f.dir == 'V' || f.dir == 'H') out << " " << num(f.angle);
        out << "\n";
    }
    if (!po.far1d.isEmpty()) {
        out << "far1dstyle = " << po.far1dStyle << "\n";
        out << "far1dcomponent = " << po.far1dComp[0] << " "
            << po.far1dComp[1] << " " << po.far1dComp[2] << "\n";
        out << "far1ddb = " << int(po.far1dDb) << "\n";
        if (po.far1dNorm) out << "far1dnorm = 1\n";
        if (po.far1dUserScale)
            out << "far1dscale = " << num(po.far1dMin) << " "
                << num(po.far1dMax) << " " << po.far1dDiv << "\n";
    }

    if (po.far2d) {
        out << "plotfar2d = " << po.far2dDivTheta << " " << po.far2dDivPhi << "\n";
        out << "far2dcomponent =";
        for (int i = 0; i < 7; ++i) out << " " << po.far2dComp[i];
        out << "\n";
        out << "far2ddb = " << int(po.far2dDb) << "\n";
        if (po.far2dUserScale)
            out << "far2dscale = " << num(po.far2dMin) << " "
                << num(po.far2dMax) << "\n";
        out << "far2dobj = " << num(po.far2dObj) << "\n";
    }

    for (const Near1d &n : po.near1d)
        out << "plotnear1d = " << n.cmp << " " << n.dir << " "
            << num(n.pos1) << " " << num(n.pos2) << "\n";
    if (!po.near1d.isEmpty()) {
        out << "near1ddb = " << int(po.near1dDb) << "\n";
        if (po.near1dUserScale)
            out << "near1dscale = " << num(po.near1dMin) << " "
                << num(po.near1dMax) << " " << po.near1dDiv << "\n";
        out << "near1dnoinc = " << int(po.near1dNoinc) << "\n";
    }

    for (const Near2d &n : po.near2d)
        out << "plotnear2d = " << n.cmp << " " << n.dir << " "
            << num(n.pos) << "\n";
    if (!po.near2d.isEmpty()) {
        out << "near2ddim = " << po.near2dDim[0] << " " << po.near2dDim[1] << "\n";
        if (po.near2dFrame) out << "near2dframe = 1\n";
        out << "near2ddb = " << int(po.near2dDb) << "\n";
        if (po.near2dUserScale)
            out << "near2dscale = " << num(po.near2dMin) << " "
                << num(po.near2dMax) << "\n";
        out << "near2dcontour = " << int(po.near2dContour) << "\n";
        out << "near2dobj = " << po.near2dObj << "\n";
        out << "near2dnoinc = " << int(po.near2dNoinc) << "\n";
        if (po.near2dZoom)
            out << "near2dzoom = " << num(po.near2dHzoom[0]) << " "
                << num(po.near2dHzoom[1]) << " " << num(po.near2dVzoom[0])
                << " " << num(po.near2dVzoom[1]) << "\n";
    }

    // keys the GUI doesn't model — preserved from the loaded file
    for (const QString &line : p.extraLines())
        out << line << "\n";

    out << "end\n";
    return text;
}

bool OfdIO::save(const QString &path, const Project &project, QString *err)
{
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    QTextStream out(&f);
    out.setEncoding(QStringConverter::Utf8);
    out.setGenerateByteOrderMark(false);    // 本家は BOM 無し
    out << serialize(project);
    return true;
}

// ── Loader ──────────────────────────────────────────────────────────────────
bool OfdIO::load(const QString &path, Project &project, QString *err)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    QTextStream in(&f);
    in.setEncoding(QStringConverter::Utf8);
    return parse(in.readAll(), project, err);
}

bool OfdIO::parse(const QString &text, Project &project, QString *err)
{
    static const QRegularExpression ws("\\s+");

    GeneralOpts &g = project.general();
    PostOpts    &po = project.post();
    po.plotiter = false;   // file decides what is on
    g.hasF1 = g.hasF2 = false;

    bool header = false;
    const QStringList lines = text.split('\n');
    for (QString rawLine : lines) {
        QString line = rawLine.trimmed();
        if (line.isEmpty() || line.startsWith('#')) continue;
        if (line.startsWith("end")) break;

        if (!header) {
            const QStringList t = line.split(ws, Qt::SkipEmptyParts);
            if (t.size() < 3 || (t[0] != "OpenFDTD" && t[0] != "OpenTHFD")) {
                if (err) *err = "not OpenFDTD/OpenTHFD data";
                return false;
            }
            header = true;
            continue;
        }

        const int eq = line.indexOf('=');
        if (eq < 0) continue;
        const QString key = line.left(eq).trimmed();
        QString val = line.mid(eq + 1).trimmed();
        const QString name = (key == "material" || key == "geometry")
                             ? takeTrailingName(val) : QString();
        const QStringList t = val.split(ws, Qt::SkipEmptyParts);
        auto d = [&t](int i) { return (i < t.size()) ? t[i].toDouble() : 0.0; };
        auto n = [&t](int i) { return (i < t.size()) ? t[i].toInt()    : 0;   };

        if (key == "title") {
            g.title = val;
        }
        else if (key == "xmesh" || key == "ymesh" || key == "zmesh") {
            const int a = (key[0] == 'x') ? 0 : (key[0] == 'y') ? 1 : 2;
            MeshAxis &ax = project.mesh(a);
            ax.nodes.clear(); ax.divs.clear();
            if (t.size() < 3 || t.size() % 2 == 0) {
                if (err) *err = "invalid " + key + " data";
                return false;
            }
            ax.nodes.push_back(d(0));
            for (int i = 1; i + 1 < t.size(); i += 2) {
                ax.divs.push_back(t[i].toInt());
                ax.nodes.push_back(t[i+1].toDouble());
            }
        }
        else if (key == "material") {
            Material m;
            if (t.size() >= 5 && (t[0] == "1" || t[0] == "2")) {
                m.type = n(0);
                if (m.type == 2) {
                    m.einf = d(1); m.ae = d(2); m.be = d(3); m.ce = d(4);
                } else {
                    m.epsr = d(1); m.esgm = d(2); m.amur = d(3); m.msgm = d(4);
                }
            } else if (t.size() >= 4) {
                // old format (version < 2.2): εr σ μr σm without type
                m.type = 1;
                m.epsr = d(0); m.esgm = d(1); m.amur = d(2); m.msgm = d(3);
            }
            m.name = name;
            project.materials().push_back(m);
        }
        else if (key == "geometry") {
            Geometry ge;
            ge.materialId = n(0);
            ge.shape      = n(1);
            const int np = Geometry::paramCount(ge.shape);
            for (int i = 0; i < np && i < 8; ++i) ge.g[i] = d(2 + i);
            ge.name = name;
            project.geometries().push_back(ge);
        }
        else if (key == "name") {
            // 本家カーネルは無視する — GUI も読み飛ばす
        }
        else if (key == "feed" && t.size() >= 7) {
            Feed fd;
            fd.dir = t[0].toUpper().at(0);
            fd.x = d(1); fd.y = d(2); fd.z = d(3);
            fd.volt = d(4); fd.delay = d(5); fd.z0 = d(6);
            project.feeds().push_back(fd);
        }
        else if (key == "planewave" && t.size() >= 3) {
            project.planewave().enabled = true;
            project.planewave().theta = d(0);
            project.planewave().phi   = d(1);
            project.planewave().pol   = n(2);
        }
        else if (key == "point" && t.size() >= 4) {
            Probe pr;
            pr.dir = t[0].toUpper().at(0);
            pr.x = d(1); pr.y = d(2); pr.z = d(3);
            if (project.probes().isEmpty() && t.size() >= 5)
                pr.propagation = t[4];
            project.probes().push_back(pr);
        }
        else if (key == "load" && t.size() >= 6) {
            Load l;
            l.dir = t[0].toUpper().at(0);
            l.x = d(1); l.y = d(2); l.z = d(3);
            l.kind = t[4].toUpper().at(0);
            l.value = d(5);
            project.loads().push_back(l);
        }
        else if (key == "rfeed")       { g.rfeed = d(0); }
        else if (key == "abc") {
            if (t.size() >= 1 && t[0].startsWith('0')) {
                g.abc = 0;
            } else if (t.size() >= 4 && t[0].startsWith('1')) {
                g.abc = 1;
                g.pmlL = n(1); g.pmlM = d(2); g.pmlR0 = d(3);
            }
        }
        else if (key == "pbc" && t.size() >= 3) {
            g.pbcX = n(0); g.pbcY = n(1); g.pbcZ = n(2);
        }
        else if (key == "frequency1" && t.size() >= 3) {
            g.hasF1 = true;
            g.f1min = d(0); g.f1max = d(1); g.f1div = n(2);
        }
        else if ((key == "frequency2" || key == "frequency") && t.size() >= 3) {
            g.hasF2 = true;
            g.f2min = d(0); g.f2max = d(1); g.f2div = n(2);
        }
        else if (key == "solver" && t.size() >= 3) {
            g.maxiter = n(0); g.nout = n(1); g.converg = d(2);
        }
        else if (key == "timestep")    { g.dt = d(0); }
        else if (key == "pulsewidth")  { g.tw = d(0); }
        else if (key == "plot3dgeom")  { g.plot3dgeom = n(0); }

        // ── post-processing keys ────────────────────────────────────────
        else if (key == "matchingloss") { po.matchingloss = n(0); }
        else if (key == "plotiter")     { po.plotiter  = n(0); }
        else if (key == "plotfeed")     { po.plotfeed  = n(0); }
        else if (key == "plotpoint")    { po.plotpoint = n(0); }
        else if (key == "plotsmith")    { po.plotsmith = n(0); }
        else if (key == "plotzin" || key == "plotyin" || key == "plotref" ||
                 key == "plotspara" || key == "plotcoupling") {
            FreqPlot &fp = (key == "plotzin")   ? po.zin
                         : (key == "plotyin")   ? po.yin
                         : (key == "plotref")   ? po.ref
                         : (key == "plotspara") ? po.spara : po.coupling;
            if (t.size() >= 1 && t[0] == "1") {
                fp.enabled = true; fp.userScale = false;
            } else if (t.size() >= 4 && t[0] == "2") {
                fp.enabled = true; fp.userScale = true;
                fp.min = d(1); fp.max = d(2); fp.div = n(3);
            }
        }
        else if (key == "freqdiv")      { po.freqdiv = n(0); }
        else if (key == "plotfar0d" && t.size() >= 3) {
            po.far0d = true;
            po.far0dTheta = d(0); po.far0dPhi = d(1);
            if (t[2] == "2" && t.size() >= 6) {
                po.far0dUserScale = true;
                po.far0dMin = d(3); po.far0dMax = d(4); po.far0dDiv = n(5);
            }
        }
        else if (key == "plotfar1d" && t.size() >= 2) {
            Far1d f;
            f.dir = t[0].toUpper().at(0);
            f.div = n(1);
            if ((f.dir == 'V' || f.dir == 'H') && t.size() >= 3) f.angle = d(2);
            po.far1d.push_back(f);
        }
        else if (key == "plotfar2d" && t.size() >= 2) {
            po.far2d = true;
            po.far2dDivTheta = n(0); po.far2dDivPhi = n(1);
        }
        else if (key == "plotnear1d" && t.size() >= 4) {
            Near1d nr;
            nr.cmp = t[0]; nr.dir = t[1].toUpper().at(0);
            nr.pos1 = d(2); nr.pos2 = d(3);
            po.near1d.push_back(nr);
        }
        else if (key == "plotnear2d" && t.size() >= 3) {
            Near2d nr;
            nr.cmp = t[0]; nr.dir = t[1].toUpper().at(0);
            nr.pos = d(2);
            po.near2d.push_back(nr);
        }
        else if (key == "far1dcomponent" && t.size() >= 3) {
            for (int i = 0; i < 3; ++i) po.far1dComp[i] = n(i);
        }
        else if (key == "far1dstyle")   { po.far1dStyle = n(0); }
        else if (key == "far1ddb")      { po.far1dDb = n(0); }
        else if (key == "far1dnorm")    { po.far1dNorm = n(0); }
        else if (key == "far1dscale" && t.size() >= 3) {
            po.far1dUserScale = true;
            po.far1dMin = d(0); po.far1dMax = d(1); po.far1dDiv = n(2);
        }
        else if (key == "far2dcomponent" && t.size() >= 7) {
            for (int i = 0; i < 7; ++i) po.far2dComp[i] = n(i);
        }
        else if (key == "far2ddb")      { po.far2dDb = n(0); }
        else if (key == "far2dscale" && t.size() >= 2) {
            po.far2dUserScale = true;
            po.far2dMin = d(0); po.far2dMax = d(1);
        }
        else if (key == "far2dobj")     { po.far2dObj = d(0); }
        else if (key == "near1ddb")     { po.near1dDb = n(0); }
        else if (key == "near1dscale" && t.size() >= 3) {
            po.near1dUserScale = true;
            po.near1dMin = d(0); po.near1dMax = d(1); po.near1dDiv = n(2);
        }
        else if (key == "near1dnoinc")  { po.near1dNoinc = n(0); }
        else if (key == "near2ddim" && t.size() >= 2) {
            po.near2dDim[0] = n(0); po.near2dDim[1] = n(1);
        }
        else if (key == "near2dframe")  { po.near2dFrame = n(0); }
        else if (key == "near2ddb")     { po.near2dDb = n(0); }
        else if (key == "near2dscale" && t.size() >= 2) {
            po.near2dUserScale = true;
            po.near2dMin = d(0); po.near2dMax = d(1);
        }
        else if (key == "near2dcontour"){ po.near2dContour = n(0); }
        else if (key == "near2dobj")    { po.near2dObj = n(0); }
        else if (key == "near2dnoinc")  { po.near2dNoinc = n(0); }
        else if (key == "near2dzoom" && t.size() >= 4) {
            po.near2dZoom = true;
            po.near2dHzoom[0] = qMin(d(0), d(1));
            po.near2dHzoom[1] = qMax(d(0), d(1));
            po.near2dVzoom[0] = qMin(d(2), d(3));
            po.near2dVzoom[1] = qMax(d(2), d(3));
        }
        else {
            // unknown key — keep verbatim for round-trip safety
            project.extraLines().push_back(line);
        }
    }

    if (!header) {
        if (err) *err = "empty file";
        return false;
    }
    return true;
}

// ── .ofdx (JSON sidecar) ────────────────────────────────────────────────────
bool OfdxIO::save(const QString &path, const Project &p, QString *err)
{
    QJsonObject root;
    root["schemaVersion"] = "1.0";
    root["domain"] = domainKey(p.activeDomain());

    {
        const OpticalOpts &o = p.optical();
        QJsonObject opt;
        opt["solver"] = int(o.solver);
        opt["mode"]   = int(o.mode);
        opt["wavelength"] = QJsonObject{
            {"min_nm", o.lambdaMin}, {"max_nm", o.lambdaMax}, {"div", o.lambdaDiv} };
        opt["rcwa"] = QJsonObject{
            {"nx", o.rcwaNx}, {"ny", o.rcwaNy},
            {"period_x_nm", o.rcwaPeriodX}, {"period_y_nm", o.rcwaPeriodY},
            {"layers", o.rcwaLayers} };
        opt["bpm"] = QJsonObject{
            {"algorithm", o.bpmAlgorithm}, {"dz_nm", o.bpmDz},
            {"ref_index", o.bpmRefIndex}, {"input_mode", o.bpmInputMode} };
        opt["fmm"] = QJsonObject{
            {"harmonics", o.fmmHarmonics}, {"li_rules", o.fmmLiRules} };
        opt["bpf"] = QJsonObject{
            {"band_nm", QJsonArray{ o.bpfBandMin, o.bpfBandMax }}, {"Q", o.bpfQ} };
        opt["ring"] = QJsonObject{
            {"radius_um", o.ringRadius_um}, {"gap_nm", o.ringGap_nm} };
        root["optical"] = opt;
    }
    {
        const AcousticOpts &a = p.acoustic();
        root["acoustic"] = QJsonObject{
            {"rt60", a.rt60}, {"c80", a.c80}, {"d50", a.d50},
            {"sti", a.sti}, {"edt", a.edt},
            {"impulse_response", a.impulseResponse},
            {"auralization", a.auralization},
            {"sample_rate", a.sampleRate},
            {"src_directivity", a.srcDirectivity},
            {"src_spl_db", a.srcSPL_dB},
            {"mic_count", a.micCount} };
    }
    {
        const UnderwaterOpts &u = p.underwater();
        QJsonArray ssp;
        for (const SSPPoint &pt : u.ssp)
            ssp.append(QJsonObject{ {"depth_m", pt.depth_m}, {"c_mps", pt.c_mps} });
        root["underwater"] = QJsonObject{
            {"temp_C", u.waterTemp_C}, {"salinity_psu", u.salinity_psu},
            {"ssp", ssp}, {"sofar", u.sofar},
            {"bottom_type", u.bottomType},
            {"bottom_c_mps", u.bottomC_mps},
            {"bottom_rho_kgm3", u.bottomRho_kgm3},
            {"sonar_freq_khz", u.sonarFreq_kHz},
            {"sonar_sl_db", u.sonarSL_dB},
            {"range_max_km", u.rangeMax_km} };
    }
    {
        // API key is NOT persisted here — it lives in QSettings
        const Tidy3dOpts &t = p.tidy3d();
        root["tidy3d"] = QJsonObject{
            {"project_name", t.projectName},
            {"resolution", t.resolution},
            {"auto_pml", t.autoPml} };
    }

    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    f.write(QJsonDocument(root).toJson(QJsonDocument::Indented));
    return true;
}

bool OfdxIO::load(const QString &path, Project &p, QString *err)
{
    QFile f(path);
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    const QJsonDocument doc = QJsonDocument::fromJson(f.readAll());
    if (!doc.isObject()) {
        if (err) *err = "invalid JSON";
        return false;
    }
    const QJsonObject root = doc.object();
    p.setActiveDomain(domainFromKey(root.value("domain").toString("em")));

    if (root.contains("optical")) {
        const QJsonObject opt = root["optical"].toObject();
        OpticalOpts &o = p.optical();
        o.solver = OpticalSolver(opt.value("solver").toInt(int(o.solver)));
        o.mode   = OpticalMode(opt.value("mode").toInt(int(o.mode)));
        const QJsonObject wl = opt["wavelength"].toObject();
        o.lambdaMin = wl.value("min_nm").toDouble(o.lambdaMin);
        o.lambdaMax = wl.value("max_nm").toDouble(o.lambdaMax);
        o.lambdaDiv = wl.value("div").toInt(o.lambdaDiv);
        const QJsonObject rc = opt["rcwa"].toObject();
        o.rcwaNx = rc.value("nx").toInt(o.rcwaNx);
        o.rcwaNy = rc.value("ny").toInt(o.rcwaNy);
        o.rcwaPeriodX = rc.value("period_x_nm").toDouble(o.rcwaPeriodX);
        o.rcwaPeriodY = rc.value("period_y_nm").toDouble(o.rcwaPeriodY);
        o.rcwaLayers = rc.value("layers").toInt(o.rcwaLayers);
        const QJsonObject bp = opt["bpm"].toObject();
        o.bpmAlgorithm = bp.value("algorithm").toInt(o.bpmAlgorithm);
        o.bpmDz = bp.value("dz_nm").toDouble(o.bpmDz);
        o.bpmRefIndex = bp.value("ref_index").toDouble(o.bpmRefIndex);
        o.bpmInputMode = bp.value("input_mode").toInt(o.bpmInputMode);
        const QJsonObject fm = opt["fmm"].toObject();
        o.fmmHarmonics = fm.value("harmonics").toInt(o.fmmHarmonics);
        o.fmmLiRules = fm.value("li_rules").toBool(o.fmmLiRules);
        const QJsonObject bpf = opt["bpf"].toObject();
        const QJsonArray band = bpf["band_nm"].toArray();
        if (band.size() >= 2) {
            o.bpfBandMin = band[0].toDouble(o.bpfBandMin);
            o.bpfBandMax = band[1].toDouble(o.bpfBandMax);
        }
        o.bpfQ = bpf.value("Q").toDouble(o.bpfQ);
        const QJsonObject ring = opt["ring"].toObject();
        o.ringRadius_um = ring.value("radius_um").toDouble(o.ringRadius_um);
        o.ringGap_nm = ring.value("gap_nm").toDouble(o.ringGap_nm);
    }
    if (root.contains("acoustic")) {
        const QJsonObject ac = root["acoustic"].toObject();
        AcousticOpts &a = p.acoustic();
        a.rt60 = ac.value("rt60").toBool(a.rt60);
        a.c80  = ac.value("c80").toBool(a.c80);
        a.d50  = ac.value("d50").toBool(a.d50);
        a.sti  = ac.value("sti").toBool(a.sti);
        a.edt  = ac.value("edt").toBool(a.edt);
        a.impulseResponse = ac.value("impulse_response").toBool(a.impulseResponse);
        a.auralization = ac.value("auralization").toBool(a.auralization);
        a.sampleRate = ac.value("sample_rate").toInt(a.sampleRate);
        a.srcDirectivity = ac.value("src_directivity").toString(a.srcDirectivity);
        a.srcSPL_dB = ac.value("src_spl_db").toDouble(a.srcSPL_dB);
        a.micCount = ac.value("mic_count").toInt(a.micCount);
    }
    if (root.contains("underwater")) {
        const QJsonObject uw = root["underwater"].toObject();
        UnderwaterOpts &u = p.underwater();
        u.waterTemp_C = uw.value("temp_C").toDouble(u.waterTemp_C);
        u.salinity_psu = uw.value("salinity_psu").toDouble(u.salinity_psu);
        if (uw.contains("ssp")) {
            u.ssp.clear();
            for (const QJsonValue &v : uw["ssp"].toArray()) {
                const QJsonObject o = v.toObject();
                u.ssp.push_back({ o.value("depth_m").toDouble(),
                                  o.value("c_mps").toDouble() });
            }
        }
        u.sofar = uw.value("sofar").toBool(u.sofar);
        u.bottomType = uw.value("bottom_type").toString(u.bottomType);
        u.bottomC_mps = uw.value("bottom_c_mps").toDouble(u.bottomC_mps);
        u.bottomRho_kgm3 = uw.value("bottom_rho_kgm3").toDouble(u.bottomRho_kgm3);
        u.sonarFreq_kHz = uw.value("sonar_freq_khz").toDouble(u.sonarFreq_kHz);
        u.sonarSL_dB = uw.value("sonar_sl_db").toDouble(u.sonarSL_dB);
        u.rangeMax_km = uw.value("range_max_km").toDouble(u.rangeMax_km);
    }
    if (root.contains("tidy3d")) {
        const QJsonObject t3 = root["tidy3d"].toObject();
        Tidy3dOpts &t = p.tidy3d();
        t.projectName = t3.value("project_name").toString(t.projectName);
        t.resolution = t3.value("resolution").toString(t.resolution);
        t.autoPml = t3.value("auto_pml").toBool(t.autoPml);
    }
    return true;
}
