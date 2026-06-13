// PostOpts.h — post-processing control, 1:1 with the keys parsed by
// post/post_data.c (plotiter / plotzin / plotfar1d / plotnear2d ...).
#pragma once
#include <QString>
#include <QVector>

namespace ofd {

// plotzin/yin/ref/spara/coupling = 1            (auto scale)
//                               or 2 min max div (user scale)
struct FreqPlot {
    bool    enabled = false;
    bool    userScale = false;
    double  min = 0, max = 0;
    int     div = 10;
};

// "plotfar1d = X div" or "plotfar1d = V div angle"
struct Far1d {
    QChar   dir = 'X';       // X/Y/Z/V/H
    int     div = 72;
    double  angle = 0;       // V/H only
};

// "plotnear1d = E Z pos1 pos2"
struct Near1d {
    QString cmp = "E";       // E/Ex/.../H/Hx/...
    QChar   dir = 'Z';
    double  pos1 = 0, pos2 = 0;
};

// "plotnear2d = E X pos"
struct Near2d {
    QString cmp = "E";
    QChar   dir = 'X';
    double  pos = 0;
};

struct PostOpts {
    bool matchingloss = false;
    bool plotiter  = true;
    bool plotfeed  = false;
    bool plotpoint = false;
    bool plotsmith = false;

    FreqPlot zin, yin, ref, spara, coupling;
    int  freqdiv = 10;

    // far0d
    bool   far0d = false;
    double far0dTheta = 0, far0dPhi = 0;
    bool   far0dUserScale = false;
    double far0dMin = 0, far0dMax = 0;
    int    far0dDiv = 10;

    // far1d
    QVector<Far1d> far1d;
    int    far1dStyle = 0;
    int    far1dComp[3] = {1, 0, 0};
    bool   far1dDb = true;
    bool   far1dNorm = false;
    bool   far1dUserScale = false;
    double far1dMin = 0, far1dMax = 0;
    int    far1dDiv = 10;

    // far2d
    bool   far2d = false;
    int    far2dDivTheta = 18, far2dDivPhi = 36;
    int    far2dComp[7] = {1, 0, 0, 0, 0, 0, 0};
    bool   far2dDb = true;
    bool   far2dUserScale = false;
    double far2dMin = 0, far2dMax = 0;
    double far2dObj = 0.5;

    // near1d
    QVector<Near1d> near1d;
    bool   near1dDb = false;
    bool   near1dNoinc = false;
    bool   near1dUserScale = false;
    double near1dMin = 0, near1dMax = 0;
    int    near1dDiv = 10;

    // near2d
    QVector<Near2d> near2d;
    int    near2dDim[2] = {1, 1};
    bool   near2dFrame = false;
    bool   near2dDb = false;
    bool   near2dUserScale = false;
    double near2dMin = 0, near2dMax = 0;
    bool   near2dContour = false;
    int    near2dObj = 1;
    bool   near2dNoinc = false;
    bool   near2dZoom = false;
    double near2dHzoom[2] = {0, 0};
    double near2dVzoom[2] = {0, 0};
};

} // namespace ofd
