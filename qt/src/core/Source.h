// Source.h — feed / plane wave / observation point,
// 1:1 with the OpenFDTD "feed =", "planewave =", "point =" lines.
#pragma once
#include <QString>

namespace ofd {

// "feed = dir x y z volt delay z0"
struct Feed {
    QChar   dir = 'Z';
    double  x = 0, y = 0, z = 0;
    double  volt  = 1.0;     // amplitude [V]
    double  delay = 0.0;     // pulse delay [deg]
    double  z0    = 50.0;    // feed internal impedance [Ω]
};

// "planewave = theta phi pol"  (pol: 1=V, 2=H)
struct PlaneWave {
    bool    enabled = false;
    double  theta = 90.0;
    double  phi   = 0.0;
    int     pol   = 1;
};

// "point = dir x y z [propagation]"  (propagation only on point #1)
struct Probe {
    QChar   dir = 'Z';
    double  x = 0, y = 0, z = 0;
    QString propagation;     // "+X" "-X" "+Y" ... (port #1 only)
};

} // namespace ofd
