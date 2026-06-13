// Material.h — material record, 1:1 with the OpenFDTD "material =" line.
//
// 本家仕様 (OpenFDTD >= 2.2):
//   material = 1 εr σ μr σm            (通常媒質)
//   material = 2 ε∞ a b c              (分散性媒質, 一次Debye/Drude型)
// 材質番号: 0=真空(空気), 1=PEC は組込み。ユーザー定義は 2 から。
#pragma once
#include <QString>

namespace ofd {

struct Material {
    int     type = 1;      // 1 = normal, 2 = dispersive
    // type 1
    double  epsr = 1.0;    // relative permittivity
    double  esgm = 0.0;    // electric conductivity [S/m]
    double  amur = 1.0;    // relative permeability
    double  msgm = 0.0;    // magnetic conductivity
    // type 2 (dispersive: εr(ω) = ε∞ + a / (b + jωc))
    double  einf = 1.0;
    double  ae = 0.0, be = 0.0, ce = 0.0;

    QString name;          // GUI only (written as a trailing comment)

    // 拡張ドメイン用の物性 (acoustic / underwater) — .ofdx サイドカーに保存。
    // 音響カーネルは εr↔(ρ,c) のマッピングで同じ .ofd スキーマを使う。
    double  rho = 1.225;          // density [kg/m^3]
    double  soundSpeed = 343.0;   // [m/s]
    double  absorption = 0.02;    // absorption coefficient α [0..1]
};

// Lumped element — "load = dir x y z {R|L|C} value" line.
struct Load {
    QChar   dir = 'Z';            // X/Y/Z
    double  x = 0, y = 0, z = 0;
    QChar   kind = 'R';           // R/L/C
    double  value = 50.0;         // [Ω] / [H] / [F]
};

} // namespace ofd
