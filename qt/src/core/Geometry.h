// Geometry.h — one geometry unit, 1:1 with the OpenFDTD "geometry =" line.
//
// 本家の形状コード (sol/ingeometry.c):
//   1            直方体      g[0..5] = X1 X2 Y1 Y2 Z1 Z2
//   2            楕円体      g[0..5] = 外接直方体
//   11/12/13     円柱 X/Y/Z  g[0..5] = 外接直方体
//   31/32/33     三角柱      g[0..7]
//   41/42/43     角錐台      g[0..7]
//   51/52/53     円錐台      g[0..7]
// 重なった領域は後のユニットが優先 (ユニット番号 = リスト順)。
#pragma once
#include <QString>

namespace ofd {

enum class Axis { X, Y, Z };

struct Geometry {
    int     materialId = 2;   // index into the material table (0=air, 1=PEC)
    int     shape = 1;        // 本家 shape code
    double  g[8] = {0,0,0,0,0,0,0,0};
    QString name;             // GUI only (emitted as "name = " line, kernel ignores it)

    static int paramCount(int shape) {
        switch (shape) {
            case 1: case 2: case 11: case 12: case 13:
                return 6;
            case 31: case 32: case 33:
            case 41: case 42: case 43:
            case 51: case 52: case 53:
                return 8;
        }
        return 6;
    }
};

} // namespace ofd
