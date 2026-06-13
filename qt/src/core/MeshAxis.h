// MeshAxis.h — non-uniform mesh for one axis,
// 1:1 with the OpenFDTD "xmesh = x0 d1 x1 d2 x2 ..." line.
//
// nodes[i] は区間境界座標 (昇順)、divs[i] は nodes[i]→nodes[i+1] の分割数。
// divs.size() == nodes.size() - 1。
#pragma once
#include <QVector>
#include <algorithm>

namespace ofd {

struct MeshAxis {
    QVector<double> nodes;
    QVector<int>    divs;

    int totalCells() const {
        int n = 0;
        for (int d : divs) n += d;
        return n;
    }
    double minSpacing() const {
        double m = 1e308;
        for (int i = 0; i < divs.size(); ++i) {
            if (divs[i] > 0)
                m = std::min(m, (nodes[i+1] - nodes[i]) / divs[i]);
        }
        return m;
    }
    double min() const { return nodes.isEmpty() ? 0.0 : nodes.first(); }
    double max() const { return nodes.isEmpty() ? 0.0 : nodes.last(); }

    bool isValid() const {
        if (nodes.size() < 2 || divs.size() != nodes.size() - 1) return false;
        for (int i = 0; i < divs.size(); ++i)
            if (nodes[i] >= nodes[i+1] || divs[i] <= 0) return false;
        return true;
    }
};

} // namespace ofd
