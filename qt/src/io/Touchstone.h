// Touchstone.h — write S-parameters in Touchstone format (.s1p / .s2p).
//
// 周波数特性プロット (plotspara) の共通出力形式。光・電磁ドメイン共通で
// 他ツール (ADS / scikit-rf / Lumerical INTERCONNECT) に持ち込める。
#pragma once
#include <QString>
#include <QVector>
#include <complex>

namespace ofd {

class Touchstone {
public:
    // 1-port: S11 only.
    static bool writeS1p(const QString &path,
                         const QVector<double> &freqHz,
                         const QVector<std::complex<double>> &s11,
                         QString *err = nullptr);

    // 2-port: S11 S21 S12 S22.
    static bool writeS2p(const QString &path,
                         const QVector<double> &freqHz,
                         const QVector<std::complex<double>> &s11,
                         const QVector<std::complex<double>> &s21,
                         const QVector<std::complex<double>> &s12,
                         const QVector<std::complex<double>> &s22,
                         QString *err = nullptr);

    // Convert input impedance to S11 against reference Z0.
    static std::complex<double> zToS(std::complex<double> z, double z0 = 50.0);
};

} // namespace ofd
