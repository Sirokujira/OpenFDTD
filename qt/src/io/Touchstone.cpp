// Touchstone.cpp
#include "Touchstone.h"

#include <QFile>
#include <QTextStream>

using namespace ofd;

static void writeHeader(QTextStream &out)
{
    out << "! OpenFDTD-X export\n";
    out << "# Hz S RI R 50\n";
}

bool Touchstone::writeS1p(const QString &path,
                          const QVector<double> &freqHz,
                          const QVector<std::complex<double>> &s11,
                          QString *err)
{
    if (freqHz.size() != s11.size()) {
        if (err) *err = "frequency/data size mismatch";
        return false;
    }
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    QTextStream out(&f);
    writeHeader(out);
    for (int i = 0; i < freqHz.size(); ++i)
        out << QString::number(freqHz[i], 'e', 9) << ' '
            << QString::number(s11[i].real(), 'e', 9) << ' '
            << QString::number(s11[i].imag(), 'e', 9) << '\n';
    return true;
}

bool Touchstone::writeS2p(const QString &path,
                          const QVector<double> &freqHz,
                          const QVector<std::complex<double>> &s11,
                          const QVector<std::complex<double>> &s21,
                          const QVector<std::complex<double>> &s12,
                          const QVector<std::complex<double>> &s22,
                          QString *err)
{
    const int n = freqHz.size();
    if (s11.size() != n || s21.size() != n || s12.size() != n || s22.size() != n) {
        if (err) *err = "frequency/data size mismatch";
        return false;
    }
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (err) *err = f.errorString();
        return false;
    }
    QTextStream out(&f);
    writeHeader(out);
    for (int i = 0; i < n; ++i) {
        out << QString::number(freqHz[i], 'e', 9);
        for (const auto *s : { &s11, &s21, &s12, &s22 })
            out << ' ' << QString::number((*s)[i].real(), 'e', 9)
                << ' ' << QString::number((*s)[i].imag(), 'e', 9);
        out << '\n';
    }
    return true;
}

std::complex<double> Touchstone::zToS(std::complex<double> z, double z0)
{
    return (z - z0) / (z + z0);
}
