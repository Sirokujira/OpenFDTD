// H5Writer.cpp
#include "H5Writer.h"
#include "../core/Project.h"

#ifdef OFD_USE_HDF5
#include <hdf5.h>
#endif

using namespace ofd;

bool H5Writer::available()
{
#ifdef OFD_USE_HDF5
    return true;
#else
    return false;
#endif
}

#ifdef OFD_USE_HDF5

static bool writeDoubleArray(hid_t file, const char *name,
                             const QVector<double> &data)
{
    const hsize_t dims[1] = { hsize_t(data.size()) };
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dset = H5Dcreate2(file, name, H5T_NATIVE_DOUBLE, space,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    const herr_t st = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                               H5P_DEFAULT, data.constData());
    H5Dclose(dset);
    H5Sclose(space);
    return st >= 0;
}

bool H5Writer::write(const QString &path, const Project &project,
                     const QVector<int> &steps,
                     const QVector<double> &eAvg,
                     const QVector<double> &hAvg,
                     QString *err)
{
    hid_t file = H5Fcreate(path.toLocal8Bit().constData(),
                           H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file < 0) {
        if (err) *err = "cannot create " + path;
        return false;
    }

    bool ok = true;
    // /mesh/{x,y,z}_nodes
    H5Gcreate2(file, "/mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    static const char *names[3] = { "/mesh/x_nodes", "/mesh/y_nodes", "/mesh/z_nodes" };
    for (int a = 0; a < 3; ++a)
        ok = ok && writeDoubleArray(file, names[a], project.mesh(a).nodes);

    // /convergence/{step,e_avg,h_avg}
    H5Gcreate2(file, "/convergence", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    QVector<double> stepsD;
    stepsD.reserve(steps.size());
    for (int s : steps) stepsD.push_back(double(s));
    ok = ok && writeDoubleArray(file, "/convergence/step", stepsD);
    ok = ok && writeDoubleArray(file, "/convergence/e_avg", eAvg);
    ok = ok && writeDoubleArray(file, "/convergence/h_avg", hAvg);

    H5Fclose(file);
    if (!ok && err) *err = "HDF5 write failed";
    return ok;
}

#else

bool H5Writer::write(const QString &, const Project &,
                     const QVector<int> &, const QVector<double> &,
                     const QVector<double> &, QString *err)
{
    if (err) *err = "built without HDF5 — reconfigure with -DUSE_HDF5=ON";
    return false;
}

#endif
