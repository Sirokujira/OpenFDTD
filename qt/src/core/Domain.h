// Domain.h — physics domain enumeration shared by all tabs
#pragma once
#include <QString>

namespace ofd {

enum class Domain {
    EM,          // electromagnetic (the original OpenFDTD)
    Optical,     // photonic / nm wavelengths
    Acoustic,    // room acoustics
    Underwater   // underwater (SOFAR, sonar)
};

inline QString domainKey(Domain d) {
    switch (d) {
        case Domain::EM:         return "em";
        case Domain::Optical:    return "optical";
        case Domain::Acoustic:   return "acoustic";
        case Domain::Underwater: return "underwater";
    }
    return "em";
}

inline Domain domainFromKey(const QString &key) {
    if (key == "optical")    return Domain::Optical;
    if (key == "acoustic")   return Domain::Acoustic;
    if (key == "underwater") return Domain::Underwater;
    return Domain::EM;
}

inline QString domainLabelKey(Domain d) {
    switch (d) {
        case Domain::EM:         return "d_em";
        case Domain::Optical:    return "d_optical";
        case Domain::Acoustic:   return "d_acoustic";
        case Domain::Underwater: return "d_underwater";
    }
    return "d_em";
}

inline QString accentColor(Domain d) {
    switch (d) {
        case Domain::EM:         return "#0078D4";
        case Domain::Optical:    return "#B83280";
        case Domain::Acoustic:   return "#2E8B57";
        case Domain::Underwater: return "#1E6FBF";
    }
    return "#0078D4";
}

} // namespace ofd
