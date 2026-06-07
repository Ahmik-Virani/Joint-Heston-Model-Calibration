#pragma once
#include "qe/params.hpp"
#include "qe/path.hpp"
#include "qe/surface.hpp"
#include <vector>

namespace qe {
    using QuantLib::Rate;
    
    void print_Surfaces(const std::vector<CallGrid>& surfaces);

    std::vector<CallGrid> build_surfaces_from_path(
        const HestonQParams& Q,
        const QuantLib::Date& today,
        const QuantLib::Calendar& cal,
        const PPath& ppath,
        int SurfaceFrequency,Rate r, Rate q);
    
} // namespace qe