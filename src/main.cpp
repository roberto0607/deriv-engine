#include "deriv-engine/vol_surface.hpp"
#include <iostream>
#include <fstream>

int main() {
    auto points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    std::ofstream out("data/vol_surface.csv");
    out << "strike,time_to_expiry,moneyness,type,market_iv,solved_iv,converged\n";
    for (const auto& p : points) {
        out << p.strike << "," << p.time_to_expiry << "," << p.moneyness << ","
            << (p.type == deriv::OptionType::Call ? "C" : "P") << ","
            << p.market_iv_reported << "," << p.solved_iv << "," << p.converged << "\n";
    }

    std::cerr << "Wrote " << points.size() << " points to data/vol_surface.csv\n";
    return 0;
}