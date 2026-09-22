#pragma once

#include "dual.hpp"
#include <vector>
#include <cstddef>
#include <limits>

namespace deriv {

constexpr std::size_t NO_PARENT2 = std::numeric_limits<std::size_t>::max();

// Second-order sibling of Tape (tape.hpp): identical structure and
// identical backward()-walk logic, but every value on it is a Dual
// (dual.hpp) instead of a plain double.
//
// Why this gets its own Tape/ADouble pair instead of templating the
// original: Phase 4 deliberately deferred gamma, on the grounds that it
// needs "a tape-of-tapes structure to differentiate the backward pass
// itself" -- a real, separate piece of machinery, not a small tweak to
// the existing one. Duplicating tape.hpp/adouble.hpp here with double
// swapped for Dual keeps that new machinery fully isolated: the
// original Tape/ADouble (and every Greek computed through them —
// delta, vega, theta, rho, already tested and shipped) is completely
// untouched by this file. If Tape2 has a bug, it cannot silently
// change an already-validated first-order Greek.
//
// The technique itself is "forward-over-reverse" AD, and it's exact
// (not a finite-difference approximation): running reverse-mode
// backward() with Dual-valued partials computes, in one pass, both the
// ordinary adjoint (d(root)/d(node), a first derivative) AND that
// adjoint's own forward-seeded directional derivative (a second
// derivative) — see the .dot field on each entry in `adjoints` after
// calling backward().
class Tape2 {
public:
    struct Node {
        std::size_t parent0 = NO_PARENT2;
        Dual partial0;
        std::size_t parent1 = NO_PARENT2;
        Dual partial1;
    };

    std::vector<Node> nodes;
    std::vector<Dual> adjoints;

    std::size_t push_leaf();
    std::size_t push_unary(std::size_t parent, Dual partial);
    std::size_t push_binary(std::size_t p0, Dual d0, std::size_t p1, Dual d1);

    void backward(std::size_t root);
    void clear();
};

Tape2& get_tape2();

}  // namespace deriv
