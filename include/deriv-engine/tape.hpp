#pragma once

#include <vector>
#include <cstddef>
#include <limits>

namespace deriv {

constexpr std::size_t NO_PARENT = std::numeric_limits<std::size_t>::max();

class Tape {
public:
    struct Node {
        std::size_t parent0 = NO_PARENT;
        double partial0 = 0.0;
        std::size_t parent1 = NO_PARENT;
        double partial1 = 0.0;
    };

    std::vector<Node> nodes;
    std::vector<double> adjoints;

    std::size_t push_leaf();
    std::size_t push_unary(std::size_t parent, double partial);
    std::size_t push_binary(std::size_t p0, double d0, std::size_t p1, double d1);

    // Walks the tape backward from `root`, filling `adjoints` with
    // d(root)/d(node) for every node on the tape.
    void backward(std::size_t root);

    void clear();
};

// One tape, shared by all ADouble values in the current thread. Using
// thread_local (rather than a single global) means each thread gets its
// own independent tape — necessary groundwork for the multithreaded
// Monte Carlo work later, even though this project is single-threaded
// for now.
Tape& get_tape();

}  // namespace deriv