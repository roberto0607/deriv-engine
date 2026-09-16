#include "deriv-engine/tape.hpp"

namespace deriv {

std::size_t Tape::push_leaf() {
    nodes.push_back(Node{});
    return nodes.size() - 1;
}

std::size_t Tape::push_unary(std::size_t parent, double partial) {
    Node n;
    n.parent0 = parent;
    n.partial0 = partial;
    nodes.push_back(n);
    return nodes.size() - 1;
}

std::size_t Tape::push_binary(std::size_t p0, double d0, std::size_t p1, double d1) {
    Node n;
    n.parent0 = p0;
    n.partial0 = d0;
    n.parent1 = p1;
    n.partial1 = d1;
    nodes.push_back(n);
    return nodes.size() - 1;
}

void Tape::backward(std::size_t root) {
    adjoints.assign(nodes.size(), 0.0);
    adjoints[root] = 1.0;  // seed: d(root)/d(root) = 1

    // Walk from the last-recorded node back to the first. Because every
    // node's parents were necessarily recorded BEFORE the node itself,
    // walking backward guarantees a node's own adjoint is fully
    // accumulated before we propagate it further back.
    for (std::size_t i = nodes.size(); i-- > 0;) {
        double a = adjoints[i];
        if (a == 0.0) continue;  // nothing to propagate

        const Node& n = nodes[i];
        if (n.parent0 != NO_PARENT) {
            adjoints[n.parent0] += a * n.partial0;
        }
        if (n.parent1 != NO_PARENT) {
            adjoints[n.parent1] += a * n.partial1;
        }
    }
}

void Tape::clear() {
    nodes.clear();
    adjoints.clear();
}

Tape& get_tape() {
    thread_local Tape tape;
    return tape;
}

}  // namespace deriv