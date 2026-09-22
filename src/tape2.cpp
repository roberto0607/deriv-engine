#include "deriv-engine/tape2.hpp"

namespace deriv {

std::size_t Tape2::push_leaf() {
    nodes.push_back(Node{});
    return nodes.size() - 1;
}

std::size_t Tape2::push_unary(std::size_t parent, Dual partial) {
    Node n;
    n.parent0 = parent;
    n.partial0 = partial;
    nodes.push_back(n);
    return nodes.size() - 1;
}

std::size_t Tape2::push_binary(std::size_t p0, Dual d0, std::size_t p1, Dual d1) {
    Node n;
    n.parent0 = p0;
    n.partial0 = d0;
    n.parent1 = p1;
    n.partial1 = d1;
    nodes.push_back(n);
    return nodes.size() - 1;
}

void Tape2::backward(std::size_t root) {
    adjoints.assign(nodes.size(), Dual(0.0, 0.0));
    adjoints[root] = Dual(1.0, 0.0);  // seed: d(root)/d(root) = 1, a constant (dot=0)

    // Identical walk order to Tape::backward, for the identical reason:
    // every node's parents were recorded before the node itself, so
    // walking last-to-first guarantees each node's adjoint is fully
    // accumulated before it's propagated further back. The only
    // difference from Tape::backward is that `a * n.partial0` here is a
    // Dual*Dual multiply (dual.hpp's operator*, product rule), not a
    // plain double multiply -- that's what carries the second-
    // derivative information through the chain rule alongside the
    // first-derivative information, in the same single backward pass.
    for (std::size_t i = nodes.size(); i-- > 0;) {
        Dual a = adjoints[i];
        if (a.val == 0.0 && a.dot == 0.0) continue;

        const Node& n = nodes[i];
        if (n.parent0 != NO_PARENT2) {
            adjoints[n.parent0] = adjoints[n.parent0] + a * n.partial0;
        }
        if (n.parent1 != NO_PARENT2) {
            adjoints[n.parent1] = adjoints[n.parent1] + a * n.partial1;
        }
    }
}

void Tape2::clear() {
    nodes.clear();
    adjoints.clear();
}

Tape2& get_tape2() {
    thread_local Tape2 tape;
    return tape;
}

}  // namespace deriv
