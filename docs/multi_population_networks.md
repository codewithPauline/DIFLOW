# Multi-population migration networks

DIFLOW scales beyond a single population pair by separating two problems:

1. which population pairs should be evaluated,
2. how supported directional edges are summarized after inference.

## Candidate-pair construction

An all-vs-all analysis grows as

n(n-1)/2

pairwise comparisons.

For dozens of localities this quickly becomes expensive and biologically
implausible. DIFLOW therefore supports sparse candidate graphs based on:

- maximum geographic distance
- k-nearest geographic neighbors
- or a union of both

Future versions will add watershed, habitat, resistance, and user-supplied
adjacency constraints.

## Directed migration network

Supported DIFLOW edges become a directed graph.

Each edge can carry:

- migration magnitude
- directional support
- asymmetry index
- model weight
- evidence status

Ambiguous edges are excluded by default but can be retained for diagnostics.

## Source-like and sink-like summaries

For each population i:

M_out(i) = sum_j m(i -> j)

M_in(i) = sum_j m(j -> i)

Net(i) = M_out(i) - M_in(i)

DIFLOW labels positive Net values as source-like and negative Net values as
sink-like within the inferred migration network.

This is deliberately cautious terminology. The summary does not by itself
establish ecological source-sink demography, demographic growth rate, or
fitness differences.
