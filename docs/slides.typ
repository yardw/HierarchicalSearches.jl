#set text(size:28pt)
#set block()
== Components

=== Manifold or SCG(Simplitical Complex Graph)
A *manifold* is a topological space that locally resembles Euclidean space near each point.
Naively, one can think of a n-dimentional manifold as a n-dimentional polyhedron in the limit of infinite subdivision.

$
n "-dim Manifold" M \
equiv 
{c^n (P) times g(P)| forall P in M}\
equiv
{"neighbors of" c(P)| forall c(P) in "SCG"}
$

- `simplitical_subdivide`: 
  - perform simplitical subdivision on the given manifold or simplitical complex graph.
  - return the isomorphic complex in form of graph. 
    - `precition_goal`:
  - TODO: `triangulation` implication
    - `triangulation` is an specific algorithm of `simplitical_subdivision`
    - `enmeshment`
      - if there is even number of vertex in the SCG, then 
  - TODO:
- `nvertex`
- `nedge`

==== Complex Graph
===== Grid Complex Graph
==== Vertex
===== Grid Vertex

==== Edge
===== Grid Edge
