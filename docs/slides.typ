#set text(size:28pt)
#set block()

== Questions

- What I want to do 
  - to build a graph from a given manifold
- What graph I want to build
  - a Graph consists of simplexes as vertex and edges
  - the graph isomorphic to the given manifold
  - the capability to subdivide the graph to get a finer graph 
- What functions I want to implement
  - `search`: search the graph with a given start vertex and a given goal.
    - `goal` can be a vertex, an edge, or a subgraph.
      - the goal can be set with a series of conditions.
    - `alg` supports various algorithms including `dfs`, `bfs`, etc. 
    - `advisor` ranks given vertices(e.g. the neighbors of a vertex).
    - `hook` will be called before or after a vertex is visited.
      - `hook` can be used to implement `advisor`.
  - `subdivide`
  - `nvertex`
  - `nedge`

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
