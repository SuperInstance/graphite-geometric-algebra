# Graphite Geometric Algebra Nodes

**Geometric algebra nodes for the [Graphite editor](https://graphite.art) — rotors, conformal embedding, reflections, and projections.**

This crate provides a set of node types that bring conformal geometric algebra (Cl(3,1)) into Graphite's node-based editing graph. It uses [ga-core](https://github.com/SuperInstance/ga-core) for the geometric algebra math library.

## What This Provides

- **Rotor rotation** — clean, gimbal-lock-free 2D rotations using geometric algebra rotors. No matrices needed.
- **Reflection** — geometric reflection through arbitrary mirror lines/planes using the fundamental GA formula `v' = v - 2(n·v)n`
- **Conformal projection** — grade-lowering projection of vectors onto axes
- **Transform composition** — compose rotations/reflections via the geometric product
- **Debug tools** — inspect multivectors, compute geometric products, and compose rotors numerically
- **Rotor module** — create, normalize, compose, interpolate (slerp), and extract rotation matrices from rotors
- **Conformal module** — embed Euclidean 3D points into 5D conformal space, reflect through planes, project, compute conformal distances

## Architecture

```
src/
├── lib.rs              — crate root; re-exports ga_core + ga_nodes
├── ga_core/
│   ├── mod.rs          — module structure
│   ├── multivector.rs  — 16-component Cl(3,1) multivector with cached multiplication table
│   ├── rotor.rs        — Rotors: even-grade elements for rotation via sandwich product
│   └── conformal.rs    — Conformal GA: embedding, reflection, projection, conformal distance
└── ga_nodes.rs         — Graphite node definitions using #[node_macro::node]
```

## Usage in Graphite

To integrate this crate into the Graphite editor:

1. Clone or copy this crate into `graphite/node-graph/nodes/`:

```bash
git clone https://github.com/SuperInstance/graphite-geometric-algebra graphite/node-graph/nodes/ga
```

2. Replace `Cargo.toml` dependencies with Graphite's workspace deps (see the comments in `Cargo.toml`).

3. The nodes will register under the **"Geometric Algebra"** category in Graphite's node library.

The node definitions in `src/ga_nodes.rs` use the `#[node_macro::node]` proc-macro from Graphite's `node-macro` crate. This is what registers each function as a node in the editor graph with typed inputs, outputs, and default values.

## GA Core

The `ga_core` module implements a full 16-component multivector for Cl(3,1) with:

- Complete cached 16×16 geometric product multiplication table
- Grade extraction (scalar, vector, bivector, trivector, pseudoscalar)
- Reverse, Clifford conjugate, Hodge dual
- Norm squared (via reverse product)
- Wedge and inner products
- Exponential-map rotor construction from bivectors
- Rotor composition and sandwich-product rotation

### Metric

The spacetime algebra uses the metric:
- e₀² = +1 (timelike in this convention)
- e₁² = +1
- e₂² = +1
- e₃² = -1

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache 2.0](LICENSE-APACHE) at your option.
