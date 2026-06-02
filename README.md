# Graphite Geometric Algebra

**Cl(3,1) conformal geometric algebra nodes for the [Graphite](https://github.com/GraphiteEditor/Graphite) 2D vector editor (26K stars).**

This crate provides geometric algebra primitives — rotors, conformal embeddings, reflections, projections — as Graphite node definitions and standalone Rust functions. Built on [ga-core](https://github.com/SuperInstance/ga-core) for the Cl(3,1) spacetime algebra.

---

## Table of Contents

- [Why Geometric Algebra?](#why-geometric-algebra)
- [Quick Start](#quick-start)
- [Node Reference](#node-reference)
- [Code Examples](#code-examples)
- [Rotors vs Quaternions vs Matrices](#rotors-vs-quaternions-vs-matrices)
- [Why Cl(3,1) Conformal GA?](#why-cl31-conformal-ga)
- [Architecture](#architecture)
- [License](#license)

---

## Why Geometric Algebra?

**Geometric algebra (GA)** unifies vectors, rotations, reflections, and higher-dimensional primitives into a single algebraic system. Instead of juggling separate math for matrices, quaternions, cross products, and normals, you get one operation: the **geometric product**.

For a 2D graphics editor, GA unlocks:

| Capability | Traditional approach | GA approach |
|---|---|---|
| 2D rotation | `2×2` matrix or complex number | Rotor: 2 bivector components |
| 3D rotation | 4×4 matrix or quaternion | Rotor: 3 bivector components |
| Smooth interpolation | Quaternion slerp | Rotor slerp (identical math) |
| Reflection | Custom formula | Sandwich product |
| Points + spheres | Separate types | Conformal embedding → one algebra |

**The key insight:** In GA, rotation isn't a matrix applied to a vector. It's a *sandwich product*:

```
v' = R * v * R̃
```

where `R` is a **rotor** (a multivector living in even grades) and `R̃` is its reverse. This is simple, composable, and generalizes to any dimension.

---

## Quick Start

### Add to a standalone project

```toml
[dependencies]
graphite-ga-nodes = { git = "https://github.com/SuperInstance/graphite-geometric-algebra" }
ga-core = { git = "https://github.com/SuperInstance/ga-core" }
```

### Drop into Graphite

Clone or symlink this crate into Graphite's node graph directory:

```
graphite/
└── node-graph/
    └── nodes/
        ├── ga/                         # ← clone here
        │   ├── Cargo.toml
        │   └── src/
        │       ├── lib.rs
        │       ├── ga_nodes.rs
        │       └── ga_core/
        │           ├── mod.rs
        │           ├── multivector.rs
        │           ├── rotor.rs
        │           └── conformal.rs
        ├── artboard.rs
        ├── color.rs
        └── ...
```

When integrating into Graphite, replace the `[dependencies]` in `Cargo.toml` with workspace paths:

```toml
[dependencies]
core-types = { workspace = true }
vector-types = { workspace = true }
graphic-types = { workspace = true }
node-macro = { workspace = true }
glam = { workspace = true }
serde = { workspace = true, optional = true }
ga-core = { git = "https://github.com/SuperInstance/ga-core" }
```

### Run examples

```bash
# Rotor sandbox — compose, apply, slerp
cargo run --example rotor_sandbox

# Conformal geometry — embed, reflect, project, distance
cargo run --example conformal_geometry
```

---

## Node Reference

The `#[node_macro::node]` attributes in `src/ga_nodes.rs` register Graphite editor nodes. Each node appears in the editor's node-graph UI with typed inputs and outputs.

### `GA Rotor`

Builds a rotor from an axis-angle or bivector.

```
┌─────────────────────────────┐
│        GA Rotor             │
│                             │
│  axis_x    ────┐            │
│  axis_y    ────┤            │
│  axis_z    ────┤───► Rotor │
│  angle     ────┘            │
└─────────────────────────────┘
```

**Inputs:** `axis_x`, `axis_y`, `axis_z` (`f64`), `angle` (`f64`)  
**Output:** `Rotor` (multivector)  

Creates a rotor from axis-angle representation: `R = cos(θ/2) - sin(θ/2) * B` where `B` is the unit spatial bivector (e23, e31, e12 in Cl(3,1)).

### `GA Compose Rotors`

Composes two rotors into one.

```
┌─────────────────────────────────┐
│      GA Compose Rotors         │
│                                 │
│  rotor_a  ─────┐               │
│                 │──► Rotor_out │
│  rotor_b  ─────┘               │
└─────────────────────────────────┘
```

**Inputs:** `rotor_a`, `rotor_b` (both `Rotor`)  
**Output:** `Rotor` (single composed rotor)  

Composition: `Rout = Rb * Ra`. Applying `Rout` to a vector is equivalent to applying `Ra` then `Rb`.

### `GA Apply Rotor`

Applies a rotor to a 3D point via the sandwich product.

```
┌──────────────────────────────┐
│       GA Apply Rotor        │
│                              │
│  rotor   ───┐               │
│              │──► Point_out │
│  x  ────────┘               │
│  y  ────────────────► y_out │
│  z  ────────────────► z_out │
└──────────────────────────────┘
```

**Inputs:** `rotor` (`Rotor`), `x`, `y`, `z` (`f64`)  
**Outputs:** `x_out`, `y_out`, `z_out` (`f64`)  

Computes: `v' = R * v * R̃`.

### `GA Rotor to Matrix`

Extracts a 3×3 rotation matrix from a rotor.

```
┌──────────────────────────────┐
│     GA Rotor to Matrix      │
│                              │
│  rotor  ─────► Matrix (3×3) │
└──────────────────────────────┘
```

**Input:** `rotor` (`Rotor`)  
**Output:** `Matrix3x3`  

Converts the rotor to a standard rotation matrix by applying it to the basis vectors `{e1, e2, e3}`.

### `GA Embed Point`

Embeds a Euclidean 3D point into Cl(3,1) conformal space.

```
┌────────────────────────────────────┐
│        GA Embed Point             │
│                                    │
│  x  ──┐                            │
│  y  ──┤──► Conformal (Multivector) │
│  z  ──┘                            │
└────────────────────────────────────┘
```

**Inputs:** `x`, `y`, `z` (`f64`)  
**Output:** `Multivector` (conformal point)  

Embedding: `P = e₊ + x·e₁ + y·e₂ + z·e₃ + ½||p||²` (simplified form).

### `GA Distance`

Computes Euclidean distance between two points.

```
┌──────────────────────────────┐
│       GA Distance           │
│                              │
│  x1 ──┐                      │
│  y1 ──┤                      │
│  z1 ──┤──► Distance (f64)   │
│  x2 ──┤                      │
│  y2 ──┤                      │
│  z2 ──┘                      │
└──────────────────────────────┘
```

**Inputs:** `x1`, `y1`, `z1`, `x2`, `y2`, `z2` (`f64`)  
**Output:** `f64` (straight-line Euclidean distance)

### `GA Reflect`

Reflects a point through a plane.

```
┌─────────────────────────────────┐
│        GA Reflect              │
│                                 │
│  px  ──┐                        │
│  py  ──┤                        │
│  pz  ──┤──► x_out              │
│  nx  ──┤     y_out             │
│  ny  ──┤     z_out             │
│  nz  ──┘                        │
│  d   ──┘                        │
└─────────────────────────────────┘
```

**Inputs:** `px`, `py`, `pz` (point), `nx`, `ny`, `nz` (plane normal), `d` (plane distance)  
**Outputs:** `x_out`, `y_out`, `z_out` (`f64`)

Reflection formula: `p' = p - 2·(n·p - d)·n`.

### `GA Project`

Projects a point onto a plane.

```
┌─────────────────────────────────┐
│        GA Project              │
│                                 │
│  px  ──┐                        │
│  py  ──┤                        │
│  pz  ──┤──► x_out              │
│  nx  ──┤     y_out             │
│  ny  ──┤     z_out             │
│  nz  ──┘                        │
│  d   ──┘                        │
└─────────────────────────────────┘
```

**Inputs:** `px`, `py`, `pz` (point), `nx`, `ny`, `nz` (plane normal), `d` (plane distance)  
**Outputs:** `x_out`, `y_out`, `z_out` (`f64`)

Projection formula: `p' = p - (n·p - d)·n`.

### `GA Rotate Point`

Rotates a point using a rotor — convenience wrapper.

```
┌──────────────────────────────┐
│     GA Rotate Point         │
│                              │
│  rotor  ──┐                  │
│  x  ──────┤──► x_out        │
│  y  ──────┤     y_out       │
│  z  ──────┘     z_out       │
└──────────────────────────────┘
```

**Inputs:** `rotor` (`Rotor`), `x`, `y`, `z` (`f64`)  
**Outputs:** `x_out`, `y_out`, `z_out` (`f64`)

### `GA Translate`

Translates a point by an offset.

```
┌──────────────────────────────┐
│     GA Translate            │
│                              │
│  px  ──┐                     │
│  py  ──┤──► x_out           │
│  pz  ──┤     y_out          │
│  tx  ──┤     z_out          │
│  ty  ──┤                     │
│  tz  ──┘                     │
└──────────────────────────────┘
```

**Inputs:** `px`, `py`, `pz` (point), `tx`, `ty`, `tz` (translation)  
**Outputs:** `x_out`, `y_out`, `z_out` (`f64`)

Translation: `p' = p + t`.

---

## Code Examples

### Standalone (Rust)

```rust
use graphite_ga_nodes::ga_core::Rotor;

// Create a 90° Z-axis rotor
let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);

// Apply to a point
let rotated = r.apply([1.0, 0.0, 0.0]);
assert!(rotated[0].abs() < 0.01);  // ≈ 0
assert!((rotated[1] - 1.0).abs() < 0.01);  // ≈ 1

// Compose rotors
let r45 = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4);
let r90 = r45.compose(&r45);

// Slerp between identity and a rotation
let id = Rotor::identity();
let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.0);
let mid = id.slerp(&r, 0.5);  // half-way rotation

// Extract rotation matrix
let matrix = r90.to_rotation_matrix();
```

### Conformal geometry

```rust
use graphite_ga_nodes::ga_core::{Conformal, Rotor};

// Embed a 3D point
let p = Conformal::embed_point([3.0, 4.0, 0.0]);
let extracted = Conformal::extract_point(&p);  // back to [3, 4, 0]

// Distance between points
let d = Conformal::conformal_distance(
    &Conformal::embed_point([0.0, 0.0, 0.0]),
    &Conformal::embed_point([3.0, 4.0, 0.0]),
);
assert!((d - 5.0).abs() < 0.001);

// Reflect through a plane
let reflected = Conformal::reflect([1.0, 2.0, 3.0], [1.0, 0.0, 0.0], 0.0);
// → (-1, 2, 3) (mirror through YZ-plane)

// Project onto a plane
let projected = Conformal::project_onto_plane(
    [0.0, 0.0, 5.0],  // point
    [0.0, 0.0, 1.0],  // normal
    2.0,               // plane offset
);
// → (0, 0, 2)

// Rotate a point
let rotor = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
let rotated = Conformal::rotate([1.0, 0.0, 0.0], &rotor);
// → (0, 1, 0)

// Translate
let translated = Conformal::translate([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
// → (5, 7, 9)
```

---

## Rotors vs Quaternions vs Matrices

| Property | 3×3 Matrix | Quaternion | Rotor |
|---|---|---|---|
| **Components** | 9 | 4 | 4 (scalar + 3 bivector) |
| **Smooth interpolation** | Difficult (QR decomposition) | Slerp (native) | Slerp (identical) |
| **2D specialization** | Awkward | Needs z = 0 | `apply_2d()` |
| **Generalizes to N-D** | No (3×3 only) | No (3D only) | Yes (any dimension) |
| **Reflection** | Separate API | Separate API | Same sandwich product |
| **Geometric meaning** | Opaque | Half-angle double-cover | Exponential of bivector |
| **Memory** | 72 bytes (f64) | 32 bytes | 128 bytes (Multivector) |

**When should you use rotors?**

- In a GA-native codebase (this crate)
- When you need a unified algebra for rotations, reflections, and projections
- When working in 2D with the option to extend to 3D later
- When composing many transformations with numerical stability

**When should you stick with matrices?**  

- Interoperability with existing rendering pipelines (OpenGL, WGPU)
- Affine transforms with non-uniform scaling and shears
- Sending to GPU shaders (which expect mat3/mat4)

---

## Why Cl(3,1) Conformal GA?

Conformal geometric algebra (CGA) adds two extra dimensions to Euclidean space, creating a 5D algebra Cl(4,1) — or equivalently Cl(3,1) if we use a spacetime-like signature.

**What CGA enables:**

| Geometric object | CGA representation |
|---|---|
| Point | `P = e₊ + p + ½||p||² e₋` |
| Plane | `n + d·e₊` (normal + offset) |
| Sphere | Through 4 points |
| Circle | Intersection of two spheres |
| Line | Through 2 points |

**All with the same sandwich product.** In standard 3D GA:
- Rotors rotate (grade-2 bivectors in spatial subspace)
- Reflectors reflect through planes

In CGA:
- Rotors still rotate (spatial subspace)
- **Translators** translate (e₊ ∧ e₋ bivectors)
- **Motors** combine rotation and translation (screw motions, SE(3))
- Inversions about spheres are also simple reflections

For Graphite's 2D + layer depth use case, the full CGA machinery may be overkill, but the *principle* — one algebra, one operation, many geometries — is the foundation.

> **Note:** This crate uses a **simplified** conformal embedding. The ½||p||² term is stored as a scalar component rather than in a true null-vector basis. This enables round-trip embed/extract and Euclidean distance computation, but true CGA operations (inversion in sphere, dual spheres, intersections) require a proper Cl(4,1) implementation.

---

## Architecture

```
graphite-geometric-algebra/
├── Cargo.toml
├── .github/workflows/
│   └── ci.yml                     # CI: test + clippy on Ubuntu
├── examples/
│   ├── rotor_sandbox.rs           # 2D/3D rotor sandbox
│   └── conformal_geometry.rs      # Conformal point operations
├── src/
│   ├── lib.rs                     # Crate root, module exports
│   ├── ga_nodes.rs                # Graphite #[node] definitions
│   └── ga_core/
│       ├── mod.rs                 # Public re-exports
│       ├── multivector.rs         # Cl(3,1) Multivector (16 components)
│       │                           # - Geometric product (cached table)
│       │                           # - Wedge, inner, dual, reverse, conjugate
│       │                           # - Grade extraction
│       ├── rotor.rs                # Rotor: from-axis-angle, from-bivector
│       │                           # - Apply (sandwich), compose, slerp
│       │                           # - to_rotation_matrix, apply_2d
│       │                           # - Normalization, reflection
│       └── conformal.rs           # Conformal embedding, distance, reflect,
│                                    # project, rotate, translate, barycenter
```

### Dependencies

- **[ga-core](https://github.com/SuperInstance/ga-core)** — Cl(3,1) geometric algebra core (the product table, blade enumeration, metric)
- **[glam](https://github.com/bitshifter/glam-rs)** — Vec types for Graphite interface (f64 precision)
- **serde** (optional) — Serialization support
- **node-macro** (Graphite workspace) — When integrated into Graphite

---

## License

Licensed under the [MIT License](LICENSE).

---

*Built for [Graphite](https://github.com/GraphiteEditor/Graphite) — the open-source 2D vector editor reimagined with geometric algebra.*
