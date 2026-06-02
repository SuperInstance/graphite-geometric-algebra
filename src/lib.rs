//! Geometric algebra nodes for Graphite.
//!
//! Provides conformal geometric algebra (Cl(3,1)) operations:
//! rotors, reflections, projections, and conformal embeddings.
//!
//! # Architecture
//!
//! ```text
//! src/
//! ├── lib.rs              — Crate root; re-exports ga_core
//! ├── ga_core/
//! │   ├── mod.rs          — Module structure
//! │   ├── multivector.rs  — 16-component Cl(3,1) multivector with cached multiplication table
//! │   ├── rotor.rs        — Rotors: even-grade elements for rotation via sandwich product
//! │   └── conformal.rs    — Conformal GA: embedding, reflection, projection, conformal distance
//! └── ga_nodes.rs         — Standalone GA function implementations (feature-gated Graphite nodes)
//! ```
//!
//! # Feature Flags
//!
//! - `serde` (default): Serialization support for multivectors
//! - `graphite-nodes`: Enable Graphite `#[node_macro::node]` proc-macro annotations

pub mod ga_core;

/// Standalone geometric-algebra function implementations (not Graphite-node-specific).
///
/// When the `graphite-nodes` feature is enabled, also provides Graphite
/// `#[node_macro::node]` wrappers. See the module docs for details.
pub mod ga_nodes;

/// Re-export the core GA types for ergonomic use.
pub use ga_core::*;
