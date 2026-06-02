//! Geometric algebra nodes for Graphite.
//!
//! Provides node types for conformal geometric algebra (Cl(3,1)):
//! Rotor, Reflect, Project, and Compose operations.

mod ga_core;

pub mod ga_nodes;
pub use ga_core::*;
