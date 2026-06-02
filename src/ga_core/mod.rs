//! GA Core — Conformal geometric algebra using Cl(3,1) spacetime algebra.
//!
//! Embedded in Graphite's geometric algebra node extension.

mod conformal;
mod multivector;
mod rotor;

pub use conformal::Conformal;
pub use multivector::Multivector;
pub use rotor::Rotor;
