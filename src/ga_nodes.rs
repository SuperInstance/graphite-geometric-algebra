//! Graphite node definitions for geometric algebra operations.
//!
//! These nodes bring conformal Cl(3,1) geometric algebra into Graphite's
//! node-based editing graph, enabling rotor-based rotation, geometric
//! reflection through arbitrary planes, conformal projection, and chained
//! transformations via the geometric product.
//!
//! # Graphite Integration
//!
//! When the `graphite-nodes` feature is enabled, each function is annotated
//! with `#[node_macro::node]` to register as a Graphite node. Without it,
//! the functions compile as plain Rust functions usable in any context.

use crate::ga_core::Rotor;

// ---------------------------------------------------------------------------
// Standalone implementations — no Graphite dependency required.
// ---------------------------------------------------------------------------

/// Rotate a 2D vector using a geometric-algebra rotor.
///
/// A rotor R = cos(θ/2) - sin(θ/2)·B where B is the rotation-plane bivector.
/// Rotations compose cleanly: R_total = R2 · R1 (geometric product).
/// No gimbal lock, no matrix decomposition needed.
///
/// Returns `(x, y)` of the rotated vector.
pub fn rotor_rotate_2d(x: f64, y: f64, angle_degrees: f64, invert: bool) -> (f64, f64) {
    let angle = angle_degrees.to_radians();
    let effective_angle = if invert { -angle } else { angle };
    let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], effective_angle);
    let v = r.apply_2d([x, y]);
    (v[0], v[1])
}

/// Reflect a 2D point/vector through a line specified by its direction angle.
///
/// Reflection through a line with direction `L̂`:
/// v' = 2 (v · L̂) L̂ - v
pub fn reflect_through_line(x: f64, y: f64, line_angle_degrees: f64) -> (f64, f64) {
    let angle = line_angle_degrees.to_radians();
    let lx = angle.cos();
    let ly = angle.sin();
    let dot = x * lx + y * ly;
    (2.0 * dot * lx - x, 2.0 * dot * ly - y)
}

/// Reflect a 2D point/vector through a plane defined by its unit normal.
///
/// v' = v - 2 (v · n̂) n̂
pub fn reflect_by_normal(x: f64, y: f64, nx: f64, ny: f64) -> (f64, f64) {
    let len = (nx * nx + ny * ny).sqrt();
    if len < 1e-15 {
        return (-x, -y); // point reflection through origin
    }
    let nxi = nx / len;
    let nyi = ny / len;
    let dot = x * nxi + y * nyi;
    (x - 2.0 * dot * nxi, y - 2.0 * dot * nyi)
}

/// Project a 2D vector onto a direction axis.
///
/// v_parallel = (v · û) û
pub fn project_onto_axis(x: f64, y: f64, ax: f64, ay: f64) -> (f64, f64) {
    let len = (ax * ax + ay * ay).sqrt();
    if len < 1e-15 {
        return (0.0, 0.0);
    }
    let nxi = ax / len;
    let nyi = ay / len;
    let dot = x * nxi + y * nyi;
    (dot * nxi, dot * nyi)
}

// ---------------------------------------------------------------------------
// Graphite node macros (requires `graphite-nodes` feature + Graphite workspace)
// ---------------------------------------------------------------------------

#[cfg(feature = "graphite-nodes")]
mod graphite_wrapper {
    use core_types::Ctx;
    use glam::{DAffine2, DVec2};

    use crate::ga_core::{Multivector, Rotor};

    /// Rotate a 2D vector using a geometric-algebra rotor.
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn rotor_rotate(
        _: impl Ctx,
        vector: DVec2,
        #[default(45.0)] angle_degrees: f64,
        invert: bool,
    ) -> DVec2 {
        let (x, y) = super::rotor_rotate_2d(vector.x, vector.y, angle_degrees, invert);
        DVec2::new(x, y)
    }

    /// Reflect a 2D vector through a line specified by an angle.
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn reflector(_: impl Ctx, vector: DVec2, #[default(0.0)] line_angle_degrees: f64) -> DVec2 {
        let (x, y) = super::reflect_through_line(vector.x, vector.y, line_angle_degrees);
        DVec2::new(x, y)
    }

    /// Reflect a 2D vector through a plane defined by its normal.
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn reflect_by_normal(
        _: impl Ctx,
        vector: DVec2,
        #[default(1.0, 0.0)] normal: DVec2,
    ) -> DVec2 {
        let (x, y) = super::reflect_by_normal(vector.x, vector.y, normal.x, normal.y);
        DVec2::new(x, y)
    }

    /// Project a 2D vector onto a line (or plane in 2D: a line through origin).
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn project_onto_axis(
        _: impl Ctx,
        vector: DVec2,
        #[default(1.0, 0.0)] onto: DVec2,
    ) -> DVec2 {
        let (x, y) = super::project_onto_axis(vector.x, vector.y, onto.x, onto.y);
        DVec2::new(x, y)
    }

    /// Compose two geometric transformations using matrix multiplication.
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn compose_ga_transforms(_: impl Ctx, first: DAffine2, second: DAffine2) -> DAffine2 {
        first * second
    }

    /// Create a rotor from an axis angle (3D version for future 3D support).
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn make_rotor(
        _: impl Ctx,
        axis_x: f64,
        axis_y: f64,
        #[default(1.0)] axis_z: f64,
        angle_degrees: f64,
    ) -> String {
        let axis = [axis_x, axis_y, axis_z];
        let r = Rotor::from_axis_angle(axis, angle_degrees.to_radians());
        format!(
            "Rotor(scalar={:.4}, e01={:.4}, e02={:.4}, e12={:.4})",
            r.inner.c[0], r.inner.c[5], r.inner.c[6], r.inner.c[8]
        )
    }

    /// Compose two rotors: R_total = R1 * R2 (geometric product).
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn compose_rotors(
        _: impl Ctx,
        r1_scalar: f64,
        r1_e01: f64,
        r1_e02: f64,
        r1_e12: f64,
        r2_scalar: f64,
        r2_e01: f64,
        r2_e02: f64,
        r2_e12: f64,
    ) -> String {
        let mut m1 = crate::Multivector::zero();
        m1.c[0] = r1_scalar;
        m1.c[5] = r1_e01;
        m1.c[6] = r1_e02;
        m1.c[8] = r1_e12;
        let mut m2 = crate::Multivector::zero();
        m2.c[0] = r2_scalar;
        m2.c[5] = r2_e01;
        m2.c[6] = r2_e02;
        m2.c[8] = r2_e12;
        let r1 = Rotor::from_multivector(m1);
        let r2 = Rotor::from_multivector(m2);
        let composed = r1.compose(&r2);
        format!(
            "Composed: scalar={:.4}, e01={:.4}, e02={:.4}, e12={:.4}",
            composed.inner.c[0], composed.inner.c[5], composed.inner.c[6], composed.inner.c[8]
        )
    }

    /// Compute the geometric product of two multivectors.
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn geometric_product_debug(
        _: impl Ctx,
        a0: f64,
        a1: f64,
        a2: f64,
        a3: f64,
        a4: f64,
        a5: f64,
        a6: f64,
        a7: f64,
        a8: f64,
        a9: f64,
        a10: f64,
        a11: f64,
        a12: f64,
        a13: f64,
        a14: f64,
        a15: f64,
        b0: f64,
        b1: f64,
        b2: f64,
        b3: f64,
        b4: f64,
        b5: f64,
        b6: f64,
        b7: f64,
        b8: f64,
        b9: f64,
        b10: f64,
        b11: f64,
        b12: f64,
        b13: f64,
        b14: f64,
        b15: f64,
    ) -> String {
        let mut ma = crate::Multivector::zero();
        let mc = &mut ma.c;
        mc[0] = a0;
        mc[1] = a1;
        mc[2] = a2;
        mc[3] = a3;
        mc[4] = a4;
        mc[5] = a5;
        mc[6] = a6;
        mc[7] = a7;
        mc[8] = a8;
        mc[9] = a9;
        mc[10] = a10;
        mc[11] = a11;
        mc[12] = a12;
        mc[13] = a13;
        mc[14] = a14;
        mc[15] = a15;
        let mut mb = crate::Multivector::zero();
        let mc = &mut mb.c;
        mc[0] = b0;
        mc[1] = b1;
        mc[2] = b2;
        mc[3] = b3;
        mc[4] = b4;
        mc[5] = b5;
        mc[6] = b6;
        mc[7] = b7;
        mc[8] = b8;
        mc[9] = b9;
        mc[10] = b10;
        mc[11] = b11;
        mc[12] = b12;
        mc[13] = b13;
        mc[14] = b14;
        mc[15] = b15;
        let result = ma.geometric_product(&mb);
        format!("{:?}", result.c)
    }

    /// Create a DVec2 from x,y components (convenience node for GA operations).
    #[node_macro::node(category("Geometric Algebra"))]
    pub fn ga_vec2(_: impl Ctx, x: f64, y: f64) -> DVec2 {
        DVec2::new(x, y)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::f64::consts::FRAC_PI_2;

        #[test]
        fn test_rotor_rotate_90_degrees() {
            let result = rotor_rotate((), DVec2::new(1.0, 0.0), 90.0, false);
            assert!((result.x).abs() < 0.01, "x should be ≈0, got {}", result.x);
            assert!(
                (result.y - 1.0).abs() < 0.01,
                "y should be ≈1, got {}",
                result.y
            );
        }

        #[test]
        fn test_rotor_rotate_180_degrees() {
            let result = rotor_rotate((), DVec2::new(1.0, 0.0), 180.0, false);
            assert!(
                (result.x + 1.0).abs() < 0.01,
                "x should be ≈-1, got {}",
                result.x
            );
            assert!((result.y).abs() < 0.01, "y should be ≈0, got {}", result.y);
        }

        #[test]
        fn test_rotor_rotate_inverted() {
            let result = rotor_rotate((), DVec2::new(1.0, 0.0), 90.0, true);
            assert!((result.x).abs() < 0.01, "x should be ≈0, got {}", result.x);
            assert!(
                (result.y + 1.0).abs() < 0.01,
                "y should be ≈-1, got {}",
                result.y
            );
        }

        #[test]
        fn test_rotor_rotate_zero_angle() {
            let result = rotor_rotate((), DVec2::new(3.0, 4.0), 0.0, false);
            assert!((result.x - 3.0).abs() < 1e-10);
            assert!((result.y - 4.0).abs() < 1e-10);
        }

        #[test]
        fn test_reflector_45_degrees() {
            let result = reflector((), DVec2::new(1.0, 0.0), 22.5);
            let angle = result.y.atan2(result.x).to_degrees();
            assert!((angle - 45.0).abs() < 0.1, "expected 45°, got {angle}°");
        }

        #[test]
        fn test_reflect_by_normal_y_axis() {
            let result = reflect_by_normal((), DVec2::new(1.0, 1.0), DVec2::new(0.0, 1.0));
            assert!((result.x - 1.0).abs() < 1e-10);
            assert!((result.y + 1.0).abs() < 1e-10);
        }

        #[test]
        fn test_reflect_double_is_identity() {
            let v = DVec2::new(3.0, 4.0);
            let n = DVec2::new(1.0, 0.0);
            let r1 = reflect_by_normal((), v, n);
            let r2 = reflect_by_normal((), r1, n);
            assert!((r2.x - v.x).abs() < 1e-10);
            assert!((r2.y - v.y).abs() < 1e-10);
        }

        #[test]
        fn test_project_onto_x_axis() {
            let v = DVec2::new(3.0, 4.0);
            let result = project_onto_axis((), v, DVec2::new(1.0, 0.0));
            assert!((result.x - 3.0).abs() < 1e-10);
            assert!((result.y).abs() < 1e-10);
        }

        #[test]
        fn test_project_onto_diagonal() {
            let v = DVec2::new(1.0, 1.0);
            let result = project_onto_axis((), v, DVec2::new(1.0, 1.0));
            assert!((result.x - 1.0).abs() < 1e-10);
            assert!((result.y - 1.0).abs() < 1e-10);
        }

        #[test]
        fn test_compose_ga_transforms_identity() {
            let ident = DAffine2::IDENTITY;
            let result = compose_ga_transforms((), ident, ident);
            assert!(
                (result.transform_point2(DVec2::new(1.0, 2.0)) - DVec2::new(1.0, 2.0)).length()
                    < 1e-10
            );
        }

        #[test]
        fn test_compose_ga_transforms_translate_then_rotate() {
            let trans = DAffine2::from_translation(DVec2::new(1.0, 0.0));
            let rot = DAffine2::from_angle(FRAC_PI_2);
            let composed = compose_ga_transforms((), trans, rot);
            let result = composed.transform_point2(DVec2::new(0.0, 0.0));
            assert!((result.x - 1.0).abs() < 1e-10);
            assert!((result.y).abs() < 1e-10);
        }
    }
}
