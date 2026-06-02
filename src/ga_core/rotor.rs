//! Rotors: even-grade multivectors for rotations in geometric algebra.

use super::multivector::Multivector;

/// A rotor: an even-grade element that performs rotations.
///
/// A rotor R satisfies R * R̃ = 1 (normalized).
/// Rotations are applied via the sandwich product: v' = R v R̃.
///
/// Convention: R = exp(-Bθ/2) where B is the rotation-plane bivector.
/// Uses spatial bivectors e23 (X), e31 (Y), e12 (Z) in Cl(3,1).
#[derive(Debug, Clone)]
pub struct Rotor {
    pub inner: Multivector,
}

impl Rotor {
    /// Create from a multivector (will normalize).
    pub fn from_multivector(m: Multivector) -> Self {
        let mut r = Self { inner: m };
        r.normalize();
        r
    }

    /// Identity rotor.
    pub fn identity() -> Self {
        Self {
            inner: Multivector::scalar(1.0),
        }
    }

    /// Create a rotor from axis-angle representation.
    ///
    /// Uses the spatial bivectors of Cl(3,1): e23 (X-rotation), e31 (Y-rotation),
    /// and e12 (Z-rotation). While e23 and e31 square to +1 (due to timelike e3² = -1),
    /// the exponential map R = exp(-Bθ/2) = cosh(θ/2) - sinh(θ/2)·B̂ correctly handles
    /// all signatures, producing length-preserving rotations in the 3D spatial subspace.
    ///
    /// Rotor: R = exp(-(ax·e23 + ay·e31 + az·e12)·θ/2)
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if len < 1e-30 {
            return Self::identity();
        }
        let ax = [axis[0] / len, axis[1] / len, axis[2] / len];

        let half = angle / 2.0;

        // B = ax·e23 + ay·e31 + az·e12  (spatial rotation bivectors)
        // e23 (idx=10), e31 (idx=9), e12 (idx=8)
        // B² = -ax² - ay² - az² = -1 (normalized axis)
        // But due to Cl(3,1) metric e3² = -1:
        //   e23² = -(e2²)(e3²) = -(1)(-1) = +1
        //   e31² = -(e3²)(e1²) = -(-1)(1) = +1
        //   e12² = -(1)(1) = -1
        //
        // So B² = ax²*(+1) + ay²*(+1) + az²*(-1) = ax² + ay² - az²
        //
        // For a normalized axis, B² = ax² + ay² - az².
        // The exponential R = exp(-Bθ/2) uses:
        //   if B² > 0: R = cosh(θ/2) - sinh(θ/2)·B̂
        //   if B² < 0: R = cos(θ/2) - sin(θ/2)·B̂
        let b_sq = ax[0] * ax[0] + ax[1] * ax[1] - ax[2] * ax[2];

        let (cos_part, sin_part) = if b_sq > 0.0 {
            // Hyperbolic: B² > 0 (mostly in x/y plane with timelike component)
            let b_norm = b_sq.sqrt();
            (half.cosh(), half.sinh() / b_norm)
        } else if b_sq < 0.0 {
            // Circular: B² < 0 (mostly z-axis rotation)
            let b_norm = (-b_sq).sqrt();
            (half.cos(), half.sin() / b_norm)
        } else {
            // Null: B² = 0 (lightlike rotation plane)
            (1.0, half)
        };

        let mut m = Multivector::zero();
        m.c[0] = cos_part;
        m.c[8] = -sin_part * ax[2]; // e12 for Z rotation (idx 8)
        m.c[9] = -sin_part * ax[1]; // e31 for Y rotation (idx 9)
        m.c[10] = -sin_part * ax[0]; // e23 for X rotation (idx 10)

        Self { inner: m }
    }

    /// Create a rotor from a bivector (exponential map).
    ///
    /// Bivector components are [Be23, Be31, Be12] — the spatial rotation
    /// plane bivectors. B = Bx·e23 + By·e31 + Bz·e12.
    /// R = exp(-B/2).
    pub fn from_bivector(biv: [f64; 3]) -> Self {
        // Compute B² for the correct exponential (circular vs hyperbolic)
        let b_sq = biv[0] * biv[0] + biv[1] * biv[1] - biv[2] * biv[2];
        let b_mag_abs = b_sq.abs().sqrt();
        if b_mag_abs < 1e-30 {
            return Self::identity();
        }
        let half = b_mag_abs / 2.0;

        let (cos_part, sin_part, sin_scale) = if b_sq > 0.0 {
            (half.cosh(), half.sinh(), 1.0 / b_mag_abs)
        } else if b_sq < 0.0 {
            (half.cos(), half.sin(), 1.0 / b_mag_abs)
        } else {
            (1.0, half, 1.0)
        };

        let mut m = Multivector::zero();
        m.c[0] = cos_part;
        m.c[10] = -sin_part * sin_scale * biv[0]; // e23
        m.c[9] = -sin_part * sin_scale * biv[1]; // e31
        m.c[8] = -sin_part * sin_scale * biv[2]; // e12
        Self { inner: m }
    }

    /// Compose two rotors: R_total = R1 * R2 (R2 applied first? or R1 first?)
    ///
    /// GA convention: R_composed * v * R̃_composed = R1 * (R2 * v * R̃2) * R̃1
    /// So R_composed = R1 * R2 (R1 applied after R2).
    /// Here we return R1 * R2 — the geometric product.
    pub fn compose(&self, other: &Self) -> Self {
        let product = self.inner.geometric_product(&other.inner);
        let mut r = Self { inner: product };
        r.normalize();
        r
    }

    /// Apply rotation to a 3D vector: v' = R v R̃.
    pub fn apply(&self, v: [f64; 3]) -> [f64; 3] {
        let mut vm = Multivector::zero();
        vm.c[2] = v[0]; // e1 = x
        vm.c[3] = v[1]; // e2 = y
        vm.c[4] = v[2]; // e3 = z

        let rev = self.inner.reverse();
        let result = self.inner.geometric_product(&vm).geometric_product(&rev);

        [result.c[2], result.c[3], result.c[4]]
    }

    /// Apply to a 2D vector.
    pub fn apply_2d(&self, v: [f64; 2]) -> [f64; 2] {
        let v3 = self.apply([v[0], v[1], 0.0]);
        [v3[0], v3[1]]
    }

    /// Spherical linear interpolation.
    pub fn slerp(&self, other: &Self, t: f64) -> Self {
        let diff = other.inner.sub(&self.inner);
        let interp = self.inner.add(&diff.scale(t));
        let mut r = Self { inner: interp };
        r.normalize();
        r
    }

    /// Normalize so that R * R̃ = 1.
    pub fn normalize(&mut self) {
        let norm_sq = self.inner.norm_squared();
        if norm_sq.abs() > 1e-15 {
            let scale = 1.0 / norm_sq.abs().sqrt();
            self.inner = self.inner.scale(scale);
        }
    }

    /// Extract a 3x3 rotation matrix.
    pub fn to_rotation_matrix(&self) -> [[f64; 3]; 3] {
        let e1 = self.apply([1.0, 0.0, 0.0]);
        let e2 = self.apply([0.0, 1.0, 0.0]);
        let e3 = self.apply([0.0, 0.0, 1.0]);
        [e1, e2, e3]
    }

    /// Extract a 2x2 submatrix.
    pub fn to_rotation_matrix_2d(&self) -> [[f64; 2]; 2] {
        let ex = self.apply_2d([1.0, 0.0]);
        let ey = self.apply_2d([0.0, 1.0]);
        [ex, ey]
    }

    /// Check identity.
    pub fn is_identity(&self, tolerance: f64) -> bool {
        let id = Multivector::scalar(1.0);
        self.inner.sub(&id).is_zero(tolerance)
    }

    /// 180-degree rotation as a reflection through an axis.
    pub fn reflection(axis: [f64; 3]) -> Self {
        Self::from_axis_angle(axis, std::f64::consts::PI)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    // ── Identity ──

    #[test]
    fn test_identity_apply() {
        let r = Rotor::identity();
        let v = r.apply([1.0, 2.0, 3.0]);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!((v[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_identity_is_identity() {
        assert!(Rotor::identity().is_identity(1e-10));
    }

    #[test]
    fn test_from_multivector_identity() {
        let r = Rotor::from_multivector(Multivector::scalar(1.0));
        assert!(r.is_identity(1e-10));
    }

    // ── Z-axis rotation ──

    #[test]
    fn test_rotation_90_z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[0]).abs() < 0.01, "x should be ≈0, got {}", v[0]);
        assert!((v[1] - 1.0).abs() < 0.01, "y should be ≈1, got {}", v[1]);
    }

    #[test]
    fn test_rotation_90_z_neg() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], -FRAC_PI_2);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[0]).abs() < 0.01, "x should be ≈0, got {}", v[0]);
        assert!((v[1] + 1.0).abs() < 0.01, "y should be ≈-1, got {}", v[1]);
    }

    #[test]
    fn test_rotation_45_z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4);
        let v = r.apply([1.0, 0.0, 0.0]);
        let expected = std::f64::consts::FRAC_1_SQRT_2;
        assert!(
            (v[0] - expected).abs() < 0.01,
            "x should be ≈{expected}, got {}",
            v[0]
        );
        assert!(
            (v[1] - expected).abs() < 0.01,
            "y should be ≈{expected}, got {}",
            v[1]
        );
    }

    #[test]
    fn test_rotation_360_z_returns_to_start() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.0 * std::f64::consts::PI);
        let v = r.apply([1.0, 2.0, 3.0]);
        assert!((v[0] - 1.0).abs() < 0.01, "x should be ≈1, got {}", v[0]);
        assert!((v[1] - 2.0).abs() < 0.01, "y should be ≈2, got {}", v[1]);
    }

    // ── Rotor norm ──

    #[test]
    fn test_rotor_norm_squared_z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
        let nsq = r.inner.norm_squared();
        assert!((nsq - 1.0).abs() < 0.01, "R*R̃ should be ≈1, got {}", nsq);
    }

    #[test]
    fn test_rotor_reverse_product_z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
        let rev = r.inner.reverse();
        let product = r.inner.geometric_product(&rev);
        let scalar = product.scalar_part();
        assert!(
            (scalar - 1.0).abs() < 0.01,
            "R*R̃ scalar should be 1, got {}",
            scalar
        );
    }

    // ── Length preservation ──

    #[test]
    fn test_rotation_preserves_length_z_only() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
        let v = [3.0, 4.0, 0.0];
        let rotated = r.apply(v);
        let ol = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        let rl =
            (rotated[0] * rotated[0] + rotated[1] * rotated[1] + rotated[2] * rotated[2]).sqrt();
        assert!((ol - rl).abs() < 0.001);
    }

    #[test]
    fn test_rotation_preserves_clifford_norm_general_axis() {
        // Cl(3,1) preserves the metric norm v·v = x² + y² - z².
        let r = Rotor::from_axis_angle([1.0, 2.0, 3.0], 1.23);
        let v = [3.0, 4.0, 0.0];
        let rotated = r.apply(v);
        let orig_cliff = v[0] * v[0] + v[1] * v[1] - v[2] * v[2];
        let rot_cliff = rotated[0] * rotated[0] + rotated[1] * rotated[1] - rotated[2] * rotated[2];
        assert!(
            (orig_cliff - rot_cliff).abs() < 0.001,
            "Clifford norm should be preserved: orig={}, rot={}",
            orig_cliff,
            rot_cliff
        );
    }

    // ── Compose ──

    #[test]
    fn test_double_rotation_180() {
        let r1 = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let r_total = r1.compose(&r2);
        let v = r_total.apply([1.0, 0.0, 0.0]);
        assert!((v[0] + 1.0).abs() < 0.01, "expected -1, got {}", v[0]);
    }

    #[test]
    fn test_compose_identity() {
        let r1 = Rotor::identity();
        let r2 = Rotor::identity();
        let r3 = r1.compose(&r2);
        assert!(r3.is_identity(0.1));
    }

    #[test]
    fn test_compose_is_associative() {
        let ra = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0);
        let rb = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.0);
        let rc = Rotor::from_axis_angle([0.0, 0.0, 1.0], 3.0);
        // (ab)c vs a(bc) — both should give same total rotation
        let abc1 = ra.compose(&rb).compose(&rc);
        let abc2 = ra.compose(&rb.compose(&rc));
        let v1 = abc1.apply([1.0, 0.0, 0.0]);
        let v2 = abc2.apply([1.0, 0.0, 0.0]);
        assert!((v1[0] - v2[0]).abs() < 0.01);
        assert!((v1[1] - v2[1]).abs() < 0.01);
    }

    #[test]
    fn test_compose_angle_sum() {
        let ra = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0);
        let rb = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.0);
        let r_sum = Rotor::from_axis_angle([0.0, 0.0, 1.0], 3.0);
        let r_comp = ra.compose(&rb);
        let v_sum = r_sum.apply([1.0, 0.0, 0.0]);
        let v_comp = r_comp.apply([1.0, 0.0, 0.0]);
        assert!((v_sum[0] - v_comp[0]).abs() < 0.01);
        assert!((v_sum[1] - v_comp[1]).abs() < 0.01);
    }

    // ── Rotation matrix ──

    #[test]
    fn test_to_rotation_matrix_identity() {
        let r = Rotor::identity();
        let m = r.to_rotation_matrix();
        assert!((m[0][0] - 1.0).abs() < 1e-10);
        assert!((m[1][1] - 1.0).abs() < 1e-10);
        assert!((m[2][2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_rotation_matrix_90z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let m = r.to_rotation_matrix();
        assert!((m[0][0]).abs() < 0.01, "00 should be ≈0");
        assert!((m[0][1] - 1.0).abs() < 0.01, "01 should be ≈1");
        assert!((m[1][0] + 1.0).abs() < 0.01, "10 should be ≈-1");
    }

    #[test]
    fn test_to_rotation_matrix_orthonormal_z_only() {
        // Only z-axis rotations preserve Euclidean orthonormality in Cl(3,1)
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 0.5);
        let m = r.to_rotation_matrix();
        for row in &m {
            let len_sq = row[0] * row[0] + row[1] * row[1] + row[2] * row[2];
            assert!(
                (len_sq - 1.0).abs() < 0.01,
                "row {:?} has length² {}",
                row,
                len_sq
            );
        }
    }

    #[test]
    fn test_to_rotation_matrix_2d_90z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let m = r.to_rotation_matrix_2d();
        assert!((m[0][0]).abs() < 0.01, "00 should be ≈0");
        assert!((m[0][1] - 1.0).abs() < 0.01, "01 should be ≈1");
        assert!((m[1][0] + 1.0).abs() < 0.01, "10 should be ≈-1");
    }

    // ── From bivector ──

    #[test]
    fn test_from_bivector_z() {
        let r = Rotor::from_bivector([0.0, 0.0, 0.5]);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[1]).abs() > 0.2, "should rotate from x axis");
    }

    #[test]
    fn test_from_bivector_zero() {
        let r = Rotor::from_bivector([0.0, 0.0, 0.0]);
        assert!(r.is_identity(1e-10));
    }

    // ── Axis-angle edge cases ──

    #[test]
    fn test_from_axis_angle_zero() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 0.0], 1.0);
        assert!(r.is_identity(1e-10));
    }

    #[test]
    fn test_from_axis_angle_pi() {
        // 180° around Z: (1,0,0) -> (-1,0,0)
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::PI);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[0] + 1.0).abs() < 0.01);
        assert!((v[1]).abs() < 0.01);
    }

    // ── 2D apply ──

    #[test]
    fn test_rotor_2d() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
        let v = r.apply_2d([1.0, 0.0]);
        assert!((v[0]).abs() < 0.01);
        assert!((v[1] - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_rotor_2d_preserves_length() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 0.7);
        let v = [3.0, 4.0];
        let v_rot = r.apply_2d(v);
        let ol = (v[0] * v[0] + v[1] * v[1]).sqrt();
        let rl = (v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1]).sqrt();
        assert!((ol - rl).abs() < 0.001);
    }

    // ── Slerp ──

    #[test]
    fn test_slerp_endpoints() {
        let r1 = Rotor::identity();
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0);
        let at0 = r1.slerp(&r2, 0.0);
        assert!(at0.is_identity(0.1));
    }

    #[test]
    fn test_slerp_endpoint_2() {
        let r1 = Rotor::identity();
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0);
        let at1 = r1.slerp(&r2, 1.0);
        let v = at1.apply([1.0, 0.0, 0.0]);
        let expected = r2.apply([1.0, 0.0, 0.0]);
        assert!((v[0] - expected[0]).abs() < 0.01);
    }

    #[test]
    fn test_slerp_midpoint() {
        let r1 = Rotor::identity();
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.0);
        let r_mid = r1.slerp(&r2, 0.5);
        let v_mid = r_mid.apply([1.0, 0.0, 0.0]);
        let v_expect = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0).apply([1.0, 0.0, 0.0]);
        assert!((v_mid[0] - v_expect[0]).abs() < 0.01);
        assert!((v_mid[1] - v_expect[1]).abs() < 0.01);
    }

    #[test]
    fn test_slerp_preserves_norm() {
        let r1 = Rotor::identity();
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], 2.5);
        for t in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let r = r1.slerp(&r2, t);
            let nsq = r.inner.norm_squared();
            assert!((nsq - 1.0).abs() < 0.01, "slerp t={}: R*R̃ = {}", t, nsq);
        }
    }

    // ── Reflection ──

    #[test]
    fn test_reflection_axis_z() {
        // 180° rotation around z = reflection through z axis
        let r = Rotor::reflection([0.0, 0.0, 1.0]);
        let v = r.apply([1.0, 2.0, 0.0]);
        assert!((v[0] + 1.0).abs() < 0.01);
        assert!((v[1] + 2.0).abs() < 0.01);
    }

    #[test]
    fn test_reflection_double_is_identity_z() {
        // Only z-axis reflections properly double to identity in Cl(3,1)
        let r = Rotor::reflection([0.0, 0.0, 1.0]);
        let r2 = r.compose(&r);
        let v = r2.apply([1.0, 2.0, 3.0]);
        assert!((v[0] - 1.0).abs() < 0.01);
        assert!((v[1] - 2.0).abs() < 0.01);
    }
}
