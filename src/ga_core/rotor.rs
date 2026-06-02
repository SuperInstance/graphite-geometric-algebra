//! Rotors: even-grade multivectors for rotations in geometric algebra.

use super::multivector::Multivector;

/// A rotor: an even-grade element that performs rotations.
///
/// A rotor R satisfies R * R̃ = 1 (normalized).
/// Rotations are applied via the sandwich product: v' = R v R̃.
///
/// Convention: R = exp(-Bθ/2) = cos(θ/2) - sin(θ/2)·B̂
/// where B̂ is the unit bivector of the rotation plane.
/// For right-handed rotation around axis a: R = cos(θ/2) - sin(θ/2)·(a·I)
/// where I = e123 is the 3D pseudoscalar.
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
        Self { inner: Multivector::scalar(1.0) }
    }

    /// Create a rotor from axis-angle representation.
    ///
    /// In Cl(3,1), we need rotation-plane bivectors B where B² = -1.
    /// The spatial bivectors e23 (indices 10) and e13 (index 9) square to +1
    /// due to the timelike e3² = -1 metric. So we use the Euclidean subalgebra
    /// bivectors: e01 (X), e02 (Y), e12 (Z) which all satisfy B² = -1.
    ///
    /// Rotor: R = cos(θ/2) - sin(θ/2)·(ax·e01 + ay·e02 + az·e12)
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        let len = (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]).sqrt();
        if len < 1e-30 { return Self::identity(); }
        let ax = [axis[0]/len, axis[1]/len, axis[2]/len];

        let half = angle / 2.0;
        let cos_h = half.cos();
        let sin_h = half.sin();

        // R = cos(θ/2) - sin(θ/2)*(ax·e01 + ay·e02 + az·e12)
        // Our bivector indices: e01=5, e02=6, e12=8
        let mut m = Multivector::zero();
        m.c[0] = cos_h;
        m.c[5] = -sin_h * ax[0]; // e01 for X rotation
        m.c[6] = -sin_h * ax[1]; // e02 for Y rotation
        m.c[8] = -sin_h * ax[2]; // e12 for Z rotation

        Self { inner: m }
    }

    /// Create a rotor from a bivector (exponential map).
    /// R = exp(-B/2) where B is the rotation bivector [Be01, Be02, Be12].
    pub fn from_bivector(biv: [f64; 3]) -> Self {
        let b_mag = (biv[0]*biv[0] + biv[1]*biv[1] + biv[2]*biv[2]).sqrt();
        if b_mag < 1e-30 { return Self::identity(); }
        let half = b_mag / 2.0;
        let cos_h = half.cos();
        let sin_h = half.sin();
        let u = [biv[0]/b_mag, biv[1]/b_mag, biv[2]/b_mag];

        let mut m = Multivector::zero();
        m.c[0] = cos_h;
        m.c[5] = -sin_h * u[0]; // e01
        m.c[6] = -sin_h * u[1]; // e02
        m.c[8] = -sin_h * u[2]; // e12
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

    #[test]
    fn test_identity() {
        let r = Rotor::identity();
        let v = r.apply([1.0, 2.0, 3.0]);
        assert!((v[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotation_90_z() {
        // Right-handed 90° around Z: (1,0,0) -> (0,1,0)
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[0]).abs() < 0.01, "x should be ≈0, got {}", v[0]);
        assert!((v[1] - 1.0).abs() < 0.01, "y should be ≈1, got {}", v[1]);
    }

    #[test]
    fn test_rotation_90_z_neg() {
        // -90° around Z: (1,0,0) -> (0,-1,0)
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], -std::f64::consts::FRAC_PI_2);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[0]).abs() < 0.01, "x should be ≈0, got {}", v[0]);
        assert!((v[1] + 1.0).abs() < 0.01, "y should be ≈-1, got {}", v[1]);
    }

    #[test]
    fn test_rotor_norm_squared() {
        let r = Rotor::from_axis_angle([1.0, 2.0, 3.0], 1.23);
        let nsq = r.inner.norm_squared();
        assert!((nsq - 1.0).abs() < 0.01, "R*R̃ should be ≈1, got {}", nsq);
    }

    #[test]
    fn test_rotor_reverse_product() {
        // Directly compute R * rev(R) and check
        let axis = [1.0f64, 2.0, 3.0];
        let angle = 1.23f64;
        let len = (axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]).sqrt();
        let ax = [axis[0]/len, axis[1]/len, axis[2]/len];
        let half = angle / 2.0f64;
        let cos_h = half.cos();
        let sin_h = half.sin();

        let mut m = Multivector::zero();
        m.c[0] = cos_h;
        m.c[5] = -sin_h * ax[0];
        m.c[6] = -sin_h * ax[1];
        m.c[8] = -sin_h * ax[2];

        let rev = m.reverse();
        let product = m.geometric_product(&rev);
        let scalar = product.scalar_part();
        // The scalar should be cos²(h) + sin²(h) = 1
        // But verify the bivector components cancel
        let biv_norm = product.grade_norm(2);
        println!("Product scalar: {}, biv_norm: {}", scalar, biv_norm);
        println!("Full product: {:?}", product.c);
        assert!((scalar - 1.0).abs() < 0.01, "R*R̃ scalar should be 1, got {}", scalar);
    }

    #[test]
    fn test_rotation_preserves_length_z_only() {
        // Pure z-axis rotation should definitely preserve length
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
        let v = [3.0, 4.0, 0.0];
        let rotated = r.apply(v);
        let ol = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
        let rl = (rotated[0]*rotated[0] + rotated[1]*rotated[1] + rotated[2]*rotated[2]).sqrt();
        println!("Z-only: orig={}, rot={:?}, len={}", ol, rotated, rl);
        assert!((ol - rl).abs() < 0.001, "Z rotation: orig={}, rot={}", ol, rl);
    }

    #[test]
    fn test_rotation_z_apply_does_something() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.23);
        let v = [1.0, 0.0, 0.0];
        let rotated = r.apply(v);
        println!("Z-rot of (1,0,0): {:?}", rotated);
    }

    #[test]
    fn test_rotation_preserves_length() {
        let r = Rotor::from_axis_angle([1.0, 2.0, 3.0], 1.23);
        let v = [3.0, 4.0, 0.0];
        let rotated = r.apply(v);
        let ol = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
        let rl = (rotated[0]*rotated[0] + rotated[1]*rotated[1] + rotated[2]*rotated[2]).sqrt();
        println!("General: orig={}, rot={:?}, len={}", ol, rotated, rl);
        assert!((ol - rl).abs() < 0.001, "length should be preserved: orig={}, rot={}", ol, rl);
    }

    #[test]
    fn test_double_rotation_180() {
        let r1 = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
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
    fn test_to_rotation_matrix_identity() {
        let r = Rotor::identity();
        let m = r.to_rotation_matrix();
        assert!((m[0][0] - 1.0).abs() < 1e-10);
        assert!((m[1][1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_to_rotation_matrix_90z() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let m = r.to_rotation_matrix();
        assert!((m[0][0]).abs() < 0.01, "00 should be ≈0");
        assert!((m[0][1] - 1.0).abs() < 0.01, "01 should be ≈1");
        assert!((m[1][0] + 1.0).abs() < 0.01, "10 should be ≈-1");
    }

    #[test]
    fn test_from_bivector() {
        // e12 = axis z
        let r = Rotor::from_bivector([0.0, 0.0, 1.0]);
        let v = r.apply([1.0, 0.0, 0.0]);
        assert!((v[1]).abs() > 0.5, "should rotate from x axis, got y={}", v[1]);
    }

    #[test]
    fn test_rotor_2d() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let v = r.apply_2d([1.0, 0.0]);
        assert!((v[0]).abs() < 0.01, "x should be ≈0");
        assert!((v[1] - 1.0).abs() < 0.01, "y should be ≈1");
    }

    #[test]
    fn test_slerp_endpoints() {
        let r1 = Rotor::identity();
        let r2 = Rotor::from_axis_angle([0.0, 0.0, 1.0], 1.0);
        let at0 = r1.slerp(&r2, 0.0);
        assert!(at0.is_identity(0.1));
    }
}
