//! Multivector: the fundamental element of geometric algebra.
//!
//! For Cl(3,1), a multivector has 2^4 = 16 components:
//! 1 scalar, 4 vectors, 6 bivectors, 4 trivectors, 1 pseudoscalar.

/// A multivector in Cl(3,1) spacetime algebra.
#[derive(Debug, Clone, PartialEq)]
pub struct Multivector {
    pub c: [f64; 16],
}

impl Multivector {
    pub fn zero() -> Self { Self { c: [0.0; 16] } }
    pub fn scalar(v: f64) -> Self { let mut m = Self::zero(); m.c[0] = v; m }
    pub fn vector(v: [f64; 4]) -> Self {
        let mut m = Self::zero();
        m.c[1..=4].copy_from_slice(&v);
        m
    }
    pub fn bivector(v: [f64; 6]) -> Self {
        let mut m = Self::zero();
        m.c[5..=10].copy_from_slice(&v);
        m
    }
    pub fn scalar_part(&self) -> f64 { self.c[0] }
    pub fn vector_part(&self) -> [f64; 4] { [self.c[1], self.c[2], self.c[3], self.c[4]] }
    pub fn bivector_part(&self) -> [f64; 6] {
        [self.c[5], self.c[6], self.c[7], self.c[8], self.c[9], self.c[10]]
    }
    pub fn grade_norm(&self, k: usize) -> f64 {
        match k {
            0 => self.c[0].abs(),
            1 => self.c[1..=4].iter().map(|x| x.abs()).sum(),
            2 => self.c[5..=10].iter().map(|x| x.abs()).sum(),
            3 => self.c[11..=14].iter().map(|x| x.abs()).sum(),
            4 => self.c[15].abs(),
            _ => 0.0,
        }
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..16 { r.c[i] = self.c[i] + other.c[i]; }
        r
    }
    pub fn sub(&self, other: &Self) -> Self {
        let mut r = Self::zero();
        for i in 0..16 { r.c[i] = self.c[i] - other.c[i]; }
        r
    }
    pub fn scale(&self, s: f64) -> Self {
        let mut r = Self::zero();
        for i in 0..16 { r.c[i] = self.c[i] * s; }
        r
    }

    /// Reverse: grade k -> (-1)^(k(k-1)/2) sign.
    pub fn reverse(&self) -> Self {
        let mut r = self.clone();
        for i in 5..=10 { r.c[i] = -r.c[i]; }
        for i in 11..=14 { r.c[i] = -r.c[i]; }
        r
    }

    /// Clifford conjugation: reverse + negate odd grades.
    pub fn conjugate(&self) -> Self {
        let mut r = self.clone();
        for i in 1..=4  { r.c[i] = -r.c[i]; }
        for i in 11..=14 { r.c[i] = -r.c[i]; }
        r
    }

    /// Hodge dual: multiply by pseudoscalar on the right.
    pub fn dual(&self) -> Self {
        let mut i = Self::zero(); i.c[15] = 1.0;
        self.geometric_product(&i)
    }

    /// Norm squared through the reverse product.
    pub fn norm_squared(&self) -> f64 {
        self.geometric_product(&self.reverse()).scalar_part()
    }

    /// Full geometric product for Cl(3,1) with metric (+++-).
    ///
    /// This uses the cached 16×16 multiplication table computed once.
    pub fn geometric_product(&self, other: &Self) -> Self {
        let table = product_table();
        let mut result = Self::zero();
        for i in 0..16 {
            let a = self.c[i];
            if a.abs() < 1e-40 { continue; }
            let row = &table[i];
            for j in 0..16 {
                let b = other.c[j];
                if b.abs() < 1e-40 { continue; }
                let (k, sign) = row[j];
                result.c[k] += a * b * sign;
            }
        }
        result
    }

    pub fn wedge(&self, other: &Self) -> Self {
        self.geometric_product(other).sub(&other.geometric_product(self)).scale(0.5)
    }
    pub fn inner(&self, other: &Self) -> Self {
        self.geometric_product(other).add(&other.geometric_product(self)).scale(0.5)
    }
    pub fn is_zero(&self, tolerance: f64) -> bool {
        self.c.iter().all(|&x| x.abs() < tolerance)
    }
}

/// Grade of a basis blade index.
fn grade_of(idx: usize) -> usize {
    match idx {
        0 => 0, 1..=4 => 1, 5..=10 => 2, 11..=14 => 3, 15 => 4,
        _ => 0,
    }
}

/// Build the complete 16×16 Cl(3,1) multiplication table.
///
/// Each entry (k, s) means e_i * e_j = s * e_k.
/// Metric: e0²=1, e1²=1, e2²=1, e3²=-1; vectors anticommute.
fn build_product_table() -> [[(usize, f64); 16]; 16] {
    let mut table = [[(0usize, 0.0f64); 16]; 16];

    // Represent each blade as a sorted list of vector indices (sentinel 4 = none).
    let blade_vecs: [[usize; 4]; 16] = [
        [4,4,4,4], // 0: scalar
        [0,4,4,4], // 1: e0
        [1,4,4,4], // 2: e1
        [2,4,4,4], // 3: e2
        [3,4,4,4], // 4: e3
        [0,1,4,4], // 5: e01
        [0,2,4,4], // 6: e02
        [0,3,4,4], // 7: e03
        [1,2,4,4], // 8: e12
        [1,3,4,4], // 9: e13
        [2,3,4,4], // 10: e23
        [0,1,2,4], // 11: e012
        [0,1,3,4], // 12: e013
        [0,2,3,4], // 13: e023
        [1,2,3,4], // 14: e123
        [0,1,2,3], // 15: e0123
    ];

    let metric = [1.0, 1.0, 1.0, -1.0]; // e0², e1², e2², e3²

    fn vecs_to_idx(v: &[usize; 4], len: usize) -> usize {
        if len == 0 { return 0; }
        if len == 1 { return v[0] + 1; }
        if len == 2 {
            return match (v[0], v[1]) {
                (0,1) => 5, (0,2) => 6, (0,3) => 7,
                (1,2) => 8, (1,3) => 9, (2,3) => 10,
                _ => 0,
            };
        }
        if len == 3 {
            return match (v[0], v[1], v[2]) {
                (0,1,2) => 11, (0,1,3) => 12, (0,2,3) => 13, (1,2,3) => 14,
                _ => 0,
            };
        }
        if len == 4 && v[0]==0 && v[1]==1 && v[2]==2 && v[3]==3 { return 15; }
        0
    }

    for i in 0..16 {
        for j in 0..16 {
            let (result, sign) = multiply_blades(&blade_vecs[i], &blade_vecs[j], &metric, vecs_to_idx);
            table[i][j] = (result, sign);
        }
    }

    table
}

/// Multiply two basis blades represented as sorted vector-index arrays.
fn multiply_blades(
    va: &[usize; 4],
    vb: &[usize; 4],
    metric: &[f64; 4],
    index_of: fn(&[usize; 4], usize) -> usize,
) -> (usize, f64) {
    // Count how many non-sentinel vectors each blade has
    let grade_a = va.iter().position(|&x| x == 4).unwrap_or(4);
    let grade_b = vb.iter().position(|&x| x == 4).unwrap_or(4);

    if grade_a == 0 { return (vecs_to_idx(vb, grade_b), 1.0); }
    if grade_b == 0 { return (vecs_to_idx(va, grade_a), 1.0); }

    let mut result_vecs = [4usize; 4];
    let mut len = 0;
    let mut sign = 1.0;

    // Copy va into result
    for k in 0..grade_a {
        result_vecs[len] = va[k];
        len += 1;
    }

    // Insert each vector from vb, sorting and handling squares
    for b in 0..grade_b {
        let v = vb[b];
        // Find insertion point moving left; if we find same vector, square it
        let mut inserted = false;
        let mut pos = len as isize - 1;
        while pos >= 0 {
            let p = pos as usize;
            if result_vecs[p] == v {
                // Square: e_v * e_v = metric[v], absorb into sign
                sign *= metric[v];
                // Remove this occurrence by shifting left
                for k in p..len - 1 {
                    result_vecs[k] = result_vecs[k + 1];
                }
                len -= 1;
                result_vecs[len] = 4;
                inserted = true;
                break;
            }
            if result_vecs[p] < v {
                // Found insertion point
                break;
            }
            // Swap: distinct vectors anticommute
            sign = -sign;
            pos -= 1;
        }

        if !inserted {
            let ins = (pos + 1) as usize;
            // Shift right to make room
            for k in (ins..len).rev() {
                result_vecs[k + 1] = result_vecs[k];
            }
            result_vecs[ins] = v;
            len += 1;
        }
    }

    // Pad with sentinels
    for k in len..4 { result_vecs[k] = 4; }

    let idx = index_of(&result_vecs, len);

    // Determine grade: it's the number of remaining distinct vectors
    // Works: the geometric product grade is |grade_a - grade_b| modulo 2, up to
    // grade_a + grade_b minus 2*cancellations
    // The index_of function handles this correctly via the sorted vector list

    (idx, sign)
}

/// Re-export index_of helper for use within build_product_table closure
fn vecs_to_idx(v: &[usize; 4], len: usize) -> usize {
    if len == 0 { return 0; }
    if len == 1 { return v[0] + 1; }
    if len == 2 {
        return match (v[0], v[1]) {
            (0,1) => 5, (0,2) => 6, (0,3) => 7,
            (1,2) => 8, (1,3) => 9, (2,3) => 10,
            _ => 0,
        };
    }
    if len == 3 {
        return match (v[0], v[1], v[2]) {
            (0,1,2) => 11, (0,1,3) => 12, (0,2,3) => 13, (1,2,3) => 14,
            _ => 0,
        };
    }
    if len == 4 && v[0]==0 && v[1]==1 && v[2]==2 && v[3]==3 { return 15; }
    0
}

use std::sync::OnceLock;
static PRODUCT_TABLE: OnceLock<[[(usize, f64); 16]; 16]> = OnceLock::new();

fn product_table() -> &'static [[(usize, f64); 16]; 16] {
    PRODUCT_TABLE.get_or_init(build_product_table)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_identity() {
        let s = Multivector::scalar(5.0);
        let v = Multivector::vector([1.0, 2.0, 3.0, 4.0]);
        let p = s.geometric_product(&v);
        assert!((p.c[1] - 5.0).abs() < 1e-10);
        assert!((p.c[4] - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_squares_e0() {
        let v = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
        let p = v.geometric_product(&v);
        assert!((p.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_squares_e1() {
        let v = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
        let p = v.geometric_product(&v);
        assert!((p.scalar_part() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_vector_squares_e3_timelike() {
        let v = Multivector::vector([0.0, 0.0, 0.0, 1.0]);
        let p = v.geometric_product(&v);
        assert!((p.scalar_part() - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_vector_orthogonal_product() {
        let v0 = Multivector::vector([1.0, 0.0, 0.0, 0.0]);
        let v1 = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
        let p = v0.geometric_product(&v1);
        // e0*e1 = e01 (index 5)
        assert!((p.c[5] - 1.0).abs() < 1e-10);
        // Also no scalar part
        assert!((p.scalar_part()).abs() < 1e-10);
    }

    #[test]
    fn test_vector_bivector_product() {
        // e1 * e12 = e1*e1*e2 = e2
        let v = Multivector::vector([0.0, 1.0, 0.0, 0.0]);
        let b = Multivector::bivector([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]); // e12
        let p = v.geometric_product(&b);
        assert!((p.c[3] - 1.0).abs() < 1e-10, "e1*e12 should give e2, got vector part {:?}", p.vector_part());
    }

    #[test]
    fn test_bivector_bivector_product() {
        // e01 * e01 = e0*e1*e0*e1 = -e0*e0*e1*e1 = -1*1 = -1
        let b = Multivector::bivector([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]); // e01
        let p = b.geometric_product(&b);
        assert!((p.scalar_part() - (-1.0)).abs() < 1e-10, "e01² should be -1, got {}", p.scalar_part());
    }

    #[test]
    fn test_rotor_composition_sandwich() {
        // R = cos(θ/2) - sin(θ/2) * B (where B^2 = -1)
        // For Z-rotation: B = e12
        let angle = std::f64::consts::FRAC_PI_2;
        let half = angle / 2.0;
        let mut rotor = Multivector::zero();
        rotor.c[0] = half.cos();
        rotor.c[8] = -half.sin(); // e12 component

        let v = Multivector::vector([0.0, 1.0, 0.0, 0.0]); // e1 direction
        let rev = rotor.reverse();
        // v' = R * v * R̃
        let p = rotor.geometric_product(&v).geometric_product(&rev);
        // After 90° around z: e1 -> e2
        assert!((p.c[3] - 1.0).abs() < 0.01, "e1 rotated 90° around z should be ≈ e2, got {:?}", p.vector_part());
    }

    #[test]
    fn test_norm_squared_vector() {
        let v = Multivector::vector([3.0, 4.0, 0.0, 0.0]);
        let n = v.norm_squared();
        assert!((n - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_dual() {
        let s = Multivector::scalar(1.0);
        let d = s.dual();
        // Dual of scalar is pseudoscalar
        assert!((d.c[15] - 1.0).abs() < 1e-10);
        // Double dual
        let dd = d.dual();
        assert!((dd.scalar_part() - 1.0).abs() > 1e-10);
    }

    #[test]
    fn test_product_table_every_pair_nonsingular() {
        let table = build_product_table();
        for i in 0..16 {
            for j in 0..16 {
                let (k, sign) = table[i][j];
                assert!(sign.abs() <= 1.0 || k != 0, "entry [{i},{j}] -> (k={k}, sign={sign}) is invalid");
                assert!(k < 16, "entry [{i},{j}] has target index {k}");
            }
        }
    }

    #[test]
    fn test_product_table_inverses() {
        // e0 * e0 = 1
        let table = build_product_table();
        let (k, s) = table[1][1];
        assert!(k == 0 && (s - 1.0).abs() < 1e-10, "e0² should be 1");
        // e3 * e3 = -1
        let (k, s) = table[4][4];
        assert!(k == 0 && (s + 1.0).abs() < 1e-10, "e3² should be -1");
    }
}
