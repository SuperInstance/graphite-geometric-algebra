//! Conformal geometric algebra: embeds Euclidean 3D into 5D conformal space.

use super::multivector::Multivector;
use super::rotor::Rotor;

/// Conformal geometric algebra operations.
pub struct Conformal;

impl Conformal {
    /// Embed a Euclidean 3D point into conformal space.
    /// P = e₊ + p + ½|p|²e₋
    pub fn embed_point(p: [f64; 3]) -> Multivector {
        let norm_sq = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        let mut m = Multivector::zero();
        m.c[1] = 1.0;   // e₊ component (e0)
        m.c[2] = p[0];  // e1 (x)
        m.c[3] = p[1];  // e2 (y)
        m.c[4] = p[2];  // e3 (z)
        m.c[0] = 0.5 * norm_sq;
        m
    }

    /// Extract a Euclidean 3D point.
    pub fn extract_point(m: &Multivector) -> [f64; 3] {
        [m.c[2], m.c[3], m.c[4]]
    }

    /// Reflect a point through a plane (normal + distance from origin).
    pub fn reflect(point: [f64; 3], normal: [f64; 3], distance: f64) -> [f64; 3] {
        let dot = point[0]*normal[0] + point[1]*normal[1] + point[2]*normal[2];
        let factor = 2.0 * (dot - distance);
        [
            point[0] - factor * normal[0],
            point[1] - factor * normal[1],
            point[2] - factor * normal[2],
        ]
    }

    /// Rotate a point using a rotor.
    pub fn rotate(point: [f64; 3], rotor: &Rotor) -> [f64; 3] {
        rotor.apply(point)
    }

    /// Distance between two embedded conformal points.
    pub fn conformal_distance(p1: &Multivector, p2: &Multivector) -> f64 {
        let ip = p1.inner(p2);
        let neg_half_dsq = ip.scalar_part();
        if neg_half_dsq < 0.0 { (-2.0 * neg_half_dsq).sqrt() } else { 0.0 }
    }

    /// Create a plane multivector from normal and distance.
    pub fn plane(normal: [f64; 3], distance: f64) -> Multivector {
        let mut m = Multivector::zero();
        m.c[2] = normal[0];
        m.c[3] = normal[1];
        m.c[4] = normal[2];
        m.c[0] = -distance;
        m
    }

    /// Project a point onto a plane.
    pub fn project_onto_plane(point: [f64; 3], normal: [f64; 3], distance: f64) -> [f64; 3] {
        let dot = point[0]*normal[0] + point[1]*normal[1] + point[2]*normal[2];
        let factor = dot - distance;
        [
            point[0] - factor * normal[0],
            point[1] - factor * normal[1],
            point[2] - factor * normal[2],
        ]
    }

    /// Midpoint.
    pub fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [(a[0]+b[0])/2.0, (a[1]+b[1])/2.0, (a[2]+b[2])/2.0]
    }

    /// Barycentric combination.
    pub fn barycenter(points: &[[f64; 3]], weights: &[f64]) -> [f64; 3] {
        let total: f64 = weights.iter().sum();
        if total.abs() < 1e-15 { return [0.0, 0.0, 0.0]; }
        let mut r = [0.0; 3];
        for (p, w) in points.iter().zip(weights) {
            r[0] += p[0]*w; r[1] += p[1]*w; r[2] += p[2]*w;
        }
        [r[0]/total, r[1]/total, r[2]/total]
    }

    /// Translate a point using conformal translation rotor.
    pub fn translate(point: [f64; 3], translation: [f64; 3]) -> [f64; 3] {
        [
            point[0] + translation[0],
            point[1] + translation[1],
            point[2] + translation[2],
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_embed_extract() {
        let p = [1.0, 2.0, 3.0];
        assert_eq!(Conformal::extract_point(&Conformal::embed_point(p)), p);
    }

    #[test]
    fn test_reflect_through_origin() {
        let p = [1.0, 2.0, 3.0];
        let r = Conformal::reflect(p, [1.0, 0.0, 0.0], 0.0);
        assert!((r[0] - (-1.0)).abs() < 1e-10);
        assert!((r[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_double_is_identity() {
        let p = [1.0, 2.0, 3.0];
        let n = [0.0, 0.0, 1.0];
        let r1 = Conformal::reflect(p, n, 0.0);
        let r2 = Conformal::reflect(r1, n, 0.0);
        for i in 0..3 { assert!((r2[i] - p[i]).abs() < 1e-10); }
    }

    #[test]
    fn test_project_onto_plane() {
        let p = [1.0, 1.0, 1.0];
        let pr = Conformal::project_onto_plane(p, [0.0, 0.0, 1.0], 0.0);
        assert!((pr[2]).abs() < 1e-10);
        assert_eq!(pr[0], 1.0);
    }

    #[test]
    fn test_barycenter() {
        let pts = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let bc = Conformal::barycenter(&pts, &[1.0, 1.0, 1.0]);
        assert!((bc[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_midpoint() {
        let m = Conformal::midpoint([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert_eq!(m, [1.0, 2.0, 3.0]);
    }
}
