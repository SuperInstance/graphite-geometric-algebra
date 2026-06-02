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
        m.c[1] = 1.0; // e₊ component (e0)
        m.c[2] = p[0]; // e1 (x)
        m.c[3] = p[1]; // e2 (y)
        m.c[4] = p[2]; // e3 (z)
        m.c[0] = 0.5 * norm_sq;
        m
    }

    /// Extract a Euclidean 3D point.
    pub fn extract_point(m: &Multivector) -> [f64; 3] {
        [m.c[2], m.c[3], m.c[4]]
    }

    /// Reflect a point through a plane (normal + distance from origin).
    pub fn reflect(point: [f64; 3], normal: [f64; 3], distance: f64) -> [f64; 3] {
        let dot = point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2];
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

    /// Euclidean distance between two points (via extracted coordinates).
    ///
    /// The simplified CGA embedding in this crate stores the norm term ½||p||²
    /// in the scalar component, which doesn't form true null vectors in Cl(3,1).
    /// For proper conformal distance the standard embedding requires a 5D
    /// representation with null vectors n₊ and n∞. Here we compute the
    /// Euclidean distance directly from the positional components.
    pub fn conformal_distance(p1: &Multivector, p2: &Multivector) -> f64 {
        let a = Self::extract_point(p1);
        let b = Self::extract_point(p2);
        let dx = a[0] - b[0];
        let dy = a[1] - b[1];
        let dz = a[2] - b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
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
        let dot = point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2];
        let factor = dot - distance;
        [
            point[0] - factor * normal[0],
            point[1] - factor * normal[1],
            point[2] - factor * normal[2],
        ]
    }

    /// Midpoint.
    pub fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            (a[0] + b[0]) / 2.0,
            (a[1] + b[1]) / 2.0,
            (a[2] + b[2]) / 2.0,
        ]
    }

    /// Barycentric combination.
    pub fn barycenter(points: &[[f64; 3]], weights: &[f64]) -> [f64; 3] {
        let total: f64 = weights.iter().sum();
        if total.abs() < 1e-15 {
            return [0.0, 0.0, 0.0];
        }
        let mut r = [0.0; 3];
        for (p, w) in points.iter().zip(weights) {
            r[0] += p[0] * w;
            r[1] += p[1] * w;
            r[2] += p[2] * w;
        }
        [r[0] / total, r[1] / total, r[2] / total]
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

    // ── Embed / Extract ──

    #[test]
    fn test_embed_extract_origin() {
        let p = [0.0, 0.0, 0.0];
        assert_eq!(Conformal::extract_point(&Conformal::embed_point(p)), p);
    }

    #[test]
    fn test_embed_extract() {
        let p = [1.0, 2.0, 3.0];
        assert_eq!(Conformal::extract_point(&Conformal::embed_point(p)), p);
    }

    #[test]
    fn test_embed_scalar_part() {
        // For p = [3,4,0], norm² = 25, so scalar should be 12.5
        let m = Conformal::embed_point([3.0, 4.0, 0.0]);
        assert!((m.scalar_part() - 12.5).abs() < 1e-10);
    }

    #[test]
    fn test_embed_e0_component() {
        // The e0 component is always 1.0
        let m = Conformal::embed_point([5.0, 6.0, 7.0]);
        assert!((m.c[1] - 1.0).abs() < 1e-10);
    }

    // ── Conformal distance ──

    #[test]
    fn test_conformal_distance_self() {
        let p = [1.0, 2.0, 3.0];
        let mp = Conformal::embed_point(p);
        let d = Conformal::conformal_distance(&mp, &mp);
        assert!(
            d.abs() < 1e-10,
            "distance from point to self should be 0, got {}",
            d
        );
    }

    #[test]
    fn test_conformal_distance_known() {
        let a = Conformal::embed_point([0.0, 0.0, 0.0]);
        let b = Conformal::embed_point([1.0, 0.0, 0.0]);
        let d = Conformal::conformal_distance(&a, &b);
        assert!((d - 1.0).abs() < 0.001, "distance should be 1, got {}", d);
    }

    #[test]
    fn test_conformal_distance_3d() {
        // Distance from (0,0,0) to (3,4,0) should be 5
        let a = Conformal::embed_point([0.0, 0.0, 0.0]);
        let b = Conformal::embed_point([3.0, 4.0, 0.0]);
        let d = Conformal::conformal_distance(&a, &b);
        assert!((d - 5.0).abs() < 0.001, "distance should be 5, got {}", d);
    }

    // ── Reflection ──

    #[test]
    fn test_reflect_through_origin() {
        let p = [1.0, 2.0, 3.0];
        let r = Conformal::reflect(p, [1.0, 0.0, 0.0], 0.0);
        assert!((r[0] - (-1.0)).abs() < 1e-10);
        assert!((r[1] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_across_plane_z() {
        // Reflect (1,2,3) across XY plane at z=2
        let p = [1.0, 2.0, 3.0];
        let r = Conformal::reflect(p, [0.0, 0.0, 1.0], 2.0);
        // distance from plane = 3-2 = 1 -> mirror to z = 1
        assert!((r[2] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_reflect_double_is_identity() {
        let p = [1.0, 2.0, 3.0];
        let n = [0.0, 0.0, 1.0];
        let r1 = Conformal::reflect(p, n, 0.0);
        let r2 = Conformal::reflect(r1, n, 0.0);
        for i in 0..3 {
            assert!((r2[i] - p[i]).abs() < 1e-10);
        }
    }

    // ── Project onto plane ──

    #[test]
    fn test_project_onto_xy_plane() {
        let p = [1.0, 1.0, 1.0];
        let pr = Conformal::project_onto_plane(p, [0.0, 0.0, 1.0], 0.0);
        assert!((pr[2]).abs() < 1e-10);
        assert_eq!(pr[0], 1.0);
    }

    #[test]
    fn test_project_onto_offset_plane() {
        // Project (0,0,5) onto plane z=2, normal [0,0,1]
        let pr = Conformal::project_onto_plane([0.0, 0.0, 5.0], [0.0, 0.0, 1.0], 2.0);
        assert!((pr[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_project_onto_plane_point_on_plane() {
        // Point already on the plane should stay put
        let pr = Conformal::project_onto_plane([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 3.0);
        assert!((pr[2] - 3.0).abs() < 1e-10);
    }

    // ── Barycenter ──

    #[test]
    fn test_barycenter_triangle() {
        let pts = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]];
        let bc = Conformal::barycenter(&pts, &[1.0, 1.0, 1.0]);
        assert!((bc[0] - 1.0).abs() < 1e-10);
        assert!((bc[1] - 2.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_barycenter_weights() {
        // weighted centroid with custom weights
        let pts = [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]];
        let bc = Conformal::barycenter(&pts, &[3.0, 1.0]);
        assert!((bc[0] - 1.0).abs() < 1e-10); // (0*3 + 4*1)/4 = 1
    }

    #[test]
    fn test_barycenter_single_point() {
        let bc = Conformal::barycenter(&[[5.0, 6.0, 7.0]], &[1.0]);
        assert_eq!(bc, [5.0, 6.0, 7.0]);
    }

    // ── Midpoint ──

    #[test]
    fn test_midpoint() {
        let m = Conformal::midpoint([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        assert_eq!(m, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_midpoint_identical_points() {
        let m = Conformal::midpoint([3.0, 3.0, 3.0], [3.0, 3.0, 3.0]);
        assert_eq!(m, [3.0, 3.0, 3.0]);
    }

    // ── Rotate ──

    #[test]
    fn test_conformal_rotate() {
        let r = Rotor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let p = Conformal::rotate([1.0, 0.0, 0.0], &r);
        assert!((p[0]).abs() < 0.01);
        assert!((p[1] - 1.0).abs() < 0.01);
    }

    // ── Translate ──

    #[test]
    fn test_conformal_translate() {
        let p = Conformal::translate([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(p, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_translate_zero_is_identity() {
        let p = Conformal::translate([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        assert_eq!(p, [1.0, 2.0, 3.0]);
    }

    // ── Plane ──

    #[test]
    fn test_plane_creation() {
        let pl = Conformal::plane([1.0, 0.0, 0.0], 5.0);
        assert!((pl.c[2] - 1.0).abs() < 1e-10);
        assert!((pl.c[0] - (-5.0)).abs() < 1e-10);
    }
}
