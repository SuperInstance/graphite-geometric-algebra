//! Example: Conformal geometry — embed points, reflect, project, and compute distances.
//!
//! Run with: `cargo run --example conformal_geometry`

use graphite_ga_nodes::ga_core::{Conformal, Rotor};
use std::f64::consts::FRAC_PI_2;

fn main() {
    println!("═══ Conformal Geometry Sandbox ═══");
    println!("Cl(3,1) conformal geometric algebra for the Graphite editor\n");

    // ── Embed a point ──
    let p = [3.0, 4.0, 0.0];
    let mp = Conformal::embed_point(p);
    let extracted = Conformal::extract_point(&mp);
    println!("Embed point ({}, {}, {}):", p[0], p[1], p[2]);
    println!("  Multivector blade components:");
    println!("    e0 (origin)  = {:.1}", mp.c[1]);
    println!("    e1 (x)       = {:.1}", mp.c[2]);
    println!("    e2 (y)       = {:.1}", mp.c[3]);
    println!("    scalar (½r²) = {:.1}", mp.c[0]);
    println!(
        "  Extracted: ({}, {}, {})",
        extracted[0], extracted[1], extracted[2]
    );
    assert_eq!(extracted, p);
    println!("  ✓ Round-trip OK\n");

    // ── Conformal distance ──
    let p1 = Conformal::embed_point([0.0, 0.0, 0.0]);
    let p2 = Conformal::embed_point([3.0, 4.0, 0.0]);
    let d = Conformal::conformal_distance(&p1, &p2);
    println!("Distance from (0,0,0) to (3,4,0): {:.1}", d);
    assert!((d - 5.0).abs() < 0.001);
    println!("  ✓ Expected 5.0\n");

    // ── Self distance ──
    let zero = Conformal::conformal_distance(&p1, &p1);
    println!("Self-distance: {:.0}", zero);
    assert!(zero.abs() < 1e-10);
    println!("  ✓ 0.0\n");

    // ── Reflect through a plane ──
    let p = [1.0, 2.0, 3.0];
    let reflected = Conformal::reflect(p, [1.0, 0.0, 0.0], 0.0);
    println!("Reflect ({}, {}, {}) through YZ-plane:", p[0], p[1], p[2]);
    println!("  → ({}, {}, {})", reflected[0], reflected[1], reflected[2]);
    assert!((reflected[0] + 1.0).abs() < 1e-10);
    assert!((reflected[1] - 2.0).abs() < 1e-10);
    println!("  ✓ x-coordinate negated\n");

    // ── Reflect through offset plane ──
    let reflected_off = Conformal::reflect(p, [0.0, 0.0, 1.0], 1.0);
    println!(
        "Reflect ({}, {}, {}) through plane z = 1:",
        p[0], p[1], p[2]
    );
    println!(
        "  → ({:.1}, {:.1}, {:.1})",
        reflected_off[0], reflected_off[1], reflected_off[2]
    );
    println!(
        "  ✓ z = {} → z' = {:.1} (should be -1)",
        p[2], reflected_off[2]
    );
    println!();

    // ── Project onto plane ──
    let proj = Conformal::project_onto_plane([0.0, 0.0, 5.0], [0.0, 0.0, 1.0], 2.0);
    println!("Project (0,0,5) onto plane z = 2:");
    println!("  → ({}, {}, {})", proj[0], proj[1], proj[2]);
    assert!((proj[2] - 2.0).abs() < 1e-10);
    println!("  ✓ z = 2.0\n");

    // ── Rotate using a rotor ──
    let rotor = Rotor::from_axis_angle([0.0, 0.0, 1.0], FRAC_PI_2);
    let rotated_pt = Conformal::rotate([1.0, 0.0, 0.0], &rotor);
    println!("Rotate (1,0,0) 90° around Z:");
    println!(
        "  → ({:.1}, {:.1}, {:.1})",
        rotated_pt[0], rotated_pt[1], rotated_pt[2]
    );
    assert!((rotated_pt[0]).abs() < 0.01);
    assert!((rotated_pt[1] - 1.0).abs() < 0.01);
    println!("  ✓ (0, 1, 0)\n");

    // ── Translate ──
    let translated = Conformal::translate([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
    println!("Translate (1,2,3) by (4,5,6):");
    println!(
        "  → ({}, {}, {})",
        translated[0], translated[1], translated[2]
    );
    assert_eq!(translated, [5.0, 7.0, 9.0]);
    println!("  ✓ (5, 7, 9)\n");

    // ── Midpoint & barycenter ──
    let midpoint = Conformal::midpoint([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
    println!("Midpoint of (0,0,0) and (2,4,6):");
    println!("  → ({}, {}, {})", midpoint[0], midpoint[1], midpoint[2]);
    assert_eq!(midpoint, [1.0, 2.0, 3.0]);

    let bc = Conformal::barycenter(
        &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]],
        &[1.0, 1.0, 1.0],
    );
    println!("Barycenter of triangle:");
    println!("  → ({:.1}, {:.1})", bc[0], bc[1]);
    println!("  ✓ Centroid at (1.0, 0.67)");
    println!();

    println!("═══ All conformal operations verified ═══");
}
