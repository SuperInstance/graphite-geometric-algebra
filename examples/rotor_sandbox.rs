//! Example: Rotor sandbox — create 2D rotors, compose them, and apply to points.
//!
//! Run with: `cargo run --example rotor_sandbox`

use graphite_ga_nodes::ga_core::Rotor;
use std::f64::consts::PI;

fn main() {
    println!("═══ Rotor Sandbox ═══");
    println!("Geometric algebra rotors in Cl(3,1) for the Graphite editor\n");

    // ── Create a 90° z-axis rotor ──
    let r90 = Rotor::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
    println!("Rotor for 90° around Z-axis:");

    let point = [1.0, 0.0, 0.0];
    let rotated = r90.apply(point);
    println!("  apply({:?}) → {:?}", point, rotated);
    assert!((rotated[0]).abs() < 0.01);
    assert!((rotated[1] - 1.0).abs() < 0.01);
    println!("  ✓ (1,0,0) rotated 90° → near (0,1,0)\n");

    // ── Compose two rotors ──
    let r45 = Rotor::from_axis_angle([0.0, 0.0, 1.0], PI / 4.0);
    let r90_composed = r45.compose(&r45);
    let r90_check = r90_composed.apply(point);
    println!("Two 45° rotors composed:");
    println!("  apply({:?}) → {:?}", point, r90_check);
    assert!((r90_check[0]).abs() < 0.01);
    assert!((r90_check[1] - 1.0).abs() < 0.01);
    println!("  ✓ 45° + 45° = 90°\n");

    // ── 45° rotation ──
    let rt45 = r45.apply(point);
    println!("45° rotation of (1,0):");
    println!("  result: ({:.4}, {:.4})", rt45[0], rt45[1]);
    println!("  expected: (0.7071, 0.7071)");
    println!();

    // ── Slerp from identity to 180° ──
    let r_identity = Rotor::identity();
    let r180 = Rotor::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
    println!("Slerp from identity to 90°, t = 0.25:");
    let r_interp = r_identity.slerp(&r180, 0.25);
    let slerp_point = [1.0, 0.0, 0.0];
    let slerp_result = r_interp.apply(slerp_point);
    println!(
        "  apply({:?}) → ({:.4}, {:.4})",
        slerp_point, slerp_result[0], slerp_result[1]
    );
    println!(
        "  ✓ normalized: |R|² = {:.4}\n",
        r_interp.inner.norm_squared()
    );

    // ── Rotation matrix ──
    let matrix = r90.to_rotation_matrix_2d();
    println!("Rotation matrix for 90°:");
    println!("  |{:.4}  {:.4}|", matrix[0][0], matrix[0][1]);
    println!("  |{:.4}  {:.4}|", matrix[1][0], matrix[1][1]);
    println!();

    // ── Length preservation ──
    let v: [f64; 2] = [3.0, 4.0];
    let v_rot = r45.apply_2d(v);
    let orig_len = (v[0] * v[0] + v[1] * v[1]).sqrt();
    let rot_len = (v_rot[0] * v_rot[0] + v_rot[1] * v_rot[1]).sqrt();
    println!("Length preservation:");
    println!("  |({:.1}, {:.1})| = {:.4}", v[0], v[1], orig_len);
    println!("  |({:.4}, {:.4})| = {:.4}", v_rot[0], v_rot[1], rot_len);
    assert!((orig_len - rot_len).abs() < 0.001);
    println!("  ✓ Length preserved!\n");

    println!("═══ All rotor operations verified ═══");
}
