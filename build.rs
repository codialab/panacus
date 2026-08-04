use shadow_rs::ShadowBuilder;
use std::process::Command;

fn main() {
    // note: add error checking yourself.
    ShadowBuilder::builder()
        .deny_const(Default::default())
        .build()
        .unwrap();
    let output = Command::new("git")
        .args(["rev-parse", "--short", "HEAD"])
        .output()
        .unwrap();
    let git_hash = String::from_utf8(output.stdout).unwrap();
    println!("cargo:rustc-env=GIT_HASH={}", git_hash);
}
