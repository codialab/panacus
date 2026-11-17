use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn histgrowth_table_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("histgrowth")
        .arg("tests/test_files/t_groups.gfa")
        .arg("--hist");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("panacus\thist\tgrowth"));
    Ok(())
}
