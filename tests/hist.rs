use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn hist_table_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("hist").arg("tests/test_files/t_groups.gfa");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("panacus\thist"));
    Ok(())
}

#[test]
fn hist_html_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("report").arg("tests/test_files/hist.yaml");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("<div id=\"chart-bar-cov-hist"))
        .stdout(predicate::str::contains("panacus\thist"));
    Ok(())
}
