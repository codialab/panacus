use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn info_table_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("info").arg("tests/test_files/t_groups.gfa");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("table-name\tkey\tvalue"));
    Ok(())
}

#[test]
fn info_html_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("report").arg("tests/test_files/info.yaml");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("<th scope=\"col\">Statistic</th>"))
        .stdout(predicate::str::contains("table-name\tkey\tvalue"));
    Ok(())
}

#[test]
fn info_table_groups_get_written() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("info")
        .arg("tests/test_files/t_groups.gfa")
        .arg("-S");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("group-lengths\tx\t50"))
        .stdout(predicate::str::contains("group-lengths\ty\t50"));
    Ok(())
}
