use assert_cmd::prelude::*; // Add methods on commands
use flate2::{write::GzEncoder, Compression};
use predicates::prelude::*; // Used for writing assertions
use std::io::Write;
use std::process::Command; // Run programs
use tempfile::Builder;

#[test]
fn growth_table_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("growth")
        .arg("tests/test_files/t_groups.hist.tsv")
        .arg("--hist");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("panacus\thist\tgrowth"));
    Ok(())
}

#[test]
fn growth_reads_coverage_header_without_zero_bin() -> Result<(), Box<dyn std::error::Error>> {
    let mut file = Builder::new().suffix(".tsv").tempfile()?;
    writeln!(file, "panacus\thist\ncoverage\tnode\n1\t21604\n2\t10241")?;
    let mut cmd = Command::cargo_bin("panacus")?;
    cmd.arg("growth").arg(file.path()).arg("--hist");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("0\t0\t0\n"))
        .stdout(predicate::str::contains("2\t10241\t31845\n"));
    Ok(())
}

#[test]
fn growth_reads_compressed_histogram_and_selected_column() -> Result<(), Box<dyn std::error::Error>>
{
    let mut file = Builder::new().suffix(".tsv.gz").tempfile()?;
    let mut encoder = GzEncoder::new(&mut file, Compression::default());
    encoder.write_all(b"panacus\thist\thist\ncount\tnode\tbp\n0\t5\t20\n1\t2\t10\n2\t3\t30\n")?;
    encoder.finish()?;
    Command::cargo_bin("panacus")?
        .arg("growth")
        .arg(file.path())
        .args(["--count", "bp", "--hist"])
        .assert()
        .success()
        .stdout(predicate::str::contains("count\tbp\tbp\n"))
        .stdout(predicate::str::contains("2\t30\t40\n"));
    Ok(())
}

#[test]
fn growth_from_generated_histogram_matches_gfa() -> Result<(), Box<dyn std::error::Error>> {
    let graph = "tests/test_files/t_groups.gfa";
    let histogram = Command::cargo_bin("panacus")?
        .arg("hist")
        .args(["--count", "bp"])
        .arg(graph)
        .assert()
        .success()
        .get_output()
        .stdout
        .clone();
    let mut file = Builder::new().suffix(".tsv").tempfile()?;
    file.write_all(&histogram)?;
    let mut outputs = Vec::new();
    for input in [std::path::Path::new(graph), file.path()] {
        let output = Command::cargo_bin("panacus")?
            .arg("growth")
            .arg(input)
            .args([
                "--coverage",
                "1,2",
                "--quorum",
                "0,0.5",
                "--hist",
                "--count",
                "bp",
            ])
            .assert()
            .success()
            .get_output()
            .stdout
            .clone();
        outputs.push(
            String::from_utf8(output)?
                .lines()
                .filter(|line| !line.starts_with('#'))
                .collect::<Vec<_>>()
                .join("\n"),
        );
    }
    assert_eq!(outputs[0], outputs[1]);
    Ok(())
}

#[test]
fn growth_rejects_invalid_tsv_and_gfa_only_options() -> Result<(), Box<dyn std::error::Error>> {
    let mut file = Builder::new().suffix(".tsv").tempfile()?;
    writeln!(file, "panacus\thist\ncoverage\tnode\n1\tbad")?;
    Command::cargo_bin("panacus")?
        .arg("growth")
        .arg(file.path())
        .assert()
        .failure()
        .stderr(predicate::str::contains("Invalid histogram count"));
    Command::cargo_bin("panacus")?
        .arg("growth")
        .arg("tests/test_files/t_groups.hist.tsv")
        .arg("--groupby-sample")
        .assert()
        .failure()
        .stderr(predicate::str::contains("require GFA input"));
    Ok(())
}

#[ignore]
#[test]
fn growth_html_gets_written_successfully() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("panacus")?;

    cmd.arg("growth")
        .arg("tests/test_files/t_groups.hist.tsv")
        .arg("-o")
        .arg("html");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains(
            "<canvas id=\"chart-bar-pan-growth-node\"></canvas>",
        ))
        .stdout(predicate::str::contains("panacus\thist\tgrowth"));
    Ok(())
}
