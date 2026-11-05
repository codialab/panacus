use std::fs::File;
use std::io::{BufRead, BufReader, Read};

use anyhow::{anyhow, ensure, Context};
use flate2::read::MultiGzDecoder;

use crate::{coverage_matrix::CoverageMatrix, hist::Hist, util::CountType};

use super::FileFormatParser;

pub struct HistParser {
    hist: Hist,
}

impl HistParser {
    pub fn new(filename: &str, count_type: CountType) -> anyhow::Result<Self> {
        log::info!("loading histogram from {}", filename);
        let file = File::open(filename)
            .with_context(|| format!("Could not open histogram {}", filename))?;
        let reader: Box<dyn Read> = if filename.ends_with(".gz") {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let hist = parse_hist(BufReader::new(reader), count_type, filename)
            .with_context(|| format!("Could not parse histogram {}", filename))?;
        Ok(Self { hist })
    }
}

impl FileFormatParser for HistParser {
    fn generate_hist(self: Box<Self>) -> Hist {
        self.hist
    }

    fn generate_matrix(self: Box<Self>) -> CoverageMatrix {
        panic!("Histogram input cannot generate a coverage matrix")
    }
}

fn parse_hist<R: BufRead>(
    reader: R,
    count_type: CountType,
    filename: &str,
) -> anyhow::Result<Hist> {
    let mut lines = reader
        .lines()
        .enumerate()
        .filter_map(|(index, line)| match line {
            Ok(line) if line.trim().is_empty() || line.starts_with('#') => None,
            line => Some((index + 1, line)),
        });
    let (_, header) = lines.next().ok_or_else(|| anyhow!("Empty histogram"))?;
    let header = header?;
    let columns: Vec<_> = header.split('\t').collect();
    ensure!(
        columns.first() == Some(&"panacus") && columns.contains(&"hist"),
        "Expected a tab-separated 'panacus\\thist' header"
    );
    let (line_number, counts) = lines
        .next()
        .ok_or_else(|| anyhow!("Missing histogram count header"))?;
    let counts = counts?;
    let counts: Vec<_> = counts.split('\t').collect();
    ensure!(
        counts.len() == columns.len() && matches!(counts[0], "count" | "coverage"),
        "Invalid histogram count header at line {}",
        line_number
    );
    let feature_type = count_type.to_string();
    let column = columns
        .iter()
        .zip(&counts)
        .position(|(kind, count)| *kind == "hist" && *count == feature_type)
        .ok_or_else(|| anyhow!("No '{}' histogram column found", feature_type))?;

    let mut values = Vec::new();
    let mut total = 0usize;
    for (line_number, line) in lines {
        let line = line?;
        let fields: Vec<_> = line.split('\t').collect();
        ensure!(
            fields.len() == columns.len(),
            "Expected {} columns at line {}",
            columns.len(),
            line_number
        );
        // Growth tables can include coverage and quorum metadata above the histogram.
        if values.is_empty() && matches!(fields[0], "coverage" | "quorum") {
            continue;
        }
        let coverage = fields[0]
            .parse::<usize>()
            .with_context(|| format!("Invalid coverage at line {}", line_number))?;
        let count = fields[column]
            .parse::<usize>()
            .with_context(|| format!("Invalid histogram count at line {}", line_number))?;
        total = total
            .checked_add(count)
            .ok_or_else(|| anyhow!("Histogram count total is too large at line {}", line_number))?;
        if let Some(&(previous, _)) = values.last() {
            ensure!(
                coverage > previous,
                "Coverage must be strictly increasing at line {}",
                line_number
            );
        }
        values.push((coverage, count));
    }
    let maximum_coverage = values
        .last()
        .map(|&(coverage, _)| coverage)
        .ok_or_else(|| anyhow!("Histogram contains no data rows"))?;
    maximum_coverage
        .checked_add(1)
        .ok_or_else(|| anyhow!("Histogram coverage is too large"))?;
    let mut hist = Hist::from_maximum_coverage(
        maximum_coverage,
        feature_type,
        filename.to_string(),
        filename.to_string(),
    );
    for (coverage, count) in values {
        hist.insert_feature_of_coverage_and_length(coverage, count);
    }
    Ok(hist)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_coverage_header_with_missing_zero_and_gaps() {
        let input = b"panacus\thist\ncoverage\tnode\n1\t21604\n3\t6941\n";
        let hist = parse_hist(&input[..], CountType::Node, "test.tsv").unwrap();
        assert_eq!(hist.get_hist_values(), &[0, 21604, 0, 6941]);
        assert_eq!(hist.get_number_of_features(), 28545);
    }

    #[test]
    fn selects_histogram_column_and_skips_metadata() {
        let input = b"# panacus hist\npanacus\thist\thist\tgrowth\ncount\tnode\tbp\tnode\ncoverage\t\t\t1\nquorum\t\t\t0\n0\t5\t20\t0\n1\t2\t10\t2\n";
        let hist = parse_hist(&input[..], CountType::Bp, "test.tsv").unwrap();
        assert_eq!(hist.get_hist_values(), &[20, 10]);
        assert_eq!(hist.get_feature_type(), "bp");
    }

    #[test]
    fn rejects_invalid_histograms() {
        for input in [
            "panacus\thist\ncount\tnode\n",
            "panacus\thist\ncount\tbp\n1\t2\n",
            "panacus\thist\ncount\tnode\n1\tbad\n",
            "panacus\thist\ncount\tnode\n1\t2\n1\t3\n",
            "panacus\thist\ncount\tnode\n2\t2\n1\t3\n",
            "panacus hist\ncoverage node\n1 2\n",
        ] {
            assert!(parse_hist(input.as_bytes(), CountType::Node, "test.tsv").is_err());
        }
    }

    #[test]
    fn rejects_histogram_count_overflow() {
        let input = format!("panacus\thist\ncount\tnode\n1\t{}\n2\t1\n", usize::MAX);
        let error = parse_hist(input.as_bytes(), CountType::Node, "test.tsv").unwrap_err();
        assert!(error.to_string().contains("count total is too large"));
    }
}
