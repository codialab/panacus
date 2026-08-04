use std::cell::OnceCell;

use crate::{
    analyses::{regional_helpers::sort_values, MatrixBasedAnalysis},
    coverage_matrix::CoverageMatrix,
    hist::Hist,
    html_report::{AnalysisSection, ReportItem, Window},
    util::get_default_plot_downloads,
};

pub struct RegionalVariation {
    order: Option<String>,
    window_size: usize,
    windows: Option<String>,
    cache: OnceCell<Vec<(String, Vec<(Variation, usize, usize)>)>>,
}

impl MatrixBasedAnalysis for RegionalVariation {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        let data = self.get_cached_data(matrix);

        let mut text = format!("reference\tstart\tend\tHill0\tHill1\tHill2\n");
        for (sequence_id, sequence) in data {
            for (variation, start, end) in sequence {
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    sequence_id.to_string(),
                    start,
                    end,
                    variation.0,
                    variation.1,
                    variation.2
                );
                text.push_str(&line);
            }
        }
        Ok(text)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let data = self.get_cached_data(matrix);

        let id_prefix = matrix
            .get_run_id()
            .to_lowercase()
            .replace(&[' ', '|', '\\'], "-");
        let id_prefix = format!("regional-variation-{}", id_prefix);
        let (mut sequences, mut values): (Vec<String>, Vec<Vec<Window>>) = data
            .into_iter()
            .map(|(sequence, values)| {
                let values = values
                    .iter()
                    .map(|(v, start, end)| Window {
                        start: *start,
                        end: *end,
                        values: vec![v.0, v.1, v.2],
                    })
                    .collect();
                (sequence.clone(), values)
            })
            .unzip();

        self.sort_values(&mut sequences, &mut values);

        let chromosomal = vec![ReportItem::Chromosomal {
            id: format!("{id_prefix}-Chrom"),
            name: "PLACEHOLDER".to_string(),
            labels: vec![
                "Richness".to_string(),
                "Shannon entropy".to_string(),
                "Simpson index".to_string(),
            ],
            is_diverging: true,
            contains_outliers: false,
            sequence: sequences,
            values,
        }];

        let table_text = self.generate_table(matrix)?;
        let table_text = format!("`{}`", table_text);
        let regional_variation_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}"),
            analysis: "Regional Variation".to_string(),
            table: Some(table_text),
            run_name: matrix.get_run_name().to_string(),
            run_id: id_prefix,
            countable: matrix.get_feature_type().to_string(),
            items: chromosomal,
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(regional_variation_tabs)
    }

    fn get_type(&self) -> String {
        "RegionalVariation".to_string()
    }
}

impl RegionalVariation {
    pub fn new(window_size: usize, order: Option<String>, windows: Option<String>) -> Self {
        Self {
            order,
            window_size,
            windows,
            cache: OnceCell::new(),
        }
    }

    fn sort_values(&self, reference_names: &mut Vec<String>, windows: &mut Vec<Vec<Window>>) {
        if let Some(filename) = self.order.as_ref() {
            sort_values(filename, reference_names, windows);
        }
    }

    fn get_cached_data(
        &self,
        matrix: &CoverageMatrix,
    ) -> &Vec<(String, Vec<(Variation, usize, usize)>)> {
        self.cache.get_or_init(|| {
            matrix
                .get_regional_hists(self.window_size, self.window_size, self.windows.as_deref())
                .map(|(r, i)| {
                    (
                        r,
                        i.map(|(start, end, hist)| {
                            let richness = calc_richness(&hist);
                            let shannon = calc_shannon_entropy(&hist);
                            let simpson = calc_simpson_index(&hist);
                            (Variation(richness, shannon, simpson), start, end)
                        })
                        .collect::<Vec<(Variation, usize, usize)>>(),
                    )
                })
                .collect::<Vec<_>>()
        })
    }
}

/// Simple triple struct holding
/// Hill number 0: count of features (sum of hist), Richness
/// Hill number 1: Shannon entropy
/// Hill number 2: Simpson index
#[derive(Debug)]
struct Variation(f64, f64, f64);

/// Calculates the (species) richness (0th hill number) for a given hist
pub fn calc_richness(hist: &Hist) -> f64 {
    hist.get_hist_values().iter().skip(1).sum::<usize>() as f64
}

/// Calculates the Shannon entropy (1st hill number)
/// for a given hist
pub fn calc_shannon_entropy(hist: &Hist) -> f64 {
    let u = hist
        .get_hist_values()
        .iter()
        .enumerate()
        .skip(1)
        .map(|(j, q)| j as f64 * (*q as f64))
        .sum::<f64>();
    hist.get_hist_values()
        .iter()
        .enumerate()
        .skip(1)
        .map(|(k, &q)| {
            let k = k as f64;
            let q = q as f64;
            (-k / u) * (k / u).ln() * q
        })
        .sum::<f64>()
        .exp()
}

/// Calculates the Simpson index (2nd hill number)
/// for a given hist
pub fn calc_simpson_index(hist: &Hist) -> f64 {
    let u = hist
        .get_hist_values()
        .iter()
        .enumerate()
        .skip(1)
        .map(|(j, q)| j as f64 * (*q as f64))
        .sum::<f64>();
    1.0 / hist
        .get_hist_values()
        .iter()
        .enumerate()
        .skip(1)
        .map(|(k, &q)| {
            let k = k as f64;
            let q = q as f64;
            (k / u).powi(2) * q
        })
        .sum::<f64>()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_approx_equal(x: f64, y: f64, delta: f64) -> bool {
        (x - y).abs() < delta
    }

    #[test]
    fn test_richness_calculation() {
        let h = Hist::from_values(vec![2, 3, 2, 0, 4]);
        let r = calc_richness(&h);
        assert_eq!(r, 9.0);
    }

    #[test]
    fn test_shannon_entropy() {
        let h = Hist::from_values(vec![2, 3, 2, 0, 4]);
        let s = calc_shannon_entropy(&h);
        assert!(is_approx_equal(s, 7.77, 0.01));
    }

    #[test]
    fn test_simpson_index() {
        let h = Hist::from_values(vec![2, 3, 2, 0, 4]);
        let s = calc_simpson_index(&h);
        assert!(is_approx_equal(s, 7.05, 0.01));
    }
}
