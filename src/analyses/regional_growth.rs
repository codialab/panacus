use anyhow::anyhow;
use linreg::linear_regression;
use nalgebra::DVector;
use std::cell::OnceCell;
use varpro::prelude::*;
use varpro::{problem::SeparableProblemBuilder, solvers::levmar::LevMarSolver};

use crate::analyses::growth::calc_growth;
use crate::analyses::regional_helpers::sort_values;
use crate::util::Threshold;
use crate::{
    analyses::MatrixBasedAnalysis,
    coverage_matrix::CoverageMatrix,
    hist::Hist,
    html_report::{AnalysisSection, ReportItem, Window},
    util::get_default_plot_downloads,
};

pub struct RegionalGrowth {
    order: Option<String>,
    window_size: usize,
    windows: Option<String>,
    cache: OnceCell<Vec<(String, Vec<(Growth, usize, usize)>)>>,
}

impl MatrixBasedAnalysis for RegionalGrowth {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        let data = self.get_cached_data(matrix);

        let mut text = "reference\tstart\tend\tHill0\tHill1\tHill2\n".to_string();
        for (sequence_id, sequence) in data {
            for (variation, start, end) in sequence {
                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\n",
                    sequence_id, start, end, variation.0, variation.1, variation.2
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
            .replace([' ', '|', '\\'], "-");
        let id_prefix = format!("regional-growth-{}", id_prefix);
        let (mut sequences, mut values): (Vec<String>, Vec<Vec<Window>>) = data
            .iter()
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
                "2x increase".to_string(),
                "3x increase".to_string(),
                "5x increase".to_string(),
            ],
            is_diverging: true,
            contains_outliers: false,
            sequence: sequences,
            values,
        }];

        let table_text = self.generate_table(matrix)?;
        let table_text = format!("`{}`", table_text);
        let regional_variation_tabs = vec![AnalysisSection {
            id: id_prefix.to_string(),
            analysis: "Regional Growth".to_string(),
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
        "RegionalGrowth".to_string()
    }
}

impl RegionalGrowth {
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
    ) -> &Vec<(String, Vec<(Growth, usize, usize)>)> {
        self.cache.get_or_init(|| {
            let mut error_counter = 0;
            let result = matrix
                .get_regional_hists(self.window_size, self.window_size, self.windows.as_deref())
                .map(|(r, i)| {
                    let windows = i
                        .filter_map(|(start, end, hist)| {
                            if hist.is_empty() {
                                return None;
                            }
                            let result = calc_growth_predictions(&hist);
                            match result {
                                Ok((d2, d3, d5)) => Some((Growth(d2, d3, d5), start, end)),
                                Err(_) => {
                                    error_counter += 1;
                                    None
                                }
                            }
                        })
                        .collect::<Vec<(Growth, usize, usize)>>();
                    (r, windows)
                })
                .collect::<Vec<_>>();
            if error_counter > 0 {
                log::warn!("Had problems in the fitting of {} windows", error_counter);
            }
            result
        })
    }
}

/// Simple triple struct holding
/// Percent increase in terms of feature if haplotypes increase 2x
/// Percent increase in terms of feature if haplotypes increase 3x
/// Percent increase in terms of feature if haplotypes increase 5x
#[derive(Debug)]
struct Growth(f64, f64, f64);

pub fn calc_growth_predictions(hist: &Hist) -> anyhow::Result<(f64, f64, f64)> {
    let d = hist.get_hist_values().iter().sum::<usize>() as f64;
    let n = hist.get_hist_values().len();
    let h1_diff = get_h1_difference(hist)?;
    let gte = calc_good_toulmin(hist, 1.8);

    let mut modified_hist = hist.clone();
    modified_hist.correct_h1(h1_diff as usize);
    let (k, alpha, f_new) = fit_fnew(&modified_hist)?;

    let d_hat1: f64 = (1..n)
        .map(|x| {
            if x < n / 2 {
                f_new[x - 1]
            } else {
                k * (x as f64).powf(-alpha)
            }
        })
        .sum::<f64>()
        + 1.0 * h1_diff;

    let d_hat1_8: f64 = (1..(1.8 * (n as f64)) as usize)
        .map(|x| {
            if x < n / 2 {
                f_new[x - 1]
            } else {
                k * (x as f64).powf(-alpha)
            }
        })
        .sum::<f64>()
        + 1.8 * h1_diff;

    let d_hat2: f64 = (1..2 * n)
        .map(|x| {
            if x < n / 2 {
                f_new[x - 1]
            } else {
                k * (x as f64).powf(-alpha)
            }
        })
        .sum::<f64>()
        + 2.0 * h1_diff;

    let d_hat3: f64 = (1..3 * n)
        .map(|x| {
            if x < n / 2 {
                f_new[x - 1]
            } else {
                k * (x as f64).powf(-alpha)
            }
        })
        .sum::<f64>()
        + 3.0 * h1_diff;

    let d_hat5: f64 = (1..5 * n)
        .map(|x| {
            if x < n / 2 {
                f_new[x - 1]
            } else {
                k * (x as f64).powf(-alpha)
            }
        })
        .sum::<f64>()
        + 5.0 * h1_diff;

    if (d - d_hat1).abs() > 0.1 * d {
        return Err(anyhow!(
            "Too much difference in model ({}) and ground-truth ({})",
            d_hat1,
            d
        ));
    }

    if (gte - d_hat1_8).abs() > 0.1 * gte {
        return Err(anyhow!(
            "Too much difference in parametric ({}) and non-parametric ({}) model",
            d_hat2,
            gte
        ));
    }

    let d2_ratio = d_hat2 / d;
    let d3_ratio = d_hat3 / d;
    let d5_ratio = d_hat5 / d;

    Ok((d2_ratio, d3_ratio, d5_ratio))
}

/// Returns K and alpha
pub fn fit_fnew(corrected_hist: &Hist) -> anyhow::Result<(f64, f64, Vec<f64>)> {
    let t_coverage = Threshold::Absolute(1);
    let t_quorum = Threshold::Relative(0.0);
    let growths = calc_growth(corrected_hist, &t_coverage, &t_quorum);
    let mut f_new: Vec<f64> = (0..growths.len())
        .map(|idx| {
            if idx == 0 {
                growths[0]
            } else {
                growths[idx] - growths[idx - 1]
            }
        })
        .collect();
    let res_fnew = f_new.clone();
    // Filter out all zeroes as these will make problems when taking the ln
    let mut f_new_idx: Vec<f64> = (1..=f_new.len())
        .filter_map(|x| {
            if f_new[x - 1] > 0.0 {
                Some(x as f64)
            } else {
                None
            }
        })
        .collect();
    for idx in f_new_idx.iter_mut() {
        *idx = idx.ln();
    }
    f_new.retain(|&x| x > 0.0);
    for v in f_new.iter_mut() {
        *v = v.ln();
    }

    if f_new.is_empty() {
        return Err(anyhow!("No elements inside F_new (likely the histogram for this window only contained singletons)"));
    }
    let (alpha, k): (f64, f64) = linear_regression(&f_new_idx[1..], &f_new[1..])
        .map_err(|_| anyhow!("Could not fit linear regression"))?;
    let alpha = -alpha;
    let k = k.exp();
    Ok((k, alpha, res_fnew))
}

pub fn power_law(x: &DVector<f64>, gamma: f64) -> DVector<f64> {
    x.map(|x| x.powf(-gamma))
}

pub fn power_law_dgamma(gvec: &DVector<f64>, gamma: f64) -> DVector<f64> {
    gvec.map(|g| -g.ln() * g.powf(-gamma))
}

/// Gets the difference between true h1 and real h1 by
/// fitting a power law curve to the first half of the histogram
pub fn get_h1_difference(hist: &Hist) -> anyhow::Result<f64> {
    let y = hist.get_hist_values();
    let num_values = y.len() / 2;
    let x = hist.get_index().take(num_values).skip(2);
    let y = &y[2..num_values];
    let x = DVector::from_iterator(y.len(), x.map(|v| v as f64));
    let y = DVector::from_iterator(y.len(), y.iter().map(|v| *v as f64));
    let model = SeparableModelBuilder::<f64>::new(&["gamma"])
        .function(&["gamma"], power_law)
        .partial_deriv("gamma", power_law_dgamma)
        .independent_variable(x)
        .initial_parameters(vec![1.0])
        .build()
        .unwrap();

    let problem = SeparableProblemBuilder::new(model)
        .observations(y)
        .build()
        .unwrap();

    let fit_result = LevMarSolver::default()
        .fit(problem)
        .map_err(|_| anyhow!("Fit did not suceed"))?;
    let _gamma: f64 = fit_result.nonlinear_parameters()[0];
    let c: f64 = fit_result
        .linear_coefficients()
        .ok_or_else(|| anyhow!("Could not get linear coefficients"))?[0];
    let h1_corr = c;
    Ok(hist[1] as f64 - h1_corr)
}

pub fn calc_good_toulmin(hist: &Hist, lambda: f64) -> f64 {
    let d = hist.get_hist_values().iter().copied().sum::<usize>() as f64;
    let d_new = hist
        .get_index()
        .zip(hist.get_hist_values().iter().copied())
        .map(|(i, h)| (1.0f64 - lambda).powi(i as i32) * (h as f64))
        .sum::<f64>();

    d - d_new
}
