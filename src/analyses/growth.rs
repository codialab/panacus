use std::cmp;

use ml_helpers::linear_regression::huber_regressor::{solve, HuberRegressor};

use crate::analyses::HistBasedAnalysis;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::hist::Hist;
use crate::html_report::ReportItem;
use crate::io::write_table;
use crate::util::{get_default_plot_downloads, Threshold};
use rayon::prelude::*;

use super::AnalysisSection;

pub struct Growth {
    coverage: Option<String>,
    quorum: Option<String>,
    add_hist: bool,
    add_alpha: bool,
    growths: Option<Vec<Vec<f64>>>,
    thresholds: Option<ThresholdContainer>,
}

impl HistBasedAnalysis for Growth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(&mut self, hist: &Hist) -> anyhow::Result<String> {
        log::info!(
            "reporting hist table, {}, {}",
            self.add_hist,
            self.add_alpha
        );

        self.set_inner(hist)?;
        let growths = self.growths.as_ref().unwrap();
        let thresholds = self.thresholds.as_ref().unwrap();

        let mut res = String::new();
        res.push_str(&format!(
            "# {}\n",
            std::env::args().collect::<Vec<String>>().join(" ")
        ));
        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            "coverage".to_string(),
            "quorum".to_string(),
        ]];
        let mut output_columns: Vec<Vec<f64>> = Vec::new();

        let hist_values = hist.get_hist_values();

        if self.add_hist {
            output_columns.push(hist_values.iter().map(|x| *x as f64).collect());
            header_cols.push(vec![
                "hist".to_string(),
                hist.get_feature_type().to_string(),
                String::new(),
                String::new(),
            ])
        }

        let growths: Vec<Vec<f64>> = growths
            .iter()
            .map(|row| {
                row.iter()
                    .map(|el| if el.is_nan() { 0.0 } else { *el })
                    .collect()
            })
            .collect();
        output_columns.extend(growths);
        let m = thresholds.coverage.len();
        header_cols.extend(
            std::iter::repeat("growth")
                .take(m)
                .zip(std::iter::repeat(hist.get_feature_type()).take(m))
                .zip(&thresholds.coverage)
                .zip(&thresholds.quorum)
                .map(|(((p, t), c), q)| {
                    vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                }),
        );
        res.push_str(&write_table(&header_cols, &output_columns)?);
        Ok(res)
    }

    fn generate_report_section(&mut self, hist: &Hist) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(hist)?;
        let hist_aux = &self.thresholds.as_ref().unwrap();
        let growth_labels = (0..hist_aux.coverage.len())
            .map(|i| {
                format!(
                    "coverage ≥ {}, quorum ≥ {}%",
                    hist_aux.coverage[i].get_string(),
                    match hist_aux.quorum[i] {
                        crate::util::Threshold::Relative(x) => (x * 100.0).to_string(),
                        crate::util::Threshold::Absolute(x) => (x * 100).to_string(),
                    }
                )
            })
            .collect::<Vec<_>>();
        let table = self.generate_table(hist)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "pan-growth-{}",
            hist.get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        let k = hist.get_feature_type();
        let growths = self.growths.as_ref().unwrap();
        let growth_tabs = vec![AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Pangenome Growth".to_string(),
            run_name: hist.get_run_name().to_owned(),
            run_id: hist.get_run_id().to_owned(),
            countable: hist.get_feature_type().to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::MultiBar {
                id: format!("{id_prefix}-{k}"),
                names: growth_labels.clone(),
                x_label: "taxa".to_string(),
                y_label: format!("#{}s", k),
                labels: (0..(growths[0].len())).map(|i| i.to_string()).collect(),
                values: growths
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|el| if el.is_nan() { 0.0 } else { *el })
                            .collect()
                    })
                    .collect(),
                curve: None,
                alpha: None,
                log_toggle: false,
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(growth_tabs)
    }
}

impl Growth {
    pub fn new(
        coverage: Option<String>,
        quorum: Option<String>,
        add_hist: bool,
        add_alpha: bool,
    ) -> Self {
        Self {
            coverage,
            quorum,
            add_hist,
            add_alpha,
            growths: None,
            thresholds: None,
        }
    }

    fn set_inner(&mut self, hist: &Hist) -> anyhow::Result<()> {
        if self.growths.is_some() {
            return Ok(());
        }
        let quorum = self.quorum.to_owned().unwrap_or("0".to_string());
        let coverage = self.coverage.to_owned().unwrap_or("1".to_string());
        let hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;

        let growths = calc_all_growths(hist, &hist_aux, true);
        self.growths = Some(growths);
        self.thresholds = Some(hist_aux);
        Ok(())
    }
}

#[allow(dead_code)]
fn get_regression(hist: &[f64]) -> (f64, f64) {
    let x: Vec<f64> = (1..hist.len()).map(|x| x as f64).collect();
    let log_x: Vec<f64> = x
        .iter()
        .map(|x| x.log10())
        .map(|x| {
            if x.is_infinite() && x.is_sign_negative() {
                -100000.0
            } else {
                x
            }
        })
        .collect();
    let y = hist.iter().skip(1).copied().collect::<Vec<f64>>();
    let log_y: Vec<f64> = y
        .iter()
        .map(|y| y.log10())
        .map(|y| {
            if y.is_infinite() && y.is_sign_negative() {
                -100000.0
            } else {
                y
            }
        })
        .collect();
    let n = cmp::max(hist.len() / 2, 10);
    let log_x2: Vec<f64> = log_x.iter().skip(1).take(n).copied().collect();
    let log_y2: Vec<f64> = log_y.iter().skip(1).take(n).copied().collect();
    let huber = HuberRegressor::from(log_x2.clone(), log_y2.clone());
    let params = solve(huber);
    let alpha = 2.0 + params[0];
    (alpha, params[1])
}

pub fn calc_growth(hist: &Hist, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage();

    if n > 0 {
        let quorum = usize::max(1, t_quorum.to_absolute(n));
        if quorum == 1 {
            calc_growth_union(hist, t_coverage)
        } else if quorum >= n {
            calc_growth_core(hist, t_coverage)
        } else {
            calc_growth_quorum(hist, t_coverage, t_quorum)
        }
    } else {
        Vec::new()
    }
}

pub fn calc_all_growths(
    hist: &Hist,
    hist_aux: &ThresholdContainer,
    insert_zero: bool,
) -> Vec<Vec<f64>> {
    let mut growths: Vec<Vec<f64>> = hist_aux
        .coverage
        .par_iter()
        .zip(&hist_aux.quorum)
        .map(|(c, q)| {
            log::info!(
                "calculating growth for coverage >= {} and quorum >= {}",
                &c,
                &q
            );
            calc_growth(hist, c, q)
        })
        .collect();
    // insert empty row for 0 element
    if insert_zero {
        for g in &mut growths {
            g.insert(0, f64::NAN);
        }
    }
    growths
}

fn calc_growth_union(hist: &Hist, t_coverage: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from 0..n (both included)
    let c = usize::max(1, t_coverage.to_absolute(n));

    let mut pangrowth: Vec<f64> = vec![0.0; n];
    let mut n_fall_m: f64 = 0.0;
    let tot = hist.get_hist_values()[c..].iter().sum::<usize>() as f64;

    // perc_mult[i] contains the percentage of combinations that
    // have an item of multiplicity i
    let mut perc_mult: Vec<f64> = Vec::with_capacity(n + 1);
    perc_mult.resize(n + 1, 0.0);

    for m in 1..n + 1 {
        let mut y: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for i in c..n - m + 1 {
            perc_mult[i] += (n as f64 - m as f64 - i as f64 + 1.0).log2();
            y += ((hist.get_hist_values()[i] as f64).log2() + perc_mult[i] - n_fall_m).exp2();
        }

        pangrowth[m - 1] = tot - y;
    }

    pangrowth
}

fn calc_growth_core(hist: &Hist, t_coverage: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from 0..n (both included)
    let c = usize::max(1, t_coverage.to_absolute(n + 1));
    let mut n_fall_m: f64 = 0.0;
    let mut pangrowth: Vec<f64> = vec![0.0; n];

    // In perc_mult[i] is contained the percentage of combinations
    // that have an item of multiplicity i
    let mut perc_mult: Vec<f64> = Vec::with_capacity(n + 1);
    perc_mult.resize(n + 1, 0.0);

    for m in 1..n + 1 {
        let mut y: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for i in usize::max(m, c)..n + 1 {
            perc_mult[i] += (i as f64 - m as f64 + 1.0).log2();
            y += ((hist.get_hist_values()[i] as f64).log2() + perc_mult[i] - n_fall_m).exp2();
        }
        pangrowth[m - 1] = y;
    }

    pangrowth
}

fn calc_growth_quorum(hist: &Hist, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
    let n = hist.get_maximum_coverage(); // hist array has length n+1: from [0..n]
    let c = usize::max(1, t_coverage.to_absolute(n));
    let quorum = t_quorum.to_relative(n);
    let mut pangrowth: Vec<f64> = vec![0.0; n];

    let mut n_fall_m: f64 = 0.0;
    let mut m_fact: f64 = 0.0;

    let mut perc_mult: Vec<f64> = vec![0.0; n + 1];
    let mut q: Vec<Vec<f64>> = vec![vec![0.0; n + 1]; n + 1];

    for m in 1..n + 1 {
        m_fact += (m as f64).log2();
        let m_quorum = (m as f64 * quorum).ceil() as usize;

        //100% quorum
        let mut yl: f64 = 0.0;
        n_fall_m += (n as f64 - m as f64 + 1.0).log2();
        for i in usize::max(m, c)..n + 1 {
            perc_mult[i] += (i as f64 - m as f64 + 1.0).log2();
            yl += ((hist.get_hist_values()[i] as f64).log2() + perc_mult[i] - n_fall_m).exp2();
        }

        //[m_quorum, 100) quorum
        let mut yr: f64 = 0.0;
        for i in m_quorum..n {
            let mut sum_q = 0.0;
            let mut add = false;
            for j in usize::max(m_quorum, c)..m {
                if n + j + 1 > i + m && j <= i {
                    if q[i][j] == 0.0 {
                        q[i][j] = choose(i, j);
                    }
                    q[i][j] += (n as f64 - i as f64 - m as f64 + 1.0 + j as f64).log2();
                    q[i][j] -= (m as f64 - j as f64).log2();
                    sum_q += (q[i][j] + m_fact - n_fall_m).exp2();
                    add = true;
                }
            }
            if add {
                yr += ((hist.get_hist_values()[i] as f64).log2() + sum_q.log2()).exp2();
            }
        }
        pangrowth[m - 1] = yl + yr;
    }
    pangrowth
}

fn choose(n: usize, k: usize) -> f64 {
    let mut res: f64 = 0.0;
    if k > n {
        return 0.0;
    }

    let k = if k > n - k { n - k } else { k };

    let n = n as f64;

    for i in 0..k {
        res += (n - i as f64).log2();
        res -= (i as f64 + 1.0).log2();
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hist_calc_growth_union() {
        let hist = Hist::from_values(vec![0, 5, 3, 2]);

        let t_coverage = Threshold::Absolute(0);
        let test_growth: Vec<f64> = vec![5.666666666666667, 8.333333333333334, 10.0];
        let growth = calc_growth_union(&hist, &t_coverage);
        assert_eq!(growth, test_growth, "Wrong growth union");
    }

    #[test]
    fn test_hist_calc_growth_core() {
        let hist = Hist::from_values(vec![0, 5, 3, 2]);

        let t_coverage = Threshold::Absolute(0);
        let test_core: Vec<f64> = vec![5.666666666666666, 3.0, 2.0];
        let core = calc_growth_core(&hist, &t_coverage);
        assert_eq!(core, test_core, "Wrong growth core");
    }

    #[test]
    fn test_hist_calc_growth_quorum() {
        let hist = Hist::from_values(vec![0, 5, 3, 2, 3, 5, 0, 4, 2, 1]);

        let t_coverage = Threshold::Absolute(0);
        let t_quorum = Threshold::Relative(0.9);
        let test_growth: Vec<f64> = vec![
            11.88888888888889,
            7.027777777777777,
            4.761904761904761,
            3.4444444444444438,
            2.5952380952380953,
            2.0,
            1.5555555555555545,
            1.2222222222222217,
            1.0,
        ];
        let growth = calc_growth_quorum(&hist, &t_coverage, &t_quorum);
        assert_eq!(growth, test_growth, "Wrong growth quorum");
    }
}
