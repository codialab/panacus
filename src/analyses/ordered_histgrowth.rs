use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;

use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::file_formats::gfa_parser::ThresholdContainer;
use crate::html_report::ReportItem;
use crate::io::{write_metadata_comments, write_ordered_table};
use crate::util::get_default_plot_downloads;
use crate::util::Threshold;
use rayon::prelude::*;

use super::AnalysisSection;

type Growths = Vec<Vec<f64>>;

pub struct OrderedHistgrowth {
    coverage: Option<String>,
    quorum: Option<String>,
    path_order: Option<Vec<String>>,
    order: Option<String>,
    inner: Option<InnerOrderedGrowth>,
}

impl MatrixBasedAnalysis for OrderedHistgrowth {
    fn get_type(&self) -> String {
        "OrderedHistgrowth".to_string()
    }

    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        self.set_inner(matrix)?;

        let mut res = write_metadata_comments()?;

        let inner = self.inner.as_ref().unwrap();
        // insert empty row for 0 element
        let m = inner.hist_aux.coverage.len();
        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            "coverage".to_string(),
            "quorum".to_string(),
        ]];
        header_cols.extend(
            std::iter::repeat("ordered-growth")
                .take(m)
                .zip(std::iter::repeat(matrix.get_feature_type()).take(m))
                .zip(inner.hist_aux.coverage.iter())
                .zip(&inner.hist_aux.quorum)
                .map(|(((p, t), c), q)| {
                    vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                })
                .collect::<Vec<Vec<String>>>(),
        );
        let mut order = self.path_order.as_ref().unwrap().clone();

        // Insert zero element
        order.insert(0, "No_Path".to_string());
        let mut growths = inner.growths.clone();
        for growth in growths.iter_mut() {
            growth.insert(0, 0.0);
        }
        let table = write_ordered_table(&header_cols, &growths, &order)?;
        res.push_str(&table);
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(matrix)?;
        let hist_aux = &self.inner.as_ref().unwrap().hist_aux;
        let growth_labels = (0..hist_aux.coverage.len())
            .map(|i| {
                format!(
                    "coverage ≥ {}, quorum ≥ {}%",
                    hist_aux.coverage[i].get_string(),
                    hist_aux.quorum[i].get_string()
                )
            })
            .collect::<Vec<_>>();
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let growths = &self.inner.as_ref().unwrap().growths;
        let id_prefix = format!(
            "pan-ordered-growth-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        let labels = matrix.get_path_names().clone();
        let count = matrix.get_feature_type();
        let growth_tabs = vec![AnalysisSection {
            id: id_prefix.to_string(),
            analysis: "Ordered Growth".to_string(),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: matrix.get_feature_type().to_string(),
            table: Some(table.clone()),
            items: vec![ReportItem::MultiBar {
                id: id_prefix.to_string(),
                names: growth_labels.clone(),
                x_label: "taxa".to_string(),
                y_label: format!("{}s", count),
                //labels: (1..growths[0].len()).map(|i| i.to_string()).collect(),
                labels,
                values: growths.clone(),
                log_toggle: false,
                alpha: None,
                curve: None,
            }],
            plot_downloads: get_default_plot_downloads(),
        }];
        Ok(growth_tabs)
    }
}

impl OrderedHistgrowth {
    pub fn new(coverage: Option<String>, quorum: Option<String>, order: Option<String>) -> Self {
        Self {
            coverage,
            quorum,
            order,
            path_order: None,
            inner: None,
        }
    }
    fn set_inner(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        let quorum = self.quorum.clone().unwrap_or("0".to_string());
        let coverage = self.coverage.clone().unwrap_or("1".to_string());
        let hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;
        let order = match &self.order {
            Some(filename) => parse_list(&filename[..])?,
            None => matrix.get_path_names().clone(),
        };

        let growths: Vec<Vec<f64>> = hist_aux
            .coverage
            .par_iter()
            .zip(&hist_aux.quorum)
            .map(|(c, q)| {
                log::info!(
                    "calculating ordered growth for coverage >= {} and quorum >= {}",
                    &c,
                    &q
                );
                calc_growth(matrix, c, q, &order)
            })
            .collect();
        self.path_order = Some(order);
        self.inner = Some(InnerOrderedGrowth { growths, hist_aux });
        Ok(())
    }
}

struct InnerOrderedGrowth {
    growths: Growths,
    hist_aux: ThresholdContainer,
}

pub fn parse_list(filename: &str) -> anyhow::Result<Vec<String>> {
    let lines = read_lines(filename)?.map_while(Result::ok).collect();
    Ok(lines)
}

pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where
    P: AsRef<Path>,
{
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

pub fn calc_growth(
    matrix: &CoverageMatrix,
    t_coverage: &Threshold,
    t_quorum: &Threshold,
    order: &[String],
) -> Vec<f64> {
    let n_paths = matrix.get_path_names().len();
    let mut res = vec![0.0; n_paths];

    let c = usize::max(1, t_coverage.to_absolute(n_paths));
    let q = f64::max(0.0, t_quorum.to_relative(n_paths));
    let feature_lengths = matrix.get_feature_lengths();

    for (i, feature_length) in feature_lengths
        .iter()
        .enumerate()
        .take(matrix.get_feature_count())
    {
        // If this feature appears less often than the coverage
        // skip it.
        if matrix.get_count_of_feature(i) < c {
            continue;
        }

        let entry_indices: Vec<usize> = matrix.get_appearances_for_order(i, order);
        let mut entry_indices_index = 0;
        for (path_idx, path) in res
            .iter_mut()
            .enumerate()
            .take(n_paths)
            .skip(entry_indices[0])
        {
            if entry_indices_index < entry_indices.len() - 1
                && path_idx >= entry_indices[entry_indices_index + 1]
            {
                entry_indices_index += 1;
            }
            let num_seen_paths = path_idx + 1;
            let num_appearances = entry_indices_index + 1;
            let ratio = num_appearances as f64 / num_seen_paths as f64;
            if ratio >= q {
                *path += *feature_length as f64;
            }
        }
    }
    res
}
