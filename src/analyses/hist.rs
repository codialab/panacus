use crate::html_report::ReportItem;
use crate::{io::write_table, util::get_default_plot_downloads};

use super::{AnalysisSection, HistBasedAnalysis};

pub struct Hist {}

impl HistBasedAnalysis for Hist {
    fn get_type(&self) -> String {
        "Hist".to_string()
    }

    fn generate_table(&mut self, hist: &crate::hist::Hist) -> anyhow::Result<String> {
        log::info!("reporting hist table");
        let mut res = String::new();
        res.push_str(&crate::io::write_metadata_comments()?);

        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            String::new(),
            String::new(),
        ]];
        let mut output_columns = Vec::new();
        output_columns.push(hist.get_hist_values().iter().map(|x| *x as f64).collect());
        header_cols.push(vec![
            "hist".to_string(),
            hist.get_feature_type().to_string(),
            String::new(),
            String::new(),
        ]);
        res.push_str(&write_table(&header_cols, &output_columns)?);
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        hist: &crate::hist::Hist,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        let table = self.generate_table(hist)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "cov-hist-{}",
            hist.get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        let histogram_tabs = AnalysisSection {
            id: format!("{id_prefix}-{}", hist.get_feature_type()),
            analysis: "Coverage Histogram".to_string(),
            table: Some(table.clone()),
            run_name: hist.get_run_name().to_string(),
            run_id: hist.get_run_id().to_string(),
            countable: hist.get_feature_type().to_string(),
            items: vec![ReportItem::Bar {
                id: format!("{id_prefix}-{}", hist.get_feature_type()),
                name: hist.get_run_name().to_string(),
                x_label: "taxa".to_string(),
                y_label: format!("#{}s", hist.get_feature_type()),
                labels: (0..hist.get_hist_values().len())
                    .map(|s| s.to_string())
                    .collect(),
                values: hist.get_hist_values().iter().map(|c| *c as f64).collect(),
                log_toggle: true,
            }],
            plot_downloads: get_default_plot_downloads(),
        };
        Ok(vec![histogram_tabs])
    }
}

impl Default for Hist {
    fn default() -> Self {
        Self::new()
    }
}

impl Hist {
    pub fn new() -> Self {
        Self {}
    }
}
