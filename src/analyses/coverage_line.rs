use crate::{
    analyses::MatrixBasedAnalysis,
    coverage_matrix::CoverageMatrix,
    hist::Hist,
    html_report::{AnalysisSection, ReportItem},
    io::write_table_with_start_index,
    util::get_default_plot_downloads,
};

pub struct CoverageLine {
    reference: Option<String>,
    hist: Option<Hist>,
}

impl MatrixBasedAnalysis for CoverageLine {
    fn get_type(&self) -> String {
        "CoverageLine".to_string()
    }

    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        log::info!("reporting coverage line table");
        let mut res = String::new();
        res.push_str(&crate::io::write_metadata_comments()?);
        if self.hist.is_none() {
            self.set_inner(matrix);
        }

        let mut header_cols = vec![vec![
            "panacus".to_string(),
            "count".to_string(),
            String::new(),
            String::new(),
        ]];
        let mut output_columns = Vec::new();
        output_columns.push(
            self.hist
                .as_ref()
                .unwrap()
                .get_hist_values()
                .iter()
                .map(|x| *x as f64)
                .skip(1)
                .collect(),
        );
        header_cols.push(vec![
            "hist".to_string(),
            matrix.get_feature_type().to_string(),
            String::new(),
            String::new(),
        ]);
        res.push_str(&write_table_with_start_index(
            &header_cols,
            &output_columns,
            1,
        )?);
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<crate::html_report::AnalysisSection>> {
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "coverage-line-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        if self.hist.is_none() {
            self.set_inner(matrix);
        }
        let mut values = self.hist.as_ref().unwrap().get_hist_values().to_vec();
        while let Some(last) = values.pop() {
            if last != 0 {
                values.push(0);
                break;
            }
        }
        let values: Vec<f32> = values.into_iter().skip(1).map(|c| c as f32).collect();
        let k = matrix.get_feature_type();
        let coverage_line_tabs = AnalysisSection {
            id: format!("{id_prefix}-{k}"),
            analysis: "Coverage Line".to_string(),
            table: Some(table.clone()),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: k.to_string(),
            items: vec![ReportItem::Line {
                id: format!("{id_prefix}-{k}"),
                name: matrix.get_run_name().to_owned(),
                x_label: "Allele count".to_string(),
                y_label: format!("#{}s", k),
                x_values: (1..=values.len()).map(|s| s as f32).collect(),
                y_values: values,
                log_x: true,
                log_y: true,
            }],
            plot_downloads: get_default_plot_downloads(),
        };
        Ok(vec![coverage_line_tabs])
    }
}

impl CoverageLine {
    pub fn new(reference: Option<String>) -> Self {
        Self {
            reference,
            hist: None,
        }
    }

    fn set_inner(&mut self, matrix: &CoverageMatrix) {
        if let Some(r) = self.reference.as_ref() {
            self.hist = Some(matrix.get_hist_for_reference(r));
        } else {
            self.hist = Some(matrix.get_hist());
        }
    }
}
