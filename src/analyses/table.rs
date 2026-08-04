use crate::{
    analyses::{ordered_histgrowth::parse_list, MatrixBasedAnalysis},
    coverage_matrix::CoverageMatrix,
    io::write_metadata_comments,
};

use super::AnalysisSection;

pub struct Table {
    total: bool,
    order: Option<String>,
}

impl MatrixBasedAnalysis for Table {
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        let order = match &self.order {
            Some(filename) => parse_list(&filename[..])?,
            None => matrix.get_path_names().clone(),
        };
        let mut string = write_metadata_comments()?;
        string.push_str("panacus\ttable\n");
        string.push_str(matrix.get_feature_type());
        if self.total {
            string.push_str("\ttotal");
        } else {
            for path_name in &order {
                string.push('\t');
                string.push_str(path_name);
            }
        }
        string.push('\n');
        for i in 0..matrix.get_feature_count() {
            string.push_str(matrix.get_feature_name(i));
            if self.total {
                let count = matrix.get_count_of_feature(i);
                string.push_str(&format!("\t{}", count));
            } else {
                for value in matrix.get_counts_for_feature_in_order(i, &order) {
                    string.push('\t');
                    string.push_str(&format!("{}", value));
                }
            }
            string.push('\n');
        }
        Ok(string)
    }

    fn get_type(&self) -> String {
        "Table".to_string()
    }

    fn generate_report_section(
        &mut self,
        _matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        Ok(Vec::new())
    }
}

impl Table {
    pub fn new(total: bool, order: Option<String>) -> Self {
        Self { total, order }
    }
}
