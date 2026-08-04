use std::cmp;
use std::io::Write;
use std::process::Command;

use once_cell::sync::Lazy;
use tempfile::NamedTempFile;

use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::html_report::ReportItem;
use crate::util::CountType;

use super::AnalysisSection;

static VIRIDIS: Lazy<Vec<(f64, (f64, f64, f64))>> = Lazy::new(|| {
    vec![
        (0.0, (253.0, 231.0, 37.0)),
        (0.25, (94.0, 201.0, 98.0)),
        (0.5, (33.0, 145.0, 140.0)),
        (0.75, (59.0, 82.0, 139.0)),
        (1.0, (68.0, 1.0, 84.0)),
    ]
});

pub struct CoverageColors {
    inner: Option<InnerCoverageColors>,
}

impl MatrixBasedAnalysis for CoverageColors {
    fn get_type(&self) -> String {
        "CoverageColors".to_string()
    }
    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        self.set_inner(matrix)?;
        let mut result = String::new();
        result.push_str("Name,Color\n");
        for (id, cov) in self
            .inner
            .as_ref()
            .expect("Setting inner worked for coverage colors")
            .coverages
            .iter()
        {
            let color = convert_to_color(
                *cov,
                self.inner
                    .as_ref()
                    .expect("Setting inner worked for coverage colors")
                    .max_cov,
            );
            result.push_str(&format!("{},{}\n", id, color));
        }
        Ok(result)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(matrix)?;
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let id_prefix = format!(
            "coverage-colors-{}",
            matrix
                .get_run_id()
                .to_lowercase()
                .replace([' ', '|', '\\'], "-")
        );
        let output_svg = self.get_image(matrix);
        let growth_tabs = vec![AnalysisSection {
            analysis: "Coverage Colors".to_string(),
            run_name: matrix.get_run_name().to_owned(),
            run_id: matrix.get_run_id().to_owned(),
            countable: CountType::Node.to_string(),
            table: Some(table),
            id: id_prefix.to_string(),
            plot_downloads: vec![("svg".to_string(), "Download as svg".to_string())],
            items: match output_svg {
                Ok(svg_path) => vec![ReportItem::Svg {
                    id: format!("svg-{id_prefix}"),
                    file: svg_path,
                }],
                Err(e) => {
                    log::error!("Error running bandage: {e}");
                    Vec::new()
                }
            },
        }];
        Ok(growth_tabs)
    }
}

impl Default for CoverageColors {
    fn default() -> Self {
        Self::new()
    }
}

impl CoverageColors {
    pub fn new() -> Self {
        Self { inner: None }
    }

    fn set_inner(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        let coverages = matrix.get_feature_counts();
        let node_names = matrix.get_feature_names();

        let coverages: Vec<(String, u64)> = coverages
            .into_iter()
            .enumerate()
            .map(|(id, cov)| (node_names[id].clone(), cov as u64))
            .collect();
        let max_cov = matrix.get_path_names().len() as u64;

        let coverages = InnerCoverageColors { coverages, max_cov };
        self.inner = Some(coverages);
        Ok(())
    }

    fn get_image(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        if self.inner.is_none() {
            panic!("get_image can only be called if coverage colors was already calculated");
        }

        let mut csv = NamedTempFile::with_suffix(".csv")?;
        let table = self.generate_table(matrix)?;
        write!(csv, "{}", table)?;
        let csv_path = csv.path().to_str().unwrap();
        let svg = NamedTempFile::with_suffix(".svg")?;
        let svg_path = svg.path().to_str().unwrap();
        // TODO get file name
        let graph_path = matrix.get_run_name();
        Command::new("Bandage")
            .args(["image", graph_path, svg_path, "--colors", csv_path])
            .output()?;
        let copied_path = svg_path.to_string();
        svg.keep()?;
        Ok(copied_path)
    }
}

fn convert_to_color(cov: u64, max_cov: u64) -> String {
    // Handle case of coverage == 1
    if cov == 1 {
        return "#FF0000".to_string();
    }

    if cov == 0 {
        return "#AAAAAA".to_string();
    }

    // Trim value
    let value = cmp::max(cmp::min(cov, max_cov), 2);
    let relative_value = (value - 2) as f64 / (max_cov - 2) as f64;

    for i in 0..VIRIDIS.len() - 1 {
        if relative_value >= VIRIDIS[i].0 && relative_value <= VIRIDIS[i + 1].0 {
            let factor = (relative_value - VIRIDIS[i].0) / (VIRIDIS[i + 1].0 - VIRIDIS[i].0);
            let r = (VIRIDIS[i].1 .0 + factor * (VIRIDIS[i + 1].1 .0 - VIRIDIS[i].1 .0)) as u64;
            let g = (VIRIDIS[i].1 .1 + factor * (VIRIDIS[i + 1].1 .1 - VIRIDIS[i].1 .1)) as u64;
            let b = (VIRIDIS[i].1 .2 + factor * (VIRIDIS[i + 1].1 .2 - VIRIDIS[i].1 .2)) as u64;
            return format!("#{:02X}{:02X}{:02X}", r, g, b);
        }
    }
    "#000000".to_string()
}

struct InnerCoverageColors {
    coverages: Vec<(String, u64)>,
    max_cov: u64,
}
