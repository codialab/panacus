use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use crate::analyses::growth::calc_all_growths;
use crate::analyses::MatrixBasedAnalysis;
use crate::coverage_matrix::CoverageMatrix;
use crate::file_formats::gfa_parser::{choose, ThresholdContainer};
use crate::html_report::ReportItem;
use crate::util::{get_default_plot_downloads, Threshold};
use rayon::prelude::*;

use super::AnalysisSection;

pub struct SectionGrowth {
    sections_file: String,
    // Map from section name to group index (depends on grouping!)
    sections: HashMap<String, Vec<usize>>,
    section_order: Vec<String>,
    coverage: String,
    quorum: String,
    inner: Option<InnerGrowth>,
}

impl MatrixBasedAnalysis for SectionGrowth {
    fn get_type(&self) -> String {
        "SectionGrowth".to_string()
    }

    fn generate_table(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<String> {
        log::info!("reporting hist table");

        self.set_inner(matrix)?;
        let inner = self.inner.as_ref().unwrap();
        let mut res = String::new();

        res.push_str(&format!(
            "# {}\n",
            std::env::args().collect::<Vec<String>>().join(" ")
        ));

        res.push_str(&format!("panacus\tsection-growth\n"));
        res.push_str(&format!("count\t{}\n", matrix.get_feature_type()));

        res.push_str("coverage");
        for c in inner.thresholds.coverage.iter() {
            res.push_str(&format!("\t{}", c.get_string()));
        }
        res.push_str("\nquorum");
        for q in inner.thresholds.quorum.iter() {
            res.push_str(&format!("\t{}", q.get_string()));
        }
        res.push_str("\n");

        let mut section_index = 0;
        for idx in 0..inner.growths[0].len() {
            if section_index < inner.sections.len() - 1
                && idx >= inner.sections[section_index + 1].1
            {
                section_index += 1;
            }
            res.push_str(&format!("{}\t{}", idx, inner.sections[section_index].0,));
            for j in 0..inner.thresholds.coverage.len() {
                res.push_str(&format!("\t{}", inner.growths[j][idx]));
            }
            res.push_str("\n");
        }
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        matrix: &CoverageMatrix,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(matrix)?;
        let id = format!("section-growth-{}", matrix.get_run_id());
        let table = self.generate_table(matrix)?;
        let table = format!("`{}`", &table);
        let hist_aux = &self.inner.as_ref().unwrap().thresholds;
        let growths = &self.inner.as_ref().unwrap().growths;
        let sections = &self.inner.as_ref().unwrap().sections;
        let section_labels = sections.iter().map(|(x, _)| x.to_string()).collect();
        let section_separators = sections
            .iter()
            .skip(1)
            .map(|(_, y)| *y as f64 - 0.5) // Subtract 0.5 to place the separators in between the groups
            .collect();
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
        let analysis_section = AnalysisSection {
            analysis: self.get_type(),
            run_name: matrix.get_run_name().to_string(),
            run_id: matrix.get_run_id().to_string(),
            countable: matrix.get_feature_type().to_string(),
            id: id.clone(),
            table: Some(table),
            plot_downloads: get_default_plot_downloads(),
            items: vec![ReportItem::SectionLine {
                id: format!("{id}-{}", matrix.get_feature_type()),
                names: growth_labels.clone(),
                x_label: "taxa".to_string(),
                y_label: format!("#{}s", matrix.get_feature_type()),
                labels: (0..(growths[0].len())).map(|i| i.to_string()).collect(),
                values: growths
                    .iter()
                    .map(|row| {
                        row.iter()
                            .map(|el| if el.is_nan() { 0.0 } else { *el })
                            .collect()
                    })
                    .collect(),
                section_separators,
                section_labels,
                log_toggle: false,
            }],
        };
        Ok(vec![analysis_section])
    }
}

impl SectionGrowth {
    pub fn new(sections_file: String, coverage: Option<String>, quorum: Option<String>) -> Self {
        Self {
            sections_file,
            sections: HashMap::new(),
            section_order: Vec::new(),
            coverage: coverage.unwrap_or("1".to_string()),
            quorum: quorum.unwrap_or("0.0".to_string()),
            inner: None,
        }
    }

    fn get_growths_for_c_q(
        &self,
        c: Threshold,
        thresholds: ThresholdContainer,
        matrix: &CoverageMatrix,
    ) -> Vec<Vec<f64>> {
        let first_group = &self.sections[&self.section_order[0]];
        let (hist, _) = matrix.get_hist_and_used_features_with_coverage(None, first_group, &c);
        let mut full_growths = calc_all_growths(&hist, &thresholds, false);
        let mut path_collection: Vec<usize> = first_group.clone();

        for section in self.section_order.iter().skip(1) {
            log::info!("Handling quorum section {}", section);
            let paths = &self.sections[section];
            let (current_abacus, _) =
                matrix.get_abacus_and_used_features_with_coverage(None, paths, &c);
            // log::info!("Current abacus: {:?}", current_abacus);
            let (previous_abacus, _) =
                matrix.get_abacus_and_used_features_with_coverage(None, &path_collection, &c);
            // log::info!("Previous abacus: {:?}", previous_abacus);
            let hist = Hist3D::from_abaci(
                &previous_abacus,
                path_collection.len(),
                &current_abacus,
                paths.len(),
                matrix,
            );
            let growths = hist.calc_all_growths(&thresholds, false);
            full_growths = full_growths
                .into_iter()
                .zip(growths.into_iter())
                .map(|(mut full_growth, growth)| {
                    full_growth.extend(growth.into_iter());
                    full_growth
                })
                .collect();
            path_collection.extend(paths.into_iter());
        }
        full_growths
    }

    fn get_growths_for_coverage(
        &self,
        c: Threshold,
        q: Vec<Threshold>,
        matrix: &CoverageMatrix,
    ) -> Vec<Vec<f64>> {
        let thresholds = ThresholdContainer {
            coverage: vec![Threshold::Absolute(1); q.len()],
            quorum: q,
        };
        // let (q0_thresholds, other_thresholds) = Self::split_thresholds(&thresholds);
        let other_growths = self.get_growths_for_c_q(c, thresholds, matrix);
        let full_growths = other_growths;

        full_growths
    }

    fn get_sections(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<Vec<(String, usize)>> {
        let sections_file = File::open(&self.sections_file)?;
        let reader = BufReader::new(sections_file);
        for line in reader.lines() {
            let line = line?;
            let fields: Vec<_> = line.trim().split_whitespace().collect();
            assert_eq!(fields.len(), 2, "Sections file can only have two columns");
            let section = fields[1].to_owned();
            if !self.section_order.contains(&section) {
                self.section_order.push(section.clone());
                self.sections.insert(section.clone(), Vec::new());
            }
            let groups = matrix.get_path_names();
            self.sections.get_mut(&section).unwrap().push(
                groups
                    .iter()
                    .position(|g| g == fields[0])
                    .expect(&format!("Group {} is not a group", fields[0])),
            );
        }

        let mut section_index = vec![];
        let mut rolling_index = 0;
        for i in 0..self.section_order.len() {
            section_index.push((self.section_order[i].clone(), rolling_index));
            rolling_index += self.sections[&self.section_order[i]].len();
        }
        Ok(section_index)
    }

    fn set_inner(&mut self, matrix: &CoverageMatrix) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        let mut section_index = self.get_sections(matrix)?;
        let thresholds = ThresholdContainer::parse_params(&self.quorum, &self.coverage)?;

        let coverage_quorums: HashMap<usize, Vec<_>> = thresholds
            .coverage
            .iter()
            .zip(thresholds.quorum.iter())
            .fold(HashMap::new(), |mut acc, (c, q)| match c {
                Threshold::Absolute(c) => {
                    acc.entry(*c).or_default().push(*q);
                    acc
                }
                _ => panic!("Coverage can not be relative for section growth"),
            });
        let mut full_growths = Vec::new();
        let mut thresholds_sorted = Vec::new();
        for (c, qs) in coverage_quorums.into_iter() {
            let c = Threshold::Absolute(c);
            let new_thresholds: Vec<(Threshold, Threshold)> = qs.iter().map(|q| (c, *q)).collect();
            thresholds_sorted.extend(new_thresholds.into_iter());
            full_growths.extend(self.get_growths_for_coverage(c, qs, matrix).into_iter());
        }
        let (coverage_sorted, quorum_sorted): (Vec<_>, Vec<_>) =
            thresholds_sorted.into_iter().unzip();
        let thresholds_sorted = ThresholdContainer {
            coverage: coverage_sorted,
            quorum: quorum_sorted,
        };

        for growth in full_growths.iter_mut() {
            growth.insert(0, f64::NAN);
        }
        for val in section_index.iter_mut().skip(1) {
            let x = &mut val.1;
            *x += 1;
        }

        self.inner = Some(InnerGrowth {
            growths: full_growths,
            sections: section_index,
            thresholds: thresholds_sorted,
        });
        Ok(())
    }

    // fn split_thresholds(
    //     thresholds: &ThresholdContainer,
    // ) -> (ThresholdContainer, ThresholdContainer) {
    //     let mut first_c = Vec::new();
    //     let mut first_q = Vec::new();
    //     let mut second_c = Vec::new();
    //     let mut second_q = Vec::new();
    //     for (c, q) in thresholds.coverage.iter().zip(thresholds.quorum.iter()) {
    //         if q.to_relative(1) == 0.0 {
    //             first_c.push(*c);
    //             first_q.push(*q);
    //         } else {
    //             second_c.push(*c);
    //             second_q.push(*q);
    //         }
    //     }
    //     (
    //         ThresholdContainer {
    //             coverage: first_c,
    //             quorum: first_q,
    //         },
    //         ThresholdContainer {
    //             coverage: second_c,
    //             quorum: second_q,
    //         },
    //     )
    // }
}

struct InnerGrowth {
    growths: Vec<Vec<f64>>,
    sections: Vec<(String, usize)>,
    thresholds: ThresholdContainer,
}

#[derive(Debug, Clone)]
pub struct Hist3D {
    pub coverage: Vec<Vec<usize>>,
}

impl Hist3D {
    pub fn from_abaci(
        a: &Vec<usize>,
        n_a: usize,
        b: &Vec<usize>,
        n_b: usize,
        matrix: &CoverageMatrix,
    ) -> Self {
        assert_eq!(
            a.len(),
            b.len(),
            "Abaci of Hist3D need to have the same length!"
        );
        let coverage = Self::construct_hist_bps(a, n_a, b, n_b, matrix);
        // log::info!("Coverage vec: {:?}", coverage);
        Self { coverage }
    }

    fn construct_hist_bps(
        a: &Vec<usize>,
        n_a: usize,
        b: &Vec<usize>,
        n_b: usize,
        matrix: &CoverageMatrix,
    ) -> Vec<Vec<usize>> {
        let a_max_coverage = n_a;
        let b_max_coverage = n_b;
        let mut hist: Vec<Vec<usize>> = vec![vec![0; b_max_coverage + 1]; a_max_coverage + 1];
        for (i, a_cov) in a.iter().enumerate() {
            let b_cov = b[i];
            if i != 0 {
                if *a_cov as usize >= hist.len() {
                    log::warn!("coverage {} of item {} in Abacus a exceeds the number of groups {}, it'll be ignored in the count", a_cov, i, a_max_coverage);
                    continue;
                }
                if b_cov as usize >= hist[0].len() {
                    log::warn!("coverage {} of item {} in Abacus b exceeds the number of groups {}, it'll be ignored in the count", b_cov, i, b_max_coverage);
                    continue;
                }
            }
            hist[*a_cov as usize][b_cov as usize] += matrix.get_feature_lengths()[i] as usize;
        }

        hist
    }

    pub fn calc_all_growths(
        &self,
        hist_aux: &ThresholdContainer,
        insert_zero: bool,
    ) -> Vec<Vec<f64>> {
        // self.to_tsv();
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
                self.calc_growth(c, q)
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

    pub fn calc_growth(&self, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
        self.calc_growth_quorum(t_coverage, t_quorum)
    }

    fn calc_growth_quorum(&self, t_coverage: &Threshold, t_quorum: &Threshold) -> Vec<f64> {
        let n1 = self.coverage.len() - 1;
        let n2 = self.coverage[0].len() - 1;
        let c = usize::max(1, t_coverage.to_absolute(n2));
        let quorum = t_quorum.to_relative(n2);
        let mut pangrowth: Vec<f64> = vec![0.0; n2];

        log::info!("Calculating quorum: {}, c: {}", quorum, c);

        for m in 1..n2 + 1 {
            let mut res = 0.0;
            let m_quorum = usize::max(((m as f64 + n1 as f64) * quorum).ceil() as usize, 1);
            let n2_choose_m = choose(n2, m);

            for i in 0..=n2 {
                for j in 0..=n1 {
                    if self.coverage[j][i] == 0 {
                        continue;
                    }
                    if m_quorum <= j {
                        res += self.coverage[j][i] as f64;
                        continue;
                    }
                    let mut sub_res = 0.0f64;
                    for k in usize::max(m_quorum - j, c)..=usize::min(i, m) {
                        let first = choose(i, k);
                        let second = choose(n2 - i, m - k);
                        if i >= k && n2 - i >= m - k {
                            sub_res += (first + second - n2_choose_m).exp2();
                        }
                    }
                    res += (self.coverage[j][i] as f64) * sub_res;
                }
            }
            pangrowth[m - 1] = res;
        }

        pangrowth
    }
}
