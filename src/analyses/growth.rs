use core::{panic, str};
use std::cmp;
use std::collections::HashSet;
use std::fs::File;
use std::io::BufReader;

use ml_helpers::linear_regression::huber_regressor::{solve, HuberRegressor};
use rayon::iter::{IntoParallelRefIterator, ParallelBridge, ParallelIterator};

use crate::analysis_parameter::AnalysisParameter;
use crate::graph_broker::{GraphBroker, Hist, ThresholdContainer};
use crate::html_report::ReportItem;
use crate::{
    io::parse_hists,
    io::write_table,
    util::{get_default_plot_downloads, CountType},
};

use super::{Analysis, AnalysisSection, ConstructibleAnalysis, InputRequirement};

type Hists = Vec<Hist>;
type Growths = Vec<(CountType, Vec<Vec<f64>>)>;
type Comments = Vec<Vec<u8>>;

pub struct Growth {
    parameter: AnalysisParameter,
    inner: Option<InnerGrowth>,
}

impl Analysis for Growth {
    fn get_type(&self) -> String {
        "Growth".to_string()
    }

    fn generate_table(
        &mut self,
        dm: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<String> {
        log::info!("reporting hist table");

        self.set_inner(dm)?;
        let growths = &self.inner.as_ref().unwrap().growths;
        let hist_aux = &self.inner.as_ref().unwrap().hist_aux;
        let comments = &self.inner.as_ref().unwrap().comments;
        let heaps_curves = &self.inner.as_ref().unwrap().heaps_curves;
        let mut res = String::new();
        for c in comments {
            res.push_str(str::from_utf8(&c[..])?);
            res.push_str("\n");
        }
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

        let hists = match &self.inner.as_ref().unwrap().hists {
            Some(h) => h.iter().collect::<Vec<_>>(),
            None => dm
                .expect("Growth needs either hist file or graph")
                .get_hists()
                .values()
                .collect::<Vec<_>>(),
        };

        if let AnalysisParameter::Growth {
            add_hist,
            add_alpha,
            ..
        } = self.parameter
        {
            if add_alpha {
                if let Some(heap_curves) = heaps_curves {
                    for ((alpha, _heaps), (count, _growth)) in
                        heap_curves.iter().zip(growths.iter())
                    {
                        res.push_str(&format!("# alpha ({}): {}\n", count.to_string(), alpha));
                    }
                }
            }
            if add_hist {
                for h in hists {
                    output_columns.push(h.coverage.iter().map(|x| *x as f64).collect());
                    header_cols.push(vec![
                        "hist".to_string(),
                        h.count.to_string(),
                        String::new(),
                        String::new(),
                    ])
                }
            }
        } else {
            panic!("Growth needs growth parameter");
        }

        for (count, g) in growths {
            output_columns.extend(g.clone());
            let m = hist_aux.coverage.len();
            header_cols.extend(
                std::iter::repeat("growth")
                    .take(m)
                    .zip(std::iter::repeat(count).take(m))
                    .zip(hist_aux.coverage.iter())
                    .zip(&hist_aux.quorum)
                    .map(|(((p, t), c), q)| {
                        vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                    }),
            );
        }
        res.push_str(&write_table(&header_cols, &output_columns)?);
        Ok(res)
    }

    fn generate_report_section(
        &mut self,
        dm: Option<&crate::graph_broker::GraphBroker>,
    ) -> anyhow::Result<Vec<AnalysisSection>> {
        self.set_inner(dm)?;
        let hist_aux = &self.inner.as_ref().unwrap().hist_aux;
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
        let table = self.generate_table(dm)?;
        let table = format!("`{}`", &table);
        let growths = &self.inner.as_ref().unwrap().growths;
        let heaps_curves = &self.inner.as_ref().unwrap().heaps_curves;
        let id_prefix = format!(
            "pan-growth-{}",
            self.get_run_id(dm.expect("Growth should be called with a graph"))
                .to_lowercase()
                .replace(&[' ', '|', '\\'], "-")
        );
        let growth_tabs = growths
            .iter()
            .enumerate()
            .map(|(i, (k, v))| AnalysisSection {
                id: format!("{id_prefix}-{k}"),
                analysis: "Pangenome Growth".to_string(),
                run_name: self.get_run_name(dm.expect("Growth should be called with a graph")),
                run_id: self.get_run_id(dm.expect("Growth should be called with a graph")),
                countable: k.to_string(),
                table: Some(table.clone()),
                items: vec![ReportItem::MultiBar {
                    id: format!("{id_prefix}-{k}"),
                    names: growth_labels.clone(),
                    x_label: "taxa".to_string(),
                    y_label: format!("#{}s", k),
                    labels: (1..v[0].len()).map(|i| i.to_string()).collect(),
                    values: v
                        .iter()
                        .map(|row| {
                            row.iter()
                                .map(|el| if el.is_nan() { 0.0 } else { *el })
                                .collect()
                        })
                        .collect(),
                    curve: heaps_curves.as_ref().map(|c| c[i].1.clone()).flatten(),
                    alpha: heaps_curves.as_ref().map(|c| c[i].0.clone()),
                    log_toggle: false,
                }],
                plot_downloads: get_default_plot_downloads(),
            })
            .collect();
        Ok(growth_tabs)
    }

    // fn get_subcommand() -> Command {
    //     Command::new("growth")
    //         .about("Calculate growth curve from coverage histogram")
    //         .args(&[
    //             arg!(hist_file: <HIST_FILE> "Coverage histogram as tab-separated value (tsv) file"),
    //             arg!(-a --hist "Also include histogram in output"),
    //             Arg::new("coverage").help("Ignore all countables with a coverage lower than the specified threshold. The coverage of a countable corresponds to the number of path/walk that contain it. Repeated appearances of a countable in the same path/walk are counted as one. You can pass a comma-separated list of coverage thresholds, each one will produce a separated growth curve (e.g., --coverage 2,3). Use --quorum to set a threshold in conjunction with each coverage (e.g., --quorum 0.5,0.9)")
    //                 .short('l').long("coverage").default_value("1"),
    //             Arg::new("quorum").help("Unlike the --coverage parameter, which specifies a minimum constant number of paths for all growth point m (1 <= m <= num_paths), --quorum adjust the threshold based on m. At each m, a countable is counted in the average growth if the countable is contained in at least floor(m*quorum) paths. Example: A quorum of 0.9 requires a countable to be in 90% of paths for each subset size m. At m=10, it must appear in at least 9 paths. At m=100, it must appear in at least 90 paths. A quorum of 1 (100%) requires presence in all paths of the subset, corresponding to the core. Default: 0, a countable counts if it is present in any path at each growth point. Specify multiple quorum values with a comma-separated list (e.g., --quorum 0.5,0.9). Use --coverage to set static path thresholds in conjunction with variable quorum percentages (e.g., --coverage 5,10).")
    //                 .short('q').long("quorum").default_value("0"),
    //         ])
    // }

    fn get_graph_requirements(&self) -> HashSet<super::InputRequirement> {
        HashSet::from([InputRequirement::Hist])
    }
}

impl ConstructibleAnalysis for Growth {
    fn from_parameter(parameter: AnalysisParameter) -> Self {
        Growth {
            parameter,
            inner: None,
        }
    }
}

impl Growth {
    pub fn generate_table_from_hist(&mut self, file: &str) -> anyhow::Result<String> {
        if let AnalysisParameter::Growth {
            quorum,
            coverage,
            add_hist,
            add_alpha,
        } = &self.parameter
        {
            log::info!("reporting hist table");

            let quorum = quorum.to_owned().unwrap_or("0".to_string());
            let coverage = coverage.to_owned().unwrap_or("1".to_string());
            let hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;
            let file = File::open(file)?;
            let mut data = BufReader::new(file);
            let (coverages, comments) = parse_hists(&mut data)?;
            let hists: Hists = coverages
                .into_iter()
                .map(|(count, coverage)| Hist { count, coverage })
                .collect();
            let growths: Growths = hists
                .par_iter()
                .map(|h| (h.count, h.calc_all_growths(&hist_aux)))
                .collect();
            let mut res = String::new();
            for c in comments {
                res.push_str(str::from_utf8(&c[..])?);
                res.push_str("\n");
            }
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

            if *add_hist {
                for h in hists.iter() {
                    output_columns.push(h.coverage.iter().map(|x| *x as f64).collect());
                    header_cols.push(vec![
                        "hist".to_string(),
                        h.count.to_string(),
                        String::new(),
                        String::new(),
                    ])
                }
            }

            for ((count, g), hist) in growths.iter().zip(hists.iter()) {
                output_columns.extend(g.clone());
                let m = hist_aux.coverage.len();
                if *add_alpha {
                    let hist: Vec<f64> = hist.coverage.iter().map(|x| *x as f64).collect();
                    let (alpha, _) = get_regression(&hist);
                    res.push_str(&format!("# alpha ({}): {}\n", count.to_string(), alpha));
                }
                header_cols.extend(
                    std::iter::repeat("growth")
                        .take(m)
                        .zip(std::iter::repeat(count).take(m))
                        .zip(hist_aux.coverage.iter())
                        .zip(&hist_aux.quorum)
                        .map(|(((p, t), c), q)| {
                            vec![p.to_string(), t.to_string(), c.get_string(), q.get_string()]
                        }),
                );
            }
            res.push_str(&write_table(&header_cols, &output_columns)?);
            Ok(res)
        } else {
            panic!("Growth needs growth parameter");
        }
    }

    fn get_run_name(&self, gb: &GraphBroker) -> String {
        format!("{}", gb.get_run_name())
    }
    fn get_run_id(&self, gb: &GraphBroker) -> String {
        format!("{}-growth", gb.get_run_id())
    }

    fn set_inner(&mut self, gb: Option<&GraphBroker>) -> anyhow::Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }
        if let AnalysisParameter::Growth {
            coverage, quorum, ..
        } = &self.parameter
        {
            let quorum = quorum.to_owned().unwrap_or("0".to_string());
            let coverage = coverage.to_owned().unwrap_or("1".to_string());
            let hist_aux = ThresholdContainer::parse_params(&quorum, &coverage)?;

            if gb.is_none() {
                unimplemented!("Have not implemented growth without graph");
            } else {
                let gb = gb.unwrap();
                let growths: Growths = gb
                    .get_hists()
                    .values()
                    .par_bridge()
                    .map(|h| (h.count, h.calc_all_growths(&hist_aux)))
                    .collect();
                let hists = gb.get_hists();
                let heaps_curves = hist_aux.has_full_growth_at_idx().map(|index| {
                    log::info!("Calculating heaps law");
                    let heaps_curves: Vec<_> = growths
                        .iter()
                        .zip(hists.iter())
                        .map(|((_count_type, growth), (_count_type2, hist))| {
                            let growth = growth[index].clone();
                            let growth_len = growth.len();
                            let growth_last = *growth.last().unwrap();
                            let hist: Vec<f64> = hist.coverage.iter().map(|x| *x as f64).collect();
                            let x1 = 2.0f64;
                            let y1 = growth[1];
                            let x2 = growth.len() as f64 - 1.0;
                            let y2 = growth_last;
                            let (alpha, _offset) = get_regression(&hist);
                            if alpha >= 10.0 {
                                // TODO change 10 back to 1
                                (alpha, None)
                            } else {
                                let gamma = 1.0 - alpha;
                                let k = (y1 - y2) / (x1.powf(gamma) - x2.powf(gamma));
                                let c = y1 - k * x1.powf(gamma);
                                let curve_values = (1..=growth_len)
                                    .map(|x| (x as f64).powf(gamma) * k + c)
                                    .collect::<Vec<_>>();
                                (alpha, Some(curve_values))
                            }
                        })
                        .collect();
                    heaps_curves
                });
                self.inner = Some(InnerGrowth {
                    growths,
                    comments: Vec::new(),
                    hist_aux,
                    hists: None,
                    heaps_curves,
                });
            }
            Ok(())
        } else {
            panic!("Growth should always contain growth parameter")
        }
    }
}

fn get_regression(hist: &Vec<f64>) -> (f64, f64) {
    let x: Vec<f64> = (1..hist.len()).map(|x| (x as f64)).collect();
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

struct InnerGrowth {
    growths: Growths,
    comments: Comments,
    hist_aux: ThresholdContainer,
    hists: Option<Hists>,
    heaps_curves: Option<Vec<(f64, Option<Vec<f64>>)>>,
}
