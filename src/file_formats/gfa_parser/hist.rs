/* standard use */
use std::io::Write;
use std::io::{Error, ErrorKind};
use std::str::FromStr;

/* external crate */
use rayon::prelude::*;

/* private use */
use crate::util::{CountType, Threshold};

use super::abacus::AbacusByTotal;
use super::graph::GraphStorage;

#[derive(Debug, Clone)]
pub struct Hist3D {
    pub count: CountType,
    pub coverage: Vec<Vec<usize>>,
}

impl Hist3D {
    pub fn from_abaci(
        a: &AbacusByTotal,
        b: &AbacusByTotal,
        graph_aux: Option<&GraphStorage>,
    ) -> Self {
        assert_eq!(
            a.count, b.count,
            "Abaci of Hist3D need to have the same count type"
        );
        assert_eq!(
            a.countable.len(),
            b.countable.len(),
            "Abaci of Hist3D need to have the same length"
        );
        Self {
            count: a.count,
            coverage: match a.count {
                CountType::Node | CountType::Edge => Self::construct_hist(a, b),
                CountType::Bp => Self::construct_hist_bps(
                    a,
                    b,
                    graph_aux.expect("Graph auxiliary is needed for Bps hist"),
                ),
            },
        }
    }

    fn construct_hist(a: &AbacusByTotal, b: &AbacusByTotal) -> Vec<Vec<usize>> {
        let mut hist: Vec<Vec<usize>> = vec![vec![0; b.groups.len() + 1]; a.groups.len() + 1];
        for (i, a_cov) in a.countable.iter().enumerate() {
            let b_cov = b.countable[i];
            if i != 0 {
                if *a_cov as usize >= hist.len() {
                    log::warn!("coverage {} of item {} in Abacus a exceeds the number of groups {}, it'll be ignored in the count", a_cov, i, a.groups.len());
                    continue;
                }
                if b_cov as usize >= hist[0].len() {
                    log::warn!("coverage {} of item {} in Abacus b exceeds the number of groups {}, it'll be ignored in the count", b_cov, i, b.groups.len());
                    continue;
                }
            }
            hist[*a_cov as usize][b_cov as usize] += 1;
        }
        hist
    }

    fn construct_hist_bps(
        a: &AbacusByTotal,
        b: &AbacusByTotal,
        graph_aux: &GraphStorage,
    ) -> Vec<Vec<usize>> {
        let mut hist: Vec<Vec<usize>> = vec![vec![0; b.groups.len() + 1]; a.groups.len() + 1];
        for (i, a_cov) in a.countable.iter().enumerate() {
            let b_cov = b.countable[i];
            if i != 0 {
                if *a_cov as usize >= hist.len() {
                    log::warn!("coverage {} of item {} in Abacus a exceeds the number of groups {}, it'll be ignored in the count", a_cov, i, a.groups.len());
                    continue;
                }
                if b_cov as usize >= hist[0].len() {
                    log::warn!("coverage {} of item {} in Abacus b exceeds the number of groups {}, it'll be ignored in the count", b_cov, i, b.groups.len());
                    continue;
                }
            }
            hist[*a_cov as usize][b_cov as usize] += graph_aux.node_lens[i] as usize;
        }

        // subtract uncovered bps
        assert_eq!(
            a.uncovered_bps, b.uncovered_bps,
            "Uncovered bps of abaci are not identical"
        );
        let uncovered_bps = a.uncovered_bps.as_ref().unwrap();
        for (id, uncov) in uncovered_bps.iter() {
            hist[a.countable[*id as usize] as usize][b.countable[*id as usize] as usize] -= uncov;
            hist[0][0] += uncov;
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
        // eprintln!("hist3d: {:?}", self.coverage);

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

    pub fn to_tsv(&self) {
        eprintln!("");
        eprintln!("");
        for i in 0..self.coverage.len() {
            for j in 0..self.coverage[i].len() - 1 {
                eprint!("{}\t", self.coverage[i][j]);
            }
            eprintln!("{}", self.coverage[i][self.coverage[i].len() - 1]);
        }
        eprintln!("");
        eprintln!("");
    }
}

#[derive(Debug, Clone)]
pub struct Hist {
    pub count: CountType,
    pub coverage: Vec<usize>,
}

pub fn choose(n: usize, k: usize) -> f64 {
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

impl Hist {
    #[allow(dead_code)]
    pub fn to_tsv<W: std::io::Write>(&self, out: &mut std::io::BufWriter<W>) -> Result<(), Error> {
        writeln!(out, "hist\t{}", self.count)?;
        for (i, c) in self.coverage.iter().enumerate() {
            writeln!(out, "{}\t{}", i, c)?;
        }

        Ok(())
    }
}

pub enum RequireThreshold {
    Absolute,
    Relative,
    #[allow(dead_code)]
    Either,
}

pub fn parse_threshold_cli(
    threshold_str: &str,
    require: RequireThreshold,
) -> Result<Vec<Threshold>, Error> {
    let mut thresholds = Vec::new();

    for (i, el) in threshold_str.split(',').enumerate() {
        let rel_val = match f64::from_str(el.trim()) {
            Ok(t) => {
                if (0.0..=1.0).contains(&t) {
                    Ok(t)
                } else {
                    Err(Error::new(
                        ErrorKind::InvalidData,
                        format!(
                            "relative threshold \"{}\" ({}. element in list) must be within [0,1].",
                            &threshold_str,
                            i + 1
                        ),
                    ))
                }
            }
            Err(_) => Err(Error::new(
                ErrorKind::InvalidData,
                format!(
                    "threshold \"{}\" ({}. element in list) is required to be float, but isn't.",
                    &threshold_str,
                    i + 1
                ),
            )),
        };

        thresholds.push(
            match require {
                RequireThreshold::Absolute => Threshold::Absolute(usize::from_str(el.trim()).map_err(|_|
                    Error::new(
                            ErrorKind::InvalidData,
                            format!("threshold \"{}\" ({}. element in list) is required to be integer, but isn't.",
                    &threshold_str,
                    i + 1)))?),
            RequireThreshold::Relative => Threshold::Relative(rel_val?),
            RequireThreshold::Either =>
        if let Ok(t) = usize::from_str(el.trim()) {
            Threshold::Absolute(t)
        } else {
            Threshold::Relative(rel_val?)
            }
            }
            );
    }
    Ok(thresholds)
}

pub struct ThresholdContainer {
    pub quorum: Vec<Threshold>,
    pub coverage: Vec<Threshold>,
}

impl ThresholdContainer {
    pub fn parse_params(quorum: &str, coverage: &str) -> Result<Self, Error> {
        let mut quorum_thresholds = Vec::new();
        if !quorum.is_empty() {
            quorum_thresholds = parse_threshold_cli(quorum, RequireThreshold::Relative)?;
            log::debug!(
                "loaded {} quorum thresholds: {}",
                quorum_thresholds.len(),
                quorum_thresholds
                    .iter()
                    .map(|t| format!("{}", t))
                    .collect::<Vec<String>>()
                    .join(", ")
            );
        }
        if quorum_thresholds.is_empty() {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "quorum threshold setting requires at least one element, but none is given",
            ));
        }

        let mut coverage_thresholds = Vec::new();
        if !coverage.is_empty() {
            coverage_thresholds = parse_threshold_cli(coverage, RequireThreshold::Absolute)?;
            log::debug!(
                "loaded {} coverage thresholds: {}",
                coverage_thresholds.len(),
                coverage_thresholds
                    .iter()
                    .map(|t| format!("{}", t))
                    .collect::<Vec<String>>()
                    .join(", ")
            );
        }
        if coverage_thresholds.is_empty() {
            return Err(Error::new(
                ErrorKind::InvalidData,
                "coverage threshold setting requires at least one element, but none is given",
            ));
        }

        if quorum_thresholds.len() != coverage_thresholds.len() {
            if quorum_thresholds.len() == 1 {
                quorum_thresholds = vec![quorum_thresholds[0]; coverage_thresholds.len()];
            } else if coverage_thresholds.len() == 1 {
                coverage_thresholds = vec![coverage_thresholds[0]; quorum_thresholds.len()];
            } else {
                return Err(Error::new(ErrorKind::InvalidData,
                        "number of coverage and quorum threshold must match, or either one must have a single value"));
            }
        }

        Ok(Self {
            quorum: quorum_thresholds,
            coverage: coverage_thresholds,
        })
    }

    pub fn has_full_growth_at_idx(&self) -> Option<usize> {
        let mut index_by_coverage: Vec<(f64, usize)> = self
            .coverage
            .iter()
            .zip(self.quorum.iter())
            .enumerate()
            .filter_map(|(i, (c, q))| {
                if q.to_relative(1) == 0.0 {
                    Some((c.to_relative(1), i))
                } else {
                    None
                }
            })
            .collect();
        index_by_coverage.sort_by(|(c1, _), (c2, _)| c1.partial_cmp(c2).unwrap());
        index_by_coverage.get(0).map(|(_, i)| *i)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_almost_eq(a: f64, b: f64) {
        let epsilon = 1e-10;
        if (a - b).abs() > epsilon {
            panic!("Values are not almost equal: {} vs {}", a, b);
        }
    }

    fn factorial(n: usize) -> f64 {
        (1..=n).fold(1.0, |acc, x| acc * x as f64)
    }

    #[test]
    fn test_choose_function() {
        assert_almost_eq(choose(5, 0), 0.0);
        assert_almost_eq(choose(5, 5), 0.0);
        assert_almost_eq(choose(5, 1), (5.0f64).log2());
        assert_almost_eq(choose(5, 4), (5.0f64).log2());
        let expected = (factorial(5) / (factorial(2) * factorial(3))).log2();
        assert_almost_eq(choose(5, 2), expected);
        assert_eq!(choose(5, 6), 0.0);
    }
}
