use std::fmt;
use std::fmt::Debug;
use std::fmt::Display;
use strum_macros::{EnumIter, EnumString, EnumVariantNames};

use serde::{Deserialize, Serialize};

use crate::analyses::coverage_colors::CoverageColors;
use crate::analyses::coverage_line::CoverageLine;
use crate::analyses::growth::Growth;
use crate::analyses::hist::Hist;
use crate::analyses::info::Info;
use crate::analyses::node_distribution::NodeDistribution;
use crate::analyses::ordered_histgrowth::OrderedHistgrowth;
use crate::analyses::regional_growth::RegionalGrowth;
use crate::analyses::regional_variation::RegionalVariation;
use crate::analyses::section_growth::SectionGrowth;
use crate::analyses::similarity::Similarity;
use crate::analyses::table::Table;
use crate::analyses::Analysis;
use crate::file_formats::vcf_parser::VcfCountType;
use crate::util::CountType;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
pub enum FileRun {
    Gfa {
        graph: String,
        #[serde(default)]
        subset: String,
        #[serde(default)]
        exclude: String,
        reference: Option<String>,
        grouping: Option<Grouping>,
        #[serde(default)]
        nice: bool,
        count_type: CountType,
        analyses: Vec<AnalysisParameter>,
    },
    Vcf {
        file: String,
        #[serde(default)]
        split_haplotypes: bool,
        #[serde(default)]
        count_type: VcfCountType,
        exclude: Option<String>,
        analyses: Vec<AnalysisParameter>,
    },
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
pub enum AnalysisParameter {
    Hist {},
    Growth {
        coverage: Option<String>,
        quorum: Option<String>,
        #[serde(default)]
        add_hist: bool,
        #[serde(default)]
        add_alpha: bool,
    },
    Table {
        total: bool,
        order: Option<String>,
    },
    NodeDistribution {
        #[serde(default = "get_radius")]
        radius: u32,
        #[serde(default = "get_threshold")]
        threshold: usize,
    },
    Info,
    OrderedGrowth {
        coverage: Option<String>,
        quorum: Option<String>,
        order: Option<String>,
    },
    CoverageLine {
        reference: Option<String>,
    },
    Similarity {
        #[serde(default)]
        cluster_method: ClusterMethod,
    },
    Custom {
        name: String,
        file: String,
    },
    RegionalVariation {
        #[serde(default = "get_window_size")]
        window_size: usize,
        order: Option<String>,
        windows: Option<String>,
        #[serde(default)]
        normalize_windows: bool,
    },
    RegionalGrowth {
        #[serde(default = "get_window_size")]
        window_size: usize,
        order: Option<String>,
        windows: Option<String>,
        #[serde(default)]
        normalize_windows: bool,
    },
    SectionGrowth {
        sections: String,
        coverage: Option<String>,
        quorum: Option<String>,
    },
    CoverageColors,
}

impl AnalysisParameter {
    /// This match has to be extended when adding a new analysis. Each analysis
    /// needs to decide, whether it needs the full matrix or whether just the
    /// hist is fine. Just using the histogram is always preferred.
    pub fn to_analysis(self) -> Analysis {
        match self {
            Self::Hist {} => Analysis::HistBased(Box::new(Hist::new())),
            Self::Growth {
                coverage,
                quorum,
                add_hist,
                add_alpha,
            } => Analysis::HistBased(Box::new(Growth::new(coverage, quorum, add_hist, add_alpha))),
            Self::Table { total, order } => {
                Analysis::MatrixBased(Box::new(Table::new(total, order)))
            }
            Self::OrderedGrowth {
                coverage,
                quorum,
                order,
            } => Analysis::MatrixBased(Box::new(OrderedHistgrowth::new(coverage, quorum, order))),
            Self::NodeDistribution { radius, threshold } => {
                Analysis::MatrixBased(Box::new(NodeDistribution::new(radius, threshold)))
            }
            Self::CoverageColors {} => Analysis::MatrixBased(Box::new(CoverageColors::new())),
            Self::CoverageLine { reference } => {
                Analysis::MatrixBased(Box::new(CoverageLine::new(reference)))
            }
            Self::SectionGrowth {
                sections,
                coverage,
                quorum,
            } => Analysis::MatrixBased(Box::new(SectionGrowth::new(sections, coverage, quorum))),
            Self::Similarity { cluster_method } => {
                Analysis::MatrixBased(Box::new(Similarity::new(cluster_method)))
            }
            Self::Info => Analysis::MatrixBased(Box::new(Info::new())),
            Self::RegionalVariation {
                window_size,
                order,
                windows,
                normalize_windows,
            } => Analysis::MatrixBased(Box::new(RegionalVariation::new(
                window_size,
                order,
                windows,
                normalize_windows,
            ))),
            Self::RegionalGrowth {
                window_size,
                order,
                windows,
                normalize_windows,
            } => Analysis::MatrixBased(Box::new(RegionalGrowth::new(
                window_size,
                order,
                windows,
                normalize_windows,
            ))),
            _ => unimplemented!("Other analyses have not been yet implemented"),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash, Clone, PartialOrd, Ord)]
pub enum Grouping {
    Sample,
    Haplotype,
    Custom(String),
}

impl Display for Grouping {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Sample => write!(f, "Group By Sample"),
            Self::Haplotype => write!(f, "Group By Haplotype"),
            Self::Custom(file) => write!(f, "Group By {}", file),
        }
    }
}

fn get_radius() -> u32 {
    20
}

fn get_window_size() -> usize {
    1_000_000
}

fn get_threshold() -> usize {
    1_000
}

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    EnumString,
    EnumVariantNames,
    EnumIter,
    Hash,
    Eq,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "lowercase")]
pub enum ClusterMethod {
    Single,
    Complete,
    Average,
    Weighted,
    Ward,
    Centroid,
    Median,
}

impl Default for ClusterMethod {
    fn default() -> Self {
        Self::Centroid
    }
}

impl ClusterMethod {
    pub fn to_kodama(self) -> kodama::Method {
        match self {
            Self::Single => kodama::Method::Single,
            Self::Complete => kodama::Method::Complete,
            Self::Average => kodama::Method::Average,
            Self::Weighted => kodama::Method::Weighted,
            Self::Ward => kodama::Method::Ward,
            Self::Centroid => kodama::Method::Centroid,
            Self::Median => kodama::Method::Median,
        }
    }
}

impl fmt::Display for ClusterMethod {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(
            formatter,
            "{}",
            match self {
                Self::Single => "single",
                Self::Complete => "complete",
                Self::Average => "average",
                Self::Weighted => "weighted",
                Self::Ward => "ward",
                Self::Centroid => "centroid",
                Self::Median => "median",
            }
        )
    }
}
