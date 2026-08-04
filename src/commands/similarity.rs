use crate::clap_enum_variants_no_all;
use clap::{arg, Arg, ArgMatches, Command};
use strum::VariantNames;

use crate::analysis_parameter::{AnalysisParameter, ClusterMethod, FileRun, Grouping};
use crate::util::CountType;

pub fn get_subcommand() -> Command {
    Command::new("similarity")
        .about("Compute coverage table for count type")
        .args(&[
            arg!(gfa_file: <GFA_FILE> "graph in GFA1 format, accepts also compressed (.gz) file"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample"),
            arg!(-a --"total" "Summarize by totaling presence/absence over all groups"),
            Arg::new("count").help("Graph quantity to be counted").default_value("node").ignore_case(true).short('c').long("count").value_parser(clap_enum_variants_no_all!(CountType)),
            Arg::new("cluster_method").help("Method for clustering results").default_value("centroid").ignore_case(true).short('m').long("method").value_parser(clap_enum_variants_no_all!(ClusterMethod)),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<anyhow::Result<Vec<FileRun>>> {
    if let Some(args) = args.subcommand_matches("similarity") {
        let graph = args
            .get_one::<String>("gfa_file")
            .expect("ordered-histgrowth has gfa file")
            .to_owned();
        let count = args
            .get_one::<CountType>("count")
            .expect("hist subcommand has count type")
            .to_owned();
        let cluster_method = args
            .get_one::<ClusterMethod>("cluster_method")
            .expect("hist subcommand has count type")
            .to_owned();
        let subset = args
            .get_one::<String>("subset")
            .cloned()
            .unwrap_or_default();
        let exclude = args
            .get_one::<String>("exclude")
            .cloned()
            .unwrap_or_default();
        let grouping = args.get_one::<String>("groupby").cloned();
        let grouping = if args.get_flag("groupby-sample") {
            Some(Grouping::Sample)
        } else if args.get_flag("groupby-haplotype") {
            Some(Grouping::Haplotype)
        } else {
            grouping.map(Grouping::Custom)
        };
        let parameters = vec![FileRun::Gfa {
            graph,
            subset,
            exclude,
            grouping,
            reference: None,
            nice: false,
            count_type: count,
            analyses: vec![AnalysisParameter::Similarity { cluster_method }],
        }];
        // log::info!("{parameters:?}");
        Some(Ok(parameters))
    } else {
        None
    }
}
