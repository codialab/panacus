use clap::{arg, ArgMatches, Command};

use crate::analysis_parameter::{AnalysisParameter, FileRun, Grouping};

pub fn get_subcommand() -> Command {
    Command::new("info")
        .about("Return general graph and paths info")
        .args(&[
            arg!(gfa_file: <GFA_FILE> "graph in GFA1 format, accepts also compressed (.gz) file"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample"),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<Result<Vec<FileRun>, anyhow::Error>> {
    if let Some(args) = args.subcommand_matches("info") {
        let graph = args
            .get_one::<String>("gfa_file")
            .expect("info subcommand has gfa file")
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
            count_type: crate::util::CountType::Node,
            analyses: vec![AnalysisParameter::Info],
        }];
        log::info!("{parameters:?}");
        Some(Ok(parameters))
    } else {
        None
    }
}
