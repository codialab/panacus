use clap::{arg, Arg, ArgMatches, Command};

use crate::analysis_parameter::{AnalysisParameter, FileRun, Grouping};

pub fn get_subcommand() -> Command {
    Command::new("node-distribution")
        .about("Return the list of bins with there coverages, log10-lengths and log10-sizes. Due to this being the values for the centers of the hexagons shown in the html plot and not real values, some values might be negative.")
        .args(&[
            arg!(gfa_file: <GFA_FILE> "graph in GFA1 format, accepts also compressed (.gz) file"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file) (ONLY IN GFA MODE)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list (ONLY IN GFA MODE)"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file (ONLY IN GFA MODE)"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype (ONLY IN GFA MODE)"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample (ONLY IN GFA MODE)"),
            Arg::new("radius")
                .help("Radius of the hexagons used to bin")
                .short('r')
                .long("radius")
                .value_parser(clap::value_parser!(u32))
                .default_value("20"),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<Result<Vec<FileRun>, anyhow::Error>> {
    if let Some(args) = args.subcommand_matches("node-distribution") {
        let graph = args
            .get_one::<String>("gfa_file")
            .expect("info subcommand has gfa file")
            .to_owned();
        let radius = args
            .get_one::<u32>("radius")
            .expect("node-distribution has radius")
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
            count_type: crate::util::CountType::Bp,
            analyses: vec![AnalysisParameter::NodeDistribution {
                radius,
                threshold: 1000,
            }],
        }];
        log::info!("{parameters:?}");
        Some(Ok(parameters))
    } else {
        None
    }
}
