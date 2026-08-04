use clap::{arg, Arg, ArgMatches, Command};

use crate::{
    analysis_parameter::{AnalysisParameter, FileRun, Grouping},
    clap_enum_variants,
    util::CountType,
};

pub fn get_subcommand() -> Command {
    Command::new("growth")
        .visible_alias("histgrowth")
        .about("Calculate growth curve from coverage histogram")
        .args(&[
            arg!(file: <FILE> "EITHER graph in GFA1 format, accepts also compressed (.gz) file OR a histogram as a .tsv"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file) (ONLY IN GFA MODE)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list (ONLY IN GFA MODE)"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file (ONLY IN GFA MODE)"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype (ONLY IN GFA MODE)"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample (ONLY IN GFA MODE)"),
            arg!(-i --hist "Also include histogram in output (ONLY IN GFA MODE)"),
            arg!(-a --alpha "Include alpha value as a comment in the output"),
            Arg::new("count").help("Graph quantity to be counted").default_value("node").ignore_case(true).short('c').long("count").value_parser(clap_enum_variants!(CountType)),
            Arg::new("coverage").help("Ignore all countables with a coverage lower than the specified threshold. The coverage of a countable corresponds to the number of path/walk that contain it. Repeated appearances of a countable in the same path/walk are counted as one. You can pass a comma-separated list of coverage thresholds, each one will produce a separated growth curve (e.g., --coverage 2,3). Use --quorum to set a threshold in conjunction with each coverage (e.g., --quorum 0.5,0.9)")
            .short('l').long("coverage").default_value("1"),
            Arg::new("quorum").help("Unlike the --coverage parameter, which specifies a minimum constant number of paths for all growth point m (1 <= m <= num_paths), --quorum adjust the threshold based on m. At each m, a countable is counted in the average growth if the countable is contained in at least floor(m*quorum) paths. Example: A quorum of 0.9 requires a countable to be in 90% of paths for each subset size m. At m=10, it must appear in at least 9 paths. At m=100, it must appear in at least 90 paths. A quorum of 1 (100%) requires presence in all paths of the subset, corresponding to the core. Default: 0, a countable counts if it is present in any path at each growth point. Specify multiple quorum values with a comma-separated list (e.g., --quorum 0.5,0.9). Use --coverage to set static path thresholds in conjunction with variable quorum percentages (e.g., --coverage 5,10).")
            .short('q').long("quorum").default_value("0"),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<Result<Vec<FileRun>, anyhow::Error>> {
    if let Some(args) = args.subcommand_matches("growth") {
        // let hist = args.get_one::<String>("hist_file").expect("").to_owned();
        let coverage = args.get_one::<String>("coverage").cloned();
        let quorum = args.get_one::<String>("quorum").cloned();
        let add_hist = args.get_flag("hist");
        let add_alpha = args.get_flag("alpha");
        let count = args
            .get_one::<CountType>("count")
            .expect("hist subcommand has count type")
            .to_owned();
        let graph = args
            .get_one::<String>("file")
            .expect("growth subcommand has gfa file")
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
        Some(Ok(vec![FileRun::Gfa {
            graph,
            subset,
            exclude,
            grouping,
            reference: None,
            nice: false,
            count_type: count,
            analyses: vec![AnalysisParameter::Growth {
                coverage,
                quorum,
                add_hist,
                add_alpha,
            }],
        }]))
    } else {
        None
    }
}
