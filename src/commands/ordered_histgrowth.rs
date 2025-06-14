use crate::clap_enum_variants_no_all;
use clap::{arg, Arg, ArgMatches, Command};
use strum::VariantNames;

use crate::analysis_parameter::{AnalysisParameter, AnalysisRun, Grouping};
use crate::util::CountType;

pub fn get_subcommand() -> Command {
    Command::new("ordered-histgrowth")
        .about("Calculate growth curve based on group file order (if order is unspecified, use path order in GFA)")
        .args(&[
            arg!(gfa_file: <GFA_FILE> "graph in GFA1 format, accepts also compressed (.gz) file"),
            arg!(-s --subset <FILE> "Produce counts by subsetting the graph to a given list of paths (1-column list) or path coordinates (3- or 12-column BED file)"),
            arg!(-e --exclude <FILE> "Exclude bp/node/edge in growth count that intersect with paths (1-column list) or path coordinates (3- or 12-column BED-file) provided by the given file; all intersecting bp/node/edge will be exluded also in other paths not part of the given list"),
            arg!(-g --groupby <FILE> "Merge counts from paths by path-group mapping from given tab-separated two-column file"),
            arg!(-H --"groupby-haplotype" "Merge counts from paths belonging to same haplotype"),
            arg!(-S --"groupby-sample" "Merge counts from paths belonging to same sample"),
            arg!(-O --order <FILE> "The ordered histogram will be produced according to order of paths/groups in the supplied file (1-column list). If this option is not used, the order is determined by the rank of paths/groups in the subset list, and if that option is not used, the order is determined by the rank of paths/groups in the GFA file."),
            Arg::new("count").help("Graph quantity to be counted").default_value("node").ignore_case(true).short('c').long("count").value_parser(clap_enum_variants_no_all!(CountType)),
            Arg::new("coverage").help("Ignore all countables with a coverage lower than the specified threshold. The coverage of a countable corresponds to the number of path/walk that contain it. Repeated appearances of a countable in the same path/walk are counted as one. You can pass a comma-separated list of coverage thresholds, each one will produce a separated growth curve (e.g., --coverage 2,3). Use --quorum to set a threshold in conjunction with each coverage (e.g., --quorum 0.5,0.9)")
                .short('l').long("coverage").default_value("1"),
            Arg::new("quorum").help("Unlike the --coverage parameter, which specifies a minimum constant number of paths for all growth point m (1 <= m <= num_paths), --quorum adjust the threshold based on m. At each m, a countable is counted in the average growth if the countable is contained in at least floor(m*quorum) paths. Example: A quorum of 0.9 requires a countable to be in 90% of paths for each subset size m. At m=10, it must appear in at least 9 paths. At m=100, it must appear in at least 90 paths. A quorum of 1 (100%) requires presence in all paths of the subset, corresponding to the core. Default: 0, a countable counts if it is present in any path at each growth point. Specify multiple quorum values with a comma-separated list (e.g., --quorum 0.5,0.9). Use --coverage to set static path thresholds in conjunction with variable quorum percentages (e.g., --coverage 5,10).")
                .short('q').long("quorum").default_value("0"),
        ])
}

pub fn get_instructions(args: &ArgMatches) -> Option<anyhow::Result<Vec<AnalysisRun>>> {
    if let Some(args) = args.subcommand_matches("ordered-histgrowth") {
        let count = args
            .get_one::<CountType>("count")
            .expect("hist subcommand has count type")
            .to_owned();
        let order = args.get_one::<String>("order").cloned();
        let coverage = args.get_one::<String>("coverage").cloned();
        let quorum = args.get_one::<String>("quorum").cloned();
        let graph = args
            .get_one::<String>("gfa_file")
            .expect("hist subcommand has gfa file")
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
            grouping.map(|g| Grouping::Custom(g))
        };
        let parameters = vec![AnalysisRun::new(
            graph,
            None,
            subset,
            exclude,
            grouping,
            false,
            vec![AnalysisParameter::OrderedGrowth {
                coverage,
                quorum,
                count_type: count,
                order,
            }],
        )];
        log::info!("{parameters:?}");
        Some(Ok(parameters))
    } else {
        None
    }
}
