use std::io::BufReader;
use std::{env, process};
use std::{fs::File, io::Read, io::Write};

use anyhow::anyhow;
use clap::{Arg, ArgAction, ArgMatches, Command};

use crate::analysis_parameter::FileRun;

const TEMPLATE: &str = "- !Gfa
  graph: ../graphs/my_graph_1.gfa
  grouping: Haplotype
  count_type: Bp
#  subset: ../graphs/my_graph_1.genes.bed
#  exclude: ../graphs/my_exclude.bed
#  reference: CHM13#0
#  nice: false
  analyses:
    - !Hist
    - !Growth
      coverage: 1,1,2
      quorum: 0,0.9,0
#    - !NodeDistribution
#      radius: 20
#      threshold: 1000
#    - !Info
#    - !CoverageLine
#      reference: \"GRCh38#0\"
#    - !Similarity
#      cluster_method: Centroid
#    - !RegionalVariation
#      window_size: 100000000
#    - !RegionalGrowth
#      window_size: 100000000
#    - !SectionGrowth
#      sections: ../graphs/my_sections.tsv
#      coverage: 1,1,2
#      quorum: 0,0.9,0
#    - !Custom
#      name: Csv
#      file: ../tables/a_simple_table.csv
#    - !Custom
#      name: Pdf
#      file: ../reports/a_simple_report.pdf
#
#- !Vcf:
#  file: ../variants/my_vcf.vcf
#  split_haplotypes: false
#  count_type: AltLength
#  analyses:
#    - !Hist
#    - !Growth
#      coverage: 1,1,2
#      quorum: 0,0.9,0
#    ... (supports all of the same analyses as !Gfa)
#
#- !Hist
#  file: ../histograms/my_hist.tsv
#  count_type: Node
#  analyses:
#    - !Hist
#    - !Growth
#      coverage: 1,1,2
#      quorum: 0,0.9,0
#
#- !Config
#  output: Html";

pub fn get_subcommand() -> Command {
    Command::new("report")
        .about("Create an html report from a YAML config file")
        .args(&[Arg::new("yaml_file")
            .required(false)
            .help("Specifies yaml config")])
        .args(&[Arg::new("dry_run")
            .required(false)
            .long("dry-run")
            .short('d')
            .action(ArgAction::SetTrue)
            .help(
                "If set, no actual computation is done, only the planned computation will be shown",
            )])
        .args(&[Arg::new("json")
                .required(false)
                .long("json")
                .short('j')
                .action(ArgAction::SetTrue)
                .help(
                    "Instead of an HTML report, a json result will be delivered. These can later be combined and rendered as a single HTML.",
            )])
        .arg(Arg::new("variable").short('r').long("variable").action(ArgAction::Append).help("A variable in the format --variable MYVAR=\"test_text\". This will replace the string \"{{MYVAR}}\" in the YAML file with \"test_text\""))
}

pub fn get_instructions(args: &ArgMatches) -> Option<anyhow::Result<(Vec<FileRun>, String)>> {
    args.subcommand_matches("report").map(parse_report_args)
}

fn parse_report_args(args: &ArgMatches) -> anyhow::Result<(Vec<FileRun>, String)> {
    let yaml_file = match args.get_one::<String>("yaml_file").cloned() {
        Some(f) => f,
        None => {
            let named_tempfile = tempfile::Builder::new()
                .suffix(".yaml")
                .disable_cleanup(true)
                .tempfile()?;
            let file_path = named_tempfile
                .path()
                .to_str()
                .ok_or_else(|| anyhow!("Could not get tempfile path"))?
                .to_string();
            let mut output = File::create(&file_path)?;
            write!(output, "{}", TEMPLATE)?;
            let editor = env::var("EDITOR")?;
            process::Command::new(editor).arg(&file_path).status()?;
            file_path
        }
    };
    let f = File::open(&yaml_file)?;
    let mut reader = BufReader::new(f);
    let mut buffer = String::new();
    reader.read_to_string(&mut buffer)?;
    let variables = args
        .get_many::<String>("variable")
        .unwrap_or_default()
        .map(|v| v.as_str())
        .collect::<Vec<_>>();
    for variable in variables {
        let mut iter = variable.split("=");
        let key = iter
            .next()
            .ok_or_else(|| anyhow!("Empty variable parameter, use it like --variable KEY=value"))?;
        let value = iter.next().ok_or_else(|| {
            anyhow!("Variable parameter without '=', use it like --variable KEY=value")
        })?;
        // Turns into "{{KEY}}"
        let text = format!("{{{{{key}}}}}");
        buffer = buffer.replace(text.as_str(), value);
    }
    let contents = serde_yaml::from_str(&buffer)?;
    Ok((contents, yaml_file))
}
