use std::io::BufReader;
use std::{fs::File, io::Read};

use anyhow::anyhow;
use clap::{Arg, ArgAction, ArgMatches, Command};

use crate::analysis_parameter::FileRun;

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

pub fn get_instructions(args: &ArgMatches) -> Option<Result<Vec<FileRun>, anyhow::Error>> {
    args.subcommand_matches("report").map(parse_report_args)
}

fn parse_report_args(args: &ArgMatches) -> Result<Vec<FileRun>, anyhow::Error> {
    if let Some(yaml_file) = args.get_one::<String>("yaml_file").cloned() {
        let f = File::open(yaml_file)?;
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
            let key = iter.next().ok_or_else(|| {
                anyhow!("Empty variable parameter, use it like --variable KEY=value")
            })?;
            let value = iter.next().ok_or_else(|| {
                anyhow!("Variable parameter without '=', use it like --variable KEY=value")
            })?;
            // Turns into "{{KEY}}"
            let text = format!("{{{{{key}}}}}");
            buffer = buffer.replace(text.as_str(), value);
        }
        let contents = serde_yaml::from_str(&buffer)?;
        Ok(contents)
    } else {
        println!(
            "
# Missing YAML file!
#
# Example YAML:
# To get started copy this into a .yaml file and edit it

- graph: ../graphs/test_graph.gfa
  grouping: Haplotype
  count_type: Bp
  analyses:
    - !Hist
    - !Growth
      coverage: 1,1,2
      quorum: 0,0.9,0

# For more information see: https://github.com/codialab/panacus/wiki
                "
        );
        Err(anyhow!("Missing YAML configuration!"))
    }
}
