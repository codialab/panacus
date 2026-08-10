/* private use */
pub mod analyses;
mod analysis_parameter;
mod commands;
pub mod coverage_matrix;
pub mod file_formats;
pub mod hist;
mod html_report;
mod io;
mod util;

use env_logger::Builder;
use itertools::Itertools;
use log::LevelFilter;
use std::io::Read;
use std::{fmt::Debug, io::Write};
use thiserror::Error;

use analyses::Analysis;
use analysis_parameter::AnalysisParameter;
use clap::{Arg, ArgAction, ArgMatches, Command};
use html_report::AnalysisSection;

use std::fs::File;
use std::io::BufReader;

use shadow_rs::shadow;

use crate::analyses::{HistBasedAnalysis, MatrixBasedAnalysis};
use crate::analysis_parameter::{FileRun, Grouping};
use crate::coverage_matrix::CoverageMatrix;
use crate::file_formats::gfa_parser::{GfaParser, GraphMaskParameters};
use crate::file_formats::vcf_parser::VcfParser;
use crate::file_formats::FileFormatParser;
use crate::hist::Hist;

shadow!(build);

#[macro_export]
macro_rules! clap_enum_variants {
    // Credit: Johan Andersson (https://github.com/repi)
    // Code from https://github.com/clap-rs/clap/discussions/4264
    ($e: ty) => {{
        use clap::builder::TypedValueParser;
        use strum::VariantNames;
        clap::builder::PossibleValuesParser::new(<$e>::VARIANTS).map(|s| s.parse::<$e>().unwrap())
    }};
}

#[macro_export]
macro_rules! clap_enum_variants_no_all {
    ($e: ty) => {{
        use clap::builder::TypedValueParser;
        clap::builder::PossibleValuesParser::new(<$e>::VARIANTS.iter().filter(|&x| x != &"all"))
            .map(|s| s.parse::<$e>().unwrap())
    }};
}

#[macro_export]
macro_rules! some_or_return {
    ($x:expr, $y:expr) => {
        match $x {
            Some(v) => v,
            None => return $y,
        }
    };
}

type HistBasedAnalyses = Vec<Box<dyn HistBasedAnalysis>>;
type MatrixBasedAnalyses = Vec<Box<dyn MatrixBasedAnalysis>>;

fn set_number_of_threads(params: &ArgMatches) {
    //if num_threads is 0 then the Rayon will select
    //the number of threads to the core number automatically
    let threads = params.get_one("threads").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(*threads)
        .build_global()
        .expect("Failed to initialize global thread pool");
    log::info!(
        "running panacus on {} threads",
        rayon::current_num_threads()
    );
}

fn set_verbosity(args: &ArgMatches) {
    if args.get_flag("verbose") {
        Builder::new().filter_level(LevelFilter::Debug).init();
    } else {
        Builder::new().filter_level(LevelFilter::Info).init();
    }
}

pub fn run_cli() -> Result<(), anyhow::Error> {
    let mut out = std::io::BufWriter::new(std::io::stdout());

    // read parameters and store them in memory
    // let params = cli::read_params();
    let args = Command::new("panacus")
        .subcommand(commands::render::get_subcommand())
        .subcommand(commands::report::get_subcommand())
        .subcommand(commands::hist::get_subcommand())
        .subcommand(commands::growth::get_subcommand())
        .subcommand(commands::info::get_subcommand())
        .subcommand(commands::ordered_histgrowth::get_subcommand())
        .subcommand(commands::table::get_subcommand())
        .subcommand(commands::node_distribution::get_subcommand())
        .subcommand(commands::similarity::get_subcommand())
        .subcommand(commands::coverage_colors::get_subcommand())
        .subcommand(commands::section_growth::get_subcommand())
        .subcommand_required(true)
        .arg(
            Arg::new("threads")
                .short('t')
                .action(ArgAction::Set)
                .value_name("COUNT")
                .default_value("0")
                .value_parser(clap::value_parser!(usize))
                .global(true)
                .help("Set the number of threads used (default: use all threads)"),
        )
        .arg(
            Arg::new("verbose")
                .short('v')
                .long("verbose")
                .action(ArgAction::SetTrue)
                .global(true)
                .help("Set the number of threads used (default: use all threads)"),
        )
        .long_version(build::CLAP_LONG_VERSION)
        .get_matches();

    set_verbosity(&args);
    set_number_of_threads(&args);

    let mut instructions: Vec<FileRun> = Vec::new();
    let mut shall_write_html = false;
    let mut dry_run = false;
    let mut json = false;
    let mut config_content = "EMPTY".to_string();

    if let Some(args) = args.subcommand_matches("render") {
        let json_files: Vec<String> = args
            .get_many::<String>("json_files")
            .unwrap()
            .cloned()
            .collect();
        let mut full_report = Vec::new();
        for file_path in &json_files {
            let file = File::open(file_path)?;
            let reader = BufReader::new(file);

            // Read the JSON contents of the file as an instance of `User`.
            let report: Vec<AnalysisSection> = serde_json::from_reader(reader)?;
            full_report.extend(report);
        }
        let mut registry = handlebars::Handlebars::new();
        let report_text = AnalysisSection::generate_report(
            full_report,
            &mut registry,
            &json_files[0],
            "-- GENERATED VIA RENDER --",
        )?;
        writeln!(&mut out, "{report_text}")?;
        return Ok(());
    }

    if let Some(report) = commands::report::get_instructions(&args) {
        shall_write_html = true;
        if let Some(report_matches) = args.subcommand_matches("report") {
            dry_run = report_matches.get_flag("dry_run");
            json = report_matches.get_flag("json");
            let config = report_matches
                .get_one::<String>("yaml_file")
                .expect("Contains required yaml config")
                .to_owned();
            let f = File::open(config)?;
            let mut reader = BufReader::new(f);
            config_content = String::new();
            reader.read_to_string(&mut config_content)?;
            config_content.push_str("\n---------------------\n");
            config_content.push_str(serde_yaml::to_string(&report.as_ref().unwrap())?.as_str());
        }
        instructions.extend(report?);
    }
    if let Some(hist) = commands::hist::get_instructions(&args) {
        instructions.extend(hist?);
    }
    if let Some(growth) = commands::growth::get_instructions(&args) {
        instructions.extend(growth?);
    }
    if let Some(info) = commands::info::get_instructions(&args) {
        instructions.extend(info?);
    }
    if let Some(coverage_colors) = commands::coverage_colors::get_instructions(&args) {
        instructions.extend(coverage_colors?);
    }
    if let Some(section_growth) = commands::section_growth::get_instructions(&args) {
        instructions.extend(section_growth?);
    }
    if let Some(ordered_histgrowth) = commands::ordered_histgrowth::get_instructions(&args) {
        instructions.extend(ordered_histgrowth?);
    }
    if let Some(table) = commands::table::get_instructions(&args) {
        instructions.extend(table?);
    }
    if let Some(counts) = commands::node_distribution::get_instructions(&args) {
        instructions.extend(counts?);
    }
    if let Some(similarity) = commands::similarity::get_instructions(&args) {
        instructions.extend(similarity?);
    }

    log::info!("{:?}", instructions);

    // ride on!
    if !dry_run {
        execute_pipeline(
            instructions,
            &mut out,
            &config_content,
            shall_write_html,
            json,
        )?;
    } else {
        println!("{:#?}", instructions);
    }

    // clean up & close down
    out.flush()?;
    Ok(())
}

fn execute_pipeline<W: Write>(
    instructions: Vec<FileRun>,
    out: &mut std::io::BufWriter<W>,
    config_content: &str,
    shall_write_html: bool,
    json: bool,
) -> anyhow::Result<()> {
    for file in instructions {
        let (file_parser, analyses) = get_file_parser(file)?;
        let (hist_based, matrix_based) = split_analyses(analyses);
        log::info!("File parser generated");
        if !matrix_based.is_empty() {
            let matrix = file_parser.generate_matrix();
            log::info!("Matrix generated");
            let report = match shall_write_html {
                true => get_matrix_reports(matrix, matrix_based, hist_based, config_content, json)?,
                false => get_matrix_tables(matrix, matrix_based, hist_based)?,
            };
            writeln!(out, "{report}")?;
        } else {
            let hist = file_parser.generate_hist();
            log::info!("Hist generated");
            let report = match shall_write_html {
                true => get_hist_reports(hist, hist_based, config_content, json)?,
                false => get_hist_tables(hist, hist_based)?,
            };
            writeln!(out, "{report}")?;
        }
    }
    Ok(())
}

fn get_hist_reports(
    hist: Hist,
    hist_based: HistBasedAnalyses,
    config_content: &str,
    json: bool,
) -> anyhow::Result<String> {
    let reports: Vec<AnalysisSection> = hist_based
        .into_iter()
        .filter_map(|mut x| {
            log::info!("Doing hist-analysis: {}", x.get_type());
            let out = x.generate_report_section(&hist).ok();
            log::info!("Finished hist-analysis: {}", x.get_type());
            out
        })
        .flatten()
        .collect();
    let mut registry = handlebars::Handlebars::new();
    if json {
        let report = serde_json::to_string_pretty(&reports)?;
        Ok(report)
    } else {
        let report = AnalysisSection::generate_report(
            reports,
            &mut registry,
            "<Placeholder Filename>",
            config_content,
        )?;
        Ok(report)
    }
}

fn get_hist_tables(hist: Hist, hist_based: HistBasedAnalyses) -> anyhow::Result<String> {
    let reports: Vec<String> = hist_based
        .into_iter()
        .filter_map(|mut x| x.generate_table(&hist).ok())
        .collect();
    let report = reports.into_iter().join("\n\n\n");
    Ok(report)
}

fn get_matrix_tables(
    matrix: CoverageMatrix,
    matrix_based: MatrixBasedAnalyses,
    hist_based: HistBasedAnalyses,
) -> anyhow::Result<String> {
    let mut reports: Vec<String> = matrix_based
        .into_iter()
        // TODO remove all filter_maps and replace with correct error handling
        .filter_map(|mut x| x.generate_table(&matrix).ok())
        .collect();
    let hist = matrix.get_hist();
    let mut hist_reports: Vec<String> = hist_based
        .into_iter()
        .filter_map(|mut x| x.generate_table(&hist).ok())
        .collect();
    reports.append(&mut hist_reports);
    let report = reports.into_iter().join("\n\n\n");
    Ok(report)
}

fn get_matrix_reports(
    matrix: CoverageMatrix,
    matrix_based: MatrixBasedAnalyses,
    hist_based: HistBasedAnalyses,
    config_content: &str,
    json: bool,
) -> anyhow::Result<String> {
    let mut reports: Vec<AnalysisSection> = matrix_based
        .into_iter()
        // TODO remove all filter_maps and replace with correct error handling
        .filter_map(|mut x| {
            log::info!("Doing matrix-analysis: {}", x.get_type());
            let out = x.generate_report_section(&matrix).ok();
            log::info!("Finished matrix-analysis: {}", x.get_type());
            out
        })
        .flatten()
        .collect();
    let hist = matrix.get_hist();
    let mut hist_reports: Vec<AnalysisSection> = hist_based
        .into_iter()
        .filter_map(|mut x| {
            log::info!("Doing hist-analysis: {}", x.get_type());
            let out = x.generate_report_section(&hist).ok();
            log::info!("Finished hist-analysis: {}", x.get_type());
            out
        })
        .flatten()
        .collect();
    reports.append(&mut hist_reports);

    if json {
        let report = serde_json::to_string_pretty(&reports)?;
        Ok(report)
    } else {
        let mut registry = handlebars::Handlebars::new();
        let report = AnalysisSection::generate_report(
            reports,
            &mut registry,
            "<Placeholder Filename>",
            config_content,
        )?;
        Ok(report)
    }
}

fn split_analyses(analyses: Vec<AnalysisParameter>) -> (HistBasedAnalyses, MatrixBasedAnalyses) {
    analyses.into_iter().map(|x| x.into_analysis()).fold(
        (Vec::new(), Vec::new()),
        |(mut hists, mut matrices), item| {
            match item {
                Analysis::HistBased(h) => hists.push(h),
                Analysis::MatrixBased(m) => matrices.push(m),
            }
            (hists, matrices)
        },
    )
}

/// This match has to be extended if a new file type is added.
fn get_file_parser(
    file: FileRun,
) -> anyhow::Result<(Box<dyn FileFormatParser>, Vec<AnalysisParameter>)> {
    match file {
        FileRun::Gfa {
            graph,
            subset,
            exclude,
            grouping,
            reference,
            nice,
            count_type,
            analyses,
        } => {
            let (groupby, groupby_sample, groupby_haplotype) = match grouping {
                Some(Grouping::Sample) => (String::new(), true, false),
                Some(Grouping::Haplotype) => (String::new(), false, true),
                Some(Grouping::Custom(file)) => (file, false, false),
                None => (String::new(), false, false),
            };
            let graph_mask_parameters = GraphMaskParameters {
                positive_list: subset,
                negative_list: exclude,
                groupby,
                groupby_sample,
                groupby_haplotype,
            };
            Ok((
                Box::new(GfaParser::new(
                    &graph,
                    count_type,
                    graph_mask_parameters,
                    reference,
                    nice,
                )?),
                analyses,
            ))
        }
        FileRun::Vcf {
            file,
            split_haplotypes,
            count_type,
            exclude,
            analyses,
        } => Ok((
            Box::new(VcfParser::new(
                &file,
                count_type,
                split_haplotypes,
                exclude,
            )?),
            analyses,
        )),
    }
}

#[derive(Error, Debug)]
pub enum ConfigParseError {
    #[error("no config block with name {name} was found")]
    NameNotFound { name: String },
}
