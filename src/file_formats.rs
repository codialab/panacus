use crate::{coverage_matrix::CoverageMatrix, hist::Hist};

pub mod gfa_parser;
pub mod hist_parser;
pub mod vcf_parser;

pub trait FileFormatParser {
    fn generate_hist(self: Box<Self>) -> Hist;
    fn generate_matrix(self: Box<Self>) -> CoverageMatrix;
}
