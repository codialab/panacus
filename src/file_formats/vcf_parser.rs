use std::{
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Read},
    str::SplitWhitespace,
};

use anyhow::Result;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    analyses::info::FileInfo, coverage_matrix::CoverageMatrix, file_formats::FileFormatParser,
    hist::Hist,
};

pub struct VcfParser {
    filename: String,
    count_type: VcfCountType,
    /// The same as groupby_haplotypes, however,
    /// this requires the VCF to have constant ploidy for this sample
    split_haplotypes: bool,
    current_variant_id: usize,
    exclude: Vec<(String, usize, usize)>,
}

impl FileFormatParser for VcfParser {
    fn generate_hist(self: Box<Self>) -> crate::hist::Hist {
        let file = File::open(&self.filename).expect("Can read file");
        let buf_reader = BufReader::new(file);
        self.generate_hist_from_reader(buf_reader)
    }

    fn generate_matrix(self: Box<Self>) -> CoverageMatrix {
        let file = File::open(&self.filename).expect("Can read file");
        let buf_reader = BufReader::new(file);
        self.generate_matrix_from_reader(buf_reader)
    }
}

#[derive(Default)]
struct VcfStats {
    num_samples: usize,
    num_variants: usize,
    num_alleles: usize,
    num_singletons: usize,
}

impl VcfStats {
    fn into_file_info(self) -> FileInfo {
        let mut info = FileInfo::new("vcf");
        info.add_info("Number of samples", &self.num_samples.to_string());
        info.add_info("Number of variants", &self.num_variants.to_string());
        info.add_info("Number of alt alleles", &self.num_alleles.to_string());
        info.add_info(
            "Number of singleton alleles",
            &self.num_singletons.to_string(),
        );
        info
    }
}

impl VcfParser {
    pub fn new(
        filename: &str,
        count_type: VcfCountType,
        split_haplotypes: bool,
        exclude: Option<String>,
    ) -> Result<Self> {
        let exclude = match exclude {
            Some(filename) => {
                let f = File::open(filename).expect("Exclude file must exist");
                let mut reader = BufReader::new(f);
                let mut buf = String::new();
                reader
                    .read_to_string(&mut buf)
                    .expect("Exclude file must be readable");
                let res = buf
                    .lines()
                    .map(|l| {
                        let mut fields = l.split('\t');
                        let ref_name = fields
                            .next()
                            .expect("Exclude file should not contain empty lines")
                            .to_owned();
                        let start = fields
                            .next()
                            .expect("Three columns are required")
                            .parse::<usize>()
                            .unwrap();
                        let end = fields
                            .next()
                            .expect("Three columns are required")
                            .parse::<usize>()
                            .unwrap();
                        (ref_name, start, end)
                    })
                    .collect();
                res
            }
            None => Vec::new(),
        };

        Ok(Self {
            filename: filename.to_owned(),
            count_type,
            split_haplotypes,
            exclude,
            current_variant_id: 1,
        })
    }

    fn generate_hist_from_reader(self: Self, buf_reader: impl BufRead) -> Hist {
        let mut lines = buf_reader.lines().map(|l| l.expect("Failed to read line"));
        let mut header = VcfHeader::parse(&mut lines).expect("Failed to parse header");
        let paths = std::mem::take(&mut header.paths);
        let run_id = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let run_name = run_id.clone();
        let mut hist =
            Hist::from_maximum_coverage(paths.len(), self.count_type.to_string(), run_id, run_name);
        for line in lines {
            let variants = self
                .parse_variant_line_to_count(&line)
                .expect("Failed to parse variant");
            for variant in variants {
                let coverage = variant.value;
                if coverage == 0 {
                    log::info!("Variant has coverage 0: {}", variant.id);
                }
                match self.count_type {
                    VcfCountType::Variants => {
                        hist.insert_feature_of_coverage_and_length(coverage, 1);
                    }
                    VcfCountType::AltLength => {
                        hist.insert_feature_of_coverage_and_length(coverage, variant.alt_length);
                    }
                }
            }
        }
        hist
    }

    fn generate_matrix_from_reader(mut self: Self, buf_reader: impl BufRead) -> CoverageMatrix {
        let mut lines = buf_reader.lines().map(|l| l.expect("Failed to read line"));

        let mut header = VcfHeader::parse(&mut lines).expect("Failed to parse header");
        let paths = std::mem::take(&mut header.paths);

        let run_id = format!(
            "{}-{}-{}",
            &self.filename,
            self.count_type,
            self.get_split_haplotype_text()
        );
        let run_name = run_id.clone();

        let mut matrix = CoverageMatrix::new(
            self.count_type.to_string(),
            run_id,
            run_name,
            FileInfo::new("vcf"), // Pass blank FileInfo initially
        );

        let mut stats = VcfStats {
            num_samples: paths.len(),
            ..Default::default()
        };

        let mut are_paths_initialized = false;

        for line in lines {
            stats.num_variants += 1;

            let variants = self
                .parse_variant_line_to_allele_list(&line)
                .expect("Failed to parse variant");

            for variant in variants {
                let feature = variant.value;

                // If the paths are not yet initialized, initialize them
                if !are_paths_initialized {
                    let matrix_paths = self.generate_path_names(&paths, &feature);
                    matrix.set_path_names(matrix_paths);
                    are_paths_initialized = true;
                }

                let feature_length = match self.count_type {
                    VcfCountType::Variants => 1,
                    VcfCountType::AltLength => variant.alt_length,
                };
                let name = match variant.id {
                    "." => self.get_variant_id(),
                    id => format!("{}-{id}", self.get_variant_id()),
                };

                if feature.iter().flatten().sum::<u32>() == 1 {
                    stats.num_singletons += 1;
                }
                matrix.insert_feature(
                    name,
                    feature_length,
                    (variant.chrom, variant.pos),
                    feature.into_iter().flatten().collect(),
                );
                stats.num_alleles += 1;
            }
        }
        let file_info = stats.into_file_info();
        matrix.set_file_info(file_info);
        matrix
    }

    fn generate_path_names(&self, paths: &[String], feature: &[Vec<u32>]) -> Vec<String> {
        if !self.split_haplotypes {
            return paths.to_vec();
        }
        paths
            .iter()
            .enumerate()
            .flat_map(|(idx, name)| {
                let ploidy = feature[idx].len();
                if ploidy == 1 {
                    vec![format!("{name}#0")] // Haploid individuals' haplotypes are named sample#0 in PanSN spec
                } else {
                    (1..=ploidy)
                        .map(|i| format!("{name}#{i}"))
                        .collect::<Vec<_>>()
                }
            })
            .collect()
    }

    fn get_split_haplotype_text(&self) -> &str {
        match self.split_haplotypes {
            true => "haplotypes",
            false => "samples",
        }
    }

    fn get_variant_id(&mut self) -> String {
        let id = self.current_variant_id;
        self.current_variant_id += 1;
        id.to_string()
    }

    /// Parses a variant and returns a list noting the presence/absence of the variant for each path
    fn parse_variant_line_to_allele_list<'a>(
        &self,
        line: &'a str,
    ) -> Result<Vec<Variant<'a, Vec<Vec<u32>>>>> {
        self.parse_variant_line_base(line, |samples, num_alts, genotype_position| {
            let samples: Vec<&str> = samples.collect();
            let mut allele_lists = if self.split_haplotypes {
                vec![vec![Vec::new(); samples.len()]; num_alts]
            } else {
                vec![vec![vec![0; 1]; samples.len()]; num_alts]
            };
            samples.into_iter().enumerate().for_each(|(idx, sample)| {
                let sample_genotype = sample
                    .split(":")
                    .nth(genotype_position)
                    .expect("Sample needs a genotype");
                if !self.split_haplotypes {
                    sample_genotype
                        .split(&['|', '/'])
                        .filter_map(|x| {
                            if x == "." || x == "0" {
                                None
                            } else {
                                Some(x.parse::<usize>().expect("Allele is integer"))
                            }
                        })
                        .unique()
                        .for_each(|x| {
                            allele_lists[x - 1][idx][0] = 1;
                        });
                } else {
                    sample_genotype.split(&['|', '/']).for_each(|x| {
                        for (i, allele_list) in allele_lists.iter_mut().enumerate().take(num_alts) {
                            if x == "." {
                                allele_list[idx].push(0);
                            } else {
                                let value = x.parse::<usize>().expect("Allele is integer");
                                if (i + 1) == value {
                                    allele_list[idx].push(1);
                                } else {
                                    allele_list[idx].push(0);
                                }
                            }
                        }
                    });
                }
            });
            allele_lists
        })
    }

    /// Parses a variant line into simple counts, how often each allele appears
    fn parse_variant_line_to_count<'a>(&self, line: &'a str) -> Result<Vec<Variant<'a, usize>>> {
        self.parse_variant_line_base(line, |samples, num_alts, genotype_position| {
            let mut counts = vec![0; num_alts];
            samples.for_each(|sample| {
                let sample_genotype = sample
                    .split(":")
                    .nth(genotype_position)
                    .expect("Sample needs a genotype");
                if !self.split_haplotypes {
                    sample_genotype
                        .split(&['|', '/'])
                        .filter_map(|x| {
                            if x == "." || x == "0" {
                                None
                            } else {
                                Some(x.parse::<usize>().expect("Allele is integer"))
                            }
                        })
                        .unique()
                        .for_each(|x| counts[x - 1] += 1);
                } else {
                    sample_genotype.split(&['|', '/']).for_each(|x| {
                        if x != "." {
                            let value = x.parse::<usize>().expect("Allele is integer");
                            if value > 0 {
                                counts[value - 1] += 1;
                            }
                        }
                    });
                }
            });
            counts
        })
    }

    /// Takes a line that contains a variant (so it does not start with #)
    /// and a closure that takes an iterator over the samples, the number of alt alleles
    /// and the position of the genotype in each sample and returns a vec with a value
    /// for each alt allele
    fn parse_variant_line_base<'a, T>(
        &self,
        line: &'a str,
        process_samples: impl FnOnce(SplitWhitespace<'a>, usize, usize) -> Vec<T>,
    ) -> Result<Vec<Variant<'a, T>>> {
        let mut fields = line.split_whitespace();
        let chrom = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs chrom field!"))?;
        let pos: usize = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs pos field!"))?
            .parse()?;
        if !self.is_variant_included(chrom, pos) {
            return Ok(Vec::new());
        }
        let id = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs id field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs ref field!"))?;
        let alt_lengths: Vec<usize> = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs alt field!"))?
            .split(",")
            .map(|x| if x == "." { 0 } else { x.len() })
            .collect();
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs qual field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs filter field!"))?;
        fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs info field!"))?;
        let format = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("Variant needs format field!"))?;
        let genotype_position = format.split(":").position(|x| x == "GT").unwrap_or(0);
        let num_alts = alt_lengths.len();
        let values = process_samples(fields, num_alts, genotype_position);
        Ok(alt_lengths
            .into_iter()
            .zip(values)
            .map(|(alt_length, value)| Variant {
                chrom,
                pos,
                id,
                alt_length,
                value,
            })
            .collect())
    }

    fn is_variant_included(&self, chrom: &str, pos: usize) -> bool {
        !self
            .exclude
            .iter()
            .any(|(r, s, e)| r == chrom && pos >= *s && pos <= *e)
    }
}

struct Variant<'a, T> {
    chrom: &'a str,
    pos: usize,
    id: &'a str,
    alt_length: usize,
    value: T,
}

struct VcfHeader {
    paths: Vec<String>,
}

impl VcfHeader {
    fn parse<I, S>(lines: &mut I) -> Result<Self>
    where
        I: Iterator<Item = S>,
        S: AsRef<str>,
    {
        let mut header = Self { paths: Vec::new() };

        for line in lines {
            let line = line.as_ref().trim_end();
            if line.starts_with("##") {
                continue;
            } else if line.starts_with("#CHROM") {
                let samples: Vec<&str> = line.split("\t").skip(9).collect();
                header.paths = samples.into_iter().map(|x| x.to_string()).collect();
                break;
            }
        }
        Ok(header)
    }
}

/// What is set for the length of a feature
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum VcfCountType {
    /// 1 for every feature (just count the number of variants)
    Variants,
    /// The length of the ALT allele is used
    AltLength,
}

impl Default for VcfCountType {
    fn default() -> Self {
        Self::AltLength
    }
}

impl Display for VcfCountType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Variants => write!(f, "Variant"),
            Self::AltLength => write!(f, "AltAlleleLength"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF_HEADER: &str = "##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP membership, build 129\">
##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">
##FILTER=<ID=q10,Description=\"Quality below 10\">
##FILTER=<ID=s50,Description=\"Less than 50% of samples have data\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description=\"Haplotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003";

    const VCF_TEST: &str = "##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##INFO=<ID=CONFLICT,Number=.,Type=String,Description=\"Sample names for which there are multiple paths in the graph with conflicting alleles\">
##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the snarl tree (0=top level)\">
##INFO=<ID=PS,Number=1,Type=String,Description=\"ID of variant corresponding to parent snarl\">
##contig=<ID=chr10,length=134758134>
##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">
##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">
##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"Decomposed from a complex record using vcflib vcfwave and alignment with WFA2-lib.\">
##INFO=<ID=INV,Number=0,Type=Flag,Description=\"Inversion detected\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	GRCh38	HG00097	HG00099	HG00126	HG00128	HG00133	HG00140	HG00146	HG002	HG00232	HG00235	HG00253	HG00280	HG00290	HG00320	HG00321	HG00323	HG00329	HG00344	HG00350	HG00408	HG00423	HG00438	HG005	HG00544	HG00558	HG00597	HG00609	HG00621	HG00639	HG00642	HG00658	HG00673	HG00706	HG00733	HG00735	HG00738	HG00741	HG01071	HG01074	HG01081	HG01099	HG01106	HG01109	HG01123	HG01150	HG01167	HG01175	HG01192	HG01243	HG01252	HG01255	HG01258	HG01261	HG01346	HG01358	HG01361	HG01433	HG01496	HG01530	HG01784	HG01786	HG01884	HG01891	HG01928	HG01934	HG01940	HG01943	HG01952	HG01960	HG01969	HG01975	HG01978	HG01981	HG01993	HG02004	HG02015	HG02027	HG02040	HG02055	HG02056	HG02071	HG02074	HG02080	HG02083	HG02109	HG02129	HG02132	HG02135	HG02145	HG02148	HG02155	HG02165	HG02178	HG02257	HG02258	HG02273	HG02280	HG02293	HG02300	HG02391	HG02392	HG02451	HG02486	HG02514	HG02523	HG02559	HG02572	HG02583	HG02602	HG02615	HG02622	HG02630	HG02647	HG02668	HG02698	HG02717	HG02723	HG02735	HG02738	HG02809	HG02818	HG02841	HG02886	HG02922	HG02965	HG02976	HG02984	HG03017	HG03041	HG03050	HG03098	HG03130	HG03139	HG03195	HG03209	HG03225	HG03239	HG03270	HG03369	HG03453	HG03470	HG03471	HG03486	HG03516	HG03521	HG03540	HG03579	HG03583	HG03654	HG03669	HG03688	HG03704	HG03710	HG03742	HG03784	HG03804	HG03816	HG03831	HG03834	HG03874	HG03927	HG03942	HG04115	HG04157	HG04160	HG04184	HG04187	HG04199	HG04204	HG04228	HG06807	NA18505	NA18508	NA18522	NA18565	NA18570	NA18608	NA18620	NA18747	NA18879	NA18906	NA18940	NA18943	NA18944	NA18945	NA18948	NA18952	NA18959	NA18960	NA18967	NA18970	NA18971	NA18974	NA18976	NA18982	NA18983	NA19036	NA19043	NA19087	NA19159	NA19185	NA19240	NA19338	NA19391	NA19443	NA19468	NA19682	NA19700	NA19776	NA19835	NA19909	NA20129	NA20282	NA20346	NA20503	NA20752	NA20762	NA20799	NA20805	NA20806	NA20809	NA20827	NA20850	NA20870	NA20905	NA21093	NA21102	NA21106	NA21110	NA21144	NA21309
chr10	1	>13046>13130_1	C	CTAACCCTAAACCCCTAACCCTAACCTAACCC,CTAACCCTAACCAACCCTAACCACCCTAACCCTAACCC,CTAACCCTAACCCT,CTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACC,CTAACCCTAACCTAACCCTAACCCTAACCTAACCC	60	.	AC=1,1,1,1,1;AF=0.0212766,0.0212766,0.0212766,0.0212766,0.0212766;AN=47;NS=232;LV=0;ORIGIN=chr10:3;LEN=31,37,13,64,34;TYPE=ins,ins,ins,ins,ins	GT	.	.|.	0|.	.|0	.|.	0|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|0	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	0|0	.|.	.|.	2|.	.|.	.|.	.|4	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	5|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	3|.	.|.	.|.	.|.	.|.	1|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
chr10	2	>13041>13046_1	T	A,TCCCTCAC	60	.	AC=1,1;AF=0.0217391,0.0217391;AN=46;NS=232;LV=0;ORIGIN=chr10:1;LEN=1,7;TYPE=snp,ins	GT	.	.|.	0|.	.|0	.|.	0|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|0	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|1	.|.	0|0	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	2|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
chr10	3	>13046>13130_2;>13046>13130_4	A	AACCCT,ACACCCTCCCCCTCCCCTACCCTAACCCCCTCACCCTACCCCACCCCCTACCCTCCGACACCCCCTCACCCCGTCCGCTCCGAAACCCCTAACCACCCTAACCCTAACCCTAACCCCTGACCCTC,ACACTAACCCGTAC	60	.	AC=1,1,1;AF=0.0212766,0.0212766,0.0212766;AN=47;NS=232;LV=0;ORIGIN=chr10:3;LEN=5,124,13;TYPE=ins,ins,ins	GT	.	.|.	0|.	.|0	.|.	0|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|2	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	0|0	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|.	.|.	.|1	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	3|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
chr10	4	>13046>13130_5	A	C	60	.	AC=1;AF=0.0212766;AN=47;NS=232;LV=0;ORIGIN=chr10:3;LEN=1;TYPE=snp	GT	.	.|.	0|.	.|0	.|.	0|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|0	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	0|0	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	1|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
chr10	6	>13046>13130_3	C	CCTTCACCCTTAACCTTGAACCCTTAACCCTTAACCCTTAACGTACTTACCCTTAACCCTTAACCCTTACCTTGACTTAACCTGTACTGACCTAACCCTAACTCTACTA	60	.	AC=1;AF=0.0212766;AN=47;NS=232;LV=0;ORIGIN=chr10:3;LEN=108;TYPE=ins	GT	.	.|.	0|.	.|0	.|.	1|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|0	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	0|0	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
chr10	7	>13046>13130_6	C	G	60	.	AC=1;AF=0.0212766;AN=47;NS=232;LV=0;ORIGIN=chr10:3;LEN=1;TYPE=snp	GT	.	.|.	0|.	.|0	.|.	0|.	.|0	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|0	.|.	.|0	.|.	.|.	.|.	0|.	.|0	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	0|0	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|0	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	.|.	.|.	0|0	.|0	0|.	.|.	.|.	0|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|.	0|.	.|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|0	0|.	0|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	0|.	.|.	.|.	.|0	.|.	.|.	.|.	.|0	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	.|.	1|.	0|0	.|.	.|.	.|.	.|.	.|.	.|.
";

    #[test]
    fn test_vcf_header() {
        let mut header_lines = VCF_HEADER.lines();
        let vcf_header = VcfHeader::parse(&mut header_lines).unwrap();
        assert_eq!(vcf_header.paths.len(), 3);
    }

    #[test]
    fn test_vcf_parse_line_to_count_haplotypes() {
        let line = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
        let parser = VcfParser::new("", VcfCountType::Variants, true, None).unwrap();
        let mut variants = parser.parse_variant_line_to_count(line).unwrap();
        assert_eq!(variants.len(), 1);
        let variant = variants.remove(0);
        assert_eq!(variant.chrom, "20");
        assert_eq!(variant.pos, 14370);
        assert_eq!(variant.id, "rs6054257");
        assert_eq!(variant.alt_length, 1);
        assert_eq!(variant.value, 3);
    }

    #[test]
    fn test_vcf_parse_line_to_count_samples() {
        let line = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
        let parser = VcfParser::new("", VcfCountType::Variants, false, None).unwrap();
        let mut variants = parser.parse_variant_line_to_count(line).unwrap();
        assert_eq!(variants.len(), 1);
        let variant = variants.remove(0);
        assert_eq!(variant.chrom, "20");
        assert_eq!(variant.pos, 14370);
        assert_eq!(variant.id, "rs6054257");
        assert_eq!(variant.alt_length, 1);
        assert_eq!(variant.value, 2);
    }

    #[test]
    fn test_vcf_parse_line_to_list_samples() {
        let line = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.";
        let parser = VcfParser::new("", VcfCountType::Variants, false, None).unwrap();
        let mut variants = parser.parse_variant_line_to_allele_list(line).unwrap();
        assert_eq!(variants.len(), 1);
        let variant = variants.remove(0);
        assert_eq!(variant.chrom, "20");
        assert_eq!(variant.pos, 14370);
        assert_eq!(variant.id, "rs6054257");
        assert_eq!(variant.alt_length, 1);
        assert_eq!(variant.value, vec![vec![0], vec![1], vec![1]]);
    }

    #[test]
    fn test_vcf_parse_line_to_list_haplotypes() {
        let line = "20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t./1:43:5:.,.";
        let parser = VcfParser::new("", VcfCountType::Variants, true, None).unwrap();
        let mut variants = parser.parse_variant_line_to_allele_list(line).unwrap();
        assert_eq!(variants.len(), 1);
        let variant = variants.remove(0);
        assert_eq!(variant.chrom, "20");
        assert_eq!(variant.pos, 14370);
        assert_eq!(variant.id, "rs6054257");
        assert_eq!(variant.alt_length, 1);
        assert_eq!(variant.value, vec![vec![0, 0], vec![1, 0], vec![0, 1]]);
    }

    #[test]
    fn test_create_empty_vcf_stats() {
        let stats = VcfStats::default();
        assert_eq!(stats.num_samples, 0);
        assert_eq!(stats.num_variants, 0);
        assert_eq!(stats.num_alleles, 0);
        assert_eq!(stats.num_singletons, 0);
    }

    #[test]
    fn test_vcf_stats_to_file_info() {
        let stats = VcfStats::default();
        let info = stats.into_file_info();
        assert_eq!(info.get_filetype(), "vcf");
        assert_eq!(info.iterate_infos().count(), 4);
    }

    #[test]
    fn test_default_vcf_count_type() {
        let c = VcfCountType::default();
        assert_eq!(c, VcfCountType::AltLength);
    }

    #[test]
    fn test_print_variant_count_type() {
        let c = VcfCountType::Variants;
        assert_eq!(&c.to_string(), "Variant");
    }

    #[test]
    fn test_print_alt_length_count_type() {
        let c = VcfCountType::AltLength;
        assert_eq!(&c.to_string(), "AltAlleleLength");
    }

    #[test]
    fn test_get_split_haplotype_text() {
        let parser = VcfParser::new("", VcfCountType::default(), true, None).unwrap();
        assert_eq!(parser.get_split_haplotype_text(), "haplotypes");
    }

    #[test]
    fn test_get_split_sample_text() {
        let parser = VcfParser::new("", VcfCountType::default(), false, None).unwrap();
        assert_eq!(parser.get_split_haplotype_text(), "samples");
    }

    #[test]
    fn test_get_variant_id() {
        let mut parser = VcfParser::new("", VcfCountType::default(), false, None).unwrap();
        let id = parser.get_variant_id();
        assert!(id < parser.get_variant_id());
    }

    #[test]
    fn test_generate_path_names_for_samples() {
        let parser = VcfParser::new("", VcfCountType::default(), false, None).unwrap();
        let samples = vec!["1".to_string(), "2".to_string(), "3".to_string()];
        let feature = Vec::new();
        let path_names = parser.generate_path_names(&samples, &feature);
        assert_eq!(samples, path_names);
    }

    #[test]
    fn test_generate_path_names_for_haplotypes() {
        let parser = VcfParser::new("", VcfCountType::default(), true, None).unwrap();
        let samples = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let feature = vec![vec![1], vec![1, 1, 1], vec![1, 1]];
        let path_names = parser.generate_path_names(&samples, &feature);
        let haplotypes = vec![
            "a#0".to_string(),
            "b#1".to_string(),
            "b#2".to_string(),
            "b#3".to_string(),
            "c#1".to_string(),
            "c#2".to_string(),
        ];
        assert_eq!(haplotypes, path_names);
    }

    #[test]
    fn test_generate_hist_samples_allele_length() {
        let parser = Box::new(VcfParser::new("", VcfCountType::AltLength, false, None).unwrap());
        let hist = parser.generate_hist_from_reader(VCF_TEST.as_bytes());
        assert_eq!(hist.get_number_of_features(), 449);
        assert_eq!(hist.get_maximum_coverage(), 232);
    }

    #[test]
    fn test_generate_hist_samples_variants() {
        let parser = Box::new(VcfParser::new("", VcfCountType::Variants, false, None).unwrap());
        let hist = parser.generate_hist_from_reader(VCF_TEST.as_bytes());
        assert_eq!(hist.get_number_of_features(), 13);
        assert_eq!(hist.get_maximum_coverage(), 232);
        let mut hist_values = vec![0; 232 + 1];
        hist_values[1] = 13;
        assert_eq!(hist.get_hist_values(), &hist_values);
    }

    #[test]
    fn test_generate_matrix_samples_variants() {
        let parser = Box::new(VcfParser::new("", VcfCountType::Variants, false, None).unwrap());
        let matrix = parser.generate_matrix_from_reader(VCF_TEST.as_bytes());
        assert_eq!(matrix.get_feature_count(), 13);
        assert_eq!(matrix.get_path_names().len(), 232);
    }

    #[test]
    #[should_panic]
    fn test_hist_panic_on_emtpy_open() {
        let parser = Box::new(VcfParser::new("", VcfCountType::AltLength, false, None).unwrap());
        parser.generate_hist();
    }

    #[test]
    #[should_panic]
    fn test_matrix_panic_on_emtpy_open() {
        let parser = Box::new(VcfParser::new("", VcfCountType::AltLength, false, None).unwrap());
        parser.generate_matrix();
    }
}
