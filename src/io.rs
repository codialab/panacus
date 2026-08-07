/* standard use */
use std::io::{BufRead, BufReader, Read};
use std::io::{Error, ErrorKind};
use std::str::{self, FromStr};

/* external use */
use flate2::read::MultiGzDecoder;
use quick_csv::Csv;

/* internal use */
use crate::file_formats::gfa_parser::PathSegment;
use crate::util::*;

pub fn bufreader_from_compressed_gfa(gfa_file: &str) -> BufReader<Box<dyn Read + Send>> {
    log::info!("loading graph from {}", &gfa_file);
    if let Ok(f) = std::fs::File::open(gfa_file) {
        let reader: Box<dyn Read + Send> = if gfa_file.ends_with(".gz") {
            log::info!("assuming that {} is gzip compressed..", &gfa_file);
            Box::new(MultiGzDecoder::new(f))
        } else {
            Box::new(f)
        };
        BufReader::new(reader)
    } else {
        log::error!(
            "Could not read file {}, please make sure it exists and is readable!",
            gfa_file
        );
        panic!("Could not read input file!")
    }
}

pub fn parse_bed_to_path_segments<R: Read>(
    data: &mut BufReader<R>,
    use_block_info: bool,
) -> Vec<PathSegment> {
    // based on https://en.wikipedia.org/wiki/BED_(file_format)
    let mut segments = Vec::new();

    for (i, line) in data.lines().enumerate() {
        let line = match line {
            Ok(l) => l,
            Err(e) => {
                panic!("error reading line {}: {}", i + 1, e);
            }
        };

        let fields = {
            let mut fields: Vec<&str> = line.split('\t').collect();
            if fields.is_empty() {
                fields = vec![&line];
            }
            fields
        };
        let path_name = fields[0];

        if path_name.starts_with("browser ")
            || path_name.starts_with("track ")
            || path_name.starts_with('#')
        {
            continue;
        }

        if fields.len() == 1 {
            segments.push(PathSegment::create_from_str(path_name));
        } else if fields.len() >= 3 {
            let start = usize::from_str(fields[1]).unwrap_or_else(|_| {
                panic!("error line {}: `{}` is not an usize", i + 1, fields[1])
            });
            let end = usize::from_str(fields[2]).unwrap_or_else(|_| {
                panic!("error line {}: `{}` is not an usize", i + 1, fields[2])
            });

            if use_block_info && fields.len() == 12 {
                let block_count = fields[9].parse::<usize>().unwrap_or(0);
                let block_sizes: Vec<usize> = fields[10]
                    .split(',')
                    .filter_map(|s| usize::from_str(s.trim()).ok())
                    .collect();
                let block_starts: Vec<usize> = fields[11]
                    .split(',')
                    .filter_map(|s| usize::from_str(s.trim()).ok())
                    .collect();

                if block_count == block_sizes.len() && block_count == block_starts.len() {
                    for (size, start_offset) in block_sizes.iter().zip(block_starts.iter()) {
                        let block_start = start + start_offset;
                        let block_end = block_start + size;
                        segments.push(PathSegment::from_str_start_end(
                            path_name,
                            block_start,
                            block_end,
                        ));
                    }
                } else {
                    panic!(
                        "error in block sizes/starts in line {}: counts do not match",
                        i + 1
                    );
                }
            } else {
                segments.push(PathSegment::from_str_start_end(path_name, start, end));
            }
        } else {
            panic!(
                "error in line {}: row must have either 1, 3, or 12 columns, but has 2",
                i + 1
            );
        }
    }

    segments
}

pub fn parse_groups<R: Read>(data: &mut BufReader<R>) -> Result<Vec<(PathSegment, String)>, Error> {
    let mut res: Vec<(PathSegment, String)> = Vec::new();

    let mut i = 1;
    let mut buf = vec![];
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        //Remove new line at the end
        if let Some(&last_byte) = buf.last() {
            if last_byte == b'\n' || last_byte == b'\r' {
                buf.pop();
            }
        }
        let line = String::from_utf8(buf.clone())
            .unwrap_or_else(|_| panic!("error in line {}: some character is not UTF-8", i));
        let columns: Vec<&str> = line.split('\t').collect();

        if columns.len() != 2 {
            let msg = format!("error in line {}: table must have exactly two columns", i);
            log::error!("{}", &msg);
            return Err(Error::new(ErrorKind::InvalidData, msg));
        }

        let path_seg = PathSegment::create_from_str(columns[0]);
        res.push((path_seg, columns[1].to_string()));

        i += 1;
        buf.clear();
    }

    Ok(res)
}

#[allow(dead_code)]
pub fn parse_threshold_file<R: Read>(data: &mut BufReader<R>) -> Result<Vec<Threshold>, Error> {
    let mut res = Vec::new();

    let reader = Csv::from_reader(data)
        .delimiter(b'\t')
        .flexible(true)
        .has_header(false);
    for (i, row) in reader.enumerate() {
        let row = row.unwrap();
        let mut row_it = row.bytes_columns();
        if let Some(col) = row_it.next() {
            let threshold_str = str::from_utf8(col).unwrap();
            if let Ok(t) = usize::from_str(threshold_str) {
                res.push(Threshold::Absolute(t));
            } else if let Ok(t) = f64::from_str(threshold_str) {
                res.push(Threshold::Relative(t));
            } else {
                return Err(Error::new(
                    ErrorKind::InvalidData,
                    &format!(
                        "threshold \"{}\" (line {}) is neither an integer nor a float",
                        &threshold_str,
                        i + 1
                    )[..],
                ));
            }
        }
    }

    Ok(res)
}

pub fn write_table(headers: &[Vec<String>], columns: &[Vec<f64>]) -> Result<String, Error> {
    write_table_with_start_index(headers, columns, 0)
}

pub fn write_table_with_start_index(
    headers: &[Vec<String>],
    columns: &[Vec<f64>],
    start_index: usize,
) -> Result<String, Error> {
    let n = headers.first().unwrap_or(&Vec::new()).len();
    let mut res = String::new();
    for i in 0..n {
        for (j, value) in headers.iter().enumerate() {
            if j > 0 {
                res.push('\t');
            }
            res.push_str(&format!("{:0}", value[i]));
        }
        res.push('\n');
    }
    let n = columns.first().unwrap_or(&Vec::new()).len();
    for i in 0..n {
        res.push_str(&(i + start_index).to_string());
        for column in columns {
            res.push_str(&format!("\t{:0}", column[i].floor()));
        }
        res.push('\n');
    }
    Ok(res)
}

pub fn write_ordered_table(
    headers: &[Vec<String>],
    columns: &Vec<Vec<f64>>,
    index: &[String],
) -> anyhow::Result<String> {
    let n = headers.first().unwrap_or(&Vec::new()).len();
    let mut res = String::new();

    for i in 0..n {
        for (j, value) in headers.iter().enumerate() {
            if j > 0 {
                res.push('\t');
            }
            res.push_str(&format!("{:0}", value[i]));
        }
        res.push('\n');
    }
    let n = columns.first().unwrap_or(&Vec::new()).len();
    for i in 0..n {
        res.push_str(&index[i].to_string());
        for column in columns {
            res.push_str(&format!("\t{:0}", column[i].floor()));
        }
        res.push('\n');
    }

    Ok(res)
}

pub fn write_metadata_comments() -> anyhow::Result<String> {
    let mut res = format!(
        "# {}\n",
        std::env::args().collect::<Vec<String>>().join(" ")
    );
    let version = option_env!("GIT_HASH").unwrap_or(env!("CARGO_PKG_VERSION"));
    let version = format!("# version {}\n", version);
    res.push_str(&version);
    Ok(res)
}
