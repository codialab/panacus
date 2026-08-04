/* standard use */
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;

use serde::{Deserialize, Serialize};
/* external use */
use strum_macros::{EnumIter, EnumString, EnumVariantNames};

use crate::file_formats::gfa_parser::ItemId;

/* internal use */

// storage space for item IDs
pub type ItemIdSize = u64;
pub type CountSize = u32;
pub type GroupSize = u64;

pub struct Wrap<T>(pub *mut T);
unsafe impl Sync for Wrap<Vec<usize>> {}
unsafe impl Sync for Wrap<Vec<u64>> {}
unsafe impl Sync for Wrap<Vec<u32>> {}
unsafe impl Sync for Wrap<Vec<u16>> {}
unsafe impl Sync for Wrap<Vec<Vec<u32>>> {}
unsafe impl Sync for Wrap<Vec<Vec<u64>>> {}
// unsafe impl Sync for Wrap<[HashMap<u64, InfixEqStorage>; SIZE_T]> {}

#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    EnumString,
    EnumVariantNames,
    EnumIter,
    Hash,
    Eq,
    PartialOrd,
    Ord,
    Serialize,
    Deserialize,
)]
#[strum(serialize_all = "lowercase")]
pub enum CountType {
    Node,
    Bp,
    Edge,
}

impl Default for CountType {
    fn default() -> Self {
        CountType::Node
    }
}

impl fmt::Display for CountType {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        write!(
            formatter,
            "{}",
            match self {
                CountType::Node => "node",
                CountType::Edge => "edge",
                CountType::Bp => "bp",
            }
        )
    }
}

pub fn get_default_plot_downloads() -> Vec<(String, String)> {
    vec![
        ("png".to_string(), "Download as png".to_string()),
        ("svg".to_string(), "Download as svg".to_string()),
        ("vega-editor".to_string(), "Open in vega editor".to_string()),
    ]
}

#[derive(Debug, Clone)]
pub struct ItemTable {
    pub items: Vec<ItemIdSize>, // For each path, each feature in one long list
    pub id_prefsum: Vec<ItemIdSize>, // In the items list, where (at what index) does each new path start (+ length)
}

impl ItemTable {
    pub fn new(num_walks_paths: usize) -> Self {
        Self {
            items: vec![],
            id_prefsum: vec![0; num_walks_paths + 1],
        }
    }
}

#[derive(Debug, Clone)]
pub struct ActiveTable {
    pub items: Vec<bool>,
    // intervall container + item len vector
    annotation: Option<IntervalContainer>,
}

impl ActiveTable {
    // if you provide item_length, then it an active table with annotation
    pub fn new(size: usize, with_annotation: bool) -> Self {
        Self {
            items: vec![false; size],
            annotation: if with_annotation {
                Some(IntervalContainer::new())
            } else {
                None
            },
        }
    }

    pub fn activate(&mut self, id: &ItemId) {
        self.items[id.0 as usize] |= true;
    }

    #[allow(dead_code)]
    pub fn is_active(&self, id: &ItemId) -> bool {
        self.items[id.0 as usize]
    }

    pub fn activate_n_annotate(
        &mut self,
        id: ItemId,
        item_len: usize,
        start: usize,
        end: usize,
    ) -> Result<(), ActiveTableError> {
        match &mut self.annotation {
            None => Err(ActiveTableError::NoAnnotation),
            Some(m) => {
                // if interval completely covers item, remove it from map
                if end - start == item_len {
                    self.items[id.0 as usize] |= true;
                    m.remove(&id);
                } else {
                    if start > end {
                        log::error!(
                            "start ({}) is larger than end ({}) for node {}",
                            start,
                            end,
                            id
                        );
                    } else {
                        m.add(id, start, end);
                    }
                    if m.get(&id).unwrap()[0] == (0, item_len) {
                        m.remove(&id);
                        self.items[id.0 as usize] |= true;
                    }
                }
                Ok(())
            }
        }
    }

    pub fn get_active_intervals(&self, id: &ItemId, item_len: usize) -> Vec<(usize, usize)> {
        if self.items[id.0 as usize] {
            vec![(0, item_len)]
        } else if let Some(container) = &self.annotation {
            match container.get(id) {
                None => Vec::new(),
                Some(v) => v.to_vec(),
            }
        } else {
            Vec::new()
        }
    }

    pub fn with_annotation(&self) -> bool {
        self.annotation.is_some()
    }

    pub fn get_annotation_keys(&self) -> impl Iterator<Item = &ItemId> + '_ {
        self.annotation
            .as_ref()
            .map(|x| x.keys())
            .into_iter()
            .flatten()
    }
}

#[derive(Debug, Clone)]
pub struct IntervalContainer {
    map: HashMap<ItemId, Vec<(usize, usize)>>,
}

impl IntervalContainer {
    pub fn new() -> Self {
        IntervalContainer {
            map: HashMap::default(),
        }
    }

    pub fn add(&mut self, id: ItemId, start: usize, end: usize) {
        log::debug!("add {}:{}-{} to interval container", id, start, end);
        // produce union of intervals
        self.map
            .entry(id)
            .and_modify(|x| {
                let i = x
                    .binary_search_by_key(&start, |&(y, _)| y)
                    .unwrap_or_else(|z| z);
                log::debug!(
                    "binary search of ({}, {}) in {:?} returned {}",
                    start,
                    end,
                    x,
                    i
                );
                if i > 0 && x[i - 1].1 >= start {
                    if x[i - 1].1 < end {
                        let mut stop = end;
                        while i < x.len() && x[i].0 <= end {
                            stop = std::cmp::max(stop, x[i].1);
                            x.remove(i);
                        }
                        x[i - 1].1 = stop;
                    }
                    // else do nothing, because the new interval is fully enclosed in the previous
                    // interval
                } else if i < x.len() && x[i].1 >= start && x[i].0 <= end {
                    x[i].0 = std::cmp::min(x[i].0, start);
                    let mut stop = std::cmp::max(x[i].1, end);
                    while i + 1 < x.len() && x[i + 1].0 <= end {
                        stop = std::cmp::max(stop, x[i + 1].1);
                        x.remove(i + 1);
                    }
                    x[i].1 = stop;
                } else {
                    x.insert(i, (start, end));
                }
            })
            .or_insert(vec![(start, end)]);
    }

    pub fn get(&self, id: &ItemId) -> Option<&[(usize, usize)]> {
        self.map.get(id).map(|x| &x[..])
    }

    pub fn contains(&self, id: &ItemId) -> bool {
        self.map.contains_key(id)
    }

    pub fn remove(&mut self, id: &ItemId) -> Option<Vec<(usize, usize)>> {
        self.map.remove(id)
    }

    pub fn total_coverage(&self, id: &ItemId, exclude: &Option<Vec<(usize, usize)>>) -> usize {
        self.map
            .get(id)
            .as_ref()
            .map(|v| match exclude {
                None => v.iter().fold(0, |x, (a, b)| x + b - a),
                Some(ex) => {
                    let mut res = 0;
                    let mut i = 0;
                    for (start, end) in v.iter() {
                        // intervals have exclusive right bound, so "<=" is the right choice here
                        while i < ex.len() && &ex[i].1 <= start {
                            i += 1;
                        }
                        if i < ex.len() && &ex[i].0 < end {
                            // interval that starts with node start and ends with exclude start or
                            // node end, whichever comes first
                            //
                            // mind the (include, exclude] character of intervals!
                            res += usize::min(ex[i].0 - 1, *end) - start;

                            // interval that starts with exclude end and ends with node end
                            //
                            // mind the [include, exclude) character of intervals!
                            if &ex[i].1 < end {
                                res += end - ex[i].1 + 1;
                            }
                        } else {
                            res += end - start;
                        }
                    }
                    res
                }
            })
            .unwrap_or(0)
    }

    #[allow(dead_code)]
    pub fn iter(&self) -> impl Iterator<Item = (&ItemId, &Vec<(usize, usize)>)> + '_ {
        self.map.iter()
    }

    pub fn keys(&self) -> impl Iterator<Item = &ItemId> + '_ {
        self.map.keys()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ActiveTableError {
    NoAnnotation,
}

impl std::error::Error for ActiveTableError {}

impl fmt::Display for ActiveTableError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            ActiveTableError::NoAnnotation => write!(f, "Active Table has no annotations"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Threshold {
    Relative(f64),
    Absolute(usize),
}

impl fmt::Display for Threshold {
    fn fmt(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Threshold::Relative(c) => write!(formatter, "{}R", c)?,
            Threshold::Absolute(c) => write!(formatter, "{}A", c)?,
        }
        Ok(())
    }
}

impl Threshold {
    pub fn get_string(&self) -> String {
        match self {
            Threshold::Relative(c) => format!("{}", c),
            Threshold::Absolute(c) => format!("{}", c),
        }
    }

    pub fn to_absolute(self, n: usize) -> usize {
        match self {
            Threshold::Absolute(c) => c,
            Threshold::Relative(c) => (n as f64 * c).ceil() as usize,
        }
    }

    pub fn to_relative(self, n: usize) -> f64 {
        match self {
            Threshold::Relative(c) => c,
            Threshold::Absolute(c) => c as f64 / n as f64,
        }
    }
}

//
// helper functions
//

pub fn intersects(v: &[(usize, usize)], el: &(usize, usize)) -> bool {
    // this code assumes that intervals of v are (i) sorted (ii) non-overlapping

    v.binary_search_by(|(s, e)| {
        if s <= &el.1 && e >= &el.0 {
            Ordering::Equal
        } else if e < &el.0 {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    })
    .is_ok()
}

pub fn is_contained(v: &[(usize, usize)], el: &(usize, usize)) -> bool {
    // this code assumes that intervals of v are (i) sorted (ii) non-overlapping

    v.binary_search_by(|(s, e)| {
        if s <= &el.0 && e >= &el.1 {
            Ordering::Equal
        } else if e <= &el.1 {
            Ordering::Less
        } else {
            Ordering::Greater
        }
    })
    .is_ok()
}

pub fn averageu32(v: &[u32]) -> f32 {
    if v.is_empty() {
        return 0.0;
    }
    (v.iter().map(|x| *x as u64).sum::<u64>() as f64 / v.len() as f64) as f32
}

pub fn averagef64(v: &[f64]) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.iter().copied().sum::<f64>() / v.len() as f64
}

//pub fn averageu64 (v: &[u64]) -> f64 {
//    v.iter().sum::<u64>() as f64 / v.len() as f64
//}

pub fn median_already_sorted(v: &[u32]) -> f64 {
    //v.sort(); this has been done before
    let n = v.len();
    let mid = n / 2;
    if n % 2 == 1 {
        v[mid] as f64
    } else {
        (v[mid - 1] as f64 + v[mid] as f64) / 2.0
    }
}

pub fn n50_already_sorted(v: &[u32]) -> Option<u32> {
    //v.sort(); this has been done before
    let total_length: u32 = v.iter().sum();

    let mut running_sum = 0;
    for &len in v.iter() {
        running_sum += len;
        if running_sum * 2 >= total_length {
            return Some(len);
        }
    }

    None
}

#[allow(dead_code)]
pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    dna.iter()
        .rev() // Reverse the sequence
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't', // Handle lowercase
            b't' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            _ => panic!("Invalid nucleotide: {}", b as char),
        })
        .collect()
}

#[allow(dead_code)]
pub fn bits2kmer(kmer_bits: u64, k: usize) -> String {
    let nucleotides = ['A', 'C', 'G', 'T'];
    let mut kmer_str = String::with_capacity(k);

    for i in 0..k {
        let index = ((kmer_bits >> (2 * (k - i - 1))) & 3) as usize;
        kmer_str.push(nucleotides[index]);
    }
    kmer_str
}

const NUCLEOTIDE_BITS: [u8; 256] = {
    let mut map = [4; 256];
    map[b'A' as usize] = 0;
    map[b'C' as usize] = 1;
    map[b'G' as usize] = 2;
    map[b'T' as usize] = 3;
    map[b'a' as usize] = 0;
    map[b'c' as usize] = 1;
    map[b'g' as usize] = 2;
    map[b't' as usize] = 3;
    map
};

//let kmer = b"ACGTacgt";
//let result = kmer_u8_to_u64(kmer);
pub fn kmer_u8_to_u64(kmer: &[u8]) -> u64 {
    let mut result: u64 = 0;
    for &nucleotide in kmer {
        let bits = NUCLEOTIDE_BITS[nucleotide as usize];
        if bits < 4 {
            result = (result << 2) | bits as u64;
        } else {
            panic!("Invalid nucleotide: {}", nucleotide as char);
        }
    }
    result
}

const LOOKUP_RC: [u64; 256] = [
    0xff, 0xbf, 0x7f, 0x3f, 0xef, 0xaf, 0x6f, 0x2f, 0xdf, 0x9f, 0x5f, 0x1f, 0xcf, 0x8f, 0x4f, 0x0f,
    0xfb, 0xbb, 0x7b, 0x3b, 0xeb, 0xab, 0x6b, 0x2b, 0xdb, 0x9b, 0x5b, 0x1b, 0xcb, 0x8b, 0x4b, 0x0b,
    0xf7, 0xb7, 0x77, 0x37, 0xe7, 0xa7, 0x67, 0x27, 0xd7, 0x97, 0x57, 0x17, 0xc7, 0x87, 0x47, 0x07,
    0xf3, 0xb3, 0x73, 0x33, 0xe3, 0xa3, 0x63, 0x23, 0xd3, 0x93, 0x53, 0x13, 0xc3, 0x83, 0x43, 0x03,
    0xfe, 0xbe, 0x7e, 0x3e, 0xee, 0xae, 0x6e, 0x2e, 0xde, 0x9e, 0x5e, 0x1e, 0xce, 0x8e, 0x4e, 0x0e,
    0xfa, 0xba, 0x7a, 0x3a, 0xea, 0xaa, 0x6a, 0x2a, 0xda, 0x9a, 0x5a, 0x1a, 0xca, 0x8a, 0x4a, 0x0a,
    0xf6, 0xb6, 0x76, 0x36, 0xe6, 0xa6, 0x66, 0x26, 0xd6, 0x96, 0x56, 0x16, 0xc6, 0x86, 0x46, 0x06,
    0xf2, 0xb2, 0x72, 0x32, 0xe2, 0xa2, 0x62, 0x22, 0xd2, 0x92, 0x52, 0x12, 0xc2, 0x82, 0x42, 0x02,
    0xfd, 0xbd, 0x7d, 0x3d, 0xed, 0xad, 0x6d, 0x2d, 0xdd, 0x9d, 0x5d, 0x1d, 0xcd, 0x8d, 0x4d, 0x0d,
    0xf9, 0xb9, 0x79, 0x39, 0xe9, 0xa9, 0x69, 0x29, 0xd9, 0x99, 0x59, 0x19, 0xc9, 0x89, 0x49, 0x09,
    0xf5, 0xb5, 0x75, 0x35, 0xe5, 0xa5, 0x65, 0x25, 0xd5, 0x95, 0x55, 0x15, 0xc5, 0x85, 0x45, 0x05,
    0xf1, 0xb1, 0x71, 0x31, 0xe1, 0xa1, 0x61, 0x21, 0xd1, 0x91, 0x51, 0x11, 0xc1, 0x81, 0x41, 0x01,
    0xfc, 0xbc, 0x7c, 0x3c, 0xec, 0xac, 0x6c, 0x2c, 0xdc, 0x9c, 0x5c, 0x1c, 0xcc, 0x8c, 0x4c, 0x0c,
    0xf8, 0xb8, 0x78, 0x38, 0xe8, 0xa8, 0x68, 0x28, 0xd8, 0x98, 0x58, 0x18, 0xc8, 0x88, 0x48, 0x08,
    0xf4, 0xb4, 0x74, 0x34, 0xe4, 0xa4, 0x64, 0x24, 0xd4, 0x94, 0x54, 0x14, 0xc4, 0x84, 0x44, 0x04,
    0xf0, 0xb0, 0x70, 0x30, 0xe0, 0xa0, 0x60, 0x20, 0xd0, 0x90, 0x50, 0x10, 0xc0, 0x80, 0x40, 0x00,
];

pub fn revcmp(kmer: u64, k: usize) -> u64 {
    (LOOKUP_RC[(kmer & 0xff) as usize] << 56
        | LOOKUP_RC[((kmer >> 8) & 0xff) as usize] << 48
        | LOOKUP_RC[((kmer >> 16) & 0xff) as usize] << 40
        | LOOKUP_RC[((kmer >> 24) & 0xff) as usize] << 32
        | LOOKUP_RC[((kmer >> 32) & 0xff) as usize] << 24
        | LOOKUP_RC[((kmer >> 40) & 0xff) as usize] << 16
        | LOOKUP_RC[((kmer >> 48) & 0xff) as usize] << 8
        | LOOKUP_RC[((kmer >> 56) & 0xff) as usize])
        >> (64 - k as u64 * 2)
}

#[allow(dead_code)]
pub fn canonical(kmer_bits: u64, k: usize) -> u64 {
    let kmer_bits_rc = revcmp(kmer_bits, k);
    if kmer_bits < kmer_bits_rc {
        kmer_bits
    } else {
        kmer_bits_rc
    }
}

pub fn to_id(s: &str) -> String {
    s.to_string()
        .to_lowercase()
        .replace([' ', '|', '/', '\\', '\'', '"'], "-")
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::file_formats::gfa_parser::ItemId;

    #[test]
    fn test_interval_container() {
        let mut ic = IntervalContainer::new();
        ic.add(ItemId(0), 5, 6);
        ic.add(ItemId(0), 9, 10);
        ic.add(ItemId(0), 7, 8);
        assert_eq!(ic.map.get(&ItemId(0)), Some(&vec![(5, 6), (7, 8), (9, 10)]));
        ic.add(ItemId(0), 4, 5);
        assert_eq!(ic.map.get(&ItemId(0)), Some(&vec![(4, 6), (7, 8), (9, 10)]));
        ic.add(ItemId(0), 0, 11);
        assert_eq!(ic.map.get(&ItemId(0)), Some(&vec![(0, 11)]));
        ic.add(ItemId(0), 11, 12);
        assert_eq!(ic.map.get(&ItemId(0)), Some(&vec![(0, 12)]));
        ic.add(ItemId(0), 13, 15);
        ic.add(ItemId(0), 16, 20);
        assert_eq!(
            ic.map.get(&ItemId(0)),
            Some(&vec![(0, 12), (13, 15), (16, 20)])
        );
        ic.add(ItemId(0), 14, 17);
        assert_eq!(ic.map.get(&ItemId(0)), Some(&vec![(0, 12), (13, 20)]));
    }
}
