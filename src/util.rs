/* standard use */
use std::cmp::Ordering;
use std::collections::HashMap;
use std::fmt;

/* external use */
use strum_macros::{EnumString, EnumVariantNames, EnumIter};

/* internal use */
use crate::graph::ItemId;

//
// storage space for item IDs
//
pub type ItemIdSize = u32;
pub type CountSize = u32;
pub type GroupSize = u16;

pub const SIZE_T: usize = 1024;
pub struct Wrap<T>(pub *mut T);
unsafe impl Sync for Wrap<Vec<usize>> {}
unsafe impl Sync for Wrap<Vec<u64>> {}
unsafe impl Sync for Wrap<Vec<u32>> {}
unsafe impl Sync for Wrap<Vec<u16>> {}
unsafe impl Sync for Wrap<[Vec<u32>; SIZE_T]> {}
unsafe impl Sync for Wrap<Vec<Vec<u32>>> {}
unsafe impl Sync for Wrap<[Vec<u64>; SIZE_T]> {}
unsafe impl Sync for Wrap<Vec<Vec<u64>>> {}

#[derive(Debug, Clone, Copy, PartialEq, EnumString, EnumVariantNames, EnumIter)]
#[strum(serialize_all = "lowercase")]
pub enum CountType {
    Node,
    Bp,
    Edge,
    All,
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
                CountType::All => "all",
            }
        )
    }
}

pub struct ItemTable {
    pub items: [Vec<ItemIdSize>; SIZE_T],
    pub id_prefsum: [Vec<ItemIdSize>; SIZE_T],
}

impl ItemTable {
    pub fn new(num_walks_paths: usize) -> Self {
        Self {
            items: [(); SIZE_T].map(|_| vec![]),
            id_prefsum: [(); SIZE_T].map(|_| vec![0; num_walks_paths + 1]),
        }
    }
}

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
        // produce union of intervals
        self.map
            .entry(id)
            .and_modify(|x| {
                let i = x
                    .binary_search_by_key(&start, |&(y, _)| y)
                    .unwrap_or_else(|z| z);
                if i > 0 && x[i - 1].1 >= start && x[i - 1].1 <= end {
                    x[i - 1].1 = end;
                } else if i < x.len() && x[i].1 >= start && x[i].1 < end {
                    x[i].1 = end;
                } else if i < x.len() && x[i].0 <= end {
                    x[i].0 = start;
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
    pub fn to_string(&self) -> String {
        match self {
            Threshold::Relative(c) => format!("{}", c),
            Threshold::Absolute(c) => format!("{}", c),
        }
    }

    pub fn to_absolute(&self, n: usize) -> usize {
        match self {
            Threshold::Absolute(c) => *c,
            Threshold::Relative(c) => (n as f64 * c).ceil() as usize,
        }
    }

    pub fn to_relative(&self, n: usize) -> f64 {
        match self {
            Threshold::Relative(c) => *c,
            Threshold::Absolute(c) => *c as f64 / n as f64,
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

pub fn average (v: &[u32]) -> f32 {
    v.iter().sum::<u32>() as f32 / v.len() as f32
}

pub fn median_already_sorted(v: &[u32]) -> f32 {
    //v.sort(); this has been done before
    let n = v.len();
    let mid = n / 2;
    if n % 2 == 1 {
        v[mid] as f32
    } else {
        (v[mid - 1] as f32 + v[mid] as f32) / 2.0
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

//const NUCLEOTIDE_BITS: [u8; 256] = {
//    let mut map = [4; 256];
//    map[b'A' as usize] = 0;
//    map[b'C' as usize] = 1;
//    map[b'G' as usize] = 2;
//    map[b'T' as usize] = 3;
//    map[b'a' as usize] = 0;
//    map[b'c' as usize] = 1;
//    map[b'g' as usize] = 2;
//    map[b't' as usize] = 3;
//    map
//}

//pub fn log2_add(a: f64, b: f64) -> f64 {
//    // we assume both a and b are log2'd
//    let (a, b) = if a < b { (a, b) } else { (b, a) };
//
//    b + (1.0 + (a - b).exp2()).log2()
//}
