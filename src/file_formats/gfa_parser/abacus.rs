/* standard use */
use std::fs;
use std::io::BufReader;
use std::io::{Error, ErrorKind};
use std::iter::FromIterator;
use std::path::Path;
//use std::sync::{Arc, Mutex};

/* external crate*/
use regex::Regex;
use std::collections::{HashMap, HashSet};

use crate::file_formats::gfa_parser::grammar::Grammar;
/* private use */
use crate::io::*;
use crate::util::*;

use super::graph::{GraphStorage, PathSegment};
use super::util::parse_gfa_paths_walks;

#[derive(Debug, PartialEq, Eq, Clone, Default)]
pub struct GraphMaskParameters {
    pub positive_list: String,
    pub negative_list: String,
    pub groupby: String,
    pub groupby_sample: bool,
    pub groupby_haplotype: bool,
}

#[derive(Debug, Clone)]
pub struct GraphMask {
    pub groups: HashMap<PathSegment, String>,
    pub include_coords: Option<Vec<PathSegment>>,
    pub exclude_coords: Option<Vec<PathSegment>>,
}

impl GraphMask {
    pub fn from_datamgr(
        params: &GraphMaskParameters,
        graph_storage: &GraphStorage,
    ) -> Result<Self, Error> {
        let groups = GraphMask::load_groups(
            &params.groupby,
            params.groupby_haplotype,
            params.groupby_sample,
            graph_storage,
        )?;
        let paths = &graph_storage.path_segments;
        let include_coords = GraphMask::complement_with_group_assignments(
            GraphMask::load_coord_list(&params.positive_list, paths)?,
            &groups,
        )?;
        let exclude_coords = GraphMask::complement_with_group_assignments(
            GraphMask::load_coord_list(&params.negative_list, paths)?,
            &groups,
        )?;

        Ok(GraphMask {
            groups,
            include_coords,
            exclude_coords,
        })
    }

    pub fn complement_with_group_assignments(
        coords: Option<Vec<PathSegment>>,
        groups: &HashMap<PathSegment, String>,
    ) -> Result<Option<Vec<PathSegment>>, Error> {
        //
        // We allow coords to be defined via groups; the following code
        // 1. complements coords with path segments from group assignments
        // 2. checks that group-based coordinates don't have start/stop information
        //
        let mut group2paths: HashMap<String, Vec<PathSegment>> = HashMap::default();
        for (p, g) in groups.iter() {
            group2paths.entry(g.clone()).or_default().push(p.clone())
        }
        let path_to_group: HashMap<PathSegment, String> = groups
            .iter()
            .map(|(ps, g)| (ps.clear_coords(), g.clone()))
            .collect();

        match coords {
            None => Ok(None),
            Some(v) => {
                v.into_iter()
                    .map(|p| {
                        // check if path segment defined in coords associated with a specific path,
                        // it is not considered a group
                        if path_to_group.contains_key(&p.clear_coords()) {
                            Ok(vec![p])
                        } else if group2paths.contains_key(&p.id()) {
                            if p.coords().is_some() {
                                let msg = format!("invalid coordinate \"{}\": group identifiers are not allowed to have start/stop information!", &p);
                                log::error!("{}", &msg);
                                Err(Error::new( ErrorKind::InvalidData, msg))
                            } else {
                                let paths = group2paths.get(&p.id()).unwrap().clone();
                                log::debug!("complementing coordinate list with {} paths associted with group {}", paths.len(), p.id());
                                Ok(paths)
                            }
                        } else {
                            let msg = format!("unknown path/group {}", &p);
                            log::error!("{}", &msg);
                            // let's not be so harsh as to throw an error, ok?
                            // Err(Error::new(ErrorKind::InvalidData, msg))
                            Ok(Vec::new())
                        }
                    })
                    .collect::<Result<Vec<Vec<PathSegment>>, Error>>().map(|x| Some(x[..]
                    .concat()))
            }
        }
    }

    pub fn load_coord_list_file(file_name: &str) -> Result<Option<Vec<PathSegment>>, Error> {
        log::info!("loading coordinates from {}", file_name);
        let mut data = BufReader::new(fs::File::open(file_name)?);
        let use_block_info = true;
        let coords = parse_bed_to_path_segments(&mut data, use_block_info);
        log::debug!("loaded {} coordinates", coords.len());
        Ok(Some(coords))
    }

    pub fn load_coord_list(
        coord_text: &str,
        paths: &[PathSegment],
    ) -> Result<Option<Vec<PathSegment>>, Error> {
        Ok(if coord_text.is_empty() {
            None
        } else {
            if Path::new(coord_text).is_file() {
                Self::load_coord_list_file(coord_text)?
            } else if let Ok(re) = Regex::new(coord_text) {
                log::info!("filtering paths based on regex {}", coord_text);
                let coords = paths
                    .iter()
                    .filter(|path| re.is_match(&path.to_string()))
                    .cloned()
                    .collect::<Vec<_>>();
                if coords.is_empty() {
                    log::warn!("filtering with regex did not find any paths!");
                }
                Some(coords)
            } else {
                log::error!(
                    "string {} is not valid! Neither as a file name nor as a regex",
                    coord_text
                );
                panic!("Invalid file name");
            }
        })
    }

    fn load_groups(
        file_name: &str,
        groupby_haplotype: bool,
        groupby_sample: bool,
        graph_storage: &GraphStorage,
    ) -> Result<HashMap<PathSegment, String>, Error> {
        if groupby_haplotype {
            log::info!("grouping by haplotype");
            Ok(graph_storage
                .path_segments
                .iter()
                .map(|x| {
                    (
                        x.clear_coords(),
                        format!(
                            "{}#{}",
                            &x.sample,
                            &x.haplotype.as_ref().unwrap_or(&String::new())
                        ),
                    )
                })
                .collect())
        } else if groupby_sample {
            log::info!("grouping by sample");
            Ok(graph_storage
                .path_segments
                .iter()
                .map(|x| (x.clear_coords(), x.sample.clone()))
                .collect())
        } else if !file_name.is_empty() {
            log::info!("loading groups from {}", file_name);
            let mut data = BufReader::new(fs::File::open(file_name)?);
            let group_assignments = parse_groups(&mut data)?;
            let mut path_to_group = HashMap::default();
            for (i, (path, group)) in group_assignments.into_iter().enumerate() {
                let path: PathSegment = path;
                let path_nocoords: PathSegment = path.clear_coords();
                match path_to_group.get(&path_nocoords) {
                    Some(g) => {
                        if g != &group {
                            let msg = format!(
                                "error in line {}: path {} cannot be assigned to more than one group, but is assigned to at least two groups: {}, {}",
                                i, &path_nocoords, &g, &group
                            );
                            log::error!("{}", &msg);
                            return Err(Error::new(ErrorKind::InvalidData, msg));
                        }
                    }
                    None => {
                        path_to_group.insert(path_nocoords, group);
                    }
                }
            }
            log::debug!("loaded {} group assignments", path_to_group.len());

            // augment the group assignments with yet unassigned path segments
            graph_storage.path_segments.iter().for_each(|x| {
                let path = x.clear_coords();
                path_to_group.entry(path).or_insert_with(|| x.id());
            });
            Ok(path_to_group)
        } else {
            log::info!("no explicit grouping instruction given, group paths by their IDs (sample ID+haplotype ID+seq ID)");
            Ok(graph_storage
                .path_segments
                .iter()
                .map(|x| (x.clear_coords(), x.id()))
                .collect())
        }
    }

    pub fn get_path_order<'a>(
        &'a self,
        path_segments: &[PathSegment],
    ) -> Vec<(ItemIdSize, &'a str)> {
        // orders elements of path_segments by the order in graph_mask.include; the returned vector
        // maps indices of path_segments to the group identifier

        let mut group_to_paths: HashMap<&'a str, Vec<(ItemIdSize, &'a str)>> = HashMap::default();

        for (i, p) in path_segments.iter().enumerate() {
            let group: &'a str = self.groups.get(&p.clear_coords()).unwrap();
            group_to_paths
                .entry(group)
                .or_default()
                .push((i as ItemIdSize, group));
        }

        let order: Vec<&PathSegment> = if let Some(include) = &self.include_coords {
            include.iter().collect()
        } else {
            let exclude: HashSet<&PathSegment> = match &self.exclude_coords {
                Some(e) => e.iter().collect(),
                None => HashSet::new(),
            };
            path_segments
                .iter()
                .filter(|x| !exclude.contains(x))
                .collect::<Vec<&PathSegment>>()
        };
        order
            .into_iter()
            .map(|p| {
                group_to_paths
                    .remove(&self.groups.get(&p.clear_coords()).unwrap()[..])
                    .unwrap_or_default()
            })
            .collect::<Vec<Vec<(ItemIdSize, &'a str)>>>()
            .concat()
    }

    #[allow(dead_code)]
    pub fn count_groups(&self) -> usize {
        HashSet::<&String>::from_iter(self.groups.values()).len()
    }

    pub fn build_subpath_map(
        path_segments: &[PathSegment],
    ) -> HashMap<String, Vec<(usize, usize)>> {
        // intervals are 0-based, and [start, end), see https://en.wikipedia.org/wiki/BED_(file_format)
        let mut res: HashMap<String, HashSet<(usize, usize)>> = HashMap::default();

        path_segments.iter().for_each(|x| {
            res.entry(x.id()).or_default().insert(match x.coords() {
                None => (0, usize::MAX),
                Some((i, j)) => (i, j),
            });
        });

        HashMap::from_iter(res.into_iter().map(|(pid, coords)| {
            let mut v: Vec<(usize, usize)> = coords.into_iter().collect();
            v.sort();
            let mut i = 1;
            // remove overlaps
            while i < v.len() {
                if v[i - 1].1 >= v[i].0 {
                    let x = v.remove(i);
                    v[i - 1].1 = std::cmp::max(v[i - 1].1, x.1);
                } else {
                    i += 1
                }
            }
            (pid, v)
        }))
    }

    pub fn load_optional_subsetting(
        &self,
        graph_storage: &GraphStorage,
        count: &CountType,
    ) -> (
        Option<IntervalContainer>,
        Option<ActiveTable>,
        HashMap<String, Vec<(usize, usize)>>,
        HashMap<String, Vec<(usize, usize)>>,
    ) {
        // *only relevant for bps count in combination with subset option*
        // this table stores the number of bps of nodes that are *partially* uncovered by subset
        // coordinates
        let subset_covered_bps: Option<IntervalContainer> =
            if count == &CountType::Bp && self.include_coords.is_some() {
                Some(IntervalContainer::new())
            } else {
                None
            };

        // this table stores information about excluded nodes *if* the exclude setting is used
        let exclude_table = self.exclude_coords.as_ref().map(|_| {
            ActiveTable::new(
                graph_storage.number_of_items(count) + 1,
                count == &CountType::Bp,
            )
        });

        // build "include" lookup table
        let include_map = match &self.include_coords {
            None => HashMap::default(),
            Some(coords) => Self::build_subpath_map(coords),
        };

        // build "exclude" lookup table
        let exclude_map = match &self.exclude_coords {
            None => HashMap::default(),
            Some(coords) => Self::build_subpath_map(coords),
        };

        (subset_covered_bps, exclude_table, include_map, exclude_map)
    }

    pub fn load_optional_subsetting_multiple(
        &self,
        graph_storage: &GraphStorage,
        count_types: &[CountType],
    ) -> (
        Option<IntervalContainer>,
        Vec<Option<ActiveTable>>,
        HashMap<String, Vec<(usize, usize)>>,
        HashMap<String, Vec<(usize, usize)>>,
    ) {
        // *only relevant for bps count in combination with subset option*
        // this table stores the number of bps of nodes that are *partially* uncovered by subset
        // coordinates
        let subset_covered_bps: Option<IntervalContainer> =
            if count_types.contains(&CountType::Bp) && self.include_coords.is_some() {
                Some(IntervalContainer::new())
            } else {
                None
            };

        // this table stores information about excluded nodes *if* the exclude setting is used
        let exclude_tables: Vec<_> = count_types
            .iter()
            .map(|count| {
                self.exclude_coords.as_ref().map(|_| {
                    ActiveTable::new(
                        graph_storage.number_of_items(count) + 1,
                        count == &CountType::Bp,
                    )
                })
            })
            .collect();

        // build "include" lookup table
        let include_map = match &self.include_coords {
            None => HashMap::default(),
            Some(coords) => Self::build_subpath_map(coords),
        };

        // build "exclude" lookup table
        let exclude_map = match &self.exclude_coords {
            None => HashMap::default(),
            Some(coords) => Self::build_subpath_map(coords),
        };

        (subset_covered_bps, exclude_tables, include_map, exclude_map)
    }
}

#[derive(Debug, Clone)]
pub struct AbacusByTotal {
    pub count: CountType,                                  // Type of feature
    pub countable: Vec<CountSize>, // Each entry gives coverage of each feature (except for bp, works similar to nodes)
    pub uncovered_bps: Option<HashMap<ItemIdSize, usize>>, // Map, for each incomplete node the number of not covered bps
    pub groups: Vec<String>, // List of group names (len == max coverage)
}

impl AbacusByTotal {
    pub fn from_gfa<R: std::io::Read>(
        data: &mut BufReader<R>,
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        grammar: &Grammar,
        count_type: CountType,
    ) -> (Self, HashMap<PathSegment, (u32, u32)>) {
        let (item_table, exclude_table, subset_covered_bps, paths_len, _) = parse_gfa_paths_walks(
            data,
            graph_mask,
            graph_storage,
            grammar,
            &count_type,
            &Vec::new(),
        );
        let abacus = Self::item_table_to_abacus(
            graph_mask,
            graph_storage,
            count_type,
            item_table,
            exclude_table,
            subset_covered_bps,
        );
        (abacus, paths_len)
    }

    pub fn from_item_table(
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        count: CountType,
        item_table: ItemTable,
        exclude_table: Option<ActiveTable>,
        subset_covered_bps: Option<IntervalContainer>,
    ) -> Self {
        Self::item_table_to_abacus(
            graph_mask,
            graph_storage,
            count,
            item_table,
            exclude_table,
            subset_covered_bps,
        )
    }

    pub fn from_gfa_multiple<R: std::io::Read>(
        _data: &mut BufReader<R>,
        _graph_mask: &GraphMask,
        _graph_storage: &GraphStorage,
        _count_types: &Vec<CountType>,
    ) -> (Vec<AbacusByTotal>, HashMap<PathSegment, (u32, u32)>) {
        (Vec::new(), HashMap::new())
    }

    pub fn item_table_to_abacus(
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        count: CountType,
        item_table: ItemTable,
        exclude_table: Option<ActiveTable>,
        subset_covered_bps: Option<IntervalContainer>,
    ) -> Self {
        log::info!("counting abacus entries..");
        // first element in countable is "zero" element. It is ignored in counting
        let mut countable: Vec<CountSize> = vec![0; graph_storage.number_of_items(&count) + 1];
        // countable with ID "0" is special and should not be considered in coverage histogram
        countable[0] = CountSize::MAX;
        let mut last: Vec<ItemIdSize> =
            vec![ItemIdSize::MAX; graph_storage.number_of_items(&count) + 1];

        let mut groups = Vec::new();
        for (path_id, group_id) in graph_mask.get_path_order(&graph_storage.path_segments) {
            if groups.is_empty() || groups.last().unwrap() != group_id {
                groups.push(group_id.to_string());
            }
            AbacusByTotal::coverage(
                &mut countable,
                &mut last,
                &item_table,
                &exclude_table,
                path_id,
                groups.len() as ItemIdSize - 1,
            );
        }

        log::info!(
            "abacus has {} path groups and {} countables",
            groups.len(),
            countable.len() - 1
        );

        Self {
            count,
            countable,
            uncovered_bps: Some(quantify_uncovered_bps(
                &exclude_table,
                &subset_covered_bps,
                graph_storage,
            )),
            groups,
        }
    }

    pub fn coverage(
        countable: &mut Vec<CountSize>,
        last: &mut Vec<ItemIdSize>,
        item_table: &ItemTable,
        exclude_table: &Option<ActiveTable>,
        path_id: ItemIdSize,
        group_id: ItemIdSize,
    ) {
        let countable_ptr = Wrap(countable);
        let last_ptr = Wrap(last);
        // Parallel node counting
        let start = item_table.id_prefsum[path_id as usize] as usize;
        let end = item_table.id_prefsum[path_id as usize + 1] as usize;
        for j in start..end {
            let sid = item_table.items[j] as usize;
            unsafe {
                if last[sid] != group_id
                    && (exclude_table.is_none() || !exclude_table.as_ref().unwrap().items[sid])
                {
                    (&mut (*countable_ptr.0))[sid] += 1;
                    (&mut (*last_ptr.0))[sid] = group_id;
                }
            }
        }
    }

    pub fn construct_hist_from_set(&self, indices: &[usize]) -> Vec<usize> {
        log::debug!("constructing histogram..");
        let mut hist: Vec<usize> = vec![0; self.groups.len() + 1];

        for (i, cov) in indices.iter().map(|&idx| self.countable[idx]).enumerate() {
            if cov as usize >= hist.len() {
                if i != 0 {
                    log::warn!("coverage {} of item {} exceeds the number of groups {}, it'll be ignored in the count", cov, i, self.groups.len());
                }
            } else {
                hist[cov as usize] += 1;
            }
        }
        hist
    }

    pub fn construct_hist(&self) -> Vec<usize> {
        log::info!("constructing histogram..");
        // hist must be of size = num_groups + 1; having an index that starts
        // from 1, instead of 0, makes easier the calculation in hist2pangrowth.
        let mut hist: Vec<usize> = vec![0; self.groups.len() + 1];

        for (i, cov) in self.countable.iter().enumerate() {
            if *cov as usize >= hist.len() {
                if i != 0 {
                    log::warn!("coverage {} of item {} exceeds the number of groups {}, it'll be ignored in the count", cov, i, self.groups.len());
                }
            } else {
                hist[*cov as usize] += 1;
            }
        }
        hist
    }

    pub fn construct_hist_bps_of_subset(
        &self,
        graph_storage: &GraphStorage,
        indices: &[usize],
        uncovered_bps: HashMap<u64, usize>,
    ) -> Vec<usize> {
        log::debug!("constructing bp histogram..");
        let mut hist: Vec<usize> = vec![0; self.groups.len() + 1];
        for (id, cov) in indices.iter().map(|&idx| (idx, self.countable[idx])) {
            if cov as usize >= hist.len() {
                if id != 0 {
                    log::info!("coverage {} of item {} exceeds the number of groups {}, it'll be ignored in the count", cov, id, self.groups.len());
                }
            } else {
                hist[cov as usize] += graph_storage.node_lens[id] as usize;
            }
        }
        for (id, uncov) in uncovered_bps.iter() {
            hist[self.countable[*id as usize] as usize] -= uncov;
            // add uncovered bps to 0-coverage count
            // hist[0] += uncov;
        }
        hist
    }

    pub fn construct_hist_bps(&self, graph_storage: &GraphStorage) -> Vec<usize> {
        log::info!("constructing bp histogram..");
        // hist must be of size = num_groups + 1; having an index that starts
        // from 1, instead of 0, makes easier the calculation in hist2pangrowth.
        let mut hist: Vec<usize> = vec![0; self.groups.len() + 1];
        for (id, cov) in self.countable.iter().enumerate() {
            if *cov as usize >= hist.len() {
                if id != 0 {
                    log::info!("coverage {} of item {} exceeds the number of groups {}, it'll be ignored in the count", cov, id, self.groups.len());
                }
            } else {
                hist[*cov as usize] += graph_storage.node_lens[id] as usize;
            }
        }

        // subtract uncovered bps
        let uncovered_bps = self.uncovered_bps.as_ref().unwrap();
        for (id, uncov) in uncovered_bps.iter() {
            hist[self.countable[*id as usize] as usize] -= uncov;
            // add uncovered bps to 0-coverage count
            hist[0] += uncov;
        }
        hist
    }
}

pub fn quantify_uncovered_bps(
    exclude_table: &Option<ActiveTable>,
    subset_covered_bps: &Option<IntervalContainer>,
    graph_storage: &GraphStorage,
) -> HashMap<ItemIdSize, usize> {
    //
    // 1. if subset is specified, then the node-based coverage calculated by the coverage()
    //    function overestimates the total coverage, because even nodes that are only partially
    //    covered are counted, thus the coverage needs to be reduced by the amount of uncovered
    //    bps from partially covered nodes
    // 2. if exclude is specified, then the coverage is overestimated by the coverage()
    //    function because partially excluded nodes are not excluded in the coverage
    //    calculation, thus the bps coverage needs to be reduced by the amount of excluded bps
    //    from partially excluded nodes
    // 3. if subset AND exclude are specified, nodes that are COMPLETELY excluded have not been
    //    counted in coverage, so they should not be considered here; all other nodes that are
    //    partially excluded / subset have contributed to the overestimation of coverage, so
    //    the bps coverage needs to be reduced by the amount of excluded or not coverered by
    //    any subset interval
    let mut res = HashMap::default();

    if let Some(subset_map) = subset_covered_bps {
        for sid in subset_map.keys() {
            // ignore COMPETELY excluded nodes
            if exclude_table.is_none() || !exclude_table.as_ref().unwrap().items[sid.0 as usize] {
                let l = graph_storage.node_len(sid) as usize;
                let covered = subset_map.total_coverage(
                    sid,
                    &exclude_table
                        .as_ref()
                        .map(|ex| ex.get_active_intervals(sid, l)),
                );
                if covered > l {
                    log::error!("oops, total coverage {} is larger than node length {} for node {}, intervals: {:?}", covered, l, sid.0, subset_map.get(sid).unwrap());
                } else {
                    // report uncovered bps
                    res.insert(sid.0, l - covered);
                }
            }
        }
    }
    res
}

#[cfg(test)]
mod tests {
    use std::io::Error;
    use std::io::Write;
    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_view_params_default() {
        let expected = GraphMaskParameters {
            positive_list: String::new(),
            negative_list: String::new(),
            groupby: String::new(),
            groupby_haplotype: false,
            groupby_sample: false,
        };
        let calculated = GraphMaskParameters::default();
        assert_eq!(calculated, expected);
    }

    fn get_graph_storage_path_segments() -> GraphStorage {
        GraphStorage::from_path_segments(vec![
            PathSegment::create_from_str("s1#2#2"),
            PathSegment::create_from_str("s1#1#2"),
            PathSegment::create_from_str("s1#1#1"),
            PathSegment::create_from_str("s2#1#2"),
        ])
    }

    fn get_load_groups_expected_hashmap(groups: [&str; 4]) -> HashMap<PathSegment, String> {
        HashMap::from([
            (
                PathSegment::create_from_str("s1#1#1"),
                groups[0].to_string(),
            ),
            (
                PathSegment::create_from_str("s1#1#2"),
                groups[1].to_string(),
            ),
            (
                PathSegment::create_from_str("s1#2#2"),
                groups[2].to_string(),
            ),
            (
                PathSegment::create_from_str("s2#1#2"),
                groups[3].to_string(),
            ),
        ])
    }

    #[test]
    fn test_load_groups_haplotype() -> Result<(), Error> {
        let expected = get_load_groups_expected_hashmap(["s1#1", "s1#1", "s1#2", "s2#1"]);
        let graph_storage = get_graph_storage_path_segments();
        let calculated = GraphMask::load_groups("", true, false, &graph_storage)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_load_groups_sample() -> Result<(), Error> {
        let expected = get_load_groups_expected_hashmap(["s1", "s1", "s1", "s2"]);
        let graph_storage = get_graph_storage_path_segments();
        let calculated = GraphMask::load_groups("", false, true, &graph_storage)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    fn get_temporary_file_name_with_content(text: &str) -> Result<(NamedTempFile, String), Error> {
        let mut f = NamedTempFile::new()?;
        writeln!(f, "{}", text)?;
        let msg = "Could not get path as &str";
        let file_name = f
            .path()
            .to_str()
            .ok_or(Error::new(ErrorKind::Other, msg))?
            .to_string();
        Ok((f, file_name))
    }

    #[test]
    fn test_load_groups_file() -> Result<(), Error> {
        let expected = get_load_groups_expected_hashmap(["g1", "g2", "g1", "g2"]);
        let graph_storage = get_graph_storage_path_segments();
        let text = "s1#1#1\tg1
s1#1#2\tg2
s1#2#2\tg1
s2#1#2\tg2";
        let (_file, file_name) = get_temporary_file_name_with_content(text)?;
        let calculated = GraphMask::load_groups(&file_name, false, false, &graph_storage)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_load_groups_none() -> Result<(), Error> {
        let expected = get_load_groups_expected_hashmap(["s1#1#1", "s1#1#2", "s1#2#2", "s2#1#2"]);
        let graph_storage = get_graph_storage_path_segments();
        let calculated = GraphMask::load_groups("", false, false, &graph_storage)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_load_coord_list_none() -> Result<(), Error> {
        let expected = None;
        let calculated = GraphMask::load_coord_list("", &Vec::new())?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_load_coord_list_file() -> Result<(), Error> {
        let expected = Some(vec![
            PathSegment::create_from_str("s1#1#1:0-99"),
            PathSegment::create_from_str("s1#1#1:25-109"),
        ]);
        let text = "s1#1#1\t0\t99
s1#1#1\t25\t109";
        let (_file, file_name) = get_temporary_file_name_with_content(text)?;
        let calculated = GraphMask::load_coord_list(&file_name, &Vec::new())?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_complement_with_group_assignments_no_coords() -> Result<(), Error> {
        let expected: Option<Vec<PathSegment>> = None;
        let groups = HashMap::new();
        let calculated = GraphMask::complement_with_group_assignments(None, &groups)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    fn get_path_segment_with_coordinates(start: usize, end: usize) -> PathSegment {
        PathSegment {
            sample: "1".to_string(),
            haplotype: Some("2".to_string()),
            seqid: Some("3".to_string()),
            start: Some(start),
            end: Some(end),
        }
    }

    #[test]
    fn test_complement_with_group_assignments_coord_with_path_name() -> Result<(), Error> {
        let expected = Some(vec![get_path_segment_with_coordinates(1, 3)]);
        let coords = Some(vec![get_path_segment_with_coordinates(1, 3)]);
        let groups = HashMap::from([(get_path_segment_with_coordinates(8, 6), "g1".to_string())]);
        let calculated = GraphMask::complement_with_group_assignments(coords, &groups)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_complement_with_group_assignments_coord_with_group_name() -> Result<(), Error> {
        let expected = Some(vec![get_path_segment_with_coordinates(8, 6)]);
        let coords = Some(vec![PathSegment::create_from_str("g1")]);
        let groups = HashMap::from([(get_path_segment_with_coordinates(8, 6), "g1".to_string())]);
        let calculated = GraphMask::complement_with_group_assignments(coords, &groups)?;
        assert_eq!(calculated, expected);
        Ok(())
    }

    #[test]
    fn test_complement_with_group_assignments_coord_with_group_name_invalid() {
        let groups = HashMap::from([
            (PathSegment::create_from_str("a#0"), "g1".to_string()),
            (PathSegment::create_from_str("b#0"), "g1".to_string()),
        ]);

        let coords = Some(vec![PathSegment::create_from_str("g1:1-5")]);
        let result = GraphMask::complement_with_group_assignments(coords, &groups);
        assert!(
            result.is_err(),
            "Expected error due to invalid group identifier with start/stop information"
        );
    }

    #[test]
    fn test_complement_with_group_assignments_coord_with_invalid_name() -> Result<(), Error> {
        let expected: Option<Vec<PathSegment>> = Some(Vec::new());
        let groups = HashMap::from([
            (PathSegment::create_from_str("a#0"), "g1".to_string()),
            (PathSegment::create_from_str("b#0"), "g1".to_string()),
        ]);

        let coords = Some(vec![PathSegment::create_from_str("invalid")]);
        let calculated = GraphMask::complement_with_group_assignments(coords, &groups)?;
        assert_eq!(calculated, expected);
        Ok(())
    }
}
