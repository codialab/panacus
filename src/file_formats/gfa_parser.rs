use anyhow::Error;
use itertools::Itertools;

use crate::{
    analyses::info::FileInfo,
    coverage_matrix::{CoverageMatrix, Positions},
    file_formats::{
        gfa_parser::{grammar::Grammar, graph::GraphStorage, util::parse_gfa_paths_walks},
        FileFormatParser,
    },
    hist::Hist,
    html_report::ReportItem,
    util::{
        averagef64, averageu32, median_already_sorted, n50_already_sorted, CountType, GroupSize,
        ItemIdSize,
    },
};

use std::collections::{HashSet, VecDeque};
use std::{collections::HashMap, str};

pub use abacus::AbacusByTotal;

use crate::io::bufreader_from_compressed_gfa;
use crate::util::{ActiveTable, IntervalContainer, ItemTable};

mod abacus;
mod grammar;
mod graph;
mod hist;
mod sparse_matrix;
mod util;
mod walk_splitter;

pub use abacus::GraphMask;
pub use abacus::GraphMaskParameters;
pub use graph::Edge;
pub use graph::ItemId;
pub use graph::Orientation;
pub use graph::PathSegment;
pub use hist::choose;
pub use hist::Hist3D;
pub use hist::ThresholdContainer;
pub use sparse_matrix::SparseMatrix;

type ItemTableWithMetatdata = (
    ItemTable,
    Vec<String>,
    Vec<usize>,
    Vec<String>,
    HashMap<PathSegment, (u32, u32)>,
    HashMap<String, f64>,
    Positions,
);
type PathLens = HashMap<PathSegment, (u32, u32)>;
type PathCollection = HashMap<PathSegment, Vec<(ItemId, Orientation)>>;

pub struct GfaParser {
    // Inputs
    filename: String,
    count_type: CountType,
    graph_mask_parameters: GraphMaskParameters,
    reference: Option<String>,

    // Generated
    graph_storage: GraphStorage,
    graph_mask: GraphMask,
    grammar: Grammar,
}

impl FileFormatParser for GfaParser {
    fn generate_hist(self: Box<Self>) -> Hist {
        let number_of_groups = self.graph_mask.count_groups();
        let mut hist = Hist::from_maximum_coverage(
            number_of_groups,
            self.count_type.to_string(),
            self.get_run_id(),
            self.get_run_name(),
        );
        let (abacus, _, _) = self.get_abacus_by_total();
        // Idx should be the id inside ItemId
        for (idx, feature_coverage) in abacus.countable.iter().enumerate().skip(1) {
            // If we look at nodes or basepairs, make sure this isn't a meta-node
            if (self.count_type == CountType::Node || self.count_type == CountType::Bp)
                && self.graph_storage.node2rule_id[idx] != usize::MAX
            {
                continue;
            }
            if self.count_type == CountType::Bp {
                let length = self.graph_storage.node_lens[idx];
                hist.insert_feature_of_coverage_and_length(
                    *feature_coverage as usize,
                    length as usize,
                );
            } else {
                hist.insert_feature_of_coverage(*feature_coverage as usize);
            }
        }
        let call_count = self.graph_storage.get_call_count();
        log::info!("Called get_node_id {} times", call_count);
        hist
    }

    fn generate_matrix(self: Box<Self>) -> CoverageMatrix {
        let paths_to_collect = match self.reference.as_ref() {
            Some(r) => vec![PathSegment::create_from_str(r)],
            None => vec![],
        };
        let (
            item_table,
            path_names,
            feature_lengths,
            feature_names,
            path_lengths,
            group_lengths,
            mut positions,
        ) = self
            .get_cleaned_item_table(
                &self.graph_mask,
                &self.graph_storage,
                self.count_type,
                &paths_to_collect,
            )
            .expect("Can parse GFA file");
        positions.cleanup();
        let file_info = self.get_file_info_value(path_lengths, group_lengths);
        let mut matrix = CoverageMatrix::new(
            self.count_type.to_string(),
            self.get_run_id(),
            self.get_run_name(),
            file_info,
        );
        matrix.insert_item_table(
            path_names,
            feature_lengths,
            positions,
            feature_names,
            item_table,
        );
        matrix
    }
}

fn get_close_nodes(
    ref_nodes: &Vec<(ItemId, Orientation)>,
    neighbors: &HashMap<(ItemId, Orientation), HashSet<ItemId>>,
) -> HashMap<(ItemId, Orientation), Vec<ItemId>> {
    let mut queue: VecDeque<ItemId> = VecDeque::new();
    let mut closest_ref: HashMap<ItemId, (ItemId, Orientation)> = HashMap::new();
    let ref_node_set: HashSet<ItemId> = ref_nodes.iter().map(|(x, _)| x).copied().collect();
    for (ref_node, orientation) in ref_nodes {
        closest_ref.insert(*ref_node, (*ref_node, *orientation));
        for neighbor in neighbors
            .get(&(*ref_node, *orientation))
            .unwrap_or(&HashSet::new())
        {
            if !ref_node_set.contains(neighbor) && !closest_ref.contains_key(neighbor) {
                closest_ref.insert(*neighbor, (*ref_node, *orientation));
                queue.push_front(*neighbor);
            }
        }
        for neighbor in neighbors
            .get(&(*ref_node, orientation.flip()))
            .unwrap_or(&HashSet::new())
        {
            if !ref_node_set.contains(neighbor) && !closest_ref.contains_key(neighbor) {
                closest_ref.insert(*neighbor, (*ref_node, orientation.flip()));
                queue.push_front(*neighbor);
            }
        }
    }
    while let Some(current_node) = queue.pop_back() {
        for neighbor in neighbors
            .get(&(current_node, Orientation::Forward))
            .unwrap_or(&HashSet::new())
        {
            if !closest_ref.contains_key(neighbor) {
                let current_closest_ref = closest_ref
                    .get(&current_node)
                    .expect("Current is in closest ref");
                closest_ref.insert(*neighbor, *current_closest_ref);
                queue.push_front(*neighbor);
            }
        }
        for neighbor in neighbors
            .get(&(current_node, Orientation::Backward))
            .unwrap_or(&HashSet::new())
        {
            if !closest_ref.contains_key(neighbor) {
                let current_closest_ref = closest_ref
                    .get(&current_node)
                    .expect("Current is in closest ref");
                closest_ref.insert(*neighbor, *current_closest_ref);
                queue.push_front(*neighbor);
            }
        }
    }
    let mut flipped: HashMap<(ItemId, Orientation), Vec<ItemId>> = HashMap::new();
    for (entry_node, ref_node) in closest_ref {
        flipped.entry(ref_node).or_default().push(entry_node);
    }
    flipped
}

fn get_positions_nodes(
    collected_paths: &HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    graph_storage: &GraphStorage,
) -> Positions {
    let edge2id = graph_storage.edge2id.as_ref().expect("Edges need to exist");

    let neighbors: HashMap<(ItemId, Orientation), HashSet<ItemId>> = edge2id
        .keys()
        .flat_map(|x| [((x.0, x.1), x.2), ((x.2, x.3), x.0)])
        .fold(HashMap::new(), |mut acc, (k, v)| {
            acc.entry(k).or_default().insert(v);
            acc
        });
    let mut all_positions = Positions::with_size(graph_storage.node_count);
    for (reference, nodes) in collected_paths.iter() {
        let initial_start = reference.start.unwrap_or(0);
        // Contains the position of each reference node
        let positions: Vec<usize> = nodes
            .iter()
            .scan(initial_start, |acc, (node, _)| {
                let length = graph_storage.node_lens[node.0 as usize] as usize;
                let old_acc = *acc;
                *acc += length;
                Some(old_acc)
            })
            .collect();
        let close_nodes = get_close_nodes(nodes, &neighbors);
        for (node, position) in nodes.iter().zip(positions.iter()) {
            all_positions.set(
                node.0 .0 as usize - 1,
                reference.to_string().as_str(),
                *position,
            );
            let closest_ones = &close_nodes[node];
            for offref_node in closest_ones {
                all_positions.set(
                    offref_node.0 as usize - 1,
                    reference.to_string().as_str(),
                    *position,
                );
            }
        }
    }
    all_positions
}

fn get_positions_edges(
    collected_paths: &HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    graph_storage: &GraphStorage,
) -> Positions {
    let node_positions = get_positions_nodes(collected_paths, graph_storage);
    let mut edge_positions = Positions::with_size(graph_storage.edge_count);
    let edge2id = graph_storage
        .edge2id
        .as_ref()
        .expect("Edges have been collected");
    // TODO: does this need to be done more intellegently? Approach seems a bit naive
    for (Edge(node1, _orientation1, _node2, _orientation2), &id) in edge2id.iter() {
        let (reference, position) = node_positions.get(node1.0 as usize - 1);
        edge_positions.set(id.0 as usize - 1, reference, position);
    }
    edge_positions
}

fn get_positions(
    count_type: CountType,
    collected_paths: &HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    graph_storage: &GraphStorage,
) -> Positions {
    match count_type {
        CountType::Node | CountType::Bp => get_positions_nodes(collected_paths, graph_storage),
        CountType::Edge => get_positions_edges(collected_paths, graph_storage),
    }
}

impl GfaParser {
    pub fn new(
        filename: &str,
        count_type: CountType,
        graph_mask_parameters: GraphMaskParameters,
        reference: Option<String>,
        is_nice: bool,
    ) -> Result<Self, Error> {
        let mut grammar = Grammar::from(count_type);
        let (mut graph_storage, has_meta_node) = GraphStorage::from_gfa(filename, is_nice);
        if has_meta_node {
            grammar.parse_gfa_parallel(filename, &mut graph_storage);
            log::info!("found {} rules", grammar.len());
        }
        let graph_mask = GraphMask::from_datamgr(&graph_mask_parameters, &graph_storage)?;
        Ok(Self {
            filename: filename.to_owned(),
            count_type,
            graph_mask_parameters,
            reference,
            graph_storage,
            graph_mask,
            grammar,
        })
    }

    fn get_file_info_value(
        &self,
        paths_len: HashMap<PathSegment, (u32, u32)>,
        group_lens: HashMap<String, f64>,
    ) -> FileInfo {
        // Created here to satisfy the borrow checker
        let empty_degree = vec![0, 0];
        let degree = self.graph_storage.degree.as_ref().unwrap_or(&empty_degree);
        let mut node_lens_sorted = self.graph_storage.node_lens[1..].to_vec();
        node_lens_sorted.sort_by(|a, b| b.cmp(a)); // decreasing, for N50
        let mut components = connected_components(
            self.graph_storage
                .edge2id
                .as_ref()
                .unwrap_or(&HashMap::new()),
            &self.graph_storage.get_nodes(),
        );
        components.sort();

        let paths_bp_len: Vec<_> = paths_len.values().map(|x| x.1).collect();
        let paths_len: Vec<_> = paths_len.values().map(|x| x.0).collect();

        let mut group_lens: Vec<(String, f64)> = group_lens.into_iter().collect();
        group_lens.sort_by(|a, b| (b.1).partial_cmp(&a.1).unwrap());
        let table = group_lens.clone();
        let (group_names, group_bp_lens): (Vec<String>, Vec<f64>) = group_lens.into_iter().unzip();

        // group_count: gb.get_group_count(),
        let mut file_info = FileInfo::new("gfa");
        file_info.add_info(
            "Number of nodes",
            self.graph_storage.node_count.to_string().as_str(),
        );
        file_info.add_info(
            "Number of edges",
            self.graph_storage.edge_count.to_string().as_str(),
        );
        file_info.add_info(
            "Average degree",
            averageu32(&degree[1..]).to_string().as_str(),
        );
        file_info.add_info(
            "Max degree",
            degree[1..].iter().max().unwrap().to_string().as_str(),
        );
        file_info.add_info(
            "Min degree",
            degree[1..].iter().min().unwrap().to_string().as_str(),
        );
        file_info.add_info(
            "Number of zero-degree nodes",
            degree[1..]
                .iter()
                .filter(|&x| *x == 0)
                .count()
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Number of connected components",
            (components.len() as u32).to_string().as_str(),
        );
        file_info.add_info(
            "Largest component size",
            (*components.iter().max().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Smallest component size",
            (*components.iter().min().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Median component size",
            (median_already_sorted(&components)).to_string().as_str(),
        );
        file_info.add_info(
            "Longest node",
            (*node_lens_sorted.iter().max().unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Shortest node",
            (*node_lens_sorted.iter().min().unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Average node size",
            (averageu32(&node_lens_sorted)).to_string().as_str(),
        );
        file_info.add_info(
            "Median node size",
            (median_already_sorted(&node_lens_sorted))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "N50 nodes",
            (n50_already_sorted(&node_lens_sorted).unwrap())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Number of basepairs",
            (node_lens_sorted.iter().sum::<u32>()).to_string().as_str(),
        );
        file_info.add_info(
            "Number of groups",
            (self.graph_mask.count_groups()).to_string().as_str(),
        );
        file_info.add_info(
            "Number of paths",
            (self.graph_storage.path_segments.len())
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Longest path in nodes",
            (*paths_len.iter().max().unwrap_or(&0)).to_string().as_str(),
        );
        file_info.add_info(
            "Shortest path in nodes",
            (*paths_len.iter().min().unwrap_or(&0)).to_string().as_str(),
        );
        file_info.add_info(
            "Average path length in nodes",
            (averageu32(&paths_len)).to_string().as_str(),
        );
        file_info.add_info(
            "Longest path in bps",
            (*paths_bp_len.iter().max().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Shortest path in bps",
            (*paths_bp_len.iter().min().unwrap_or(&0))
                .to_string()
                .as_str(),
        );
        file_info.add_info(
            "Average path length in bps",
            (averageu32(&paths_bp_len)).to_string().as_str(),
        );

        file_info.add_info(
            "Longest group length in bps",
            &format!("{} ({})", group_bp_lens[0], group_names[0]),
        );
        file_info.add_info(
            "Shortest group length in bps",
            &format!(
                "{} ({})",
                group_bp_lens[group_bp_lens.len() - 1],
                group_names[group_names.len() - 1]
            ),
        );
        file_info.add_info(
            "Average group length in bps",
            (averagef64(&group_bp_lens).round() as u64)
                .to_string()
                .as_str(),
        );

        let paths_plot = ReportItem::Bar {
            id: format!("info-{}-paths-plot", self.filename),
            name: format!("info-{}-paths-plot", self.filename),
            x_label: "Group".to_string(),
            y_label: "Length in bps".to_string(),
            labels: group_names,
            values: group_bp_lens,
            log_toggle: true,
        };
        file_info.add_plot(paths_plot);

        file_info.add_table("group-lengths".to_string(), table);

        file_info
    }

    fn get_run_name(&self) -> String {
        format!(
            "{} {} {}",
            self.filename,
            self.graph_mask_parameters.positive_list,
            self.graph_mask_parameters.groupby
        )
    }

    fn get_run_id(&self) -> String {
        format!(
            "{}-{}-{}",
            self.filename,
            self.graph_mask_parameters.positive_list,
            self.graph_mask_parameters.groupby
        )
    }

    pub fn get_abacus_by_total(&self) -> (AbacusByTotal, PathLens, PathCollection) {
        let mut data = bufreader_from_compressed_gfa(&self.filename);
        let (abacus, path_lens) = AbacusByTotal::from_gfa(
            &mut data,
            &self.graph_mask,
            &self.graph_storage,
            &self.grammar,
            self.count_type,
        );
        let collected_paths = HashMap::new();
        (abacus, path_lens, collected_paths)
    }

    fn get_cleaned_item_table(
        &self,
        graph_mask: &GraphMask,
        graph_storage: &GraphStorage,
        count: CountType,
        paths_to_collect: &[PathSegment],
    ) -> anyhow::Result<ItemTableWithMetatdata> {
        log::info!("parsing path + walk sequences");
        let mut data = bufreader_from_compressed_gfa(&self.filename);
        let (item_table, exclude_table, subset_covered_bps, paths_len, collected_paths) =
            parse_gfa_paths_walks(
                &mut data,
                graph_mask,
                graph_storage,
                &self.grammar,
                &count,
                paths_to_collect,
            );

        let mut path_order: Vec<(ItemIdSize, GroupSize)> = Vec::new();
        let mut groups: Vec<String> = Vec::new();

        for (path_id, group_id) in graph_mask.get_path_order(&graph_storage.path_segments) {
            log::debug!(
                "processing path {} (group {})",
                &graph_storage.path_segments[path_id as usize],
                group_id
            );
            if groups.is_empty() || groups.last().unwrap() != group_id {
                groups.push(group_id.to_string());
            }
            //if groups.len() > 65534 {
            //    panic!("data has more than 65534 path groups, but command is not supported for more than 65534");
            //}
            path_order.push((path_id, (groups.len() - 1) as GroupSize));
        }

        log::info!("Collected {} paths", collected_paths.len());
        let mut positions = if collected_paths.is_empty() {
            Positions::with_size(graph_storage.node_count)
        } else {
            get_positions(count, &collected_paths, graph_storage)
        };

        let path_names: Vec<String> = self
            .graph_storage
            .path_segments
            .iter()
            .map(|x| x.to_string())
            .collect();
        let group_names: Vec<String> = self
            .graph_storage
            .path_segments
            .iter()
            .map(|x| {
                let x = &x.clear_coords();
                if let Some(name) = graph_mask.groups.get(x) {
                    name.clone()
                } else {
                    panic!("Could not get {} from {:?}", x, graph_mask.groups)
                }
            })
            .collect_vec();
        let mut group_lens: HashMap<String, f64> = HashMap::new();
        paths_len.iter().for_each(|(k, v)| {
            let p = &k.clear_coords();
            let g = if let Some(name) = graph_mask.groups.get(p) {
                name.clone()
            } else {
                panic!("Could not get {} from {:?}", p, graph_mask.groups)
            };
            let bp_len = v.1 as f64;
            *group_lens.entry(g).or_default() += bp_len;
        });
        let mgroup_names: HashSet<&str> = graph_mask
            .get_path_order(&self.graph_storage.path_segments)
            .into_iter()
            .map(|x| x.1)
            .collect();
        let (item_table, groups) = if group_names == path_names
            && group_names.len() == mgroup_names.len()
            && exclude_table.is_none()
        {
            (item_table, group_names)
        } else {
            Self::collapse_item_table(item_table, group_names, mgroup_names, &exclude_table)
        };
        let (feature_lengths, feature_names) = Self::get_feature_lengths(
            graph_storage,
            &exclude_table,
            &subset_covered_bps,
            count,
            &mut positions,
        );

        Ok((
            item_table,
            groups,
            feature_lengths,
            feature_names,
            paths_len,
            group_lens,
            positions,
        ))
    }

    fn get_feature_lengths(
        graph_storage: &GraphStorage,
        exclude_table: &Option<ActiveTable>,
        subset_covered_bps: &Option<IntervalContainer>,
        count_type: CountType,
        positions: &mut Positions,
    ) -> (Vec<usize>, Vec<String>) {
        let mut feature_lengths = match count_type {
            CountType::Node => vec![1; graph_storage.node_count],
            CountType::Edge => vec![1; graph_storage.edge_count],
            CountType::Bp => graph_storage
                .node_lens
                .iter()
                .skip(1)
                .map(|x| *x as usize)
                .collect_vec(),
        };
        let mut feature_names = match count_type {
            CountType::Node | CountType::Bp => {
                let rev: HashMap<&ItemId, &str> = graph_storage
                    .node2id
                    .iter()
                    .map(|(k, v)| (v, std::str::from_utf8(k).unwrap()))
                    .collect();
                let names: Vec<String> = (1..feature_lengths.len() + 1)
                    .map(|x| rev[&ItemId(x as u64)].to_string())
                    .collect();
                names
            }
            CountType::Edge => {
                let rev: HashMap<&ItemId, String> = graph_storage
                    .edge2id
                    .as_ref()
                    .unwrap()
                    .iter()
                    .map(|(k, v)| (v, k.to_string()))
                    .collect();
                let names: Vec<String> = (1..feature_lengths.len() + 1)
                    .map(|x| rev[&ItemId(x as u64)].to_string())
                    .collect();
                names
            }
        };

        if let Some(e) = exclude_table.as_ref() {
            // Shorten length
            for annotated_feature in e.get_annotation_keys() {
                let idx = annotated_feature.0 as usize - 1;
                for (start, end) in e.get_active_intervals(annotated_feature, 0) {
                    let len = end - start;
                    feature_lengths[idx] -= len;
                }
            }
            // Iterator over excluded features (this has a zero element)
            // TODO can we just remove the zero element in the exclude_table
            // definition?
            let mut iter = e.items.iter().skip(1);
            feature_lengths.retain(|_| !iter.next().unwrap());

            // Do the same for the feature names
            let mut iter = e.items.iter().skip(1);
            feature_names.retain(|_| !iter.next().unwrap());

            let flipped: Vec<bool> = e.items.iter().skip(1).map(|x| !x).collect();
            positions.apply_mask(&flipped);
        }

        if let Some(s) = subset_covered_bps.as_ref() {
            if count_type == CountType::Bp {
                for key in s.keys() {
                    let feature_idx = key.0 as usize - 1;
                    let subset_length: usize = s
                        .get(key)
                        .unwrap()
                        .iter()
                        .map(|(start, end)| end - start)
                        .sum();
                    if subset_length < feature_lengths[feature_idx] {
                        feature_lengths[feature_idx] = subset_length;
                    }
                }
            }
        }

        if count_type == CountType::Node || count_type == CountType::Bp {
            let to_keep: Vec<bool> = feature_names
                .iter()
                .map(|x| {
                    let item_id = graph_storage.get_node_id(x.as_bytes()).unwrap();
                    graph_storage.node2rule_id[item_id.0 as usize] == usize::MAX
                })
                .collect();
            let mut to_keep_lengths = to_keep.iter();
            let mut to_keep_names = to_keep.iter();
            feature_lengths.retain(|_| *to_keep_lengths.next().unwrap());
            feature_names.retain(|_| *to_keep_names.next().unwrap());
            positions.apply_mask(&to_keep);
        }

        (feature_lengths, feature_names)
    }

    fn collapse_item_table(
        item_table: ItemTable,
        groups: Vec<String>,
        masked_groups: HashSet<&str>,
        exclude_table: &Option<ActiveTable>,
    ) -> (ItemTable, Vec<String>) {
        let mut items: HashMap<&str, Vec<ItemIdSize>> = HashMap::new();
        // Translation table to compress id space (in case features are excluded)
        let translation_of_ids: Option<HashMap<usize, usize>> = exclude_table.as_ref().map(|e| {
            e.items
                .iter()
                .enumerate()
                .filter_map(|(idx, active)| if *active { None } else { Some(idx) })
                .enumerate()
                .map(|(new_idx, old_idx)| (old_idx, new_idx))
                .collect()
        });
        //log::info!("Translation table: {:?}", translation_of_ids);
        for (idx, group) in groups.iter().enumerate() {
            let (start, end) = (
                item_table.id_prefsum[idx] as usize,
                item_table.id_prefsum[idx + 1] as usize,
            );
            items
                .entry(group)
                .or_default()
                .extend(&item_table.items[start..end]);
        }
        let mut new_items = Vec::new();
        let mut new_idprefsum = vec![0];
        let mut keys = Vec::new();
        for (key, mut value) in items.into_iter() {
            // If path is excluded, skip it
            if !masked_groups.contains(key) {
                continue;
            }
            if let Some(e) = exclude_table.as_ref() {
                // If the feature should be removed, remove it
                value.retain(|x| !e.is_active(&ItemId(*x)));
                // Compress the id space (excluded features should take no index)
                value.iter_mut().for_each(|x| {
                    let new_idx = translation_of_ids.as_ref().unwrap()[&(*x as usize)];
                    *x = new_idx as u64;
                });
            }
            let len = value.len();
            new_items.append(&mut value);
            new_idprefsum.push(new_idprefsum.last().unwrap() + len as u64);
            keys.push(key.to_string());
        }
        (
            ItemTable {
                items: new_items,
                id_prefsum: new_idprefsum,
            },
            keys,
        )
    }

    // fn parse_paths_walks(
    //     &self,
    //     count_type: &CountType,
    //     paths_to_collect: &Vec<PathSegment>,
    // ) -> (
    //     ItemTable,
    //     Option<ActiveTable>,
    //     Option<IntervalContainer>,
    //     HashMap<PathSegment, (u32, u32)>,
    //     HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
    // ) {
    //     log::info!("parsing path + walk sequences");
    //     log::info!("collecting: {:?}", paths_to_collect);
    //     let mut data = bufreader_from_compressed_gfa(&self.filename);
    //     let graph_storage = &self.graph_storage;
    //     let graph_mask = &self.graph_mask;
    //     let mut item_table = ItemTable::new(graph_storage.path_segments.len());

    //     let (mut subset_covered_bps, mut exclude_table, include_map, exclude_map) =
    //         graph_mask.load_optional_subsetting(graph_storage, count_type);

    //     let mut collected_paths: HashMap<PathSegment, Vec<(ItemId, Orientation)>> = HashMap::new();
    //     let mut num_path = 0;
    //     let complete: Vec<(usize, usize)> = vec![(0, usize::MAX)];
    //     let mut paths_len: HashMap<PathSegment, (u32, u32)> = HashMap::new();

    //     // Used to check whether a node was already fully included in a previous path
    //     let mut fully_included = vec![false; graph_storage.node_count + 1];

    //     let mut buf = vec![];
    //     let timer = Instant::now();
    //     while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
    //         if buf[0] == b'P' || buf[0] == b'W' {
    //             let (path_seg, buf_path_seg) = match buf[0] {
    //                 b'P' => parse_path_identifier(&buf),
    //                 b'W' => parse_walk_identifier(&buf),
    //                 _ => unreachable!(),
    //             };

    //             log::debug!("processing path {}", &path_seg);

    //             let include_coords =
    //                 get_include_coords(graph_mask, &complete, &include_map, &path_seg);
    //             let exclude_coords = get_exclude_maps(graph_mask, &exclude_map, &path_seg);

    //             let (start, end) = path_seg.coords().unwrap_or((0, usize::MAX));

    //             // do not process the path sequence if path is neither part of subset nor exclude
    //             if graph_mask.include_coords.is_some()
    //                 && !intersects(include_coords, &(start, end))
    //                 && !intersects(exclude_coords, &(start, end))
    //             {
    //                 log::debug!("path {} does not intersect with subset coordinates {:?} nor with exclude coordinates {:?} and therefore is skipped from processing", &path_seg, &include_coords, &exclude_coords);
    //                 skip_path(&mut item_table, &mut num_path, &mut buf);
    //                 continue;
    //             }

    //             if *count_type != CountType::Edge
    //                 && (graph_mask.include_coords.is_none()
    //                     || is_contained(include_coords, &(start, end)))
    //                 && (graph_mask.exclude_coords.is_none()
    //                     || is_contained(exclude_coords, &(start, end)))
    //             {
    //                 log::debug!("path {} is fully contained within subset coordinates {:?} and is eligible for full parallel processing", path_seg, include_coords);
    //                 let mut none = None;
    //                 let ex: &mut Option<ActiveTable> = if exclude_coords.is_empty() {
    //                     &mut none
    //                 } else {
    //                     &mut exclude_table
    //                 };
    //                 let (num_added_nodes, bp_len) = match buf[0] {
    //                     b'P' => parse_path_seq_update_tables(
    //                         buf_path_seg,
    //                         graph_storage,
    //                         &mut item_table,
    //                         ex.as_mut(),
    //                         num_path,
    //                     ),
    //                     b'W' => parse_walk_seq_update_tables(
    //                         buf_path_seg,
    //                         graph_storage,
    //                         &self.grammar,
    //                         &mut item_table,
    //                         ex.as_mut(),
    //                         num_path,
    //                     ),
    //                     _ => unreachable!(),
    //                 };
    //                 paths_len.insert(path_seg.clone(), (num_added_nodes, bp_len));
    //             } else {
    //                 let sids = match buf[0] {
    //                     b'P' => parse_path_seq_to_item_vec(buf_path_seg, graph_storage),
    //                     b'W' => {
    //                         parse_walk_seq_to_item_vec(buf_path_seg, graph_storage, &self.grammar)
    //                     }
    //                     _ => unreachable!(),
    //                 };
    //                 if paths_to_collect.iter().any(|p| path_seg.is_part_of(p)) {
    //                     collected_paths.insert(path_seg.clone(), sids.clone());
    //                 }
    //                 match count_type {
    //                     CountType::Node | CountType::Bp => {
    //                         let (node_len, bp_len) = update_tables(
    //                             &mut item_table,
    //                             &mut subset_covered_bps.as_mut(),
    //                             &mut fully_included,
    //                             &mut exclude_table.as_mut(),
    //                             num_path,
    //                             graph_storage,
    //                             sids,
    //                             include_coords,
    //                             exclude_coords,
    //                             start,
    //                         );
    //                         paths_len.insert(path_seg.clone(), (node_len as u32, bp_len as u32));
    //                     }
    //                     CountType::Edge => update_tables_edgecount(
    //                         &mut item_table,
    //                         &mut exclude_table.as_mut(),
    //                         num_path,
    //                         graph_storage,
    //                         sids,
    //                         include_coords,
    //                         exclude_coords,
    //                         start,
    //                     ),
    //                 };
    //             }
    //             num_path += 1;
    //         }
    //         buf.clear();
    //     }
    //     let duration = timer.elapsed();
    //     log::info!(
    //         "func done; count: {:?}; time elapsed: {:?}",
    //         count_type,
    //         duration
    //     );

    //     (
    //         item_table,
    //         exclude_table,
    //         subset_covered_bps,
    //         paths_len,
    //         collected_paths,
    //     )
    // }
}

fn connected_components(edge2id: &HashMap<Edge, ItemId>, nodes: &Vec<ItemId>) -> Vec<u32> {
    let mut component_lengths = Vec::new();
    let mut visited: HashSet<ItemId> = HashSet::new();
    let edges: HashMap<ItemId, Vec<ItemId>> = edge2id
        .keys()
        .map(|x| (x.0, x.2))
        .chain(edge2id.keys().map(|x| (x.2, x.0)))
        .fold(HashMap::new(), |mut acc, (k, v)| {
            acc.entry(k).and_modify(|x| x.push(v)).or_insert(vec![v]);
            acc
        });
    for node in nodes {
        if !visited.contains(node) {
            component_lengths.push(dfs(&edges, *node, &mut visited));
        }
    }
    component_lengths
}

fn dfs(edges: &HashMap<ItemId, Vec<ItemId>>, node: ItemId, visited: &mut HashSet<ItemId>) -> u32 {
    let mut s = Vec::new();
    let mut length = 0;
    s.push(node);
    while let Some(v) = s.pop() {
        if visited.contains(&v) {
            continue;
        }
        visited.insert(v);
        length += 1;
        if !edges.contains_key(&v) {
            continue;
        }
        for neigh in &edges[&v] {
            if !visited.contains(neigh) {
                s.push(*neigh);
            }
        }
    }
    length
}
