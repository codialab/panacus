use crossbeam_channel::bounded;
use lazy_static::lazy_static;
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::sync::OnceLock;
use std::{collections::HashSet, io::BufRead, thread, time::Instant};

use regex::Regex;

use crate::util::CountType;
use crate::{
    file_formats::gfa_parser::{
        graph::GraphStorage,
        walk_splitter::{WalkByteSplitter, WalkSplitter},
        ItemId, Orientation,
    },
    io::bufreader_from_compressed_gfa,
};

lazy_static! {
    static ref RE_WALK: Regex = Regex::new(r"([><])([!-;=?-~]+)").unwrap();
}

const MAX_CACHE_LEN: usize = 64;

pub struct Grammar {
    rules: Vec<usize>,
    child_nodes: Vec<ItemId>,
    // This is only used if the count_type is Edge
    child_orientations: Vec<Orientation>,
    rule_node_cache: Vec<OnceLock<Box<[ItemId]>>>,
    rule_orientation_cache: Vec<OnceLock<Box<[Orientation]>>>,
    count_type: CountType,
}

impl Grammar {
    pub fn from(count_type: CountType) -> Self {
        Self {
            rules: vec![0],
            child_nodes: Vec::new(),
            child_orientations: Vec::new(),
            rule_node_cache: Vec::new(),
            rule_orientation_cache: Vec::new(),
            count_type,
        }
    }
}

impl Grammar {
    pub fn insert(&mut self, rule: usize, rule_text: Vec<(ItemId, Orientation)>) {
        if rule != self.rules.len() - 1 {
            log::error!("Grammar rules should be inserted in order");
        }
        self.rules
            .push(self.rules.last().unwrap() + rule_text.len());
        let (child_nodes, child_orientations): (Vec<ItemId>, Vec<Orientation>) =
            rule_text.into_iter().unzip();
        self.child_nodes.extend(child_nodes);
        if self.count_type == CountType::Edge {
            self.child_orientations.extend(child_orientations);
            self.rule_orientation_cache.push(OnceLock::new());
        }
        self.rule_node_cache.push(OnceLock::new());
    }

    fn get_topological_sort(&self, node2rule_id: &[usize]) -> Vec<ItemId> {
        log::info!("Calculating Topological Sort of Grammar");
        let number_of_rules = self.rules.len() - 1;
        let mut l: Vec<ItemId> = Vec::with_capacity(number_of_rules);
        let mut marked = vec![false; number_of_rules];
        for (idx, &rule_id) in node2rule_id.iter().enumerate() {
            if rule_id == usize::MAX {
                continue;
            }
            if !marked[rule_id] {
                let item_id = ItemId(idx as u64);
                self.visit(rule_id, item_id, &mut marked, &mut l, node2rule_id);
            }
        }
        l.reverse();
        l
    }

    fn visit(
        &self,
        rule: usize,
        item_id: ItemId,
        marked: &mut [bool],
        l: &mut Vec<ItemId>,
        node2rule_id: &[usize],
    ) {
        if marked[rule] {
            return;
        }
        marked[rule] = true;

        let (start, end) = (self.rules[rule], self.rules[rule + 1]);
        for child in &self.child_nodes[start..end] {
            if node2rule_id[child.0 as usize] != usize::MAX {
                self.visit(
                    node2rule_id[child.0 as usize],
                    *child,
                    marked,
                    l,
                    node2rule_id,
                );
            }
        }
        l.push(item_id)
    }

    pub fn push_down_values(&self, values: &mut [u32], node2rule_id: &[usize]) {
        let sorting = self.get_topological_sort(node2rule_id);
        // check topological sorting:
        let mut seen: HashSet<ItemId> = HashSet::new();
        for rule in sorting.iter().rev() {
            seen.insert(*rule);
            let rule_id = node2rule_id[rule.0 as usize];
            let (start, end) = (self.rules[rule_id], self.rules[rule_id + 1]);
            for child in &self.child_nodes[start..end] {
                let child_id = child.0 as usize;
                if node2rule_id[child_id] != usize::MAX && !seen.contains(child) {
                    log::error!(
                        "This non-terminal {:?} was not seen before, topological sort is wrong!",
                        child
                    );
                    panic!();
                }
            }
        }
        for rule in sorting {
            let value = values[rule.0 as usize];
            let rule_id = node2rule_id[rule.0 as usize];
            values[rule.0 as usize] = 0;
            let (start, end) = (self.rules[rule_id], self.rules[rule_id + 1]);
            for child in &self.child_nodes[start..end] {
                let child_id = child.0 as usize;
                values[child_id] += value;
            }
        }
    }

    pub fn get_nodes(&self, rule: usize) -> &[ItemId] {
        let (start, end) = (self.rules[rule], self.rules[rule + 1]);
        &self.child_nodes[start..end]
    }

    pub fn get_orientations(&self, rule: usize) -> &[Orientation] {
        if self.count_type != CountType::Edge {
            panic!("Orientations are only set when CountType is Edge")
        }
        let (start, end) = (self.rules[rule], self.rules[rule + 1]);
        &self.child_orientations[start..end]
    }

    pub fn parse_gfa(&mut self, filename: &str, graph_storage: &GraphStorage) {
        let buf = bufreader_from_compressed_gfa(filename);
        self.parse_buffer(buf, graph_storage);
    }

    pub fn parse_gfa_parallel(&mut self, filename: &str, graph_storage: &mut GraphStorage) {
        let buf = bufreader_from_compressed_gfa(filename);
        self.parse_buffer_in_parallel(buf, graph_storage);
    }

    pub fn is_empty(&self) -> bool {
        self.child_nodes.is_empty()
    }

    pub fn len(&self) -> usize {
        self.rules.len() - 1
    }

    pub fn decompress_unordered(
        &self,
        segment_id: ItemId,
        graph_storage: &GraphStorage,
    ) -> Vec<ItemId> {
        let mut sids = Vec::new();
        if self.is_empty() {
            sids.push(segment_id);
        } else {
            let mut stack: Vec<DecompOp<ItemId>> = vec![DecompOp::Element(segment_id)];

            while let Some(op) = stack.pop() {
                match op {
                    DecompOp::Element(current) => {
                        let rule_idx = graph_storage.node2rule_id[current.0 as usize];
                        if rule_idx != usize::MAX {
                            // We have cached the rule already
                            if let Some(cached) = self.rule_node_cache[rule_idx].get() {
                                sids.extend_from_slice(cached);
                                continue;
                            }
                            stack.push(DecompOp::RuleEnd {
                                rule_id: rule_idx,
                                is_forward: true,
                                output_start: sids.len(),
                            });
                            let children = self.get_nodes(rule_idx);
                            stack.extend(children.iter().rev().map(|&x| DecompOp::Element(x)));
                        } else {
                            sids.push(current);
                        }
                    }
                    DecompOp::RuleEnd {
                        rule_id,
                        is_forward: _,
                        output_start,
                    } => {
                        let rule_len = sids.len() - output_start;

                        if rule_len > 0 && rule_len <= MAX_CACHE_LEN {
                            let expanded_slice: Box<[ItemId]> = sids[output_start..].into();
                            let _ = self.rule_node_cache[rule_id].set(expanded_slice);
                        }
                    }
                }
            }
        }
        sids
    }

    pub fn decompress_ordered(
        &self,
        segment_id: ItemId,
        orientation: Orientation,
        graph_storage: &GraphStorage,
    ) -> Vec<(ItemId, Orientation)> {
        let mut sids = Vec::new();
        if self.is_empty() {
            sids.push((segment_id, orientation))
        } else {
            let mut stack = vec![DecompOp::Element((segment_id, orientation))];
            while let Some(op) = stack.pop() {
                match op {
                    DecompOp::Element((current_node, current_orientation)) => {
                        let current_id = graph_storage.node2rule_id[current_node.0 as usize];
                        if current_id != usize::MAX {
                            if let (Some(cached_nodes), Some(cached_orientations)) = (
                                self.rule_node_cache[current_id].get(),
                                self.rule_orientation_cache[current_id].get(),
                            ) {
                                if current_orientation == Orientation::Forward {
                                    let iter = cached_nodes
                                        .iter()
                                        .copied()
                                        .zip(cached_orientations.iter().copied());
                                    sids.extend(iter);
                                } else {
                                    // If we need the backward direction add reverse complement of cache
                                    let iter = cached_nodes
                                        .iter()
                                        .copied()
                                        .zip(cached_orientations.iter().map(|o| o.flip()))
                                        .rev();
                                    sids.extend(iter);
                                }
                                continue;
                            }
                            stack.push(DecompOp::RuleEnd {
                                rule_id: current_id,
                                is_forward: current_orientation == Orientation::Forward,
                                output_start: sids.len(),
                            });
                            match current_orientation {
                                Orientation::Forward => {
                                    stack.extend(
                                        self.get_nodes(current_id)
                                            .iter()
                                            .copied()
                                            .zip(self.get_orientations(current_id).iter().copied())
                                            .map(DecompOp::Element)
                                            .rev(),
                                    );
                                }
                                Orientation::Backward => {
                                    stack.extend(
                                        self.get_nodes(current_id)
                                            .iter()
                                            .zip(self.get_orientations(current_id).iter())
                                            .map(|(i, o)| DecompOp::Element((*i, o.flip()))),
                                    );
                                }
                            }
                        } else {
                            sids.push((current_node, current_orientation));
                        }
                    }
                    DecompOp::RuleEnd {
                        rule_id,
                        is_forward,
                        output_start,
                    } => {
                        let rule_len = sids.len() - output_start;

                        if rule_len > 0 && rule_len <= MAX_CACHE_LEN {
                            let (node_slice, orientation_slice): (Vec<ItemId>, Vec<Orientation>) =
                                match is_forward {
                                    true => sids[output_start..].iter().copied().unzip(),
                                    // If we put the reverse rule into sids, put the forward rule in the cache
                                    false => sids[output_start..]
                                        .iter()
                                        .rev()
                                        .map(|&(n, o)| (n, o.flip()))
                                        .unzip(),
                                };
                            let node_slice = node_slice.into_boxed_slice();
                            let orientation_slice = orientation_slice.into_boxed_slice();
                            let _ = self.rule_node_cache[rule_id].set(node_slice);
                            let _ = self.rule_orientation_cache[rule_id].set(orientation_slice);
                        }
                    }
                }
            }
        }
        sids
    }

    fn parse_buffer(&mut self, buf: impl BufRead, graph_storage: &GraphStorage) {
        let lines = buf.lines().map(|l| l.expect("Failed to read line"));

        let timer = Instant::now();
        for line in lines {
            if line.starts_with("Q\t") {
                let mut fields = line.split("\t").skip(1);
                let name: &[u8] = fields.next().expect("Q line has segment name").as_bytes();
                let node_id = graph_storage.get_node_id(name).unwrap();
                let rule_id = graph_storage.node2rule_id[node_id.0 as usize];
                let content = fields.next().expect("Q line has content");
                let splitter = WalkSplitter::new(content);
                let nodes: Vec<(ItemId, Orientation)> = splitter
                    .map(|m| {
                        let (orientation, node) = m.split_at(1);
                        let orientation = Orientation::from_lg(orientation.as_bytes()[0]);
                        let node = graph_storage.get_node_id(node.as_bytes()).unwrap();
                        (node, orientation)
                    })
                    .collect();
                self.insert(rule_id, nodes);
            }
        }
        let duration = timer.elapsed();
        log::info!(
            "grammar parsing done; count: {:?}; time elapsed: {:?}",
            self.rules.len() - 1,
            duration
        );
    }

    fn parse_buffer_in_parallel<R>(&mut self, mut buf: R, graph_storage: &mut GraphStorage)
    where
        R: BufRead + Send + 'static,
    {
        let (raw_tx, raw_rx) = bounded::<Vec<u8>>(16);
        let (writer_tx, writer_rx) = bounded::<(ItemId, Vec<ItemId>, Vec<Orientation>)>(128);
        let lines_per_chunk = 100_000;

        // Read file into chunks and put them into raw channel
        let io_thread = thread::spawn(move || {
            while let Some(chunk) = read_lines_to_chunk(&mut buf, lines_per_chunk) {
                if raw_tx.send(chunk).is_err() {
                    break;
                }
            }
        });

        let writer_thread = thread::spawn(move || {
            let mut rules = vec![0];
            let mut child_nodes: Vec<ItemId> = Vec::new();
            let mut child_orientations: Vec<Orientation> = Vec::new();
            let mut order: Vec<ItemId> = Vec::new();

            for (rule, rule_children, rule_orientations) in writer_rx {
                child_nodes.extend_from_slice(&rule_children);
                child_orientations.extend_from_slice(&rule_orientations);
                rules.push(child_nodes.len());
                order.push(rule);
            }
            (rules, child_nodes, child_orientations, order)
        });

        raw_rx.into_iter().par_bridge().for_each(|chunk_bytes| {
            for line in chunk_bytes.split(|&b| b == b'\n') {
                if !line.starts_with(b"Q") {
                    continue;
                }
                let mut fields = line.split(|&c| c == b'\t').skip(1);
                let name: &[u8] = fields.next().expect("Q line has segment name");
                let node_id = graph_storage.get_node_id(name).unwrap();
                let content = fields.next().expect("Q line has content");
                let splitter = WalkByteSplitter::new(content);
                let (nodes, orientations): (Vec<ItemId>, Vec<Orientation>) = splitter
                    .map(|m| {
                        let (orientation, node) = m.split_at(1);
                        let orientation = Orientation::from_lg(orientation[0]);
                        let node = graph_storage.get_node_id(node).unwrap();
                        (node, orientation)
                    })
                    .unzip();
                writer_tx.send((node_id, nodes, orientations)).unwrap();
            }
        });

        drop(writer_tx);
        io_thread.join().unwrap();
        let (rules, child_nodes, child_orientations, order) = writer_thread.join().unwrap();
        self.rules = rules;
        self.child_nodes = child_nodes;
        if self.count_type == CountType::Edge {
            self.child_orientations = child_orientations;
            self.rule_orientation_cache = std::iter::repeat_with(OnceLock::new)
                .take(self.child_orientations.len())
                .collect();
        }
        self.rule_node_cache = std::iter::repeat_with(OnceLock::new)
            .take(self.child_nodes.len())
            .collect();
        for (idx, rule) in order.iter().enumerate() {
            graph_storage.node2rule_id[rule.0 as usize] = idx;
        }
    }
}

fn read_lines_to_chunk<R: BufRead>(reader: &mut R, lines_per_chunk: usize) -> Option<Vec<u8>> {
    let mut chunk = Vec::with_capacity(lines_per_chunk * 128);
    for _ in 0..lines_per_chunk {
        let bytes_read = reader
            .read_until(b'\n', &mut chunk)
            .expect("Failed to read disk");

        if bytes_read == 0 {
            break;
        }
    }
    if chunk.is_empty() {
        None
    } else {
        Some(chunk)
    }
}

#[derive(Debug)]
enum DecompOp<T> {
    Element(T),
    RuleEnd {
        rule_id: usize,
        is_forward: bool,
        output_start: usize,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let grammar = Grammar::from(CountType::Node);
        assert_eq!(grammar.rules.len(), 1);
        assert_eq!(grammar.child_nodes.len(), 0);
    }

    #[test]
    fn test_insert_rule() {
        let mut grammar = Grammar::from(CountType::Node);
        grammar.insert(
            0,
            vec![
                (ItemId(1), Orientation::Forward),
                (ItemId(2), Orientation::Backward),
            ],
        );
        assert_eq!(grammar.rules.len(), 2);
        assert_eq!(grammar.child_nodes.len(), 2);
    }

    #[test]
    fn test_lookup_rule() {
        let mut grammar = Grammar::from(CountType::Node);
        let text = vec![
            (ItemId(1), Orientation::Forward),
            (ItemId(2), Orientation::Backward),
        ];
        grammar.insert(0, text.clone());
        let result = grammar.get_nodes(0);
        let expected = vec![ItemId(1), ItemId(2)];
        assert_eq!(result, &expected);
    }
}
