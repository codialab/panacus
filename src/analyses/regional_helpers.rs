use std::{
    collections::{HashMap, HashSet, VecDeque},
    fs::File,
    io::{BufReader, Read},
    iter::{self},
};

use crate::{
    file_formats::gfa_parser::{Edge, ItemId, Orientation, PathSegment},
    html_report,
};

pub fn sort_values(
    file: &str,
    reference_names: &mut Vec<String>,
    windows: &mut Vec<Vec<html_report::Window>>,
) {
    let f = File::open(file).expect("Order file must exist");
    let mut reader = BufReader::new(f);
    let mut buf = String::new();
    reader
        .read_to_string(&mut buf)
        .expect("Order file must be readable");
    let order_names: Vec<_> = buf.lines().collect();
    let lookup: Vec<usize> = reference_names
        .iter()
        .map(|r| {
            order_names
                .iter()
                .position(|&o| o == r)
                .expect("Ordering contains all references")
        })
        .collect();
    reorder_in_place(reference_names, lookup.clone());
    reorder_in_place(windows, lookup.clone());
}

fn reorder_in_place<T>(x: &mut [T], mut indices: Vec<usize>) {
    for i in 0..x.len() {
        while indices[i] != i {
            let target = indices[i];

            x.swap(i, target);
            indices.swap(i, target);
        }
    }
}

pub fn get_close_nodes(
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

// pub fn get_close_nodes_with_edges(
//     ref_nodes: &Vec<ItemId>,
//     neighbors: &HashMap<ItemId, Vec<ItemId>>,
//     edges: &HashMap<Edge, ItemId>,
// ) -> HashMap<ItemId, Vec<ItemId>> {
//     let mut queue: VecDeque<ItemId> = VecDeque::new();
//     let mut closest_ref: HashMap<ItemId, (ItemId, ItemId)> = HashMap::new();
//     let ref_node_set: HashSet<ItemId> = ref_nodes.iter().copied().collect();
//     for ref_node in ref_nodes {
//         closest_ref.insert(*ref_node, *ref_node);
//         for neighbor in neighbors.get(ref_node).unwrap_or(&Vec::new()) {
//             if !ref_node_set.contains(neighbor) && !closest_ref.contains_key(neighbor) {
//                 closest_ref.insert(*neighbor, *ref_node);
//                 queue.push_front(*neighbor);
//             }
//         }
//     }
//     while let Some(current_node) = queue.pop_back() {
//         for neighbor in neighbors.get(&current_node).unwrap_or(&Vec::new()) {
//             if !closest_ref.contains_key(neighbor) {
//                 closest_ref.insert(*neighbor, current_node);
//                 queue.push_front(*neighbor);
//             }
//         }
//     }
//     let mut flipped: HashMap<ItemId, Vec<ItemId>> = HashMap::new();
//     for (entry_node, ref_node) in closest_ref {
//         flipped.entry(ref_node).or_default().push(entry_node);
//     }
//     flipped
// }

pub fn get_ref_length(ref_nodes: &Vec<(ItemId, Orientation)>, node_lens: &Vec<u32>) -> u32 {
    let length = ref_nodes
        .iter()
        .map(|(node, _)| node_lens[node.0 as usize])
        .sum();
    length
}

pub fn split_ref_paths(
    ref_paths: HashMap<PathSegment, &Vec<(ItemId, Orientation)>>,
) -> HashMap<PathSegment, HashMap<PathSegment, &Vec<(ItemId, Orientation)>>> {
    let mut result: HashMap<PathSegment, HashMap<PathSegment, &Vec<(ItemId, Orientation)>>> =
        HashMap::new();
    for (path, path_content) in ref_paths {
        let coord_less_path = path.clear_coords();
        result
            .entry(coord_less_path)
            .or_default()
            .insert(path, path_content);
    }
    result
}

fn get_window_grid(
    window_size: usize,
    contig_start: usize,
    ref_length: usize,
) -> Vec<(usize, usize)> {
    let no_of_windows = (contig_start + ref_length) / window_size + 1;
    let window_starts = (0..no_of_windows).map(|x| x * window_size);
    let mut windows: Vec<(usize, usize)> = window_starts
        .into_iter()
        .map(|s| (s, s + window_size))
        .filter(|(s, e)| *e > contig_start && *s <= contig_start + ref_length) // Filter out windows before and after the contig
        .collect();

    // Set the start of the first window to the start of the contig
    windows[0] = (contig_start, windows[0].1);
    // Set the end of the last window to the end of the contig
    let last_index = windows.len() - 1;
    let start_of_last = windows[last_index].0;
    windows[last_index] = (start_of_last, contig_start + ref_length);
    windows
}

fn get_window_ranges(
    window_size: usize,
    contig_start: usize,
    ref_length: usize,
    allowed_segments: &[PathSegment],
) -> Vec<(usize, usize)> {
    let windows = get_window_grid(window_size, contig_start, ref_length);

    if allowed_segments.is_empty() {
        return windows;
    }

    // Overlap prospective windows with allowed_segments
    let mut allowed_segments = allowed_segments.to_vec();
    allowed_segments.sort();
    let allowed_segments: Vec<(usize, usize)> = allowed_segments
        .into_iter()
        .map(|s| {
            (
                s.start.expect("Ref subset has start"),
                s.end.expect("Ref subset has end"),
            )
        })
        .collect();

    let mut wi = 0usize;
    let mut si = 0usize;
    let mut new_windows: Vec<(usize, usize)> = Vec::new();
    while wi < windows.len() && si < allowed_segments.len() {
        let window = windows[wi];
        let segment = allowed_segments[si];

        // No overlaps
        if window.1 < segment.0 {
            wi += 1;
            continue;
        } else if segment.1 < window.0 {
            si += 1;
            continue;
        }

        // Case:
        // --win--         or      ---win---
        //     --seg--              -seg-
        if window.0 <= segment.0 {
            if segment.1 > window.1 {
                new_windows.push((segment.0, window.1));
                wi += 1;
            } else {
                new_windows.push((segment.0, segment.1));
                si += 1;
            }
        }
        // Case:
        //     --win--      or        -win-
        //  --seg--                ---seg---
        else {
            if window.1 > segment.1 {
                new_windows.push((window.0, segment.1));
                si += 1;
            } else {
                new_windows.push((window.0, window.1));
                wi += 1;
            }
        }
    }
    new_windows
}

fn cover_windows_with_reference_nodes(
    window_ranges: Vec<(usize, usize)>,
    ref_nodes: &[(ItemId, Orientation)],
    node_lens: &[u32],
    contig_start: usize,
) -> Vec<Window> {
    let mut current_window_index = 0;
    let mut current_node_index = 0;
    let mut already_used_bps_of_node = 0;
    let mut ref_windows: Vec<Window> = window_ranges
        .into_iter()
        .map(|(s, e)| (Vec::new(), s, e))
        .collect();

    let mut bp_pos = contig_start;

    // While last node is not yet used
    while current_node_index < ref_nodes.len() {
        // Quit if there are still nodes left, but no window
        if current_window_index >= ref_windows.len() {
            break;
        }

        // Check if end of current node is before the start of current window
        if (bp_pos + node_lens[ref_nodes[current_node_index].0 .0 as usize] as usize)
            < ref_windows[current_window_index].1
        {
            bp_pos += node_lens[ref_nodes[current_node_index].0 .0 as usize] as usize;
            current_node_index += 1;
            already_used_bps_of_node = 0;
            continue;
        }

        // Length the window should have
        let window_length =
            ref_windows[current_window_index].2 - ref_windows[current_window_index].1;

        // Length the window currently has
        let current_window_length: usize = ref_windows[current_window_index]
            .0
            .iter()
            .map(|(_, _, l, _)| l)
            .sum();

        // If this is the first node of the window, all the bps that are before the start of the
        // window have already been used
        if current_window_length == 0 {
            already_used_bps_of_node = ref_windows[current_window_index].1 - bp_pos;
        }

        let remaining_bps_in_window = window_length - current_window_length;
        let remaining_node_length = node_lens[ref_nodes[current_node_index].0 .0 as usize] as usize
            - already_used_bps_of_node;

        // Remainder of node fits inside remainder of window
        if remaining_node_length <= remaining_bps_in_window {
            ref_windows[current_window_index].0.push((
                ref_nodes[current_node_index].0,
                ref_nodes[current_node_index].1,
                remaining_node_length,
                // either both or back
                if current_window_length == 0 && already_used_bps_of_node > 0 {
                    IncludedEnds::Back
                } else {
                    IncludedEnds::Both
                },
            ));
            bp_pos += node_lens[ref_nodes[current_node_index].0 .0 as usize] as usize;
            current_node_index += 1;
            already_used_bps_of_node = 0;
            if remaining_node_length == remaining_bps_in_window {
                current_window_index += 1;
            }
        } else {
            let cut_node_length = remaining_bps_in_window;
            ref_windows[current_window_index].0.push((
                ref_nodes[current_node_index].0,
                ref_nodes[current_node_index].1,
                cut_node_length, // either none or front
                if current_window_length == 0 && already_used_bps_of_node > 0 {
                    IncludedEnds::None
                } else {
                    IncludedEnds::Front
                },
            ));
            current_window_index += 1;
            already_used_bps_of_node += cut_node_length;
        };
    }

    ref_windows
}

fn augment_windows_with_non_reference_nodes(
    ref_nodes: &Vec<(ItemId, Orientation)>,
    ref_windows: Vec<Window>,
    neighbors: &HashMap<(ItemId, Orientation), HashSet<ItemId>>,
) -> Vec<(Vec<(ItemId, usize)>, usize, usize)> {
    let close_nodes = get_close_nodes(ref_nodes, neighbors);
    let windows: Vec<(Vec<(ItemId, usize)>, usize, usize)> = ref_windows
        .into_iter()
        .map(|(window, start, end)| {
            (
                window
                    .into_iter()
                    .flat_map(|(ref_node, ref_orientation, ref_length, included_ends)| {
                        match included_ends {
                            IncludedEnds::None => vec![(ref_node, ref_length)],
                            IncludedEnds::Back => close_nodes
                                .get(&(ref_node, ref_orientation))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, 0))
                                .chain(iter::once((ref_node, ref_length)))
                                .collect(),
                            IncludedEnds::Front => close_nodes
                                .get(&(ref_node, ref_orientation.flip()))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, 0))
                                .chain(iter::once((ref_node, ref_length)))
                                .collect(),
                            IncludedEnds::Both => close_nodes
                                .get(&(ref_node, ref_orientation))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, 0))
                                .chain(
                                    close_nodes
                                        .get(&(ref_node, ref_orientation.flip()))
                                        .unwrap_or(&Vec::new())
                                        .iter()
                                        .map(|node| (*node, 0)),
                                )
                                .chain(iter::once((ref_node, ref_length)))
                                .collect(),
                        }
                    })
                    .collect::<Vec<_>>(),
                start,
                end,
            )
        })
        .collect();
    windows
}

type Window = (
    Vec<(ItemId, Orientation, usize, IncludedEnds)>,
    usize,
    usize,
);

fn merge_small_windows(
    ref_windows: Vec<Window>,
    window_size: usize,
    contig_start: usize,
    ref_length: usize,
) -> Vec<Window> {
    let window_grid = get_window_grid(window_size, contig_start, ref_length);
    let mut window_grid: Vec<Window> = window_grid
        .into_iter()
        .map(|(s, e)| (Vec::new(), s, e))
        .collect();
    for window in ref_windows {
        for grid_window in window_grid.iter_mut() {
            if window.1 >= grid_window.1 && window.2 <= grid_window.2 {
                grid_window.0.extend(window.0);
                break;
            }
        }
    }
    window_grid
}

pub fn get_windows(
    ref_nodes: &Vec<(ItemId, Orientation)>,
    node_lens: &Vec<u32>,
    window_size: usize,
    neighbors: &HashMap<(ItemId, Orientation), HashSet<ItemId>>,
    contig_start: usize,
    _log_windows: bool,
    allowed_segments: &[PathSegment],
    should_merge_small_windows: bool,
) -> Vec<(Vec<(ItemId, usize)>, usize, usize)> {
    let ref_length = get_ref_length(ref_nodes, node_lens) as usize;

    // Get window ranges, i.e., (start, end)-coordinates
    let window_ranges: Vec<(usize, usize)> =
        get_window_ranges(window_size, contig_start, ref_length, allowed_segments);

    // Cover the windows with reference nodes
    let ref_windows: Vec<Window> =
        cover_windows_with_reference_nodes(window_ranges, ref_nodes, node_lens, contig_start);

    // If the flag is set merge windows that are inside the same window
    let ref_windows = if should_merge_small_windows {
        merge_small_windows(ref_windows, window_size, contig_start, ref_length)
    } else {
        ref_windows
    };

    // Add close non-reference nodes to windows

    augment_windows_with_non_reference_nodes(ref_nodes, ref_windows, neighbors)
}

type DirectedNodeBucket = (Vec<(ItemId, IncludedEnds)>, usize, usize);

pub fn get_edge_windows(
    ref_nodes: &Vec<(ItemId, Orientation)>,
    node_lens: &Vec<u32>,
    window_size: usize,
    neighbors: &HashMap<(ItemId, Orientation), HashSet<(ItemId, Orientation)>>,
    edge2id: &HashMap<Edge, ItemId>,
    contig_start: usize,
    allowed_segments: &[PathSegment],
    should_merge_small_windows: bool,
) -> Vec<(Vec<ItemId>, usize, usize)> {
    let close_nodes = get_close_nodes(
        ref_nodes,
        &neighbors
            .iter()
            .map(|(idx, others)| (*idx, others.iter().map(|(a, _)| *a).collect()))
            .collect(),
    );
    let ref_length = get_ref_length(ref_nodes, node_lens) as usize;

    // Get window ranges, i.e., (start, end)-coordinates
    let window_ranges: Vec<(usize, usize)> =
        get_window_ranges(window_size, contig_start, ref_length, allowed_segments);

    // Cover the windows with reference nodes
    let ref_windows: Vec<Window> =
        cover_windows_with_reference_nodes(window_ranges, ref_nodes, node_lens, contig_start);

    // if true {
    //     ref_windows.iter().for_each(|(v, _s, _e)| {
    //         eprintln!("Nodes: {}", v.len());
    //         //if v.len() == 1 {
    //         //    eprintln!("S {}-{} contains only a single reference node", s, e);
    //         //}
    //     });
    // }

    // If the flag is set merge windows that are inside the same window
    let ref_windows = if should_merge_small_windows {
        merge_small_windows(ref_windows, window_size, contig_start, ref_length)
    } else {
        ref_windows
    };

    let full_windows_with_end: Vec<DirectedNodeBucket> = ref_windows
        .into_iter()
        .map(|(window, start, end)| {
            (
                window
                    .into_iter()
                    .flat_map(|(ref_node, ref_orientation, _ref_length, included_ends)| {
                        match included_ends {
                            IncludedEnds::None => vec![(ref_node, IncludedEnds::None)],
                            IncludedEnds::Back => close_nodes
                                .get(&(ref_node, ref_orientation))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, IncludedEnds::Both))
                                .chain(iter::once((
                                    ref_node,
                                    IncludedEnds::get_back_from_orientation(ref_orientation),
                                )))
                                .collect(),
                            IncludedEnds::Front => close_nodes
                                .get(&(ref_node, ref_orientation.flip()))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, IncludedEnds::Both))
                                .chain(iter::once((
                                    ref_node,
                                    IncludedEnds::get_front_from_orientation(ref_orientation),
                                )))
                                .collect(),
                            IncludedEnds::Both => close_nodes
                                .get(&(ref_node, ref_orientation))
                                .unwrap_or(&Vec::new())
                                .iter()
                                .map(|node| (*node, IncludedEnds::Both))
                                .chain(
                                    close_nodes
                                        .get(&(ref_node, ref_orientation.flip()))
                                        .unwrap_or(&Vec::new())
                                        .iter()
                                        .map(|node| (*node, IncludedEnds::Both)),
                                )
                                .chain(iter::once((ref_node, IncludedEnds::Both)))
                                .collect(),
                        }
                    })
                    .collect::<Vec<_>>(),
                start,
                end,
            )
        })
        .collect::<Vec<_>>();
    let windows: Vec<(Vec<ItemId>, usize, usize)> = full_windows_with_end
        .into_iter()
        .map(|(window, start, end)| {
            let nodes_in_window = window.iter().map(|(node, _)| *node).collect::<HashSet<_>>();
            (
                window
                    .into_iter()
                    .flat_map(|(node, included_ends)| {
                        let empty = HashSet::new();
                        let neighbor_iter = match included_ends {
                            IncludedEnds::None => {
                                Box::new(empty.iter()) as Box<dyn Iterator<Item = &_>>
                            }
                            IncludedEnds::Back => Box::new(
                                neighbors
                                    .get(&(node, Orientation::Forward))
                                    .unwrap_or(&empty)
                                    .iter(),
                            ),
                            IncludedEnds::Front => Box::new(
                                neighbors
                                    .get(&(node, Orientation::Backward))
                                    .unwrap_or(&empty)
                                    .iter(),
                            ),
                            IncludedEnds::Both => Box::new(
                                neighbors
                                    .get(&(node, Orientation::Forward))
                                    .unwrap_or(&empty)
                                    .iter()
                                    .chain(
                                        neighbors
                                            .get(&(node, Orientation::Backward))
                                            .unwrap_or(&empty)
                                            .iter(),
                                    ),
                            ),
                        };
                        neighbor_iter
                            .flat_map(|node| {
                                neighbors
                                    .get(node)
                                    .unwrap_or(&HashSet::new())
                                    .iter()
                                    .filter_map(|other| {
                                        if nodes_in_window.contains(&other.0) {
                                            edge2id
                                                .get(&Edge(node.0, node.1, other.0, other.1))
                                                .copied()
                                        } else {
                                            None
                                        }
                                    })
                                    .collect::<Vec<ItemId>>()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>(),
                start,
                end,
            )
        })
        .collect();
    windows
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum IncludedEnds {
    None,
    Front,
    Back,
    Both,
}

impl IncludedEnds {
    fn get_back_from_orientation(o: Orientation) -> Self {
        match o {
            Orientation::Forward => IncludedEnds::Back,
            Orientation::Backward => IncludedEnds::Front,
        }
    }

    fn get_front_from_orientation(o: Orientation) -> Self {
        match o {
            Orientation::Forward => IncludedEnds::Front,
            Orientation::Backward => IncludedEnds::Back,
        }
    }
}
