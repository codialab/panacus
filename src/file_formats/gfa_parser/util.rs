use itertools::Itertools;
use memchr::{memchr, memchr2};
use std::str::{self, FromStr};
use std::time::Instant;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Read},
};

use rayon::prelude::*;

use crate::file_formats::gfa_parser::grammar::Grammar;
use crate::{
    file_formats::gfa_parser::Edge,
    util::{intersects, is_contained, ActiveTable, CountType, IntervalContainer, ItemTable},
};

use super::{abacus::GraphMask, graph::GraphStorage, ItemId, Orientation, PathSegment};

const CHUNK_SIZE: usize = 8192;

pub fn parse_gfa_paths_walks<R: Read>(
    data: &mut BufReader<R>,
    graph_mask: &GraphMask,
    graph_storage: &GraphStorage,
    grammar: &Grammar,
    count: &CountType,
    paths_to_collect: &Vec<PathSegment>,
) -> (
    ItemTable,
    Option<ActiveTable>,
    Option<IntervalContainer>,
    HashMap<PathSegment, (u32, u32)>,
    HashMap<PathSegment, Vec<(ItemId, Orientation)>>,
) {
    log::info!("parsing path + walk sequences");
    let mut item_table = ItemTable::new(graph_storage.path_segments.len());

    let (mut subset_covered_bps, mut exclude_table, include_map, exclude_map) =
        graph_mask.load_optional_subsetting(graph_storage, count);

    let mut collected_paths: HashMap<PathSegment, Vec<(ItemId, Orientation)>> = HashMap::new();
    let mut num_path = 0;
    let complete: Vec<(usize, usize)> = vec![(0, usize::MAX)];
    let mut paths_len: HashMap<PathSegment, (u32, u32)> = HashMap::new();

    // Used to check whether a node was already fully included in a previous path
    let mut fully_included = vec![false; graph_storage.node_count + 1];

    let mut buf = vec![];
    let timer = Instant::now();
    while data.read_until(b'\n', &mut buf).unwrap_or(0) > 0 {
        if buf[0] == b'P' || buf[0] == b'W' {
            let (path_seg, buf_path_seg) = match buf[0] {
                b'P' => parse_path_identifier(&buf),
                b'W' => parse_walk_identifier(&buf),
                _ => unreachable!(),
            };

            log::debug!("processing path {:?}", &path_seg);

            let include_coords = if graph_mask.include_coords.is_none() {
                &complete[..]
            } else {
                match include_map.get(&path_seg.id()) {
                    None => &[],
                    Some(coords) => {
                        log::debug!(
                            "found include coords {:?} for path segment {}",
                            &coords[..],
                            &path_seg.id()
                        );
                        &coords[..]
                    }
                }
            };
            let exclude_coords = if graph_mask.exclude_coords.is_none() {
                &[]
            } else {
                match exclude_map.get(&path_seg.id()) {
                    None => &[],
                    Some(coords) => {
                        log::debug!(
                            "found exclude coords {:?} for path segment {}",
                            &coords[..],
                            &path_seg.id()
                        );
                        &coords[..]
                    }
                }
            };

            let (start, end) = path_seg.coords().unwrap_or((0, usize::MAX));

            // do not process the path sequence if path is neither part of subset nor exclude
            if graph_mask.include_coords.is_some()
                && !intersects(include_coords, &(start, end))
                && !intersects(exclude_coords, &(start, end))
            {
                log::debug!("path {} does not intersect with subset coordinates {:?} nor with exclude coordinates {:?} and therefore is skipped from processing", &path_seg, &include_coords, &exclude_coords);

                // update prefix sum
                item_table.id_prefsum[num_path + 1] += item_table.id_prefsum[num_path];

                num_path += 1;
                buf.clear();
                continue;
            }

            if count != &CountType::Edge
                && (graph_mask.include_coords.is_none()
                    || is_contained(include_coords, &(start, end)))
                && (graph_mask.exclude_coords.is_none()
                    || is_contained(exclude_coords, &(start, end)))
            {
                log::debug!("path {} is fully contained within subset coordinates {:?} and is eligible for full parallel processing", path_seg, include_coords);
                let ex = if exclude_coords.is_empty() {
                    None
                } else {
                    exclude_table.as_mut()
                };
                let should_collect_path = paths_to_collect.iter().any(|p| path_seg.is_part_of(p));
                let (num_added_nodes, bp_len, path) = match buf[0] {
                    b'P' => parse_path_seq_update_tables(
                        buf_path_seg,
                        graph_storage,
                        &mut item_table,
                        ex,
                        num_path,
                        should_collect_path,
                    ),
                    b'W' => parse_walk_seq_update_tables(
                        buf_path_seg,
                        graph_storage,
                        grammar,
                        &mut item_table,
                        ex,
                        num_path,
                        should_collect_path,
                    ),
                    _ => unreachable!(),
                };
                if should_collect_path {
                    log::info!("Full collected: {:?}", path_seg);
                    collected_paths.insert(path_seg.clone(), path);
                }
                paths_len.insert(path_seg, (num_added_nodes, bp_len));
            } else {
                let sids = match buf[0] {
                    b'P' => parse_path_seq_to_item_vec(buf_path_seg, graph_storage),
                    b'W' => parse_walk_seq_to_item_vec(buf_path_seg, graph_storage, grammar),
                    _ => unreachable!(),
                };

                if paths_to_collect.iter().any(|p| path_seg.is_part_of(p)) {
                    log::info!("Collected: {:?}", path_seg);
                    collected_paths.insert(path_seg.clone(), sids.clone());
                }
                match count {
                    CountType::Node | CountType::Bp => {
                        let (node_len, bp_len) = update_tables(
                            &mut item_table,
                            &mut subset_covered_bps.as_mut(),
                            &mut fully_included,
                            &mut exclude_table.as_mut(),
                            num_path,
                            graph_storage,
                            sids,
                            include_coords,
                            exclude_coords,
                            start,
                        );
                        paths_len.insert(path_seg, (node_len as u32, bp_len as u32));
                    }
                    CountType::Edge => update_tables_edgecount(
                        &mut item_table,
                        &mut exclude_table.as_mut(),
                        num_path,
                        graph_storage,
                        sids,
                        include_coords,
                        exclude_coords,
                        start,
                    ),
                };
            }
            num_path += 1;
        }
        buf.clear();
    }
    let duration = timer.elapsed();
    log::info!(
        "func done; count: {:?}; time elapsed: {:?}",
        count,
        duration
    );
    (
        item_table,
        exclude_table,
        subset_covered_bps,
        paths_len,
        collected_paths,
    )
}

pub fn parse_path_seq_update_tables(
    data: &[u8],
    graph_storage: &GraphStorage,
    item_table: &mut ItemTable,
    exclude_table: Option<&mut ActiveTable>,
    num_path: usize,
    should_collect_path: bool,
) -> (u32, u32, Vec<(ItemId, Orientation)>) {
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    log::debug!("parsing path sequences of size {} bytes..", end);

    let (segment_ids, bp_len) = get_path_segment_ids(data, graph_storage, end, CHUNK_SIZE);

    // This function is only called with nodes, so orientation shouldn't matter
    let path = if should_collect_path {
        segment_ids
            .iter()
            .map(|x| (*x, Orientation::Forward))
            .collect()
    } else {
        Vec::new()
    };

    segment_ids.into_iter().for_each(|segment_id| {
        item_table.items.push(segment_id.0);
        item_table.id_prefsum[num_path + 1] += 1;
    });

    // compute prefix sum
    let mut num_nodes_path = 0;
    num_nodes_path += item_table.id_prefsum[num_path + 1];
    item_table.id_prefsum[num_path + 1] += item_table.id_prefsum[num_path];

    // is exclude table is given, we assume that all nodes of the path are excluded
    if let Some(ex) = exclude_table {
        log::debug!("flagging nodes of path as excluded");
        for j in (item_table.id_prefsum[num_path] as usize)
            ..(item_table.id_prefsum[num_path + 1] as usize)
        {
            ex.items[item_table.items[j] as usize] |= true;
        }
    }

    log::debug!("..done");
    (num_nodes_path as u32, bp_len, path)
}

fn get_path_segment_ids(
    data: &[u8],
    graph_storage: &GraphStorage,
    end: usize,
    chunk_size: usize,
) -> (Vec<ItemId>, u32) {
    let (segment_ids, bp_lens): (Vec<_>, Vec<_>) = (0..end)
        .into_par_iter()
        .step_by(chunk_size)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + chunk_size].iter().min().unwrap();
            let mut bp_len: u32 = 0;

            // sits after first comma in chunk
            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr(b',', &data[chunk_start..chunk_end])
                    .map(|v| v + 1)                     // move *after* comma
                    .unwrap_or(chunk_size + 3)          // add enough to chunk_size, so that while
                                                        // loop does not run, if no comma in chunk
                    + x
                } // add offset back
            };
            let mut segment_ids = Vec::new();
            while curr_pos - chunk_start < chunk_size + 1 {
                // sits on comma at the end of the current segment
                let segment_end = match memchr(b',', &data[curr_pos..]) {
                    None => end,
                    Some(idx) => curr_pos + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id = get_segment_id(&data[curr_pos..segment_end], graph_storage);
                bp_len += graph_storage.node_len(&segment_id);
                segment_ids.push(segment_id);
                // move curr_pos forward (after next comma)
                curr_pos = segment_end + 1;
            }
            (segment_ids, bp_len)
        })
        .unzip();

    let segment_ids = segment_ids.into_iter().concat();
    let bp_len = bp_lens.into_iter().sum();

    (segment_ids, bp_len)
}

pub fn parse_walk_seq_update_tables(
    data: &[u8],
    graph_storage: &GraphStorage,
    grammar: &Grammar,
    item_table: &mut ItemTable,
    exclude_table: Option<&mut ActiveTable>,
    num_path: usize,
    should_collect_path: bool,
) -> (u32, u32, Vec<(ItemId, Orientation)>) {
    // later codes assumes that data is non-empty...
    if data.is_empty() {
        return (0, 0, Vec::new());
    }

    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap();

    log::debug!("parsing walk sequences of size {}..", end);

    let (segment_ids, bp_len) = get_walk_segment_ids(data, graph_storage, grammar, end, CHUNK_SIZE);

    let path = if should_collect_path {
        segment_ids
            .iter()
            .map(|x| (*x, Orientation::Forward))
            .collect()
    } else {
        Vec::new()
    };

    let initial_len = item_table.items.len();
    item_table
        .items
        .extend(segment_ids.into_iter().map(|s| s.0));
    let items_added = item_table.items.len() - initial_len;
    item_table.id_prefsum[num_path + 1] += items_added as u64;

    // compute prefix sum
    let mut num_nodes_path = 0;
    num_nodes_path += item_table.id_prefsum[num_path + 1];
    item_table.id_prefsum[num_path + 1] += item_table.id_prefsum[num_path];

    // is exclude table is given, we assume that all nodes of the path are excluded
    if let Some(ex) = exclude_table {
        log::error!("flagging nodes of path as excluded");
        for j in (item_table.id_prefsum[num_path] as usize)
            ..(item_table.id_prefsum[num_path + 1] as usize)
        {
            ex.items[item_table.items[j] as usize] |= true;
        }
    }

    log::debug!("..done");
    (num_nodes_path as u32, bp_len, path)
}

fn get_walk_segment_ids(
    data: &[u8],
    graph_storage: &GraphStorage,
    grammar: &Grammar,
    end: usize,
    chunk_size: usize,
) -> (Vec<ItemId>, u32) {
    let (segment_ids, bp_lens): (Vec<_>, Vec<_>) = (0..end)
        .into_par_iter()
        .step_by(chunk_size)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + chunk_size].iter().min().unwrap();
            let mut bp_len: u32 = 0;

            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr2(b'<', b'>', &data[chunk_start..chunk_end]).unwrap_or(chunk_size + 3) + x
                }
            };

            let mut segment_ids = Vec::new();

            while curr_pos - chunk_start < chunk_size {
                let segment_end = match memchr2(b'<', b'>', &data[curr_pos + 1..]) {
                    None => end,
                    Some(idx) => curr_pos + 1 + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id = get_walk_segment_id_grammar(
                    &data[curr_pos..segment_end],
                    graph_storage,
                    grammar,
                );
                for segment in &segment_id {
                    bp_len += graph_storage.node_len(segment);
                }
                segment_ids.extend(segment_id);
                // move curr_pos forward (after next comma)
                curr_pos = segment_end;
            }
            (segment_ids, bp_len)
        })
        .unzip();

    let segment_ids = segment_ids.into_iter().concat();
    let bp_len = bp_lens.into_iter().sum();
    (segment_ids, bp_len)
}

pub fn parse_walk_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut six_col: Vec<&str> = Vec::with_capacity(6);

    let mut it = data.iter();
    let mut i = 0;
    for _ in 0..6 {
        let j = it.position(|x| x == &b'\t').unwrap();
        six_col.push(str::from_utf8(&data[i..i + j]).unwrap());
        i += j + 1;
    }

    let seq_start = match six_col[4] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let seq_end = match six_col[5] {
        "*" => None,
        a => Some(usize::from_str(a).unwrap()),
    };

    let path_seg = PathSegment::new(
        six_col[1].to_string(),
        six_col[2].to_string(),
        six_col[3].to_string(),
        seq_start,
        seq_end,
    );

    (path_seg, &data[i..])
}

pub fn parse_path_identifier(data: &[u8]) -> (PathSegment, &[u8]) {
    let mut iter = data.iter();

    let start = iter.position(|&x| x == b'\t').unwrap() + 1;
    let offset = iter.position(|&x| x == b'\t').unwrap();
    let path_name = str::from_utf8(&data[start..start + offset]).unwrap();
    (
        PathSegment::from_str(path_name),
        &data[start + offset + 1..],
    )
}

pub fn update_tables(
    item_table: &mut ItemTable,
    subset_covered_bps: &mut Option<&mut IntervalContainer>,
    fully_included: &mut Vec<bool>,
    exclude_table: &mut Option<&mut ActiveTable>,
    num_path: usize,
    graph_storage: &GraphStorage,
    path: Vec<(ItemId, Orientation)>,
    include_coords: &[(usize, usize)],
    exclude_coords: &[(usize, usize)],
    offset: usize,
) -> (usize, usize) {
    let mut i = 0;
    let mut j = 0;
    let mut p = offset;

    let mut included = 0;
    let mut included_bp = 0;
    let mut excluded = 0;

    log::debug!(
        "checking inclusion/exclusion criteria on {} nodes..",
        path.len()
    );
    if path.is_empty() {
        return (included, included_bp);
    }

    for (sid, o) in &path {
        let l = graph_storage.node_len(sid) as usize;

        // this implementation of include coords for bps is *not exact* as illustrated by the
        // following scenario:
        //
        //   subset intervals:           ____________________________
        //                ______________|_____________________________
        //               |
        //      ___________________________________________     ____
        //  ---|                some node                  |---|
        //      -------------------------------------------     ----
        //
        //
        //   what the following code does:
        //                ___________________________________________
        //               |
        //               |             coverage count
        //      ___________________________________________     ____
        //  ---|                some node                  |---|
        //      -------------------------------------------     ----
        //
        //
        // node count handling: node is only counted if *completely* covered by subset
        //
        //
        // update current pointer in include_coords list

        // end is not inclusive, so if end <= p (=offset) then advance to the next interval
        let mut stop_here = false;
        while i < include_coords.len() && include_coords[i].0 < p + l && !stop_here {
            if include_coords[i].1 > p {
                let mut a = include_coords[i].0.saturating_sub(p);
                let mut b = if include_coords[i].1 < p + l {
                    // advance to the next interval
                    i += 1;
                    include_coords[i - 1].1 - p
                } else {
                    stop_here = true;
                    l
                };

                // reverse coverage interval in case of backward orientation
                if o == &Orientation::Backward {
                    (a, b) = (l - b, l - a);
                }

                item_table.items.push(sid.0);
                item_table.id_prefsum[num_path + 1] += 1;
                if let Some(int) = subset_covered_bps.as_mut() {
                    // if fully covered, we do not need to store anything in the map
                    if b - a == l {
                        if int.contains(sid) {
                            int.remove(sid);
                        }
                        fully_included[sid.0 as usize] = true;
                    } else {
                        if !fully_included[sid.0 as usize] {
                            int.add(*sid, a, b);
                        }
                    }
                }
                included += 1;
                included_bp += b - a;
            } else {
                // advance to the next interval
                i += 1;
            }
        }

        let mut stop_here = false;
        while j < exclude_coords.len() && exclude_coords[j].0 < p + l && !stop_here {
            if exclude_coords[j].1 > p {
                let mut a = exclude_coords[j].0.saturating_sub(p);
                let mut b = if exclude_coords[j].1 < p + l {
                    // advance to the next interval for the next iteration
                    j += 1;
                    exclude_coords[j - 1].1 - p
                } else {
                    stop_here = true;
                    l
                };

                // reverse coverage interval in case of backward orientation
                if o == &Orientation::Backward {
                    (a, b) = (l - b, l - a);
                }

                if let Some(map) = exclude_table {
                    if map.with_annotation() {
                        map.activate_n_annotate(*sid, l, a, b)
                            .expect("this error should never occur");
                    } else {
                        map.activate(sid);
                    }
                    excluded += 1;
                }
            } else {
                j += 1;
            }
        }

        if i >= include_coords.len() && j >= exclude_coords.len() {
            // terminate parse if all "include" and "exclude" coords are processed
            break;
        }
        p += l;
    }

    log::debug!(
        "found {} included nodes ({} included bps) and {} excluded nodes, and discarded the rest",
        included,
        included_bp,
        excluded,
    );

    // Compute prefix sum
    item_table.id_prefsum[num_path + 1] += item_table.id_prefsum[num_path];
    log::debug!("..done");
    (included, included_bp)
}

pub fn update_tables_edgecount(
    item_table: &mut ItemTable,
    exclude_table: &mut Option<&mut ActiveTable>,
    num_path: usize,
    graph_storage: &GraphStorage,
    path: Vec<(ItemId, Orientation)>,
    include_coords: &[(usize, usize)],
    exclude_coords: &[(usize, usize)],
    offset: usize,
) {
    let mut i = 0;
    let mut j = 0;
    let mut p = offset;

    // edges are positioned between nodes, offset by the first node
    if !path.is_empty() {
        p += graph_storage.node_len(&path[0].0) as usize;
    }

    log::debug!("checking inclusion/exclusion criteria on {} nodes, inserting successful candidates to corresponding data structures..", path.len());

    for ((sid1, o1), (sid2, o2)) in path.into_iter().tuple_windows() {
        // update current pointer in include_coords list
        while i < include_coords.len() && include_coords[i].1 <= p {
            i += 1;
        }

        // update current pointer in exclude_coords list
        while j < exclude_coords.len() && exclude_coords[j].1 <= p {
            j += 1;
        }

        let l = graph_storage.node_len(&sid2) as usize;

        let e = Edge::canonical(sid1, o1, sid2, o2);
        let eid = graph_storage
            .edge2id
            .as_ref()
            .expect("update_tables_edgecount requires edge2id map in GraphStorage")
            .get(&e)
            .unwrap_or_else(|| {
                panic!(
                    "unknown edge {}. Is flipped edge known? {}",
                    &e,
                    if graph_storage
                        .edge2id
                        .as_ref()
                        .unwrap()
                        .contains_key(&e.flip())
                    {
                        "Yes"
                    } else {
                        "No"
                    }
                )
            });
        // check if the current position fits within active segment
        if i < include_coords.len() && include_coords[i].0 < p + l {
            item_table.items.push(eid.0);
            item_table.id_prefsum[num_path + 1] += 1;
        }
        if exclude_table.is_some() && j < exclude_coords.len() && exclude_coords[j].0 < p + l {
            exclude_table.as_mut().unwrap().activate(eid);
        } else if i >= include_coords.len() && j >= exclude_coords.len() {
            // terminate parse if all "include" and "exclude" coords are processed
            break;
        }
        p += l;
    }
    // Compute prefix sum
    item_table.id_prefsum[num_path + 1] += item_table.id_prefsum[num_path];
    log::debug!("..done");
}

pub fn parse_walk_seq_to_item_vec(
    data: &[u8],
    graph_storage: &GraphStorage,
    grammar: &Grammar,
) -> Vec<(ItemId, Orientation)> {
    // later codes assumes that data is non-empty...
    if data.is_empty() {
        return Vec::new();
    }

    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap_or_else(|| it.len());

    log::debug!("parsing walk sequences of size {}..", end);

    let segment_ids: Vec<_> = (0..end)
        .into_par_iter()
        .step_by(CHUNK_SIZE)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + CHUNK_SIZE].iter().min().unwrap();

            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr2(b'<', b'>', &data[chunk_start..chunk_end]).unwrap_or(CHUNK_SIZE + 3) + x
                }
            };

            let mut segment_ids = Vec::new();

            while curr_pos - chunk_start < CHUNK_SIZE {
                let segment_end = match memchr2(b'<', b'>', &data[curr_pos + 1..]) {
                    None => end,
                    Some(idx) => curr_pos + 1 + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id = get_walk_segment_id(&data[curr_pos..segment_end], graph_storage);
                let orientation = Orientation::from_lg(data[curr_pos]);

                let sids = grammar.decompress_ordered(segment_id, orientation, graph_storage);
                segment_ids.extend(sids);
                // move curr_pos forward (after next comma)
                curr_pos = segment_end;
            }
            segment_ids
        })
        .collect();

    let sids = segment_ids.into_iter().concat();

    log::debug!("..done");
    sids
}

pub fn parse_path_seq_to_item_vec(
    data: &[u8],
    graph_storage: &GraphStorage,
) -> Vec<(ItemId, Orientation)> {
    let mut it = data.iter();
    let end = it
        .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
        .unwrap_or_else(|| it.len());

    log::debug!("parsing path sequences of size {}..", end);

    let segment_ids: Vec<_> = (0..end)
        .into_par_iter()
        .step_by(CHUNK_SIZE)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + CHUNK_SIZE].iter().min().unwrap();

            // sits after first comma in chunk
            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr(b',', &data[chunk_start..chunk_end])
                    .map(|v| v + 1)                     // move *after* comma
                    .unwrap_or(CHUNK_SIZE + 3)          // add enough to chunk_size, so that while
                                                        // loop does not run, if no comma in chunk
                    + x
                } // add offset back
            };
            let mut segment_ids = Vec::new();
            while curr_pos - chunk_start < CHUNK_SIZE + 1 {
                // sits on comma at the end of the current segment
                let segment_end = match memchr(b',', &data[curr_pos..]) {
                    None => end,
                    Some(idx) => curr_pos + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id = get_segment_id(&data[curr_pos..segment_end], graph_storage);
                let orientation = Orientation::from_pm(data[segment_end - 1]);
                segment_ids.push((segment_id, orientation));
                // move curr_pos forward (after next comma)
                curr_pos = segment_end + 1;
            }
            segment_ids
        })
        .collect();

    log::debug!("..done");

    segment_ids.into_iter().concat()
}

fn get_segment_id(node: &[u8], graph_storage: &GraphStorage) -> ItemId {
    let segment_id = graph_storage
        .get_node_id(&node[0..node.len() - 1])
        .unwrap_or_else(|| panic!("unknown node {}", str::from_utf8(node).unwrap()));
    // TODO: Is orientation really necessary?
    let orientation = node[node.len() - 1];
    assert!(
        orientation == b'-' || orientation == b'+',
        "unknown orientation of segment {}",
        str::from_utf8(node).unwrap()
    );
    //plus_strands[rayon::current_thread_index().unwrap()] += (orientation == b'+') as u32;
    segment_id
}

fn get_walk_segment_id(node: &[u8], graph_storage: &GraphStorage) -> ItemId {
    let segment_id = graph_storage
        .get_node_id(&node[1..node.len()])
        .unwrap_or_else(|| panic!("unknown node {}", str::from_utf8(node).unwrap()));
    // TODO: Is orientation really necessary?
    let orientation = node[0];
    assert!(
        orientation == b'<' || orientation == b'>',
        "unknown orientation of segment {}",
        str::from_utf8(node).unwrap()
    );
    //plus_strands[rayon::current_thread_index().unwrap()] += (orientation == b'+') as u32;
    segment_id
}

fn get_walk_segment_id_grammar(
    node: &[u8],
    graph_storage: &GraphStorage,
    grammar: &Grammar,
) -> Vec<ItemId> {
    let segment_id = graph_storage
        .get_node_id(&node[1..node.len()])
        .unwrap_or_else(|| panic!("unknown node {}", str::from_utf8(node).unwrap()));
    // TODO: Is orientation really necessary?
    let orientation = node[0];
    assert!(
        orientation == b'<' || orientation == b'>',
        "unknown orientation of segment {}",
        str::from_utf8(node).unwrap()
    );

    let sids = grammar.decompress_unordered(segment_id, graph_storage);
    //plus_strands[rayon::current_thread_index().unwrap()] += (orientation == b'+') as u32;
    sids
}

#[cfg(test)]
fn get_walk_segment_id_with_orientation(
    node: &[u8],
    graph_storage: &GraphStorage,
) -> (ItemId, Orientation) {
    let segment_id = graph_storage
        .get_node_id(&node[1..node.len()])
        .unwrap_or_else(|| panic!("unknown node {}", str::from_utf8(node).unwrap()));
    // TODO: Is orientation really necessary?
    let orientation = node[0];
    assert!(
        orientation == b'<' || orientation == b'>',
        "unknown orientation of segment {}",
        str::from_utf8(node).unwrap()
    );
    //plus_strands[rayon::current_thread_index().unwrap()] += (orientation == b'+') as u32;
    (segment_id, Orientation::from_lg(orientation))
}

#[cfg(test)]
fn get_walk_segment_ids_with_orientation(
    data: &[u8],
    graph_storage: &GraphStorage,
    end: usize,
    chunk_size: usize,
) -> (Vec<(ItemId, Orientation)>, u32) {
    let (segment_ids, bp_lens): (Vec<_>, Vec<_>) = (0..end)
        .step_by(chunk_size)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + chunk_size].iter().min().unwrap();
            let mut bp_len: u32 = 0;

            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr2(b'<', b'>', &data[chunk_start..chunk_end]).unwrap_or(chunk_size + 3) + x
                }
            };

            let mut segment_ids = Vec::new();

            while curr_pos - chunk_start < chunk_size {
                let segment_end = match memchr2(b'<', b'>', &data[curr_pos + 1..]) {
                    None => end,
                    Some(idx) => curr_pos + 1 + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id = get_walk_segment_id_with_orientation(
                    &data[curr_pos..segment_end],
                    graph_storage,
                );
                bp_len += graph_storage.node_len(&segment_id.0);
                segment_ids.push(segment_id);
                // move curr_pos forward (after next comma)
                curr_pos = segment_end;
            }
            (segment_ids, bp_len)
        })
        .unzip();

    let segment_ids = segment_ids.into_iter().concat();
    let bp_len = bp_lens.into_iter().sum();
    (segment_ids, bp_len)
}

#[cfg(test)]
fn get_segment_id_with_orientation(
    node: &[u8],
    graph_storage: &GraphStorage,
) -> (ItemId, Orientation) {
    let segment_id = graph_storage
        .get_node_id(&node[0..node.len() - 1])
        .unwrap_or_else(|| panic!("unknown node {}", str::from_utf8(node).unwrap()));
    // TODO: Is orientation really necessary?
    let orientation = node[node.len() - 1];
    assert!(
        orientation == b'-' || orientation == b'+',
        "unknown orientation of segment {}",
        str::from_utf8(node).unwrap()
    );
    let orientation = Orientation::from_pm(orientation);
    //plus_strands[rayon::current_thread_index().unwrap()] += (orientation == b'+') as u32;
    (segment_id, orientation)
}

#[cfg(test)]
fn get_path_segment_ids_with_orientation(
    data: &[u8],
    graph_storage: &GraphStorage,
    end: usize,
    chunk_size: usize,
) -> (Vec<(ItemId, Orientation)>, u32) {
    let (segment_ids, bp_lens): (Vec<_>, Vec<_>) = (0..end)
        .into_par_iter()
        .step_by(chunk_size)
        .map(|chunk_start| {
            let chunk_end = *[end, chunk_start + chunk_size].iter().min().unwrap();
            let mut bp_len: u32 = 0;

            // sits after first comma in chunk
            let mut curr_pos = match chunk_start {
                0 => 0,
                x => {
                    memchr(b',', &data[chunk_start..chunk_end])
                    .map(|v| v + 1)                     // move *after* comma
                    .unwrap_or(chunk_size + 3)          // add enough to chunk_size, so that while
                                                        // loop does not run, if no comma in chunk
                    + x
                } // add offset back
            };
            let mut segment_ids = Vec::new();
            while curr_pos - chunk_start < chunk_size + 1 {
                // sits on comma at the end of the current segment
                let segment_end = match memchr(b',', &data[curr_pos..]) {
                    None => end,
                    Some(idx) => curr_pos + idx,
                };
                let segment_end = if segment_end < end { segment_end } else { end };
                if curr_pos >= segment_end {
                    break;
                }
                let segment_id =
                    get_segment_id_with_orientation(&data[curr_pos..segment_end], graph_storage);
                bp_len += graph_storage.node_len(&segment_id.0);
                segment_ids.push(segment_id);
                // move curr_pos forward (after next comma)
                curr_pos = segment_end + 1;
            }
            (segment_ids, bp_len)
        })
        .unzip();

    let segment_ids = segment_ids.into_iter().concat();
    let bp_len = bp_lens.into_iter().sum();

    (segment_ids, bp_len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_path_chunk_sizes() {
        let data = "1+,3+,5+,6+,8+,9+,11+,12+,14+,15+\t8M,1M,1M,3M,1M,19M,1M,4M,1M,11M".as_bytes();
        let mut it = data.iter();
        let end = it
            .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
            .unwrap();
        let (graph_storage, _) = GraphStorage::from_gfa("tests/test_files/t_groups.gfa", true);
        let exp = vec![
            (ItemId(1), Orientation::Forward),
            (ItemId(3), Orientation::Forward),
            (ItemId(5), Orientation::Forward),
            (ItemId(6), Orientation::Forward),
            (ItemId(8), Orientation::Forward),
            (ItemId(9), Orientation::Forward),
            (ItemId(11), Orientation::Forward),
            (ItemId(12), Orientation::Forward),
            (ItemId(14), Orientation::Forward),
            (ItemId(15), Orientation::Forward),
        ];
        for i in 1..35 {
            let (res, _) = get_path_segment_ids_with_orientation(data, &graph_storage, end, i);
            assert_eq!(res, exp);
        }
    }

    #[test]
    fn test_walk_chunk_sizes() {
        let data = ">1>3>5>6>8>9>11>12>14>15\n".as_bytes();
        let mut it = data.iter();
        let end = it
            .position(|x| x == &b'\t' || x == &b'\n' || x == &b'\r')
            .unwrap();
        let (graph_storage, _) = GraphStorage::from_gfa("tests/test_files/t_groups.gfa", true);
        let exp = vec![
            (ItemId(1), Orientation::Forward),
            (ItemId(3), Orientation::Forward),
            (ItemId(5), Orientation::Forward),
            (ItemId(6), Orientation::Forward),
            (ItemId(8), Orientation::Forward),
            (ItemId(9), Orientation::Forward),
            (ItemId(11), Orientation::Forward),
            (ItemId(12), Orientation::Forward),
            (ItemId(14), Orientation::Forward),
            (ItemId(15), Orientation::Forward),
        ];
        for i in 1..35 {
            let (res, _) = get_walk_segment_ids_with_orientation(data, &graph_storage, end, i);
            assert_eq!(res, exp);
        }
    }
}
