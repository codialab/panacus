use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use itertools::{Itertools, MinMaxResult};

use crate::{
    analyses::info::FileInfo,
    file_formats::gfa_parser::SparseMatrix,
    hist::Hist,
    util::{ItemTable, Threshold},
};

type Buckets = HashMap<usize, Vec<(usize, usize)>>;

/// (start_coord, end_coord, list with feature ids)
type ReferenceBucket = Vec<(usize, usize, Vec<usize>)>;

#[derive(Debug)]
pub struct CoverageMatrix {
    count_of_features: usize,
    feature_lengths: Vec<usize>,
    feature_positions: Positions,
    feature_names: Vec<String>,
    path_names: Vec<String>,
    matrix: SparseMatrix,
    feature_type: String,
    run_id: String,
    run_name: String,
    file_info: FileInfo,
}

impl CoverageMatrix {
    // TODO: Add FILE info object
    pub fn new(
        feature_type: String,
        run_id: String,
        run_name: String,
        file_info: FileInfo,
    ) -> Self {
        Self {
            count_of_features: 0,
            feature_lengths: Vec::new(),
            feature_positions: Positions::new(),
            feature_names: Vec::new(),
            path_names: Vec::new(),
            matrix: SparseMatrix::new(),
            feature_type,
            run_id,
            run_name,
            file_info,
        }
    }

    pub fn set_path_names(&mut self, path_names: Vec<String>) {
        self.path_names = path_names;
    }

    pub fn set_file_info(&mut self, file_info: FileInfo) {
        self.file_info = file_info;
    }

    /// Calculates the histogram from the matrix. This should
    /// only be used if the CoverageMatrix is needed anyways.
    /// If it isn't needed generate the hist directly from the
    /// file (using `FileFormatParser`'s `generate_hist()`).
    pub fn get_hist(&self) -> Hist {
        let mut hist = Hist::from_maximum_coverage(
            self.path_names.len(),
            self.feature_type.clone(),
            self.run_id.clone(),
            self.run_name.clone(),
        );
        self.get_feature_counts()
            .into_iter()
            .zip(self.feature_lengths.iter())
            .for_each(|(c, l)| hist.insert_feature_of_coverage_and_length(c, *l));
        hist
    }

    pub fn get_hist_for_reference(&self, reference: &str) -> Hist {
        let mut hist = Hist::from_maximum_coverage(
            self.path_names.len() - 1,
            self.feature_type.clone(),
            self.run_id.clone(),
            self.run_name.clone(),
        );
        let r_idx = self
            .path_names
            .iter()
            .enumerate()
            .filter_map(|(i, p)| if p == reference { Some(i) } else { None })
            .next()
            .expect("Reference is part of paths");
        self.get_feature_counts()
            .into_iter()
            .enumerate()
            .zip(self.feature_lengths.iter())
            .for_each(|((i, c), l)| {
                // Insert the feature only if it is not part of the reference
                if !self.matrix.contains(i, r_idx as u64) {
                    hist.insert_feature_of_coverage_and_length(c, *l)
                }
            });
        hist
    }

    fn get_number_of_windows(
        min: usize,
        max: usize,
        window_size: usize,
        slide_step: usize,
    ) -> usize {
        if max == min + 1 {
            1
        } else {
            ((max - window_size - min + 1) as f64 / slide_step as f64).ceil() as usize + 1
        }
    }

    fn get_feature_buckets(
        &self,
        ref_id: usize,
        window_size: usize,
        slide_step: usize,
    ) -> (usize, Vec<Vec<usize>>) {
        let (min, mut max) = match self.feature_positions.get_pos_iter_ref_id(ref_id).minmax() {
            MinMaxResult::NoElements => unimplemented!("Should return empty iterator"),
            MinMaxResult::OneElement(x) => (x, x + 1),
            MinMaxResult::MinMax(min, max) => (min, max),
        };
        if max == min {
            max += 1;
        }
        let number_of_windows = Self::get_number_of_windows(min, max, window_size, slide_step);

        let mut feature_buckets: Vec<Vec<usize>> = vec![Vec::new(); number_of_windows];
        for (feature_idx, position) in self.feature_positions.get_idpos_iter_ref_id(ref_id) {
            let start_window =
                ((position - window_size - min + 1) as f64 / slide_step as f64).ceil() as usize;
            let end_window = ((position - min) as f64 / slide_step as f64).floor() as usize;
            for window in feature_buckets
                .iter_mut()
                .take(end_window + 1)
                .skip(start_window)
            {
                window.push(feature_idx);
            }
        }
        (min, feature_buckets)
    }

    // Looks at all the features and sorts them into buckets based on the buckets data structure
    fn get_feature_buckets_for_given_buckets(
        &self,
        ref_id: usize,
        buckets: &Buckets,
    ) -> Option<(usize, Vec<Vec<usize>>)> {
        let bucket = buckets.get(&ref_id)?;
        let mut feature_buckets: Vec<Vec<usize>> = vec![Vec::new(); bucket.len()];
        for (feature_idx, position) in self.feature_positions.get_idpos_iter_ref_id(ref_id) {
            // This is the first window to definitely not contain the position, so:
            // [0, last_window_idx) is our valid interval
            let last_window_idx = bucket.partition_point(|w| w.0 <= position);

            let mut first_window_idx = if last_window_idx > 0 {
                last_window_idx - 1
            } else {
                last_window_idx
            };

            // Step back until the end coordinate is before our position
            while first_window_idx > 0 && bucket[first_window_idx].1 >= position {
                first_window_idx -= 1;
            }
            // If we went outside, we need to step back in
            if bucket[first_window_idx].1 < position {
                first_window_idx += 1;
            }

            // Our valid interval is now [first_window_idx, last_window_idx)
            for window in feature_buckets
                .iter_mut()
                .take(last_window_idx)
                .skip(first_window_idx)
            {
                window.push(feature_idx);
            }
        }
        Some((0, feature_buckets))
    }

    // Opens a bed file and reads it into the Buckets data structure
    fn get_buckets_from_bucket_file(&self, file_path: &str) -> Buckets {
        let file = match File::open(file_path) {
            Ok(f) => f,
            Err(e) => {
                log::error!("Could not open file {} ({})", file_path, e);
                return HashMap::new();
            }
        };
        let reader = BufReader::new(file);
        let reference_lookup: HashMap<&str, usize> = self
            .feature_positions
            .references
            .iter()
            .enumerate()
            .map(|(id, ref_name)| (&ref_name[..], id))
            .collect();
        let mut buckets: Buckets = HashMap::new();
        for line in reader.lines() {
            let line = match line {
                Ok(l) => l,
                Err(e) => {
                    log::error!("Could not read line in file {} ({})", file_path, e);
                    continue;
                }
            };
            let mut iter = line.split("\t");
            let (Some(reference), Some(start), Some(end), None) =
                (iter.next(), iter.next(), iter.next(), iter.next())
            else {
                log::error!(
                    "Line {} of file {} does not contain exactly three tab-separated fields",
                    line,
                    file_path
                );
                continue;
            };
            let Some(ref_id) = reference_lookup.get(reference) else {
                log::error!(
                    "Reference {} used inside file {} is not part of the reference of the input file",
                    reference,
                    file_path
                );
                continue;
            };
            let (Ok(start), Ok(end)) = (start.parse::<usize>(), end.parse::<usize>()) else {
                log::error!("Start and/or end coordinates ({}, {}) in line {} in file {} are not positive integers", start, end, line, file_path);
                continue;
            };
            buckets.entry(*ref_id).or_default().push((start, end));
        }

        // Sort the windows
        // Since sorting is implemented for tuples, Rust correctly uses start
        // as the first priority and then end
        for (_, windows) in buckets.iter_mut() {
            windows.sort_unstable();
        }
        buckets
    }

    /// Takes a window size and by what number that window should step forward (often the same number)
    /// Returns an iterator over references by name and all of the starts, ends and hists for the reference
    pub fn get_regional_hists(
        &self,
        window_size: usize,
        slide_step: usize,
        window_file: Option<&str>,
    ) -> impl Iterator<Item = (String, impl Iterator<Item = (usize, usize, Hist)> + '_)> + '_ {
        let bed_buckets: Option<HashMap<usize, Vec<(usize, usize)>>> =
            window_file.map(|file_path| self.get_buckets_from_bucket_file(file_path));
        let mut all_references_buckets: Vec<(String, ReferenceBucket)> = Vec::new();
        for (ref_id, reference) in self.feature_positions.references.iter().enumerate() {
            let (min, feature_buckets) = match bed_buckets.as_ref() {
                Some(b) => {
                    if let Some(res) = self.get_feature_buckets_for_given_buckets(ref_id, b) {
                        res
                    } else {
                        continue;
                    }
                }
                None => self.get_feature_buckets(ref_id, window_size, slide_step),
            };
            let annotated_feature_buckets: Vec<(usize, usize, Vec<usize>)> = feature_buckets
                .into_iter()
                .enumerate()
                .map(|(idx, bucket)| {
                    let (start, end) = if let Some(bed_buckets) = bed_buckets.as_ref() {
                        bed_buckets[&ref_id][idx]
                    } else {
                        (min + idx * slide_step, min + idx * slide_step + window_size)
                    };
                    (start, end, bucket)
                })
                .collect();
            all_references_buckets.push((reference.to_string(), annotated_feature_buckets));
        }
        all_references_buckets
            .into_iter()
            .map(move |(reference, buckets)| {
                (
                    reference,
                    buckets.into_iter().map(move |(start, end, bucket)| {
                        let hist = self.get_hist_for_features(&bucket);
                        (start, end, hist)
                    }),
                )
            })
    }

    pub fn get_hist_for_features(&self, features: &[usize]) -> Hist {
        let mut hist = Hist::from_maximum_coverage(
            self.path_names.len(),
            self.get_feature_type().to_owned(),
            self.get_run_id().to_owned(),
            self.get_run_name().to_owned(),
        );
        for &feature in features {
            // TODO: this might be inefficient if windows are highly overlapping
            // (if this function is called from get_regional_hists)
            let count = self.get_count_of_feature(feature);
            hist.insert_feature_of_coverage_and_length(count, self.feature_lengths[feature]);
        }
        hist
    }

    pub fn insert_feature(
        &mut self,
        feature_name: String,
        feature_length: usize,
        feature_position: (&str, usize),
        feature: Vec<u32>,
    ) {
        // Check if we have the same number of entries as we have paths
        assert_eq!(feature.len(), self.path_names.len());
        self.count_of_features += 1;
        self.feature_names.push(feature_name);
        self.feature_lengths.push(feature_length);
        self.feature_positions
            .insert(feature_position.0, feature_position.1);
        self.matrix.insert_row(feature);
    }

    pub fn insert_item_table(
        &mut self,
        path_names: Vec<String>,
        feature_lengths: Vec<usize>,
        feature_positions: Positions,
        feature_names: Vec<String>,
        item_table: ItemTable,
    ) {
        self.path_names = path_names;
        self.feature_lengths = feature_lengths;
        self.count_of_features = self.feature_lengths.len();
        self.feature_names = feature_names;
        self.feature_positions = feature_positions;
        self.matrix
            .insert_item_table(self.count_of_features, item_table);
    }

    /// Returns a histogram and the features
    /// that really appear in the histogram (i.e. are non-zero).
    /// If features is Some(list), the list is used to subset the features
    pub fn get_hist_and_used_features(
        &self,
        features: Option<&Vec<usize>>,
        paths: &[usize],
    ) -> (Hist, Vec<usize>) {
        let c = Threshold::Absolute(1);
        self.get_hist_and_used_features_with_coverage(features, paths, &c)
    }

    pub fn get_hist_and_used_features_with_coverage(
        &self,
        features: Option<&Vec<usize>>,
        paths: &[usize],
        c: &Threshold,
    ) -> (Hist, Vec<usize>) {
        let (abacus, used_features) =
            self.get_abacus_and_used_features_with_coverage(features, paths, c);
        let mut hist = Hist::from_maximum_coverage(
            paths.len(),
            self.feature_type.to_string(),
            self.run_id.to_string(),
            self.run_name.to_string(),
        );
        for (i, c) in abacus.iter().enumerate() {
            hist.insert_feature_of_coverage_and_length(*c, self.feature_lengths[i]);
        }
        (hist, used_features)
    }

    /// Gets a CSC (compressed-sparse-column) table, i.e. a table
    /// containing the occuring features by path (as opposed to by
    /// feature)
    pub fn get_csc(&self) -> (Vec<usize>, Vec<usize>) {
        self.matrix.get_csc(self.path_names.len())
    }

    pub fn get_csr(&self) -> (Vec<usize>, Vec<usize>) {
        self.matrix.get_csr()
    }

    pub fn get_abacus_and_used_features_with_coverage(
        &self,
        features: Option<&Vec<usize>>,
        paths: &[usize],
        c: &Threshold,
    ) -> (Vec<usize>, Vec<usize>) {
        let all_features = (0..self.count_of_features).collect();
        let features = match features {
            Some(f) => f,
            None => &all_features,
        };
        let mut paths = paths.to_vec();
        paths.sort();
        let mut non_zeroes = Vec::new();
        let mut abacus = vec![0; self.count_of_features];
        for feature in features {
            let mut path_idx = 0;
            let mut count = 0;
            let occurrences = self.matrix.get_occurrences(*feature);
            if occurrences.len() < c.to_absolute(self.count_of_features) {
                continue;
            }
            // log::info!(
            //     "Looking at feature: {} | {}",
            //     *feature,
            //     occurrences.len(),
            // );
            for o in occurrences {
                let o = *o as usize;
                while path_idx < paths.len() && paths[path_idx] < o {
                    path_idx += 1;
                }
                if path_idx >= paths.len() {
                    break;
                } else if o < paths[path_idx] {
                    continue;
                } else if o == paths[path_idx] {
                    count += 1;
                    path_idx += 1;
                }
            }
            // log::info!("Count: {}", count);
            // if count >= c.to_absolute(self.count_of_features) {
            non_zeroes.push(*feature);
            abacus[*feature] = count;
            // }
        }
        (abacus, non_zeroes)
    }

    pub fn get_feature_name(&self, feature: usize) -> &str {
        &self.feature_names[feature]
    }

    /// Creates an iterator over indices in the order (in order),
    /// containing only the indices for which there exists an element
    pub fn get_appearances_for_order(&self, feature: usize, order: &[String]) -> Vec<usize> {
        let translation_table: HashMap<&String, usize> = self
            .path_names
            .iter()
            .enumerate()
            .map(|(idx, path_name)| (path_name, idx))
            .collect();
        order
            .iter()
            .enumerate()
            .filter_map(|(idx, path_name)| {
                let translated_idx = translation_table[path_name];
                if self.matrix.contains(feature, translated_idx as u64) {
                    Some(idx)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn get_feature_counts(&self) -> Vec<usize> {
        (0..self.feature_lengths.len())
            .map(|i| self.matrix.get_feature_occurrence_count(i))
            .collect()
    }

    pub fn get_feature_names(&self) -> Vec<String> {
        self.feature_names.clone()
    }

    pub fn get_feature_lengths(&self) -> &Vec<usize> {
        &self.feature_lengths
    }

    pub fn get_count_of_feature(&self, feature: usize) -> usize {
        self.matrix.get_feature_occurrence_count(feature)
    }

    pub fn get_counts_for_feature(&self, id: usize) -> Vec<usize> {
        self.matrix
            .get_counts_for_feature(id, self.path_names.len())
    }

    pub fn get_counts_for_feature_in_order(&self, id: usize, order: &[String]) -> Vec<usize> {
        let translation_table: HashMap<&String, usize> = self
            .path_names
            .iter()
            .enumerate()
            .map(|(idx, path_name)| (path_name, idx))
            .collect();
        order
            .iter()
            .map(|path_name| {
                let translated_idx = translation_table[path_name];
                if self.matrix.contains(id, translated_idx as u64) {
                    1
                } else {
                    0
                }
            })
            .collect()
    }

    pub fn get_path_names(&self) -> &Vec<String> {
        &self.path_names
    }

    pub fn get_feature_count(&self) -> usize {
        self.count_of_features
    }

    pub fn get_feature_type(&self) -> &str {
        &self.feature_type
    }

    pub fn get_run_id(&self) -> &str {
        &self.run_id
    }

    pub fn get_run_name(&self) -> &str {
        &self.run_name
    }

    pub fn get_file_info(&self) -> &FileInfo {
        &self.file_info
    }
}

#[derive(Debug)]
pub struct Positions {
    references: Vec<String>,
    reference_lookup: HashMap<String, u8>,
    feature_refs: Vec<u8>,
    feature_positions: Vec<usize>,
}

impl Default for Positions {
    fn default() -> Self {
        Self::new()
    }
}

impl Positions {
    pub fn new() -> Self {
        Self {
            references: Vec::new(),
            reference_lookup: HashMap::new(),
            feature_refs: Vec::new(),
            feature_positions: Vec::new(),
        }
    }

    pub fn with_size(size: usize) -> Self {
        Self {
            references: Vec::new(),
            reference_lookup: HashMap::new(),
            feature_refs: vec![u8::MAX; size],
            feature_positions: vec![0; size],
        }
    }

    pub fn set(&mut self, idx: usize, reference: &str, position: usize) {
        if let Some(&id) = self.reference_lookup.get(reference) {
            self.feature_refs[idx] = id;
            self.feature_positions[idx] = position;
        } else {
            self.references.push(reference.to_string());
            let id = (self.references.len() - 1) as u8;
            self.reference_lookup.insert(reference.to_string(), id);
            self.feature_refs[idx] = id;
            self.feature_positions[idx] = position;
        }
    }

    pub fn get(&self, id: usize) -> (&str, usize) {
        let ref_id = self.feature_refs[id];
        let reference = if ref_id != u8::MAX {
            &self.references[ref_id as usize]
        } else {
            "DEFAULT"
        };
        let pos = self.feature_positions[id];
        (reference, pos)
    }

    pub fn insert(&mut self, reference: &str, position: usize) {
        if let Some(&id) = self.reference_lookup.get(reference) {
            self.feature_refs.push(id);
            self.feature_positions.push(position);
        } else {
            self.references.push(reference.to_string());
            let id = (self.references.len() - 1) as u8;
            self.reference_lookup.insert(reference.to_string(), id);
            self.feature_refs.push(id);
            self.feature_positions.push(position);
        }
    }

    pub fn apply_mask(&mut self, mask: &[bool]) {
        let mut iter = mask.iter();
        self.feature_refs.retain(|_| *iter.next().unwrap());
        let mut iter = mask.iter();
        self.feature_positions.retain(|_| *iter.next().unwrap());
    }

    pub fn get_pos_iter_ref_id(&self, ref_id: usize) -> impl Iterator<Item = usize> + '_ {
        self.feature_positions
            .iter()
            .enumerate()
            .filter_map(move |(idx, &pos)| {
                if self.feature_refs[idx] == ref_id as u8 {
                    Some(pos)
                } else {
                    None
                }
            })
    }

    /// This function cleans up invalid (i.e. unset) positions
    pub fn cleanup(&mut self) {
        if self.feature_refs.contains(&u8::MAX) {
            self.references.push("DEFAULT".to_string());
            let id = (self.references.len() - 1) as u8;
            self.reference_lookup.insert("DEFAULT".to_string(), id);
            self.feature_refs.iter_mut().for_each(|x| {
                if *x == u8::MAX {
                    *x = id
                }
            });
        }
    }

    /// Gets an iterator over tuples with feature index and position
    /// for all features that are part of the reference with the given
    /// id.
    pub fn get_idpos_iter_ref_id(
        &self,
        ref_id: usize,
    ) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.feature_positions
            .iter()
            .enumerate()
            .filter_map(move |(idx, &pos)| {
                if self.feature_refs[idx] == ref_id as u8 {
                    Some((idx, pos))
                } else {
                    None
                }
            })
    }
}
