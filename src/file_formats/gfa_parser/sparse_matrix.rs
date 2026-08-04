use itertools::Itertools;

use crate::util::{CountSize, GroupSize, ItemTable};

#[derive(Debug, Clone)]
pub struct SparseMatrix {
    r: Vec<usize>,
    v: Option<Vec<CountSize>>,
    c: Vec<GroupSize>,
}

impl Default for SparseMatrix {
    fn default() -> Self {
        Self::new()
    }
}

impl SparseMatrix {
    pub fn new() -> Self {
        Self {
            r: vec![0],
            c: Vec::new(),
            v: Some(Vec::new()),
        }
    }

    pub fn get_value(&self, feature: usize, path: GroupSize) -> CountSize {
        let start = self.r[feature];
        let end = self.r[feature + 1];
        for i in start..end {
            if path == self.c[i] {
                return self.v.as_ref().expect("Has value")[i];
            }
        }
        0
    }

    pub fn get_occurrences(&self, feature: usize) -> &[u64] {
        let start = self.r[feature];
        let end = self.r[feature + 1];
        &self.c[start..end]
    }

    pub fn get_feature_count(&self, feature: usize) -> CountSize {
        let mut result = 0;
        let start = self.r[feature];
        let end = self.r[feature + 1];
        for i in start..end {
            result += self.v.as_ref().expect("Has value")[i];
        }
        result
    }

    /// Transforms the contents into sparse column format
    /// Used if column-wise access is necessary
    pub fn get_csc(&self, n_paths: usize) -> (Vec<usize>, Vec<usize>) {
        let non_zeroes = self.c.len();
        let mut c = vec![0; n_paths + 1];
        let mut r = vec![0; non_zeroes];

        // Count the number of non-zero entries for each column
        // (stored with an offset of one for the prefix sum)
        for &col in self.c.iter() {
            c[col as usize + 1] += 1;
        }

        // Prefix sum
        for i in 0..n_paths {
            c[i + 1] += c[i];
        }

        // Scatter-stage
        // Stores where we are currently inserting elements for a path
        // (is one element longer than it has to be, but removing it would
        // be more runtime than the memory cost of having it)
        let mut current_col_pointers = c.clone();
        for (row, (&start, &end)) in self.r.iter().tuple_windows().enumerate() {
            for &path in &self.c[start..end] {
                let destination = current_col_pointers[path as usize];
                r[destination] = row;
                // Increase the position we write to for this path
                current_col_pointers[path as usize] += 1;
            }
        }

        (c, r)
    }

    pub fn get_csr(&self) -> (Vec<usize>, Vec<usize>) {
        (self.r.clone(), self.c.iter().map(|&x| x as usize).collect())
    }

    pub fn get_feature_occurrence_count(&self, feature: usize) -> usize {
        self.r[feature + 1] - self.r[feature]
    }

    pub fn contains(&self, feature: usize, path: GroupSize) -> bool {
        for i in self.r[feature]..self.r[feature + 1] {
            if self.c[i] == path {
                return true;
            }
        }
        false
    }

    pub fn insert_item_table(&mut self, number_of_features: usize, item_table: ItemTable) {
        log::info!(
            "Inserting item table with {} features, max feature: {}",
            number_of_features,
            item_table.items.iter().max().unwrap()
        );
        self.r = compute_row_storage_space(&item_table, number_of_features);
        (self.v, self.c) = compute_column_values(&item_table, &self.r, true);
    }

    pub fn get_counts_for_feature(&self, id: usize, num_paths: usize) -> Vec<usize> {
        let (start, end) = (self.r[id], self.r[id + 1]);
        let mut res = vec![0; num_paths];
        for i in start..end {
            let idx = self.c[i] as usize;
            let value = match &self.v {
                Some(v) => v[i] as usize,
                None => 1,
            };
            res[idx] = value;
        }
        res
    }

    pub fn insert_row(&mut self, row: Vec<u32>) {
        let mut number_of_inserted = 0;
        for (pos, value) in row.into_iter().enumerate() {
            if value > 0 {
                self.c.push(pos as GroupSize);
                self.v.as_mut().expect("has value").push(value);
                number_of_inserted += 1;
            }
        }
        self.r
            .push(number_of_inserted + self.r.last().expect("r has at least an initial zero"));
    }
}

fn compute_column_values(
    item_table: &ItemTable,
    r: &[usize],
    report_values: bool,
) -> (Option<Vec<CountSize>>, Vec<GroupSize>) {
    let n = { *r.last().unwrap() };
    log::info!("allocating storage for group-based coverage table..");
    let mut v = if report_values {
        vec![0; n]
    } else {
        // we produce a dummy
        vec![0; 1]
    };
    let mut c: Vec<GroupSize> = vec![GroupSize::MAX; n];
    log::info!("done");

    log::info!("computing group-based coverage..");

    // group id is monotone increasing from 0 to #groups
    for path_id in 0..(item_table.id_prefsum.len() - 1) {
        let start = item_table.id_prefsum[path_id] as usize;
        let end = item_table.id_prefsum[path_id + 1] as usize;
        for j in start..end {
            let sid = item_table.items[j] as usize - 1;
            let cv_start = r[sid];
            let mut cv_end = r[sid + 1];
            if cv_end != cv_start {
                // look up storage location for node cur_sid: we use the last position
                // of interval cv_start..cv_end, which is associated to coverage counts
                // of the current node (sid), in the "c" array as pointer to the
                // current column (group) / value (coverage) position. If the current group
                // id does not match the one associated with the current position, we move
                // on to the next. If cv_start + p == cv_end - 1, this means that we are
                // currently writing the last element in that interval, and we need to make
                // sure that we are no longer using it as pointer.
                if cv_end - 1 > c.len() {
                    log::error!(
                        "oops, cv_end-1 is larger than the length of c for sid={}",
                        sid
                    );
                    cv_end = c.len() - 1;
                }

                let mut p = c[cv_end - 1] as usize;
                // we  look at an untouched interval, so let's get the pointer game
                // started...
                if c[cv_end - 1] == GroupSize::MAX {
                    c[cv_start] = path_id as GroupSize;
                    // if it's just a single value in this interval, the pointer game
                    // ends before it started
                    if cv_start < cv_end - 1 {
                        c[cv_end - 1] = 0;
                    }
                    if report_values {
                        v[cv_start] += 1;
                    }
                } else if cv_start + p < cv_end - 1 {
                    // if group id of current slot does not match current group id
                    // (remember group id's are strictly monotically increasing), then
                    // move on to the next slot
                    if c[cv_start + p] < path_id as GroupSize {
                        // move on to the next slot
                        c[cv_end - 1] += 1;
                        // update local pointer
                        p += 1;
                        c[cv_start + p] = path_id as GroupSize
                    }
                    if report_values {
                        v[cv_start + p] += 1;
                    }
                } else if report_values {
                    // make sure it points to the last element and not beyond
                    v[cv_end - 1] += 1;
                }
            }
        }
    }
    log::info!("done");
    (if report_values { Some(v) } else { None }, c)
}

fn compute_row_storage_space(item_table: &ItemTable, n_items: usize) -> Vec<usize> {
    log::info!("computing space allocating storage for group-based coverage table:");
    let mut last: Vec<GroupSize> = vec![GroupSize::MAX; n_items + 1];

    let mut r: Vec<usize> = vec![0; n_items + 1];
    for path_id in 0..(item_table.id_prefsum.len() - 1) {
        let start = item_table.id_prefsum[path_id] as usize;
        let end = item_table.id_prefsum[path_id + 1] as usize;
        for j in start..end {
            let sid = item_table.items[j] as usize - 1;
            if last[sid] != path_id as GroupSize {
                r[sid] += 1;
                last[sid] = path_id as GroupSize;
            }
        }
    }
    log::info!(" ++ assigning storage locations");
    let mut c = 0;
    // can this be simplified?
    for item in &mut r {
        let tmp = *item;
        *item = c;
        c += tmp;
    }
    log::info!(
        " ++ group-aware table has {} non-zero elements",
        r.last().unwrap()
    );
    r
}
