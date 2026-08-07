use std::ops::Index;

#[derive(Debug, Clone, Default)]
pub struct Hist {
    count_of_features: usize,
    hist: Vec<usize>,
    feature_type: String,
    run_id: String,
    run_name: String,
}

impl Hist {
    /// Creates a histogram with a given maximum coverage.
    pub fn from_maximum_coverage(
        maximum_coverage: usize,
        feature_type: String,
        run_id: String,
        run_name: String,
    ) -> Self {
        Self {
            count_of_features: 0,
            hist: vec![0; maximum_coverage + 1],
            feature_type,
            run_id,
            run_name,
        }
    }

    /// Hold out setting maximum coverage for file types
    /// where we only know the maximum coverage at the first
    /// element (e.g. VCF files where we need to parse the
    /// first variant to know the max coverage)
    pub fn without_maximum_coverage(
        feature_type: String,
        run_id: String,
        run_name: String,
    ) -> Self {
        Self {
            count_of_features: 0,
            hist: Vec::new(),
            feature_type,
            run_id,
            run_name,
        }
    }

    /// Increase the size of the internal histogram,
    /// to be used in combination with the without_maximum_coverage()
    /// constructor.
    pub fn set_maximum_coverage(&mut self, maximum_coverage: usize) {
        if maximum_coverage + 1 < self.hist.len() {
            log::error!("Truncating histogram!");
        }
        self.hist.resize(maximum_coverage + 1, 0);
    }

    #[cfg(test)]
    pub fn from_values(values: Vec<usize>) -> Self {
        Self {
            count_of_features: 0,
            hist: values,
            feature_type: "Test".to_string(),
            run_id: "".to_string(),
            run_name: "".to_string(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.hist.is_empty() || self.hist.iter().all(|&v| v == 0)
    }

    /// Subtracts correction from h[1], used to correct
    /// for an increase in singletons often seen in Histograms
    /// that can be a problem when fitting curves to the Growth
    /// curve / F_new curve
    pub fn correct_h1(&mut self, correction: usize) {
        self.hist[1] -= correction;
    }

    pub fn get_index(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.hist.len()
    }

    /// Inserts a single feature with a given coverage into the histogram.
    /// Will panic if the coverage is bigger than the maximum coverage.
    pub fn insert_feature_of_coverage(&mut self, coverage: usize) {
        self.hist[coverage] += 1;
        self.count_of_features += 1;
    }

    pub fn insert_feature_of_coverage_and_length(&mut self, coverage: usize, length: usize) {
        self.hist[coverage] += length;
        self.count_of_features += length;
    }

    /// Gets the number of features that are currently inserted into
    /// the histogram (this is multiplied by the length of the feature)
    pub fn get_number_of_features(&self) -> usize {
        self.count_of_features
    }

    /// Gets the maximum allowed coverage for this histogram.
    #[inline]
    pub fn get_maximum_coverage(&self) -> usize {
        self.hist.len() - 1
    }

    pub fn get_hist_values(&self) -> &[usize] {
        &self.hist
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
}

impl Index<usize> for Hist {
    type Output = usize;

    fn index(&self, index: usize) -> &Self::Output {
        &self.hist[index]
    }
}
