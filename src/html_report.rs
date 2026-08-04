use base64::engine::general_purpose::STANDARD;
use std::ffi::OsStr;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::{collections::HashMap, str::from_utf8};
use std::{f64, fmt};

use base64::{engine::general_purpose, Engine};
use handlebars::{to_json, Handlebars, RenderError};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use time::{macros::format_description, OffsetDateTime};

use crate::util::{get_default_plot_downloads, to_id};
use shadow_rs::shadow;

shadow!(build);

type JsVars = HashMap<String, HashMap<String, String>>;
type RenderedHTML = Result<(String, JsVars), RenderError>;

pub const BOOTSTRAP_COLOR_MODES_JS: &[u8] = include_bytes!("../etc/color-modes.min.js");
pub const BOOTSTRAP_CSS: &[u8] = include_bytes!("../etc/bootstrap.min.css");
pub const BOOTSTRAP_JS: &[u8] = include_bytes!("../etc/bootstrap.bundle.min.js");
pub const CUSTOM_CSS: &[u8] = include_bytes!("../etc/custom.css");
pub const CUSTOM_LIB_JS: &[u8] = include_bytes!("../etc/lib.js");
pub const HOOK_AFTER_JS: &[u8] = include_bytes!("../etc/hook_after.js");
pub const PANACUS_LOGO: &[u8] = include_bytes!("../etc/panacus-illustration-small.png");
pub const SYMBOLS_SVG: &[u8] = include_bytes!("../etc/symbols.svg");
pub const VEGA: &[u8] = include_bytes!("../etc/vega@6.0.0.min.js");
pub const VEGA_EMBED: &[u8] = include_bytes!("../etc/vega-embed@6.29.0.min.js");
pub const VEGA_LITE: &[u8] = include_bytes!("../etc/vega-lite@6.1.0.min.js");

pub const REPORT_HBS: &[u8] = include_bytes!("../hbs/report.hbs");
pub const BAR_HBS: &[u8] = include_bytes!("../hbs/bar.hbs");
pub const TREE_HBS: &[u8] = include_bytes!("../hbs/tree.hbs");
pub const TABLE_HBS: &[u8] = include_bytes!("../hbs/table.hbs");
pub const HEATMAP_HBS: &[u8] = include_bytes!("../hbs/heatmap.hbs");
pub const CHROMOSOMAL_HBS: &[u8] = include_bytes!("../hbs/chromosomal.hbs");
pub const ANALYSIS_TAB_HBS: &[u8] = include_bytes!("../hbs/analysis_tab.hbs");
pub const REPORT_CONTENT_HBS: &[u8] = include_bytes!("../hbs/report_content.hbs");
pub const HEXBIN_HBS: &[u8] = include_bytes!("../hbs/hexbin.hbs");
pub const LINE_HBS: &[u8] = include_bytes!("../hbs/line.hbs");
pub const PNG_HBS: &[u8] = include_bytes!("../hbs/png.hbs");
pub const SVG_HBS: &[u8] = include_bytes!("../hbs/svg.hbs");
pub const PDF_HBS: &[u8] = include_bytes!("../hbs/pdf.hbs");

fn combine_vars(mut a: JsVars, b: JsVars) -> JsVars {
    for (k, v) in b {
        if let Some(x) = a.get_mut(&k) {
            x.extend(v);
        }
    }
    a
}

#[derive(Serialize, Deserialize, Debug)]
pub struct AnalysisSection {
    pub analysis: String,
    pub run_name: String,
    pub run_id: String,
    pub countable: String,
    pub items: Vec<ReportItem>,
    pub id: String,
    pub table: Option<String>,
    pub plot_downloads: Vec<(String, String)>,
}

impl AnalysisSection {
    fn into_html(self, registry: &mut Handlebars) -> RenderedHTML {
        if !registry.has_template("analysis_tab") {
            registry
                .register_template_string("analysis_tab", from_utf8(ANALYSIS_TAB_HBS).unwrap())?;
        }
        let plots = if self.items.len() > 1 {
            self.items
                .iter()
                .map(|item| HashMap::from([("id", item.get_id()), ("name", item.get_name())]))
                .collect()
        } else if self.items.is_empty() {
            vec![]
        } else {
            vec![HashMap::from([
                ("id", self.items[0].get_id()),
                ("name", "".to_string()),
            ])]
        };
        let items = self
            .items
            .into_iter()
            .map(|x| x.into_html(registry))
            .collect::<Result<Vec<_>, _>>()?;
        let (items, mut js_objects): (Vec<_>, Vec<_>) = items.into_iter().unzip();
        if let Some(table) = &self.table {
            if !js_objects.is_empty() {
                js_objects[0].insert(
                    "tables".to_string(),
                    HashMap::from([(self.id.clone(), table.clone())]),
                );
            }
        }
        let js_objects = js_objects
            .into_iter()
            .reduce(combine_vars)
            .unwrap_or_default();
        let plot_downloads: Vec<HashMap<&str, String>> = self
            .plot_downloads
            .iter()
            .map(|(format, text)| {
                HashMap::from([("type", format.to_owned()), ("text", text.to_owned())])
            })
            .collect();
        let vars = HashMap::from([
            ("id", to_json(&self.id)),
            ("analysis", to_json(&self.analysis)),
            ("run_name", to_json(&self.run_name)),
            ("run_id", to_json(&self.run_id)),
            ("countable", to_json(&self.countable)),
            ("has_table", to_json(self.table.is_some())),
            ("has_graph", to_json(!self.plot_downloads.is_empty())),
            (
                "has_multiple_plot_types",
                to_json(self.plot_downloads.len() > 1),
            ),
            ("plot_type", to_json(plot_downloads)),
            ("plot", to_json(plots)),
            ("items", to_json(items)),
        ]);
        Ok((registry.render("analysis_tab", &vars)?, js_objects))
    }

    pub fn generate_custom_section(name: String, file: String) -> anyhow::Result<Vec<Self>> {
        let id = name.to_lowercase().replace([' ', '|', '\\'], "-");
        let id = format!("custom-{id}");
        let mut table: Option<String> = None;
        let mut plot_downloads = Vec::new();
        let report_item = match get_extension_from_filename(&file) {
            Some("svg") => {
                plot_downloads = vec![("svg".to_string(), "Download as svg".to_string())];
                ReportItem::Svg {
                    id: format!("svg-{id}"),
                    file,
                }
            }
            Some("png") => {
                plot_downloads = vec![("png".to_string(), "Download as png".to_string())];
                ReportItem::Png {
                    id: format!("png-{id}"),
                    file,
                }
            }
            Some("json") => {
                plot_downloads = get_default_plot_downloads();
                ReportItem::Json {
                    id: format!("json-{id}"),
                    file,
                }
            }
            Some(t @ "csv") | Some(t @ "tsv") => {
                let f = File::open(&file)?;
                let mut reader = BufReader::new(f);
                let mut buffer = String::new();
                reader.read_to_string(&mut buffer)?;
                table = Some(format!("`{}`", buffer));
                let split_char = if t == "csv" { "," } else { "\t" };
                let mut lines = buffer.lines();
                let header = lines
                    .next()
                    .unwrap_or_else(|| {
                        panic!("{} file {} should contain at least one line", t, file)
                    })
                    .split(split_char)
                    .map(|x| x.trim().to_owned())
                    .collect();
                let values = lines
                    .map(|l| {
                        l.split(split_char)
                            .map(|x| x.trim().to_owned())
                            .collect::<Vec<String>>()
                    })
                    .collect();
                ReportItem::Table {
                    id: format!("{t}-{id}"),
                    header,
                    values,
                }
            }
            Some("pdf") => ReportItem::Pdf {
                id: format!("pdf-{id}"),
                file,
            },
            _ => unimplemented!("Other formats have not been implemented yet"),
        };
        Ok(vec![AnalysisSection {
            id,
            analysis: "Custom".to_string(),
            run_name: "My Custom Run".to_string(),
            run_id: "MyCustomRun".to_string(),
            countable: name,
            table,
            items: vec![report_item],
            plot_downloads,
        }])
    }
}

fn get_extension_from_filename(filename: &str) -> Option<&str> {
    Path::new(filename).extension().and_then(OsStr::to_str)
}

fn get_js_objects_string(objects: JsVars) -> String {
    let mut res = String::from("{");
    for (k, v) in objects {
        res.push('"');
        res.push_str(&k);
        res.push_str("\": {");
        for (subkey, subvalue) in v {
            res.push('"');
            res.push_str(&subkey);
            res.push_str("\": ");
            res.push_str(&subvalue);
            res.push_str(", ");
        }
        res.push_str("}, ");
    }
    res.push('}');
    res
}

impl AnalysisSection {
    pub fn generate_report(
        sections: Vec<Self>,
        registry: &mut Handlebars,
        filename: &str,
        config: &str,
    ) -> Result<String, RenderError> {
        if !registry.has_template("report") {
            registry.register_template_string("report", from_utf8(REPORT_HBS).unwrap())?;
        }

        let tree = Self::get_tree(&sections, registry)?;

        let (content, js_objects) = Self::generate_report_content(sections, registry, config)?;
        let mut vars = Self::get_variables();
        vars.insert("content", content);
        vars.insert("data_hook", get_js_objects_string(js_objects));
        vars.insert("fname", filename.to_string());
        vars.insert("tree", tree);
        registry.render("report", &vars)
    }

    fn get_tree(sections: &Vec<Self>, registry: &mut Handlebars) -> Result<String, RenderError> {
        let analysis_names = sections.iter().map(|x| x.analysis.clone()).unique();
        let mut analyses = Vec::new();
        for analysis_name in analysis_names {
            let run_ids = sections
                .iter()
                .filter(|x| x.analysis == analysis_name)
                .map(|x| (x.run_id.clone(), x.run_name.clone()))
                .unique();
            let analysis_sections = sections
                .iter()
                .filter(|x| x.analysis == analysis_name)
                .collect::<Vec<_>>();
            let mut runs = Vec::new();
            for (run_id, run_name) in run_ids {
                let run_sections = analysis_sections
                    .iter()
                    .filter(|x| x.run_id == run_id)
                    .collect::<Vec<_>>();
                if run_sections.is_empty() {
                    continue;
                }
                let mut countables = Vec::new();
                for section in &run_sections {
                    let content = HashMap::from([
                        ("title", to_json(&section.countable)),
                        ("id", to_json(to_id(&section.countable))),
                        ("href", to_json(&section.id)),
                    ]);
                    countables.push(to_json(content));
                }
                let run_id = run_sections
                    .first()
                    .expect("Run section has at least one run")
                    .run_id
                    .clone();
                let content = HashMap::from([
                    ("title", to_json(&run_name)),
                    ("id", to_json(to_id(&run_id))),
                    ("countables", to_json(countables)),
                ]);
                runs.push(to_json(content));
            }
            let content = HashMap::from([
                ("title", to_json(&analysis_name)),
                ("id", to_json(to_id(&analysis_name))),
                ("icon", to_json("icon-id")),
                ("runs", to_json(runs)),
            ]);
            analyses.push(to_json(content));
        }

        let mut vars = HashMap::from([("analyses", to_json(analyses))]);
        //let hash = option_env!("GIT_HASH").unwrap_or("nogit");
        // let hash = build::COMMIT_HASH;
        // let version = env!("CARGO_PKG_VERSION");
        let version_text = build::VERSION;
        // let version_text = format!("v{version}-{hash}");
        vars.insert("version", to_json(version_text));
        let now = OffsetDateTime::now_utc();
        vars.insert(
            "timestamp",
            to_json(
                now.format(&format_description!(
                    "[year]-[month]-[day]T[hour]:[minute]:[second]Z"
                ))
                .unwrap(),
            ),
        );
        if !registry.has_template("tree") {
            registry.register_template_string("tree", from_utf8(TREE_HBS).unwrap())?;
        }
        let tree = registry.render("tree", &vars)?;
        Ok(tree)
    }

    fn get_variables() -> HashMap<&'static str, String> {
        let mut vars = HashMap::new();
        vars.insert(
            "bootstrap_color_modes_js",
            String::from_utf8_lossy(BOOTSTRAP_COLOR_MODES_JS).into_owned(),
        );
        vars.insert(
            "bootstrap_css",
            String::from_utf8_lossy(BOOTSTRAP_CSS).into_owned(),
        );
        vars.insert(
            "bootstrap_js",
            String::from_utf8_lossy(BOOTSTRAP_JS).into_owned(),
        );
        vars.insert("vega", String::from_utf8_lossy(VEGA).into_owned());
        vars.insert(
            "vega_embed",
            String::from_utf8_lossy(VEGA_EMBED).into_owned(),
        );
        vars.insert("vega_lite", String::from_utf8_lossy(VEGA_LITE).into_owned());
        vars.insert(
            "custom_css",
            String::from_utf8_lossy(CUSTOM_CSS).into_owned(),
        );
        vars.insert(
            "custom_lib_js",
            String::from_utf8_lossy(CUSTOM_LIB_JS).into_owned(),
        );
        vars.insert(
            "hook_after_js",
            String::from_utf8_lossy(HOOK_AFTER_JS).into_owned(),
        );
        vars.insert(
            "panacus_logo",
            general_purpose::STANDARD_NO_PAD.encode(PANACUS_LOGO),
        );
        vars.insert(
            "symbols_svg",
            String::from_utf8_lossy(SYMBOLS_SVG).into_owned(),
        );
        vars
    }

    fn generate_report_content(
        sections: Vec<Self>,
        registry: &mut Handlebars,
        config: &str,
    ) -> RenderedHTML {
        if !registry.has_template("report_content") {
            registry.register_template_string(
                "report_content",
                from_utf8(REPORT_CONTENT_HBS).unwrap(),
            )?;
        }
        let mut js_objects = Vec::new();
        let sections = sections
            .into_iter()
            .map(|s| {
                let (content, js_object) = s.into_html(registry).unwrap();
                js_objects.push(js_object);
                content
            })
            .collect::<Vec<String>>();
        let text = registry.render("report_content", &sections)?;
        let mut js_objects = js_objects
            .into_iter()
            .reduce(combine_vars)
            .expect("Report needs to contain at least one item");
        let config_content = format!("`{}`", config);
        js_objects.insert(
            "config".to_string(),
            HashMap::from([("first".to_string(), config_content)]),
        );
        Ok((text, js_objects))
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum ReportItem {
    Bar {
        id: String,
        name: String,
        x_label: String,
        y_label: String,
        labels: Vec<String>,
        values: Vec<f64>,
        log_toggle: bool,
    },
    MultiBar {
        id: String,
        names: Vec<String>,
        x_label: String,
        y_label: String,
        labels: Vec<String>,
        values: Vec<Vec<f64>>,
        log_toggle: bool,
        curve: Option<Vec<f64>>,
        alpha: Option<f64>,
    },
    SectionLine {
        id: String,
        names: Vec<String>,
        x_label: String,
        y_label: String,
        labels: Vec<String>,
        values: Vec<Vec<f64>>,
        section_separators: Vec<f64>,
        section_labels: Vec<String>,
        log_toggle: bool,
    },
    Table {
        id: String,
        header: Vec<String>,
        values: Vec<Vec<String>>,
    },
    Hexbin {
        id: String,
        bins: Vec<Bin>,
        threshold: usize,
    },
    Heatmap {
        id: String,
        name: String,
        x_labels: Vec<String>,
        y_labels: Vec<String>,
        values: Vec<Vec<f32>>,
    },
    Line {
        id: String,
        name: String,
        x_label: String,
        y_label: String,
        x_values: Vec<f32>,
        y_values: Vec<f32>,
        log_x: bool,
        log_y: bool,
    },
    Png {
        id: String,
        file: String,
    },
    Svg {
        id: String,
        file: String,
    },
    Json {
        id: String,
        file: String,
    },
    Pdf {
        id: String,
        file: String,
    },
    Chromosomal {
        id: String,
        name: String,
        labels: Vec<String>,
        is_diverging: bool,
        contains_outliers: bool,
        // Each chromosomal plot might contain multiple chromosomes
        sequence: Vec<String>,
        values: Vec<Vec<Window>>,
    },
}

impl ReportItem {
    fn get_id(&self) -> String {
        match self {
            Self::Bar { id, .. } => id.to_string(),
            Self::MultiBar { id, .. } => id.to_string(),
            Self::SectionLine { id, .. } => id.to_string(),
            Self::Table { id, .. } => id.to_string(),
            Self::Heatmap { id, .. } => id.to_string(),
            Self::Hexbin { id, .. } => id.to_string(),
            Self::Line { id, .. } => id.to_string(),
            Self::Png { id, .. } => id.to_string(),
            Self::Svg { id, .. } => id.to_string(),
            Self::Json { id, .. } => id.to_string(),
            Self::Pdf { id, .. } => id.to_string(),
            Self::Chromosomal { id, .. } => id.to_string(),
        }
    }

    fn get_name(&self) -> String {
        match self {
            Self::Bar { name, .. } => name.to_string(),
            Self::MultiBar { .. } => "MultiBar".to_string(),
            Self::SectionLine { .. } => "SectionLine".to_string(),
            Self::Table { .. } => "Table".to_string(),
            Self::Heatmap { name, .. } => name.to_string(),
            Self::Hexbin { .. } => "Hexbin".to_string(),
            Self::Line { name, .. } => name.to_string(),
            Self::Png { .. } => "Png".to_string(),
            Self::Svg { .. } => "Svg".to_string(),
            Self::Json { .. } => "Json".to_string(),
            Self::Pdf { .. } => "Pdf".to_string(),
            Self::Chromosomal { .. } => "Chromosomal".to_string(),
        }
    }

    fn into_html(self, registry: &mut Handlebars) -> RenderedHTML {
        match self {
            Self::Table { id, header, values } => {
                if !registry.has_template("table") {
                    registry.register_template_string("table", from_utf8(TABLE_HBS).unwrap())?;
                }
                let data = HashMap::from([
                    ("id".to_string(), to_json(id)),
                    ("header".to_string(), to_json(header)),
                    ("values".to_string(), to_json(values)),
                ]);
                Ok((
                    registry.render("table", &data)?,
                    HashMap::from([("datasets".to_string(), HashMap::new())]),
                ))
            }
            Self::Heatmap {
                id,
                name,
                x_labels,
                y_labels,
                values,
            } => {
                if !registry.has_template("heatmap") {
                    registry
                        .register_template_string("heatmap", from_utf8(HEATMAP_HBS).unwrap())?;
                }
                let mut data_set = "{ 'values': [".to_string();
                for (row_i, row) in values.iter().enumerate() {
                    for (col_i, cell) in row.iter().enumerate() {
                        data_set.push_str(&format!(
                            "{{ x: '{}', y: '{}', value: {} }},",
                            x_labels[row_i], y_labels[col_i], cell
                        ));
                    }
                }
                data_set.push_str("]}");
                let js_object = format!("new Heatmap('{}', '{}', {})", id, name, data_set,);
                let max_scale = format!(
                    "{:.2}",
                    values
                        .iter()
                        .flatten()
                        .fold(f32::INFINITY, |a, &b| a.min(b))
                );
                let data = HashMap::from([
                    ("id".to_string(), to_json(&id)),
                    ("max".to_string(), to_json(max_scale)),
                ]);
                Ok((
                    registry.render("heatmap", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Chromosomal {
                id,
                name,
                labels,
                is_diverging,
                contains_outliers,
                sequence,
                values,
            } => {
                if !registry.has_template("chromosomal") {
                    registry.register_template_string(
                        "chromosomal",
                        from_utf8(CHROMOSOMAL_HBS).unwrap(),
                    )?;
                }

                let mut values = values;
                for j in 0..values.len() {
                    for i in 0..values[j].len() {
                        if i == values[j].len() - 1 {
                            continue;
                        }

                        // if two following values are continuous
                        // introduce overlap in plot (avoids gaps due to rounding errors)
                        if values[j][i].end == values[j][i + 1].start {
                            values[j].get_mut(i).expect("values has value").end =
                                values[j][i + 1].end;
                        }
                    }
                }

                let data: Vec<String> = values
                    .into_iter()
                    .map(|r| {
                        let inner = r
                            .into_iter()
                            .map(|w| {
                                let mut text = format!("{{'x': {}, 'x2': {}", w.start, w.end,);
                                for (i, label) in labels.iter().enumerate() {
                                    text.push_str(&format!(", '{}': {}", label, w.values[i]));
                                }
                                text.push('}');
                                text
                            })
                            .collect::<Vec<String>>()
                            .join(",");
                        let text = format!("[{inner}]");
                        text
                    })
                    .collect();
                let mut data_text = "{'values': [".to_string();
                for datum in data {
                    data_text.push_str(&datum);
                    data_text.push_str(", ");
                }
                data_text.push_str("]}");
                let sequences = sequence
                    .into_iter()
                    .map(|x| format!("'{x}'"))
                    .collect::<Vec<String>>()
                    .join(", ");
                let sequences = format!("[{sequences}]");
                let js_object = format!(
                    "new Chromosomal('{id}', '{name}', {:?}, {is_diverging}, {contains_outliers}, {sequences}, {data_text})", labels
                );
                let data = HashMap::from([("id".to_string(), to_json(&id))]);
                Ok((
                    registry.render("chromosomal", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Bar {
                id,
                name,
                x_label,
                y_label,
                labels,
                values,
                log_toggle,
            } => {
                if !registry.has_template("bar") {
                    registry.register_template_string("bar", from_utf8(BAR_HBS).unwrap())?;
                }
                let ordinal = labels.iter().all(|l| l.parse::<f64>().is_ok());
                let data: Vec<String> = labels
                    .into_iter()
                    .zip(values)
                    .map(|(l, v)| format!("{{ 'label': '{}', 'value': {} }}", l, v))
                    .collect();
                let mut data_text = "{'values': [".to_string();
                for datum in data {
                    data_text.push_str(&datum);
                    data_text.push_str(", ");
                }
                data_text.push_str("]}");
                let js_object = format!(
                    "new Bar('{}', '{}', '{}', '{}', {}, {}, {})",
                    id, name, x_label, y_label, data_text, log_toggle, ordinal
                );
                let data = HashMap::from([
                    ("id".to_string(), to_json(&id)),
                    ("log_toggle".to_string(), to_json(log_toggle)),
                ]);
                Ok((
                    registry.render("bar", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::MultiBar {
                id,
                names,
                x_label,
                y_label,
                labels,
                values,
                log_toggle,
                curve,
                alpha,
            } => {
                if !registry.has_template("bar") {
                    registry.register_template_string("bar", from_utf8(BAR_HBS).unwrap())?;
                }
                let ordinal = labels.iter().all(|l| l.parse::<f64>().is_ok());
                let data_text = (0..labels.len())
                    .cartesian_product(0..names.len())
                    .map(|(l, n)| {
                        format!(
                            "{{'label': '{}', 'name': '{}', 'value': {}}}",
                            labels[l], names[n], values[n][l]
                        )
                    })
                    .join(",");
                let data_text = format!("{{'values': [{}]}}", data_text);
                let curve_text = match curve {
                    Some(c) => {
                        let curve_text = labels
                            .iter()
                            .zip(c)
                            .map(|(l, c)| format!("{{'label': '{}', 'value': {}}}", l, c))
                            .join(",");
                        format!("{{'values': [{}]}}", curve_text)
                    }
                    None => "{}".to_string(),
                };
                let js_object = format!(
                    "new MultiBar('{}', '{}', '{}', {}, {}, {}, {}, {})",
                    id,
                    x_label,
                    y_label,
                    log_toggle,
                    data_text,
                    ordinal,
                    alpha.unwrap_or(f64::NAN),
                    curve_text,
                );
                let data = HashMap::from([
                    ("id".to_string(), to_json(&id)),
                    ("log_toggle".to_string(), to_json(log_toggle)),
                ]);
                Ok((
                    registry.render("bar", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::SectionLine {
                id,
                names,
                x_label,
                y_label,
                labels,
                values,
                section_separators,
                section_labels,
                log_toggle,
            } => {
                if !registry.has_template("bar") {
                    registry.register_template_string("bar", from_utf8(BAR_HBS).unwrap())?;
                }
                let ordinal = labels.iter().all(|l| l.parse::<f64>().is_ok());
                let data_text = (0..labels.len())
                    .cartesian_product(0..names.len())
                    .map(|(l, n)| {
                        format!(
                            "{{'label': '{}', 'name': '{}', 'value': {}}}",
                            labels[l], names[n], values[n][l]
                        )
                    })
                    .join(",");
                let data_text = format!("{{'values': [{}]}}", data_text);
                let section_separator_text =
                    section_separators.iter().map(|s| s.to_string()).join(",");
                let section_separator_text = format!("{{'values': [{}]}}", section_separator_text);
                let section_labels_text =
                    section_labels.iter().map(|s| format!("'{}'", s)).join(",");
                let section_labels_text = format!("{{'values': [{}]}}", section_labels_text);
                let js_object = format!(
                    "new SectionLine('{}', '{}', '{}', {}, {}, {}, {}, {})",
                    id,
                    x_label,
                    y_label,
                    log_toggle,
                    data_text,
                    ordinal,
                    section_separator_text,
                    section_labels_text
                );
                let data = HashMap::from([
                    ("id".to_string(), to_json(&id)),
                    ("log_toggle".to_string(), to_json(log_toggle)),
                ]);
                Ok((
                    registry.render("bar", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Hexbin {
                id,
                bins,
                threshold,
            } => {
                if !registry.has_template("hexbin") {
                    registry.register_template_string("hexbin", from_utf8(HEXBIN_HBS).unwrap())?;
                }
                let mut js_object = format!("new Hexbin('{}', {{'values': [", id,);
                for bin in bins.iter() {
                    js_object.push_str(&format!(
                        "{{ coverage: {}, length: {}, size: {} }}, ",
                        bin.x, bin.y, bin.size,
                    ));
                }
                js_object.push_str("]}, [");
                for bin in bins.into_iter() {
                    js_object.push_str("[");
                    for node in bin.content.iter().take(threshold) {
                        js_object.push_str(&format!("'{}',", node,));
                    }
                    if threshold < bin.content.len() {
                        js_object.push_str("-1, ");
                    }
                    js_object.push_str("],");
                }
                js_object.push_str("])");
                let data = HashMap::from([("id".to_string(), to_json(&id))]);
                Ok((
                    registry.render("hexbin", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Line {
                id,
                name,
                x_values,
                y_values,
                log_x,
                log_y,
                x_label,
                y_label,
            } => {
                if !registry.has_template("line") {
                    registry.register_template_string("line", from_utf8(LINE_HBS).unwrap())?;
                }

                let data: Vec<String> = x_values
                    .into_iter()
                    .zip(y_values)
                    .map(|(l, v)| format!("{{ 'x': '{}', 'y': {} }}", l, v))
                    .collect();
                let mut data_text = "{'values': [".to_string();
                for datum in data {
                    data_text.push_str(&datum);
                    data_text.push_str(", ");
                }
                data_text.push_str("]}");
                let js_object = format!(
                    "new Line('{}', '{}', '{}', '{}', {}, {}, {})",
                    id, name, x_label, y_label, log_x, log_y, data_text
                );

                let data = HashMap::from([("id".to_string(), to_json(&id))]);
                Ok((
                    registry.render("line", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Png { id, file } => {
                if !registry.has_template("png") {
                    registry.register_template_string("png", from_utf8(PNG_HBS).unwrap())?;
                }
                let f = File::open(file)?;
                let mut reader = BufReader::new(f);
                let mut buffer = Vec::new();
                reader.read_to_end(&mut buffer)?;
                let base64_text = STANDARD.encode(buffer);
                let data = HashMap::from([("base64", &base64_text), ("id", &id)]);
                let js_object = format!("new DownloadHelper('{}', 'png')", id,);
                Ok((
                    registry.render("png", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Svg { id, file } => {
                if !registry.has_template("svg") {
                    registry.register_template_string("svg", from_utf8(SVG_HBS).unwrap())?;
                }
                let f = File::open(file)?;
                let mut reader = BufReader::new(f);
                let mut buffer = String::new();
                reader.read_to_string(&mut buffer)?;
                let svg_content = buffer;
                let data = HashMap::from([("svg_content", &svg_content), ("id", &id)]);
                let js_object = format!("new DownloadHelper('{}', 'svg')", id,);
                Ok((
                    registry.render("svg", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Json { id, file } => {
                if !registry.has_template("line") {
                    registry.register_template_string("line", from_utf8(LINE_HBS).unwrap())?;
                }

                let f = File::open(file)?;
                let mut reader = BufReader::new(f);
                let mut buffer = String::new();
                reader.read_to_string(&mut buffer)?;
                let json_content = buffer;
                let js_object = format!("new VegaPlot('{}', {})", id, json_content);

                let data = HashMap::from([("id".to_string(), to_json(&id))]);
                Ok((
                    registry.render("line", &data)?,
                    HashMap::from([(
                        "datasets".to_string(),
                        HashMap::from([(id.clone(), js_object)]),
                    )]),
                ))
            }
            Self::Pdf { id: _id, file } => {
                if !registry.has_template("pdf") {
                    registry.register_template_string("pdf", from_utf8(PDF_HBS).unwrap())?;
                }
                let f = File::open(file)?;
                let mut reader = BufReader::new(f);
                let mut buffer = Vec::new();
                reader.read_to_end(&mut buffer)?;
                let base64_text = STANDARD.encode(buffer);
                let data = HashMap::from([("base64", &base64_text)]);
                Ok((
                    registry.render("pdf", &data)?,
                    HashMap::from([("datasets".to_string(), HashMap::new())]),
                ))
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bin {
    pub size: u64,
    pub x: f64,
    pub y: f64,
    pub content: Vec<String>,
}

impl fmt::Display for Bin {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{x: {}, y: {}, size: {:?} }}",
            self.x, self.y, self.size
        )
    }
}

impl Bin {
    pub fn hexbin(points: &Vec<(String, u32, f64)>, nx: u32, ny: u32) -> Vec<Self> {
        let max_coverage = points
            .iter()
            .map(|(_i, c, _l)| *c)
            .max()
            .expect("At least one point");
        let max_length = points.iter().map(|(_i, _c, l)| *l).fold(0.04, f64::max);
        let dx = max_coverage as f64 / (nx - 1) as f64;
        let _t = dx / 3f64.sqrt();
        let dy = max_length / (ny - 1) as f64;
        let mut bins: HashMap<(bool, i64, i64), Self> = HashMap::new();
        for point in points {
            // Calculate positions in both grids
            let mut black_x = (point.1 as f64 / dx).floor() * dx;
            let mut black_y = (point.2 / dy).floor() * dy;
            let mut green_x = ((point.1 as f64 - dx / 2.0) / dx).floor() * dx + dx / 2.0;
            let mut green_y = ((point.2 - dy / 2.0) / dy).floor() * dy + dy / 2.0;

            if black_x < green_x {
                black_x += dx;
            } else {
                green_x += dx;
            }

            if black_y < green_y {
                black_y += dy;
            } else {
                green_y += dy;
            }

            if Self::distance(point.1 as f64, point.2, black_x, black_y)
                < Self::distance(point.1 as f64, point.2, green_x, green_y)
            {
                bins.entry((false, (black_x / dx) as i64, (black_y / dy) as i64))
                    .or_insert(Self {
                        x: black_x,
                        y: black_y,
                        size: 0,
                        content: Vec::new(),
                    })
                    .content
                    .push(point.0.clone());
            } else {
                bins.entry((
                    true,
                    ((green_x - dx / 2.0) / dx) as i64,
                    ((green_y - dy / 2.0) / dy) as i64,
                ))
                .or_insert(Self {
                    x: green_x,
                    y: green_y,
                    size: 0,
                    content: Vec::new(),
                })
                .content
                .push(point.0.clone());
            }
        }
        let mut bins: Vec<Bin> = bins.into_values().collect();
        for bin in &mut bins {
            bin.size = bin.content.len() as u64;
        }
        bins
    }

    fn distance(x1: f64, y1: f64, x2: f64, y2: f64) -> f64 {
        ((x1 - x2).powf(2.0) + (y1 - y2).powf(2.0)).sqrt()
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Window {
    pub start: usize,
    pub end: usize,
    pub values: Vec<f64>,
}
