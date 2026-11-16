/*!
  * Panacus JS library
  */

const PCOLORS = ['#f77189', '#bb9832', '#50b131', '#36ada4', '#3ba3ec', '#e866f4'];

class Bar {
    constructor(id, name, x_label, y_label, data, log_toggle, ordinal) {
        this.id = id;
        this.name = name;
        this.x_label = x_label;
        this.y_label = y_label;
        this.data = data;
        this.log_toggle = log_toggle;
        this.ordinal = ordinal;
    }
}

class MultiBar {
    constructor(id, x_label, y_label, log_toggle, data, ordinal, alpha, heaps_curve) {
        this.id = id;
        this.x_label = x_label;
        this.y_label = y_label;
        this.log_toggle = log_toggle;
        this.data = data;
        this.ordinal = ordinal;
        this.alpha = alpha;
        this.heaps_curve = heaps_curve;
    }
}

class SectionLine {
    constructor(id, x_label, y_label, log_toggle, data, ordinal, separators, section_labels) {
        this.id = id;
        this.x_label = x_label;
        this.y_label = y_label;
        this.log_toggle = log_toggle;
        this.data = data;
        this.ordinal = ordinal;
        this.separators = separators;
        this.section_labels = section_labels;
    }
}

class Hexbin {
    constructor(id, bins, bin_content) {
        this.id = id;
        this.bins = bins;
        this.bin_content = bin_content;
    }
}

class Heatmap {
    constructor(id, name, data_set, dendrogram = null) {
        this.id = id;
        this.name = name;
        this.data_set = data_set;
        this.dendrogram = dendrogram;
    }
}

class Chromosomal {
    constructor(id, name, labels, diverging, contains_outliers, sequence, data) {
        this.id = id;
        this.name = name;
        this.labels = labels;
        this.diverging = diverging;
        this.contains_outliers = contains_outliers;
        this.sequence = sequence;
        this.data = data;
    }
}

class VegaPlot {
    constructor(id, jsonContent) {
        this.id = id;
        this.jsonContent = jsonContent;
    }
}

class Line {
    constructor(id, name, x_label, y_label, log_x, log_y, data) {
        this.id = id;
        this.name = name;
        this.x_label = x_label;
        this.y_label = y_label;
        this.log_x = log_x;
        this.log_y = log_y;
        this.data = data;
    }
}

class DownloadHelper {
    constructor(id, type) {
        this.id = id;
        this.type = type;
    }
}

function buildPlotDownload(chart, obj, prefix) {
    document.getElementById('btn-download-plot-' + obj).onclick = function() {
        var a = document.createElement('a');
        a.href = chart.toBase64Image();
        a.download = prefix + '_' + obj + '.png';
        a.click();
    };
}

function buildTableDownload(table, id, prefix) {
    document.getElementById('btn-download-table-' + id).onclick = function() {
        let blob = new Blob([table], {type: 'text/plain'});
        var a = document.createElement('a');
        a.href = URL.createObjectURL(blob);
        a.download = prefix + '_table.tsv';
        a.click();
    };
}

function buildLogToggle(chart, name) {
    document.getElementById('btn-logscale-plot-' + name).addEventListener('change', function(event) {
        if (event.currentTarget.checked) {
            chart.options.scales.y.type = 'logarithmic';
        } else {
            chart.options.scales.y.type = 'linear';
        }
        chart.update();
    });
}

function getColor(value, zero) {
    const corrected = (value - zero) / (1.0 - zero);
    const flipped = 1.0 - corrected;
    const r = 255.0 - flipped * (255.0 - 179.0);
    const g = 255.0 - flipped * 255.0;
    const b = 255.0 - flipped * 255.0;
    const text = 'rgb(' + r + ',' + g + ',' + b + ')';
    return text;
}

function buildColorSlider(chart, name) {
    document.getElementById('btn-colorscale-' + name).addEventListener('change', function(event) {
        chart.data.datasets[0].backgroundColor = (context) => {
            const value = context.dataset.data[context.dataIndex].v;
            return getColor(value, event.currentTarget.value)
        };
        chart.options.plugins.tooltip.callbacks = {
            title() {
                return '';
            },
            label(context) {
                const v = context.dataset.data[context.dataIndex];
                return [v.x_label + ' - ', v.y_label + ':', v.v];
            }
        }
        chart.update();
    });
}
