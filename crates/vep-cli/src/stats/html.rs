// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

//! The HTML statistics report (`Stats.pm` `dump_html` and `stats_html_head`):
//! one self-contained page with Ensembl VEP's sections, element ids and order,
//! whose charts are inline SVG drawn here, so the page loads no script, no
//! remote chart library and no image. Under each chart's heading a comment
//! holds the statement Ensembl VEP's page hands its chart library for that
//! chart (`Stats.pm:955-963`), so a program that reads a chart's counts from
//! that line of Ensembl VEP's page reads them from this one.

use std::f64::consts::{FRAC_PI_2, TAU};
use std::fmt::Write as _;

use super::{Chart, Report};

/// Every chart is this wide (`dump_html`'s `width: 800px`).
const WIDTH: f64 = 800.0;
/// A pie chart's height when `generate_chart_data` names none.
const PIE_HEIGHT: u32 = 400;
/// The height `generate_chart_data` names for a pie, and every bar and area
/// chart's.
const CHART_HEIGHT: u32 = 200;

/// The margins of a bar or area chart's plot: the count labels sit in the left
/// one, the key labels and the axis title below the baseline.
const PLOT_LEFT: f64 = 50.0;
const PLOT_RIGHT: f64 = 10.0;
const PLOT_TOP: f64 = 10.0;
const PLOT_WIDTH: f64 = WIDTH - PLOT_LEFT - PLOT_RIGHT;

/// The colours of a chart whose keys name none (the default series colours of
/// the chart library Ensembl VEP's page loads), taken in key order.
const PALETTE: [&str; 20] = [
    "#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477", "#66aa00",
    "#b82e2e", "#316395", "#994499", "#22aa99", "#aaaa11", "#6633cc", "#e67300", "#8b0707",
    "#651067", "#329262", "#5574a6", "#3b3eac",
];

/// `COLOUR_KEYS` (`Constants.pm:299-346`), keyed as the chart rows are.
const CONSEQUENCE_COLOURS: [(&str, &str); 32] = [
    ("intergenic_variant", "gray"),
    ("intron_variant", "#02599c"),
    ("upstream_gene_variant", "#a2b5cd"),
    ("downstream_gene_variant", "#a2b5cd"),
    ("5_prime_utr_variant", "#7ac5cd"),
    ("3_prime_utr_variant", "#7ac5cd"),
    ("splice_region_variant", "#ff7f50"),
    ("splice_donor_variant", "#ff7f50"),
    ("splice_acceptor_variant", "#ff7f50"),
    ("frameshift_variant", "#ff69b4"),
    ("transcript_ablation", "#ff0000"),
    ("transcript_amplification", "#ff69b4"),
    ("inframe_insertion", "#ff69b4"),
    ("inframe_deletion", "#ff69b4"),
    ("synonymous_variant", "#76ee00"),
    ("stop_retained_variant", "#76ee00"),
    ("missense_variant", "#ffd700"),
    ("initiator_codon_variant", "#ffd700"),
    ("stop_gained", "#ff0000"),
    ("stop_lost", "#ff0000"),
    ("mature_mirna_variant", "#458b00"),
    ("non_coding_exon_variant", "#32cd32"),
    ("nc_transcript_variant", "#32cd32"),
    ("incomplete_terminal_codon_variant", "#ff00ff"),
    ("nmd_transcript_variant", "#ff4500"),
    ("coding_sequence_variant", "#458b00"),
    ("tfbs_ablation", "brown"),
    ("tfbs_amplification", "brown"),
    ("tf_binding_site_variant", "brown"),
    ("regulatory_region_variant", "brown"),
    ("regulatory_region_ablation", "brown"),
    ("regulatory_region_amplification", "brown"),
];
const SIFT_COLOURS: [(&str, &str); 2] = [("tolerated", "green"), ("deleterious", "red")];
const POLYPHEN_COLOURS: [(&str, &str); 4] = [
    ("unknown", "blue"),
    ("benign", "green"),
    ("possibly damaging", "orange"),
    ("probably damaging", "red"),
];

/// The chart types `generate_chart_data` assigns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Pie,
    Bar,
    Area,
}

impl Kind {
    /// The suffix of the chart element's id (`<id>_<type>`).
    fn name(self) -> &'static str {
        match self {
            Kind::Pie => "pie",
            Kind::Bar => "bar",
            Kind::Area => "area",
        }
    }

    /// The function a chart's statement calls: `draw` and the type with its
    /// first letter in upper case (`ucfirst`, `Stats.pm:940`).
    fn drawer(self) -> &'static str {
        match self {
            Kind::Pie => "drawPie",
            Kind::Bar => "drawBar",
            Kind::Area => "drawArea",
        }
    }
}

/// How `dump_html` draws one chart (`generate_chart_data`, `Stats.pm:505-609`):
/// its type, the two column headers of its table, its height, whether the table
/// prints (`no_table` unset), whether the side menu links to it (`no_link`
/// unset), the colours its keys name (`colours`) and the chart library options
/// its statement carries (`options`). No chart names both colours and options.
struct Shape {
    kind: Kind,
    header: [&'static str; 2],
    height: u32,
    table: bool,
    linked: bool,
    colours: &'static [(&'static str, &'static str)],
    options: Option<&'static str>,
}

/// The shape of the chart with element id `id`.
fn shape(id: &str) -> Shape {
    let pie = |header, height, colours| Shape {
        kind: Kind::Pie,
        header: [header, "Count"],
        height,
        table: true,
        linked: true,
        colours,
        options: None,
    };
    match id {
        "var_class" => pie("Variant class", CHART_HEIGHT, &[]),
        "var_cons" | "consequences" | "coding" => {
            pie("Consequence type", PIE_HEIGHT, &CONSEQUENCE_COLOURS)
        }
        "sift" => pie("Prediction", CHART_HEIGHT, &SIFT_COLOURS),
        "polyphen" => pie("Prediction", CHART_HEIGHT, &POLYPHEN_COLOURS),
        "chr" => Shape {
            kind: Kind::Bar,
            header: ["Chromosome", "Count"],
            height: CHART_HEIGHT,
            table: true,
            linked: true,
            colours: &[],
            options: Some(r#"{legend: {position: "none"}}"#),
        },
        "protein" => Shape {
            kind: Kind::Bar,
            header: ["Position in protein (percentile)", "Count"],
            height: CHART_HEIGHT,
            table: false,
            linked: true,
            colours: &[],
            options: Some(
                r#"{hAxis: {title: "Position in protein (percentile)", textStyle: {fontSize: 10}}, legend: {position: "none"}}"#,
            ),
        },
        // Every other id is a `Distribution of variants on chromosome` chart
        // (`chr_<name>`, `Stats.pm:580-595`).
        _ => Shape {
            kind: Kind::Area,
            header: ["Position (mb)", "Count"],
            height: CHART_HEIGHT,
            table: false,
            linked: false,
            colours: &[],
            options: Some(
                r#"{hAxis: {title: "Position (mb)", textStyle: {fontSize: 8}}, legend: {position: "none"}}"#,
            ),
        },
    }
}

/// The side menu (`dump_html`, `Stats.pm:855-871`): the four fixed sections,
/// then every chart whose `no_link` is unset, as `id` and title pairs.
fn side_menu_entries(report: &Report) -> Vec<(String, String)> {
    let mut entries = vec![
        ("masthead".to_string(), "Top of page".to_string()),
        ("run_stats".to_string(), "VEP run statistics".to_string()),
        ("data_version".to_string(), "Data version".to_string()),
        ("gen_stats".to_string(), "General statistics".to_string()),
    ];
    entries.extend(
        report
            .charts
            .iter()
            .filter(|chart| shape(&chart.id).linked)
            .map(|chart| (chart.id.clone(), chart.title.clone())),
    );
    entries
}

/// The page (`dump_html`, `Stats.pm:846-915`, under the head and body
/// skeleton of `stats_html_head`, `:998-1148`).
pub(super) fn page(report: &Report) -> String {
    let mut out = String::new();
    out.push_str(HEAD);
    out.push_str("<div class=\"sidemenu\"><div class=\"sidemenu_head\">Links</div>");
    out.push_str("<div class=\"sidemenu_body\"><ul>");
    for (id, title) in side_menu_entries(report) {
        let _ = write!(out, "<li><a href=\"#{id}\">{}</a></li>", escape(&title));
    }
    out.push_str("</ul></div></div>");

    out.push_str("<div class='main_content'>");
    out.push_str("<h3 id=\"run_stats\">VEP run statistics</h3>");
    stats_table(
        &mut out,
        report
            .run_stats
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str())),
    );
    out.push_str("<h3 id=\"data_version\">Data version</h3>");
    if report.data_version.is_empty() {
        out.push_str("No data version information available.");
    } else {
        stats_table(
            &mut out,
            report
                .data_version
                .iter()
                .map(|(k, v)| (k.as_str(), v.as_str())),
        );
    }
    out.push_str("<h3 id=\"gen_stats\">General statistics</h3>");
    stats_table(
        &mut out,
        report.general_stats.iter().map(|(label, value)| {
            // `($_ || 0)` (`Stats.pm:895`): an undefined or empty value prints `0`.
            let value = match value {
                Some(value) if !value.is_empty() => value.as_str(),
                _ => "0",
            };
            (label.as_str(), value)
        }),
    );

    for chart in &report.charts {
        chart_section(&mut out, chart);
    }
    out.push_str("</div>");
    out.push_str("\n</div></body>\n</html>\n");
    out
}

/// A `stats_table` of two-cell rows, each cell's text printed as the report
/// holds it (`Stats.pm:878, 885, 895`): the run statistics carry their own
/// `<pre>` and `<kbd>` elements.
fn stats_table<'a>(out: &mut String, rows: impl Iterator<Item = (&'a str, &'a str)>) {
    out.push_str("<table class=\"stats_table\">");
    for (label, value) in rows {
        let _ = write!(out, "<tr><td>{label}</td><td>{value}</td></tr>");
    }
    out.push_str("</table>");
}

/// One chart (`Stats.pm:898-909`): a rule, the titled heading, the chart's
/// statement in a comment on a line of its own, the chart element holding the
/// drawing, and the table of its rows unless the chart has none. The statement,
/// the drawing and the table are written from the same rows.
fn chart_section(out: &mut String, chart: &Chart) {
    let shape = shape(&chart.id);
    let _ = write!(
        out,
        "<hr/><h3 id=\"{}\">{}</h3>",
        chart.id,
        escape(&chart.title)
    );
    // A program that reads a chart's counts takes them from the first line
    // that holds the chart's title before the statement, so the statement
    // shares its line with no heading and no other chart's statement.
    let _ = write!(out, "\n<!-- {} -->\n", chart_statement(chart, &shape));
    let _ = write!(
        out,
        "<div id=\"{}_{}\" style=\"width: {}px; height: {}px\">",
        chart.id,
        shape.kind.name(),
        WIDTH,
        shape.height
    );
    out.push_str(&chart_svg(chart, &shape));
    out.push_str("</div>");
    if shape.table {
        let _ = write!(
            out,
            "<table id=\"{}_table\" class=\"stats_table\" style=\"width: {}px\"><tr><th>{}</th><th>{}</th></tr>",
            chart.id, WIDTH, shape.header[0], shape.header[1]
        );
        for (key, count) in &chart.rows {
            let _ = write!(out, "<tr><td>{}</td><td>{count}</td></tr>", escape(key));
        }
        out.push_str("</table>");
    }
}

/// The drawing of `chart` as an inline `<svg>` element of its shape's size.
fn chart_svg(chart: &Chart, shape: &Shape) -> String {
    let height = f64::from(shape.height);
    let mut svg = String::new();
    let _ = write!(
        svg,
        "<svg width=\"{WIDTH}\" height=\"{height}\" viewBox=\"0 0 {WIDTH} {height}\" font-family=\"arial, sans-serif\">"
    );
    match shape.kind {
        Kind::Pie => pie(&mut svg, &chart.rows, height, shape.colours),
        Kind::Bar => bar(&mut svg, &chart.rows, height),
        Kind::Area => area(&mut svg, &chart.rows, height, shape.header[0]),
    }
    svg.push_str("</svg>");
    svg
}

/// The statement `stats_html_head` writes for one chart (`Stats.pm:955-963`):
/// the chart element's id, the function of its type, its title, the two column
/// headers and every row as `['key',count]`, then the chart's options. Nothing
/// in the page runs it: it stands where a program that collects reports finds
/// a chart's counts in Ensembl VEP's page.
fn chart_statement(chart: &Chart, shape: &Shape) -> String {
    let element = format!("{}_{}", script_text(&chart.id), shape.kind.name());
    let mut statement = format!(
        "var {element} = {}('{element}', '{}', google.visualization.arrayToDataTable([['{}','{}'],",
        shape.kind.drawer(),
        script_text(&chart.title),
        shape.header[0],
        shape.header[1]
    );
    for (index, (key, count)) in chart.rows.iter().enumerate() {
        if index > 0 {
            statement.push(',');
        }
        let _ = write!(statement, "['{}',{count}]", script_text(key));
    }
    let _ = write!(statement, "]), {});", statement_options(chart, shape));
    statement
}

/// The last argument of a chart's statement (`Stats.pm:942-952, 962`): a
/// `slices` list with one entry per row for a chart whose keys name colours,
/// the colour where the chart names the row's key and `{}` where it does not;
/// else the chart's own options, or `null` for a chart with neither.
fn statement_options(chart: &Chart, shape: &Shape) -> String {
    if shape.colours.is_empty() {
        return shape.options.unwrap_or("null").to_string();
    }
    let slices: Vec<String> = chart
        .rows
        .iter()
        .map(|(key, _)| match named_colour(key, shape.colours) {
            Some(colour) => format!("{{color: \"{colour}\"}}"),
            None => "{}".to_string(),
        })
        .collect();
    format!("{{slices: [{}]}}", slices.join(", "))
}

/// `text` as it stands in a chart's statement: a backslash and a single quote
/// take a backslash, and `<`, `>` and a control character are written as
/// `\xNN`, which a single-quoted string reads back as the character. So no
/// text ends the string, the line or the comment its statement stands in.
fn script_text(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\'' => out.push_str("\\'"),
            c if matches!(c, '<' | '>') || c.is_control() => {
                let _ = write!(out, "\\x{:02x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

/// The colour the chart's `colours` give `key`, for a key they name.
fn named_colour(key: &str, colours: &[(&'static str, &'static str)]) -> Option<&'static str> {
    colours
        .iter()
        .find(|(name, _)| *name == key)
        .map(|(_, colour)| *colour)
}

/// The colour of the key at `index`: the chart's own for a key it names, else
/// the palette's, taken in key order.
fn colour_of(key: &str, index: usize, colours: &[(&'static str, &'static str)]) -> &'static str {
    named_colour(key, colours).unwrap_or(PALETTE[index % PALETTE.len()])
}

/// A pie of the rows' shares, drawn clockwise from the top in row order, with a
/// legend of `key: count` entries beside it in as many columns as fit.
fn pie(
    svg: &mut String,
    rows: &[(String, u64)],
    height: f64,
    colours: &[(&'static str, &'static str)],
) {
    let total: u64 = rows.iter().map(|(_, count)| count).sum();
    let radius = height / 2.0 - 10.0;
    let centre = height / 2.0;
    let mut start = -FRAC_PI_2;
    for (index, (key, count)) in rows.iter().enumerate() {
        if total == 0 || *count == 0 {
            continue;
        }
        let colour = colour_of(key, index, colours);
        let share = *count as f64 / total as f64;
        let title = format!("{}: {count} ({:.1}%)", escape(key), 100.0 * share);
        if *count == total {
            let _ = write!(
                svg,
                "<circle cx=\"{}\" cy=\"{}\" r=\"{}\" fill=\"{colour}\"><title>{title}</title></circle>",
                num(centre),
                num(centre),
                num(radius)
            );
            continue;
        }
        let end = start + TAU * share;
        let (x1, y1) = (centre + radius * start.cos(), centre + radius * start.sin());
        let (x2, y2) = (centre + radius * end.cos(), centre + radius * end.sin());
        let large_arc = u8::from(share > 0.5);
        let _ = write!(
            svg,
            "<path d=\"M {c} {c} L {x1} {y1} A {r} {r} 0 {large_arc} 1 {x2} {y2} Z\" fill=\"{colour}\" stroke=\"white\" stroke-width=\"1\"><title>{title}</title></path>",
            c = num(centre),
            x1 = num(x1),
            y1 = num(y1),
            r = num(radius),
            x2 = num(x2),
            y2 = num(y2),
        );
        start = end;
    }

    const LINE: f64 = 14.0;
    const COLUMN: f64 = 190.0;
    let left = height + 20.0;
    let lines_per_column = ((height - 20.0) / LINE).max(1.0) as usize;
    let columns = ((WIDTH - left) / COLUMN).max(1.0) as usize;
    for (index, (key, count)) in rows.iter().enumerate().take(lines_per_column * columns) {
        let colour = colour_of(key, index, colours);
        let x = left + COLUMN * (index / lines_per_column) as f64;
        let y = 10.0 + LINE * (index % lines_per_column) as f64;
        let _ = write!(
            svg,
            "<rect x=\"{}\" y=\"{}\" width=\"10\" height=\"10\" fill=\"{colour}\"/><text x=\"{}\" y=\"{}\" font-size=\"11\">{}: {count}</text>",
            num(x),
            num(y),
            num(x + 14.0),
            num(y + 9.0),
            escape(key)
        );
    }
}

/// A column chart: one rectangle per key, scaled to the largest count, its key
/// beneath it (turned when the keys are too many to sit upright) and the
/// largest count at the top of the axis.
fn bar(svg: &mut String, rows: &[(String, u64)], height: f64) {
    const UPRIGHT_LABELS: usize = 24;
    let turned = rows.len() > UPRIGHT_LABELS;
    let bottom = if turned { 60.0 } else { 30.0 };
    let plot_height = height - PLOT_TOP - bottom;
    let baseline = PLOT_TOP + plot_height;
    let max = largest_count(rows);
    axis(svg, baseline, max);
    let slot = PLOT_WIDTH / rows.len().max(1) as f64;
    for (index, (key, count)) in rows.iter().enumerate() {
        let x = PLOT_LEFT + slot * index as f64;
        let bar_height = plot_height * (*count as f64 / max as f64);
        let _ = write!(
            svg,
            "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\"><title>{}: {count}</title></rect>",
            num(x + slot * 0.15),
            num(baseline - bar_height),
            num(slot * 0.7),
            num(bar_height),
            PALETTE[0],
            escape(key)
        );
        let label_x = x + slot / 2.0;
        if turned {
            let _ = write!(
                svg,
                "<text x=\"{x}\" y=\"{y}\" font-size=\"10\" text-anchor=\"end\" transform=\"rotate(-60 {x} {y})\">{}</text>",
                escape(key),
                x = num(label_x),
                y = num(baseline + 12.0),
            );
        } else {
            let _ = write!(
                svg,
                "<text x=\"{}\" y=\"{}\" font-size=\"10\" text-anchor=\"middle\">{}</text>",
                num(label_x),
                num(baseline + 12.0),
                escape(key)
            );
        }
    }
}

/// An area chart over consecutive bins: one polygon from the first bin to the
/// last, scaled to the largest count, over a position axis whose ticks name
/// every `step`th bin and carry the axis title.
fn area(svg: &mut String, rows: &[(String, u64)], height: f64, axis_title: &str) {
    const BOTTOM: f64 = 40.0;
    const TICKS: usize = 12;
    let plot_height = height - PLOT_TOP - BOTTOM;
    let baseline = PLOT_TOP + plot_height;
    let max = largest_count(rows);
    axis(svg, baseline, max);
    let x_of = |index: usize| -> f64 {
        if rows.len() < 2 {
            PLOT_LEFT
        } else {
            PLOT_LEFT + PLOT_WIDTH * index as f64 / (rows.len() - 1) as f64
        }
    };
    if !rows.is_empty() {
        let mut points = format!("{},{}", num(PLOT_LEFT), num(baseline));
        for (index, (_, count)) in rows.iter().enumerate() {
            let y = baseline - plot_height * (*count as f64 / max as f64);
            let _ = write!(points, " {},{}", num(x_of(index)), num(y));
        }
        let _ = write!(points, " {},{}", num(x_of(rows.len() - 1)), num(baseline));
        let _ = write!(
            svg,
            "<polygon points=\"{points}\" fill=\"{0}\" fill-opacity=\"0.4\" stroke=\"{0}\" stroke-width=\"1\"/>",
            PALETTE[0]
        );
    }
    let step = [1, 2, 5, 10, 20, 50, 100, 200, 500]
        .into_iter()
        .find(|step| rows.len() / step <= TICKS)
        .unwrap_or(1000);
    for (index, (key, _)) in rows.iter().enumerate().step_by(step) {
        let _ = write!(
            svg,
            "<text x=\"{}\" y=\"{}\" font-size=\"8\" text-anchor=\"middle\">{}</text>",
            num(x_of(index)),
            num(baseline + 12.0),
            escape(key)
        );
    }
    let _ = write!(
        svg,
        "<text x=\"{}\" y=\"{}\" font-size=\"10\" text-anchor=\"middle\">{}</text>",
        num(PLOT_LEFT + PLOT_WIDTH / 2.0),
        num(baseline + 30.0),
        escape(axis_title)
    );
}

/// The largest count of `rows`, or 1 when no count is above 0, so the scale
/// of a bar or area chart is never 0.
fn largest_count(rows: &[(String, u64)]) -> u64 {
    rows.iter()
        .map(|(_, count)| *count)
        .max()
        .unwrap_or(0)
        .max(1)
}

/// The two axis lines of a bar or area chart, the largest count at the top of
/// the vertical one and `0` at its foot.
fn axis(svg: &mut String, baseline: f64, max: u64) {
    let _ = write!(
        svg,
        "<line x1=\"{left}\" y1=\"{top}\" x2=\"{left}\" y2=\"{base}\" stroke=\"#333\"/><line x1=\"{left}\" y1=\"{base}\" x2=\"{right}\" y2=\"{base}\" stroke=\"#333\"/>",
        left = num(PLOT_LEFT),
        top = num(PLOT_TOP),
        base = num(baseline),
        right = num(WIDTH - PLOT_RIGHT),
    );
    let _ = write!(
        svg,
        "<text x=\"{}\" y=\"{}\" font-size=\"10\" text-anchor=\"end\">{max}</text><text x=\"{}\" y=\"{}\" font-size=\"10\" text-anchor=\"end\">0</text>",
        num(PLOT_LEFT - 4.0),
        num(PLOT_TOP + 4.0),
        num(PLOT_LEFT - 4.0),
        num(baseline)
    );
}

/// A coordinate or length with one decimal, which is how the computed
/// positions and sizes of the drawn elements are written.
fn num(value: f64) -> String {
    format!("{value:.1}")
}

/// `text` as HTML character data: the four characters with a meaning there
/// replaced by their references.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            c => out.push(c),
        }
    }
    out
}

/// The head and the body's opening (`stats_html_head`, `Stats.pm:998-1148`)
/// with the chart library and its loader left out, the style sheet (`:1041-1133`)
/// kept, and a text masthead in place of the two linked images (`:1137-1143`).
const HEAD: &str = r#"<html>
<head>
  <title>VEP summary</title>
  <style type="text/css">
    body {
      font-family: arial, sans-serif;
      margin: 0px;
      padding: 0px;
    }

    a {color: #36b;}
    a.visited {color: #006;}

    .stats_table {
      margin: 5px;
      table-layout: fixed;
      width: 100%;
    }

    pre {
      white-space: normal;
      word-break: keep-all;
    }

    td {
      padding: 5px;
      word-break: break-word;
    }

    td:nth-child(2) {
      width: 80%;
    }

    th {
      text-align: left;
      padding: 5px;
    }

    .stats_table tr:nth-child(odd) {
      background-color: #f0f0f0;
    }

    h3 {
      color: #666;
    }

    .masthead {
      background-color: white;
      color: rgb(204, 221, 255);
      height: 80px;
      width: 100%;
      padding: 0px;
    }

    .masthead .title {
      color: white;
      font-size: 28px;
      line-height: 80px;
      padding-left: 20px;
    }

    .main {
      padding: 10px;
    }

    .gradient {
      background: #333366;
      background: linear-gradient(to right, #333366 0%, #ffffff 100%);
      padding: 0px;
      height: 80px;
      width: 700px;
    }

    .main_content {
      margin-left: 300px;
    }

    .sidemenu {
      width: 260px;
      position: fixed;
      border-style: solid;
      border-width: 2px;
      border-color: rgb(51, 51, 102);
    }

    .sidemenu_head {
      width: 250px;
      background-color: rgb(51, 51, 102);
      color: rgb(204, 221, 255);
      padding: 5px;
    }

    .sidemenu_body {
      width: 250px;
      padding: 5px;
    }
  </style>
</head>
<body>
<div id="masthead" class="masthead">
  <div class="gradient"><span class="title">Ensembl Variant Effect Predictor</span></div>
</div>
<div class="main">
"#;

#[cfg(test)]
mod tests {
    use super::super::tests::context;
    use super::super::{Report, ReportContext, RunStats, ThreadStats};
    use super::*;
    use vep_core::transcript::Transcript;

    use std::sync::Arc;

    use smallvec::smallvec;
    use vep_core::consequence::{Consequence, FeatureType, TranscriptConsequence};
    use vep_core::prediction::PredictionTerm;
    use vep_core::variant::InputVariant;
    use vep_io::output::fields::Row;

    /// The ids of every chart `generate_chart_data` builds for a run with a
    /// FASTA, both prediction tools and records on chromosome 21.
    const CHART_IDS: [&str; 9] = [
        "var_class",
        "var_cons",
        "consequences",
        "coding",
        "sift",
        "polyphen",
        "chr",
        "chr_21",
        "protein",
    ];

    /// A record of one allele at `start` on chromosome 21.
    fn record(start: u64, reference: &str, alternate: &str) -> InputVariant {
        let mut v = InputVariant::new(
            "21".to_string(),
            start,
            start + reference.len() as u64 - 1,
            reference.as_bytes().to_vec(),
            alternate.as_bytes().to_vec(),
        );
        v.display_consequence = Some(Consequence::MissenseVariant);
        v
    }

    /// A transcript row with a SIFT and a PolyPhen term in the middle of its
    /// protein.
    fn transcript_row(id: &str) -> TranscriptConsequence {
        let mut consequence = TranscriptConsequence::default();
        consequence.transcript = Some(Arc::new({
            let mut transcript = Transcript::default();
            transcript.stable_id = Arc::from(id);
            transcript.gene_stable_id = Arc::from("ENSG1");
            transcript
        }));
        consequence.consequences = smallvec![Consequence::MissenseVariant];
        consequence.feature_type = FeatureType::Transcript;
        consequence.sift_term = Some(PredictionTerm::Deleterious);
        consequence.polyphen_term = Some(PredictionTerm::Benign);
        consequence.protein_start = Some(50);
        consequence.protein_total = Some(100);
        consequence
    }

    /// The statistics of three records with two transcript rows each, under
    /// `--sift b --polyphen b` and a FASTA of chromosome 21, as [`Report::new`]
    /// reads them.
    fn full_report() -> Report {
        let mut thread = ThreadStats::new();
        let rows = [transcript_row("ENST1"), transcript_row("ENST2")];
        for start in [100, 1_500_000, 2_000_000] {
            let v = record(start, "A", "G");
            let record_rows: Vec<Row<'_>> = rows.iter().map(|tc| (&v, Some(tc))).collect();
            thread.log_record(&[&v], &record_rows);
        }
        let mut stats = RunStats::new();
        stats.merge(thread);
        stats.set_lines_read(3);
        let ctx = ReportContext {
            sift: Some("b".to_string()),
            polyphen: Some("b".to_string()),
            chr_lengths: vec![("21".to_string(), 46_709_983)],
            ..context()
        };
        Report::new(&stats, &ctx)
    }

    /// A report of fixed rows, for a drawing that must not change.
    fn fixed_report() -> Report {
        let chart = |id: &str, title: &str, rows: &[(&str, u64)]| Chart {
            id: id.to_string(),
            title: title.to_string(),
            rows: rows.iter().map(|(k, n)| (k.to_string(), *n)).collect(),
        };
        Report {
            run_stats: vec![("Species".to_string(), "homo_sapiens".to_string())],
            data_version: Vec::new(),
            general_stats: vec![("Lines of input read".to_string(), None)],
            charts: vec![
                chart(
                    "var_cons",
                    "Consequences (most severe)",
                    &[("missense_variant", 3), ("intron_variant", 5), ("other", 2)],
                ),
                chart("chr", "Variants by chromosome", &[("21", 7), ("X", 3)]),
                chart(
                    "chr_21",
                    "Distribution of variants on chromosome 21",
                    &[("0", 1), ("1", 0), ("2", 4), ("3", 2)],
                ),
            ],
        }
    }

    /// The chart table a statement hands the chart library, read back: the two
    /// column headers, then each row's key and count. A string is single-quoted
    /// and reads `\\`, `\'` and `\xNN` as one character each.
    fn read_table(literal: &str) -> ([String; 2], Vec<(String, u64)>) {
        fn string(chars: &mut std::iter::Peekable<std::str::Chars<'_>>) -> String {
            assert_eq!(chars.next(), Some('\''));
            let mut text = String::new();
            loop {
                match chars.next().expect("a string ends at its quote") {
                    '\'' => return text,
                    '\\' => match chars.next().expect("an escape names a character") {
                        'x' => {
                            let code: String = chars.by_ref().take(2).collect();
                            text.push(char::from(u8::from_str_radix(&code, 16).unwrap()));
                        }
                        c => text.push(c),
                    },
                    c => text.push(c),
                }
            }
        }
        let mut chars = literal.chars().peekable();
        assert_eq!(chars.next(), Some('['));
        assert_eq!(chars.next(), Some('['));
        let first = string(&mut chars);
        assert_eq!(chars.next(), Some(','));
        let second = string(&mut chars);
        assert_eq!(chars.next(), Some(']'));
        assert_eq!(chars.next(), Some(','));
        let mut rows = Vec::new();
        while chars.next_if_eq(&'[').is_some() {
            let key = string(&mut chars);
            assert_eq!(chars.next(), Some(','));
            let mut count = String::new();
            while let Some(digit) = chars.next_if(char::is_ascii_digit) {
                count.push(digit);
            }
            assert_eq!(chars.next(), Some(']'));
            rows.push((key, count.parse().unwrap()));
            chars.next_if_eq(&',');
        }
        assert_eq!(chars.next(), Some(']'));
        assert_eq!(chars.next(), None);
        ([first, second], rows)
    }

    /// What a program that collects reports takes for the chart titled `title`:
    /// the first line of the page that holds the title, then
    /// `google.visualization`, then a `;`, and in that stretch the text from the
    /// first `[[` to the last `]]`.
    fn collected<'a>(html: &'a str, title: &str) -> Option<&'a str> {
        html.lines().find_map(|line| {
            let titled = &line[line.find(title)?..];
            let call = &titled[..titled.rfind(';')?];
            call[title.len()..].find("google.visualization")?;
            Some(&call[call.find("[[")?..call.rfind("]]")? + 2])
        })
    }

    /// The rows of the chart's table as the page prints them.
    fn table_rows(html: &str, id: &str) -> Vec<(String, u64)> {
        let table = &html[html.find(&format!("<table id=\"{id}_table\"")).unwrap()..];
        let table = &table[..table.find("</table>").unwrap()];
        table
            .split("<tr><td>")
            .skip(1)
            .map(|row| {
                let (key, rest) = row.split_once("</td><td>").unwrap();
                let count = rest.strip_suffix("</td></tr>").unwrap();
                (key.to_string(), count.parse().unwrap())
            })
            .collect()
    }

    /// `stats_html_head` (`Stats.pm:942-963`): a chart's statement names its
    /// element twice, the function of its type, its title, its two headers and
    /// its rows in order, then its options: a `slices` entry per row for a
    /// chart whose keys name colours, the chart's own options, or `null`.
    #[test]
    fn chart_statements_are_the_lines_stats_html_head_writes() {
        let statement = |chart: &Chart| chart_statement(chart, &shape(&chart.id));
        let report = fixed_report();
        let statements: Vec<String> = report.charts.iter().map(statement).collect();
        assert_eq!(
            statements,
            [
                r##"var var_cons_pie = drawPie('var_cons_pie', 'Consequences (most severe)', google.visualization.arrayToDataTable([['Consequence type','Count'],['missense_variant',3],['intron_variant',5],['other',2]]), {slices: [{color: "#ffd700"}, {color: "#02599c"}, {}]});"##,
                r##"var chr_bar = drawBar('chr_bar', 'Variants by chromosome', google.visualization.arrayToDataTable([['Chromosome','Count'],['21',7],['X',3]]), {legend: {position: "none"}});"##,
                r##"var chr_21_area = drawArea('chr_21_area', 'Distribution of variants on chromosome 21', google.visualization.arrayToDataTable([['Position (mb)','Count'],['0',1],['1',0],['2',4],['3',2]]), {hAxis: {title: "Position (mb)", textStyle: {fontSize: 8}}, legend: {position: "none"}});"##,
            ]
        );
        let chart = |id: &str, title: &str, rows: &[(&str, u64)]| Chart {
            id: id.to_string(),
            title: title.to_string(),
            rows: rows.iter().map(|(k, n)| (k.to_string(), *n)).collect(),
        };
        assert_eq!(
            statement(&chart(
                "var_class",
                "Variant classes",
                &[("insertion", 2), ("SNV", 9)]
            )),
            "var var_class_pie = drawPie('var_class_pie', 'Variant classes', google.visualization.arrayToDataTable([['Variant class','Count'],['insertion',2],['SNV',9]]), null);"
        );
        assert_eq!(
            statement(&chart(
                "sift",
                "SIFT summary",
                &[
                    ("tolerated_low_confidence", 1),
                    ("tolerated", 4),
                    ("deleterious", 6)
                ]
            )),
            r##"var sift_pie = drawPie('sift_pie', 'SIFT summary', google.visualization.arrayToDataTable([['Prediction','Count'],['tolerated_low_confidence',1],['tolerated',4],['deleterious',6]]), {slices: [{}, {color: "green"}, {color: "red"}]});"##
        );
        assert_eq!(
            statement(&chart(
                "protein",
                "Position in protein",
                &[("00-10%", 4), ("10-20%", 0)]
            )),
            r##"var protein_bar = drawBar('protein_bar', 'Position in protein', google.visualization.arrayToDataTable([['Position in protein (percentile)','Count'],['00-10%',4],['10-20%',0]]), {hAxis: {title: "Position in protein (percentile)", textStyle: {fontSize: 10}}, legend: {position: "none"}});"##
        );
        // A chart without rows keeps its headers and an empty `slices` list.
        assert_eq!(
            statement(&chart("coding", "Coding consequences", &[])),
            "var coding_pie = drawPie('coding_pie', 'Coding consequences', google.visualization.arrayToDataTable([['Consequence type','Count'],]), {slices: []});"
        );
    }

    /// Each chart's heading is followed by its statement in a comment on a line
    /// of its own, then by its chart element. What a program that collects
    /// reports takes from that line reads back as the chart's headers and rows,
    /// which are the rows of its table; and with the statement lines taken out
    /// the body is the one line of sections it is without them.
    #[test]
    fn each_chart_title_is_followed_by_its_counts_as_a_collector_reads_them() {
        for report in [full_report(), fixed_report()] {
            let html = report.html();
            let mut without_statements = html.clone();
            for chart in &report.charts {
                let shape = shape(&chart.id);
                let statement = chart_statement(chart, &shape);
                assert!(
                    html.contains(&format!(
                        "<h3 id=\"{}\">{}</h3>\n<!-- {statement} -->\n<div id=\"{}_{}\"",
                        chart.id,
                        escape(&chart.title),
                        chart.id,
                        shape.kind.name()
                    )),
                    "{}",
                    chart.id
                );
                let (header, rows) = read_table(
                    collected(&html, &chart.title)
                        .unwrap_or_else(|| panic!("no statement is read for {}", chart.id)),
                );
                assert_eq!(header, shape.header, "{}", chart.id);
                assert_eq!(rows, chart.rows, "{}", chart.id);
                if shape.table {
                    assert_eq!(rows, table_rows(&html, &chart.id), "{}", chart.id);
                }
                without_statements =
                    without_statements.replace(&format!("\n<!-- {statement} -->\n"), "");
            }
            let naming_the_library: Vec<&str> = html
                .lines()
                .filter(|line| line.contains("google.visualization"))
                .collect();
            assert_eq!(naming_the_library.len(), report.charts.len());
            for line in naming_the_library {
                assert!(
                    line.starts_with("<!-- var ") && line.ends_with("); -->"),
                    "{line}"
                );
            }
            assert!(!without_statements.contains("<!--"));
            let body = without_statements
                .strip_prefix(HEAD)
                .and_then(|rest| rest.strip_suffix("\n</div></body>\n</html>\n"))
                .unwrap();
            assert!(!body.contains('\n'));
        }
    }

    /// A key that holds a quote, a backslash, a line break or the marks that
    /// open and close a comment stands in its statement as one string on one
    /// line and reads back as itself; the page still holds one comment a chart.
    #[test]
    fn statement_text_ends_neither_its_string_nor_its_comment() {
        let keys = ["it's", "a\\b", "x-->y", "<!--", "two\nlines", "a]]b"];
        let report = Report {
            charts: vec![Chart {
                id: "chr".to_string(),
                title: "Variants by chromosome".to_string(),
                rows: keys
                    .iter()
                    .zip(1..)
                    .map(|(key, count)| (key.to_string(), count))
                    .collect(),
            }],
            ..fixed_report()
        };
        let html = report.html();
        assert_eq!(html.matches("<!--").count(), 1);
        assert_eq!(html.matches("-->").count(), 1);
        let (_, rows) = read_table(collected(&html, "Variants by chromosome").unwrap());
        assert_eq!(rows, report.charts[0].rows);
        assert_eq!(
            script_text("it's <b> \\ \u{7f}\r"),
            r"it\'s \x3cb\x3e \\ \x7f\x0d"
        );
    }

    /// The page loads nothing: no script, no remote resource; its sections and
    /// charts carry `Stats.pm`'s element ids, the prediction pies only under
    /// their flags, a distribution chart per chromosome with a FASTA.
    #[test]
    fn html_report_is_self_contained() {
        let html = full_report().html();
        assert!(html.starts_with("<html>\n<head>\n  <title>VEP summary</title>"));
        assert!(html.ends_with("\n</div></body>\n</html>\n"));
        assert!(!html.contains("<script"), "{html}");
        assert!(!html.contains("http://"), "{html}");
        assert!(!html.contains("https://"), "{html}");
        assert!(!html.contains("<img"), "{html}");
        for id in ["masthead", "run_stats", "data_version", "gen_stats"] {
            assert!(html.contains(&format!(" id=\"{id}\"")), "{id}");
        }
        for id in CHART_IDS {
            assert!(html.contains(&format!("<h3 id=\"{id}\">")), "{id}");
            let kind = shape(id).kind.name();
            assert!(
                html.contains(&format!(
                    "<div id=\"{id}_{kind}\" style=\"width: 800px; height: {}px\"><svg ",
                    shape(id).height
                )),
                "{id}"
            );
            assert_eq!(
                html.contains(&format!("<table id=\"{id}_table\"")),
                shape(id).table,
                "{id}"
            );
        }
        assert!(html.contains("<svg width=\"800\" height=\"400\""));
        assert!(html.contains("<td>Species</td><td>homo_sapiens</td>"));
        assert!(
            html.contains("<td>Command line options</td><td><pre>vep --cache --offline</pre></td>")
        );
        assert!(html.contains("<td>Lines of input read</td><td>3</td>"));
        assert!(html.contains("No data version information available."));
        // One key holding every count is a full circle, not an arc.
        assert!(
            html.contains("<circle cx=\"200.0\" cy=\"200.0\" r=\"190.0\" fill=\"#ffd700\"><title>missense_variant: 3 (100.0%)</title></circle>"),
            "{html}"
        );
        assert!(html.contains("<circle cx=\"100.0\" cy=\"100.0\" r=\"90.0\" fill=\"red\"><title>deleterious: 6 (100.0%)</title></circle>"));
        assert!(html.contains("<circle cx=\"100.0\" cy=\"100.0\" r=\"90.0\" fill=\"green\"><title>benign: 6 (100.0%)</title></circle>"));

        let mut without_flags = full_report();
        without_flags
            .charts
            .retain(|chart| chart.id != "sift" && chart.id != "polyphen" && chart.id != "chr_21");
        let html = without_flags.html();
        for id in ["sift", "polyphen", "chr_21"] {
            assert!(!html.contains(&format!("id=\"{id}\"")), "{id}");
            assert!(!html.contains(&format!("id=\"{id}_")), "{id}");
        }
        let plain = Report::new(&RunStats::new(), &context());
        let html = plain.html();
        assert_eq!(
            plain
                .charts
                .iter()
                .map(|c| c.id.as_str())
                .collect::<Vec<_>>(),
            ["var_class", "var_cons", "consequences", "coding"]
        );
        assert!(!html.contains("<h3 id=\"chr\">"));
        assert!(!html.contains("id=\"sift"));
        assert!(html.contains("<td>Lines of input read</td><td>0</td>"));
        assert!(html.contains("<td>Variants processed</td><td>0</td>"));
    }

    /// A chart's drawing is a function of its rows alone: two renders of one
    /// fixture are equal byte for byte, and the fixture's three kinds draw the
    /// shapes expected of them.
    #[test]
    fn pie_bar_and_area_svg_are_deterministic_for_fixed_counts() {
        let report = fixed_report();
        let first = report.html();
        assert_eq!(first, fixed_report().html());
        let [pie, bar, area] =
            [0, 1, 2].map(|i| chart_svg(&report.charts[i], &shape(&report.charts[i].id)));
        assert_eq!(pie, chart_svg(&report.charts[0], &shape("var_cons")));
        assert_eq!(pie.matches("<path ").count(), 3);
        assert_eq!(
            pie.matches("<rect ").count(),
            3,
            "one legend swatch per key"
        );
        assert!(
            pie.contains("fill=\"#ffd700\""),
            "missense_variant's colour"
        );
        assert!(pie.contains("fill=\"#02599c\""), "intron_variant's colour");
        assert!(
            pie.contains(&format!("fill=\"{}\"", PALETTE[2])),
            "the third palette colour for `other`"
        );
        assert!(pie.contains("<title>intron_variant: 5 (50.0%)</title>"));
        // The first slice, 3 of 10, runs clockwise from the top through 0.6 of a
        // half-turn: its far point is (200 + 190 cos 0.1 pi, 200 + 190 sin 0.1 pi).
        assert!(
            pie.contains(
                "<path d=\"M 200.0 200.0 L 200.0 10.0 A 190.0 190.0 0 0 1 380.7 258.7 Z\" fill=\"#ffd700\""
            ),
            "{pie}"
        );
        // The second, 5 of 10, is the half-turn that follows; at exactly half the
        // arc is not the large one.
        assert!(
            pie.contains("A 190.0 190.0 0 0 1 19.3 141.3 Z\" fill=\"#02599c\""),
            "{pie}"
        );
        // Two bars share 740 px: a 370 px slot each, the bar 70 percent of it,
        // the taller filling the 160 px plot and the other 3/7 of it.
        assert_eq!(bar.matches("<rect ").count(), 2);
        assert!(
            bar.contains("<rect x=\"105.5\" y=\"10.0\" width=\"259.0\" height=\"160.0\""),
            "{bar}"
        );
        assert!(
            bar.contains("<rect x=\"475.5\" y=\"101.4\" width=\"259.0\" height=\"68.6\""),
            "{bar}"
        );
        assert!(bar.contains(
            "<text x=\"235.0\" y=\"182.0\" font-size=\"10\" text-anchor=\"middle\">21</text>"
        ));
        assert!(!bar.contains("rotate("));
        assert_eq!(area.matches("<polygon ").count(), 1);
        assert!(
            area.contains(
                "points=\"50.0,160.0 50.0,122.5 296.7,160.0 543.3,10.0 790.0,85.0 790.0,160.0\""
            ),
            "{area}"
        );
        assert!(area.contains(">Position (mb)</text>"));
        for key in ["0", "1", "2", "3"] {
            assert!(
                area.contains(&format!("text-anchor=\"middle\">{key}</text>")),
                "{key}"
            );
        }
        assert!(!area.contains("<table"));
    }

    /// `Stats.pm:855-871`: the four fixed entries, then each chart in order except
    /// the distribution charts (`no_link`, `:593`).
    #[test]
    fn side_menu_lists_every_linked_chart_id() {
        let report = full_report();
        let entries = side_menu_entries(&report);
        let ids: Vec<&str> = entries.iter().map(|(id, _)| id.as_str()).collect();
        assert_eq!(
            ids,
            [
                "masthead",
                "run_stats",
                "data_version",
                "gen_stats",
                "var_class",
                "var_cons",
                "consequences",
                "coding",
                "sift",
                "polyphen",
                "chr",
                "protein",
            ]
        );
        let html = report.html();
        assert!(html.contains("<li><a href=\"#masthead\">Top of page</a></li><li><a href=\"#run_stats\">VEP run statistics</a></li>"));
        assert!(html.contains("<li><a href=\"#protein\">Position in protein</a></li></ul>"));
        assert!(!html.contains("href=\"#chr_21\""));
        assert!(html.contains("<h3 id=\"chr_21\">Distribution of variants on chromosome 21</h3>"));
    }

    /// Every element id `dump_html` writes is on the page, in its order: the
    /// masthead, the three tables, then per chart its heading, its chart element
    /// (`<id>_<type>`) and its table (`<id>_table`) unless `no_table`.
    #[test]
    fn html_has_every_section_id_of_stats_pm() {
        let report = full_report();
        let html = report.html();
        let mut expected = vec![
            "masthead".to_string(),
            "run_stats".to_string(),
            "data_version".to_string(),
            "gen_stats".to_string(),
        ];
        for chart in &report.charts {
            let shape = shape(&chart.id);
            expected.push(chart.id.clone());
            expected.push(format!("{}_{}", chart.id, shape.kind.name()));
            if shape.table {
                expected.push(format!("{}_table", chart.id));
            }
        }
        let ids: Vec<&str> = html
            .match_indices(" id=\"")
            .map(|(at, _)| {
                let rest = &html[at + 5..];
                &rest[..rest.find('"').unwrap()]
            })
            .collect();
        assert_eq!(ids, expected);
        assert_eq!(
            expected,
            [
                "masthead",
                "run_stats",
                "data_version",
                "gen_stats",
                "var_class",
                "var_class_pie",
                "var_class_table",
                "var_cons",
                "var_cons_pie",
                "var_cons_table",
                "consequences",
                "consequences_pie",
                "consequences_table",
                "coding",
                "coding_pie",
                "coding_table",
                "sift",
                "sift_pie",
                "sift_table",
                "polyphen",
                "polyphen_pie",
                "polyphen_table",
                "chr",
                "chr_bar",
                "chr_table",
                "chr_21",
                "chr_21_area",
                "protein",
                "protein_bar",
            ]
        );
    }

    /// `generate_chart_data` (`Stats.pm:505-609`): each chart's type, headers,
    /// height and flags; a height of 200 where the hash names one, 400 for the
    /// other pies, 200 for every bar and area chart.
    #[test]
    fn chart_shapes_follow_generate_chart_data() {
        let var_class = shape("var_class");
        assert_eq!((var_class.kind, var_class.height), (Kind::Pie, 200));
        assert_eq!(var_class.header, ["Variant class", "Count"]);
        assert!(var_class.colours.is_empty());
        for id in ["var_cons", "consequences", "coding"] {
            let s = shape(id);
            assert_eq!(
                (s.kind, s.height, s.table, s.linked),
                (Kind::Pie, 400, true, true)
            );
            assert_eq!(s.header, ["Consequence type", "Count"]);
            assert_eq!(colour_of("stop_gained", 0, s.colours), "#ff0000");
            assert_eq!(colour_of("unlisted", 1, s.colours), PALETTE[1]);
        }
        assert_eq!(shape("sift").height, 200);
        assert_eq!(colour_of("tolerated", 0, shape("sift").colours), "green");
        assert_eq!(
            colour_of("possibly damaging", 0, shape("polyphen").colours),
            "orange"
        );
        assert_eq!(
            colour_of("possibly_damaging", 0, shape("polyphen").colours),
            PALETTE[0]
        );
        let chr = shape("chr");
        assert_eq!((chr.kind, chr.table, chr.linked), (Kind::Bar, true, true));
        assert_eq!(chr.header, ["Chromosome", "Count"]);
        let protein = shape("protein");
        assert_eq!(
            (protein.kind, protein.table, protein.linked),
            (Kind::Bar, false, true)
        );
        assert_eq!(protein.header[0], "Position in protein (percentile)");
        let area = shape("chr_X_1");
        assert_eq!(
            (area.kind, area.table, area.linked),
            (Kind::Area, false, false)
        );
        assert_eq!(area.header, ["Position (mb)", "Count"]);
    }

    /// A key is written as character data wherever the page prints it; a bar
    /// chart of more than 24 keys turns its labels.
    #[test]
    fn keys_are_escaped_and_many_bar_labels_turn() {
        let report = Report {
            charts: vec![Chart {
                id: "chr".to_string(),
                title: "Variants by chromosome".to_string(),
                rows: (0..30).map(|i| (format!("c<{i}>&"), 1)).collect(),
            }],
            ..fixed_report()
        };
        let html = report.html();
        assert!(html.contains("<td>c&lt;0&gt;&amp;</td><td>1</td>"));
        assert!(html.contains("c&lt;29&gt;&amp;</text>"));
        assert!(!html.contains("c<0>"));
        assert_eq!(html.matches("rotate(-60 ").count(), 30);
        assert_eq!(escape("a \"b\" & <c>"), "a &quot;b&quot; &amp; &lt;c&gt;");
    }
}
