// =============================================================================
// Supplementary Information for
// Evaluating Low-Pass Whole Genome Sequencing as a Cost-Effective Method
// for Copy Number Variant Detection
//
// Three parts, in the order journals expect them: Supplementary Notes
// (self-contained technical analyses cited from the main text), Supplementary
// Figures, Supplementary Tables. Each part starts on a new page, so the file can
// be split into three documents later by cutting at the part headings.
//
// Numbering. Notes are numbered (S1, S2, ...) because they are a new series.
// Figures and tables keep the named placeholders the main text uses
// ("Supplementary Figure Breakpoint Precedence") until the numbering pass in
// PLAN.md assigns every supplemental a number at once; auto-numbering them here
// would collide with the Supplementary Table S1--S3 the Methods already cite.
//
// The preamble mirrors cnv-benchmark-paper.typ; keep the two in step.
// Compile from the repository root: typst compile --root . manuscript/supplementary.typ
// Prose uses semantic line breaks: one sentence per line, no column limit.
// =============================================================================

#set document(
  title: "Supplementary Information: Evaluating Low-Pass Whole Genome Sequencing as a Cost-Effective Method for Copy Number Variant Detection",
  author: ("Lionel Sequeira", "Thomas V Fernandez", "Gary A Heiman", "Jinchuan Xing"),
)

// Supplementary pages are numbered S1, S2, ... so they cannot be confused with
// the main text's pages when both are cited.
#set page(paper: "us-letter", margin: 1in, numbering: "S1")
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.72em, spacing: 1.1em)
#set heading(numbering: none)

#show heading.where(level: 1): it => block(above: 1.8em, below: 0.9em)[
  #set text(size: 15pt, weight: "bold")
  #it.body
]
#show heading.where(level: 2): it => block(above: 1.4em, below: 0.7em)[
  #set text(size: 12.5pt, weight: "bold")
  #it.body
]
#show heading.where(level: 3): it => block(above: 1.2em, below: 0.6em)[
  #set text(size: 11pt, weight: "bold", style: "italic")
  #it.body
]
#show link: set text(fill: rgb("#1a4f8a"))
#show table.cell.where(y: 0): strong
#set table(stroke: (x, y) => (
  top: if y <= 1 { 0.6pt } else { 0pt },
  bottom: 0.6pt,
))

// -----------------------------------------------------------------------------
// Helpers, as in the main text
// -----------------------------------------------------------------------------

#let c(body) = [#box[\[#body\]]]

#let cap(label, body) = block(width: 100%, above: 0.7em)[
  #set text(size: 9.5pt)
  #strong[#emph[#label]] #emph[#body]
]

// -----------------------------------------------------------------------------
// Title block
// -----------------------------------------------------------------------------

#align(center)[
  #text(size: 13pt, weight: "bold")[Supplementary Information for]

  #v(0.3em)

  #block(width: 100%)[
    #text(size: 15pt, weight: "bold")[
      Evaluating Low-Pass Whole Genome Sequencing as a Cost-Effective
      Method for Copy Number Variant Detection
    ]
  ]

  #v(0.6em)

  Lionel Sequeira, Thomas V Fernandez, Gary A Heiman, Jinchuan Xing
]

#v(1.2em)

#outline(title: [Contents], depth: 2, indent: 1.2em)

// =============================================================================
#pagebreak()
= Supplementary Notes

== Supplementary Note S1: Breakpoint precedence under union merging <note-s1>

A consensus call spans the minimum start and the maximum end of its member calls, so each of its boundaries comes from whichever caller reached furthest in that direction.
This note identifies which caller that is.

=== Approach

We used the 4,836 components of the 30x 2-of-3 consensus call set, each carrying at least two callers, at the adopted 50% reciprocal overlap.
Each caller was represented in a component by the minimum start and maximum end of its member calls, and for each boundary we recorded whether a single caller reaches it alone or two callers tie for it.
The analysis was repeated at 6x, 4x, and 2x.

=== Results

CNVpytor and GATK-gCNV both place breakpoints on the 1 kb bin grid, and they reported identical coordinates at both ends in 61.1% of the 3,339 components they share; Delly coincided exactly with neither (Supplementary Figure Breakpoint Precedence A).
Delly was the wider call in 73.9% of the components it shares with GATK-gCNV, extending a median of 263 bp beyond it at the start and 251 bp at the end, but was close to symmetric with CNVpytor, wider in 56.2% of the components they share.

Because Delly almost never ties with another caller, these small differences decide the boundary.
Delly alone reached the union start in 57.6% and the union end in 51.9% of the components it appears in, against 21.3% and 28.8% for CNVpytor and 8.6% and 6.7% for GATK-gCNV (Supplementary Table Breakpoint Precedence).
The precedence is strongest on deletions, where 64.1% of components carry at least one boundary from Delly alone.
On duplications CNVpytor alone sets the start and end in 54.4% and 54.1% of the components it appears in, against 31.4% and 34.4% for Delly, although 50.4% of duplication components still carry one boundary from Delly alone (Supplementary Figure Breakpoint Precedence B).
Delly holds the largest share of sole boundaries at every coverage.

=== Conclusion

Delly takes precedence under a union merge not because its calls are larger, but because the two read-depth callers share a bin grid and often agree exactly, which leaves Delly's off-grid breakpoints to decide the boundary.

// =============================================================================
#pagebreak()
= Supplementary Figures

#figure(
  image("/results/breakpoint_precedence/breakpoint_precedence.png", width: 100%)
)

#cap("Supplementary Figure Breakpoint Precedence:")[
  Breakpoint precedence among the three callers in the 30x 2-of-3 consensus call set (Supplementary Note S1).
  (A) Outward extension of the first-named caller over the second at the start and end of every component carrying both, positive where the first-named caller extends further; boxes span the interquartile range and whiskers the 5th to 95th percentiles.
  The fraction of shared components in which the two callers report identical coordinates at both ends is given beneath each pair.
  (B) Of the deletion and duplication components each caller appears in, the percentage in which it alone reaches the union start or the union end; boundaries tied between two callers are not counted.
]

// =============================================================================
#pagebreak()
= Supplementary Tables

#block(width: 100%)[
#set text(hyphenate: false, size: 9.5pt)
#show table.cell.where(y: 1): strong
#table(
  columns: (auto, 1fr, auto, auto, auto, auto, auto, auto),
  align: (left, left, right, right, right, right, right, right),
  stroke: none,
  table.hline(stroke: 0.6pt),
  table.header(
    table.cell(rowspan: 2, align: left + bottom)[Coverage],
    table.cell(rowspan: 2, align: left + bottom)[Caller],
    table.cell(rowspan: 2, align: right + bottom)[Components],
    table.cell(colspan: 2, align: center)[Holds (%)],
    table.cell(colspan: 2, align: center)[Holds alone (%)],
    table.cell(rowspan: 2, align: right + bottom)[Span Ratio],
    table.hline(start: 3, end: 5, stroke: 0.3pt),
    table.hline(start: 5, end: 7, stroke: 0.3pt),
    [Start], [End], [Start], [End],
  ),
  table.hline(stroke: 0.6pt),
  [30x], [CNVpytor], [4,207], [53.6], [60.6], [21.3], [28.8], [0.967],
  [], [Delly], [3,897], [57.6], [52.0], [57.6], [51.9], [0.967],
  [], [GATK-gCNV], [3,968], [42.7], [40.4], [8.6], [6.7], [0.890],
  table.hline(stroke: 0.3pt),
  [6x], [CNVpytor], [1,443], [55.9], [60.7], [27.7], [31.8], [0.972],
  [], [Delly], [1,261], [58.2], [55.8], [58.2], [55.8], [0.979],
  [], [GATK-gCNV], [1,400], [41.1], [38.9], [12.0], [9.1], [0.875],
  table.hline(stroke: 0.3pt),
  [4x], [CNVpytor], [944], [53.8], [56.6], [37.7], [38.2], [0.967],
  [], [Delly], [889], [57.3], [53.4], [57.1], [53.4], [0.975],
  [], [GATK-gCNV], [713], [35.8], [39.7], [14.3], [15.4], [0.868],
  table.hline(stroke: 0.3pt),
  [2x], [CNVpytor], [600], [59.0], [62.5], [43.7], [42.8], [0.978],
  [], [Delly], [508], [53.3], [47.2], [53.3], [47.2], [0.975],
  [], [GATK-gCNV], [511], [36.8], [43.8], [18.8], [20.7], [0.867],
  table.hline(stroke: 0.6pt),
)
]

#cap("Supplementary Table Breakpoint Precedence:")[
  Breakpoint precedence of each caller in the 2-of-3 consensus call set at each coverage, at the adopted 50% reciprocal overlap (Supplementary Note S1).
  Components is the number of consensus components carrying at least two callers in which the caller appears.
  Holds is the percentage of those components in which the caller's call reaches the union start or end, counting boundaries tied with another caller; Holds alone counts only boundaries no other caller reaches.
  Span Ratio is the median ratio of the caller's own extent within a component to the component's union span.
]
