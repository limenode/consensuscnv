// =============================================================================
// Supplementary Information for
// Evaluating Low-Pass Whole Genome Sequencing as a Cost-Effective Method
// for Copy Number Variant Detection
//
// Three parts: Supplementary Notes, Supplementary Figures, Supplementary Tables.
// Each part starts on a new page, so the file can be split into three documents
// later by cutting at the part headings.
// =============================================================================

#set document(
  title: "Supplementary Information: Evaluating Low-Pass Whole Genome Sequencing as a Cost-Effective Method for Copy Number Variant Detection",
  author: ("Lionel Sequeira", "Thomas V Fernandez", "Gary A Heiman", "Jinchuan Xing"),
)

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
#show heading.where(level: 4): it => block(above: 1.0em, below: 0.5em)[
  #set text(size: 10.5pt, weight: "regular", style: "italic")
  #it.body
]
#show link: set text(fill: rgb("#1a4f8a"))
#show table.cell.where(y: 0): strong
#set table(stroke: (x, y) => (
  top: if y <= 1 { 0.6pt } else { 0pt },
  bottom: 0.6pt,
))

// -----------------------------------------------------------------------------
// Helpers
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

== Supplementary Note S1: Parameterizing the comparison <note-s1>

Four user-defined parameters govern the comparison: the padding applied to the benchmark records before the three benchmark sets are merged, the size floor below which no event is scored on either side, the reciprocal-overlap threshold at which calls from different callers are merged into a consensus call, and the reciprocal-overlap threshold at which a query call is credited against a benchmark interval.
Each was bounded on physical grounds in the Methods of the main text.
What follows reports what each one does to the comparison across its range, taking them in the order the pipeline applies them and holding the other three at their adopted values throughout: zero padding, a 1 kb floor, and a reciprocal overlap of 0.5 for both consensus construction and classification.

=== Benchmark Padding

The merged benchmark is built from three sets produced by different technologies and assemblies, so their breakpoints for the same variant do not coincide exactly, and padding is a standard remedy for that.
The merged benchmark has a median interval of 132 bp (IQR \[71--338\] bp), so padding of even a few hundred base pairs is comparable in size to the intervals it is meant to reconcile.
We swept padding from 0 to 100 kb and re-evaluated the comparison at every point (Supplementary Figure Benchmark Padding).

The truth set moves in two directions at once (Supplementary Figure Benchmark PaddingA).
Counted over all sizes the merged benchmark shrinks, from 180,064 intervals unpadded to 164,065 at 1 kb of padding and 96,743 at 100 kb, as records are absorbed into their neighbors.
Counted at or above the 1 kb floor it grows, from 23,193 intervals to 25,562 and then to 37,626 across the same range.
The second movement is the consequential one, since it is the intervals that clear the floor which form the recall denominator, and the growth is not a discovery of additional variants.
It is runs of sub-kilobase records fused into single intervals long enough to clear the 1 kb floor.

Those intervals can be counted directly.
An interval was labeled manufactured when the longest benchmark record inside it is itself shorter than 1 kb, so that the interval exists in the truth set only because padding joined the run.
Manufactured intervals are 413 of the 23,193 intervals in the unpadded truth set (1.8%), 3,084 of 25,562 at 1 kb of padding (12.1%), and 19,716 of 37,626 at 100 kb (52.4%).
They are also almost never recovered.
At 1 kb of padding the 30x 2-of-3 consensus finds 19.4% of the native intervals and 0.13% of the manufactured ones, a separation of more than two orders of magnitude that holds across the whole sweep (Supplementary Figure Benchmark PaddingB).
Padding therefore adds to the recall denominator a population that is by construction outside the resolution of every caller in this study.

The effect on the metrics follows from that arithmetic alone (Supplementary Figure Benchmark PaddingD).
Between zero padding and the 1 kb cap, the number of benchmark intervals the 2-of-3 consensus recovers rises from 4,344 to 4,360, an increase of 0.4%, while the denominator rises from 23,193 to 25,562, an increase of 10.2%.
Recall falls from 0.187 to 0.171 and F1 from 0.310 to 0.287 as a result, with no change in detection behind either.
Recall and F1 sit at their maxima at the bottom of the range for all six 30x call sets and decline monotonically above roughly 100 bp of padding (Supplemental Figure Benchmark Padding Profiles).

Precision behaves differently, and it is the one place in the sweep where padding does what it is intended to do.
It improves slightly over the first two kilobases, reaching 0.902 at 1,823 bp of padding against 0.898 unpadded for the 2-of-3 consensus, because fusing benchmark fragments into single intervals converts a small number of boundary near-misses into matches.
The maxima for the other five call sets fall between 290 bp and 2,154 bp of padding, so the effect is real and consistently located, but it is worth less than one percentage point of precision.
Above a few kilobases precision falls with everything else, reaching 0.485 at 100 kb.

The matching itself is unaffected throughout.
At the 0.5 classification threshold used in this study, no query call matched more than one benchmark interval and no benchmark interval was split across more than one query call at any padding in the sweep, as the geometry of the threshold requires.
The matching is therefore strictly one-to-one over the entire range, and none of the movement above is an artifact of a single call being credited against several intervals at once.
Structure appears only in the permissive regime, and there it reverses direction (Supplementary Figure Benchmark PaddingC).
At a 0.1 threshold the number of query calls spanning more than one benchmark interval falls from 42 to none as fragments fuse into single partners, while the number of benchmark intervals split across more than one query call rises from 17 to a maximum of 41, both against totals of between 3,130 and 4,535 matched pairs.

Padding was therefore set to zero rather than capped.
There is no interior optimum available to select: recall and F1 are maximal at the bottom of the range, and the sub-percentage-point gain in precision available at a kilobase is bought with a 1.7 percentage point loss of recall.
Zero padding is not the same as disabling the operation, since it still bridges intervals that touch exactly; that difference amounts to 67 intervals out of 180,131, and it is retained because two benchmark records abutting at a shared coordinate describe one variant rather than two.
Padding is nonetheless carried through the joint parameter grid over \[0, 1 kb\], so that its interaction with the other three parameters is measured rather than assumed.

=== Size Floor

With padding fixed at zero, the merged benchmark holds 180,064 intervals across the thirteen samples, of which only 12.9% reach 1 kb and 2.7% reach 10 kb.
The preceding subsection took the 1 kb floor as given in order to isolate the effect of padding; this one asks whether that value is the right one.
The six 30x query call sets are between 7 and 75 times smaller than the benchmark and, because every caller was run at a 1 kb bin size, hold almost nothing below that width (Supplementary Figure Size Floor).
Scoring them against the benchmark without a size restriction therefore primarily measures the mismatch in resolution between the benchmark and the callers rather than an effect of sequencing depth, which is what this study set out to isolate.
We swept the floor symmetrically over both sides and examined how the comparison behaves as it rises (Supplementary Figure Size Floor).

Three features of this sweep together identify the usable domain.

First, the benchmark loses intervals far faster than any query call set.
At the bottom of the range the truth set outnumbers the largest query set more than sevenfold, but it falls below the 1-of-3 consensus set above a 943 bp floor and below CNVpytor above 1,954 bp, while the other four call sets remain smaller than the truth set throughout the sweep (Supplementary Figure Size FloorA).
Above those points the comparison has inverted for the two largest sets: they report more CNVs than the benchmark contains, and their precision is bounded by the size of the truth set rather than by the accuracy of the calls.

Second, recall is constrained by a ceiling that is a property of the two call set sizes rather than of detection (Supplementary Figure Size FloorB).
At an unrestricted floor (size floor = 0 bp) the maximum attainable recall is 0.134 for the 1-of-3 set and 0.013 for the 3-of-3 set, so even a caller that matched a distinct benchmark interval with every single one of its calls could not exceed those values.
Raising the floor lifts the ceiling for every call set, but it does so unevenly: the 1-of-3 and CNVpytor sets saturate at 1.0 above floors of 943 bp and 1,954 bp respectively, while the other four never do, with maxima of 0.60 for Delly, 0.38 for GATK-gCNV, 0.31 for the 2-of-3 set, and 0.18 for the 3-of-3 set.
Recall is therefore interpretable as a detection measurement only below the point at which a given call set saturates.

Third, precision is flat from 1 bp to approximately 1 kb for every call set and declines above it (Supplementary Figure Size FloorC).
The flat region is the direct consequence of the bin size: over that range the floor removes benchmark intervals almost exclusively, because the callers had produced essentially nothing there to remove, so the query sets and their precision are unchanged while the truth set falls from 180,064 intervals to 23,193.
Above 1 kb the floor begins to remove query calls as well, and precision falls for every set, most steeply for Delly (0.642 to 0.208 between the unrestricted case and a 100 kb floor), by about half for CNVpytor and the 1-of-3 consensus (0.419 to 0.217 and 0.399 to 0.183), and most gradually for GATK-gCNV, which is the most size-stable of the callers (0.603 to 0.418 over the same range, the upper end of which lies beyond the 64.6 kb floor at which its curve is cut for low counts).

Taken together, these place the usable domain immediately above the bin size.
We fixed the floor at 1 kb, chosen on the physical grounds that no caller in this study can resolve a CNV narrower than its bin.
The sweep supports that choice: F1 reaches its maximum between 1,689 bp and 3,027 bp for all six call sets, with the 1-of-3 consensus highest at 0.418 and the 2-of-3 consensus at 0.395 (Supplementary Figure Size FloorD).
A 1 kb floor therefore sits just below the empirical optimum for every call set simultaneously.

This floor is applied to every call set and to the benchmark for all analyses that follow, at all four coverages, and it is the value at which the floor is held while the two overlap thresholds are profiled below.
It was selected using 30x data only, which is the arm with the finest resolution and therefore the most permissive: a floor set there admits calls at 2x that fall below the resolution attainable at that depth, biasing the comparison against the low-coverage hypothesis rather than in its favor.

=== Consensus reciprocal overlap threshold
A consensus component is a union of all pairwise overlaps that clear the reciprocal overlap threshold.
We swept the threshold from 0.05 to 0.95 in steps of 0.05 (Supplementary Figure Consensus Overlap).

The behavior for the number of calls greater than 1kb in length differs across consensus stringencies (Supplementary Figure Consensus OverlapA).
Raising the threshold splits components between callers that fail to meet threshold requirements. Each split adds a call to lower stringency sets while removing the agreement that contributed to the counts in the levels above it.
Between 0.05 and 0.95 the 1-of-3 set grows from 21,847 to 28,413 calls above the floor, while the 2-of-3 set falls from 5,623 to 2,369 and the 3-of-3 set from 2,666 to 369.

Taking the ratio of a component's span to the span of the longest single call inside it, the median is 1.00 at every threshold and for every level, the 95th percentile never exceeds 1.16, and the largest ratio anywhere in the sweep is 2.83.
Components at the permissive end are therefore sets of calls that agree, not loci collapsed into intervals no caller reported.

What the threshold does instead is trade the two sides of the comparison against each other (Supplementary Figure Consensus OverlapB, C).
For the 2-of-3 consensus, precision rises from 0.827 to 0.936 across the range while recall falls from 0.201 to 0.096.
For the 3-of-3 consensus precision is already near its ceiling and barely moves, from 0.945 to a maximum of 0.983 at 0.65, while recall falls from 0.109 to 0.015.
The 1-of-3 set is the exception: its recall is roughly flat, 0.385 to 0.399, throughout; its precision falls to a minimum of 0.391 at 0.45 before rising to 0.470.

Over the lower half of the range those two movements nearly cancel (Supplementary Figure Consensus OverlapD).
Between 0.05 and 0.50 the F1 of the 2-of-3 consensus varies by 0.013, from 0.323 at 0.05 to 0.310 at 0.50, while its precision gains 7.1 percentage points.
The F1 maxima themselves are shallow and disagree between levels, falling at 0.05 for the 2-of-3 and 3-of-3 sets and at the top of the range for the 1-of-3 set.

We chose to adopt a threshold of 0.5 for subsequent analyses.
It is the conventional reciprocal criterion, it is the smallest threshold at which neither member of a merged pair can be more than twice the size of the other, and it costs the 2-of-3 consensus 4.0% of its attainable F1 while gaining 7.1 of the 10.9 percentage points of precision available across the whole range.

=== Classification reciprocal overlap threshold

Unlike the other three parameters this one changes neither call set, and instead strictly influences binary classification calculations.
We swept the classification reciprocal overlap threshold, which defines the overlap required to register matches between query and truth sets as valid, from 0 to 0.99 in steps of 0.01 (Supplementary Figure Classification Overlap).

All three metrics decline monotonically and the parameter has no interior optimum (Supplementary Figure Classification OverlapA).
For the 2-of-3 consensus, precision falls from 0.934 at a threshold of zero to 0.898 at 0.5 and 0.440 at 0.9, and F1 from 0.326 to 0.310 to 0.152.
The choice is therefore not between values that perform differently but between definitions of what a match is required to mean, and the informative quantity is the match topology rather than the metrics.

At a threshold of zero a single shared base pair is a match, and one query call is credited against as many as 15 benchmark intervals, depending on the call set.
The number of query calls with more than one partner falls to zero between 0.43 and 0.46 depending on the call set, ahead of the 0.5 at which the geometry requires one-to-one matching (Supplementary Figure Classification OverlapB).
The merged benchmark is internally disjoint within each sample, chromosome, and variant type, so no query call can cover half of two of its intervals; the query sets carry no such guarantee, because components built at a reciprocal overlap of 0.5 may still overlap one another below that threshold.
The 1-of-3 consensus is the case in which that matters: it still splits 127 benchmark intervals across two query calls at a threshold of 0.5, and the count reaches zero only at 0.70 (Supplementary Figure Classification OverlapC).
For the other five call sets the matching is strictly one-to-one at 0.5, and the query-side and truth-side true-positive counts coincide exactly; for the 2-of-3 consensus they are both 4,344 at 0.5.

According to F1 scores, the 1-of-3 consensus is the highest-scoring set at every threshold up to 0.82 and CNVpytor from 0.83 to 0.98, above which all six sets lie within 0.02 of one another; the 3-of-3 set is the lowest until 0.91 (Supplementary Figure Classification OverlapD).
We adopted 0.5, the conventional reciprocal-overlap criterion, which is also the smallest threshold at which no query call can be credited against two benchmark intervals at once.

== Supplementary Note S2: Variance-based sensitivity analysis and Pareto front <note-s2>

=== Methods

The profiles of Supplementary Note S1 vary one parameter at a time.
They describe the pipeline completely only if the effect of moving one parameter does not depend on where the other three are held, which is an assumption about the shape of the metric field rather than something the profiles themselves can show.
All four parameters were therefore also varied jointly, over the full factorial grid of the ranges fixed in the Methods of the main text, with padding taking \[0, 10, 25, 50, 100, 200, 400, 700, 1000\], the size floor taking \[0, 250, 500, 1000, 2000, 5000, 10000\], and both reciprocal-overlap thresholds taking \[0.05--0.95\] in steps of 0.05, and precision, recall, and F1 were evaluated at every combination together with the number of query and truth intervals entering it.
What follows asks how much of the variation in each metric the one-at-a-time reading accounts for, and where a joint reading is required instead.

==== Sensitivity Indices

The contribution of each parameter was quantified by variance-based sensitivity analysis (Sobol', 1993).
Writing $Y$ for a metric (precision, recall, f1) and $x_1, ..., x_4$ for the parameters, the first order Sobol index $S_i$ of parameter $i$ is the fraction of the metric's variance attributable to that parameter acting alone.

$ S_i = ("Var" (EE [Y | X_i])) / ("Var" (Y)) $

The total-order index $S_(T i)$ additionally includes every interaction involving $i$.

$ S_(T i) = 1 - ("Var" (EE [Y | bold(X)_(tilde i)])) / ("Var" (Y)) $

The difference $S_(T i) - S_i$ is the variance a parameter contributes only jointly with others.
Because the design is a complete factorial grid, each conditional expectation is a marginal mean over the grid and both indices were computed exactly.
The complete decomposition, comprising every term up to fourth order, sums to one, and this was verified numerically.

==== Additivity

The same decomposition expresses the metric field as a grand mean plus one univariate function per parameter:

$ Y approx mu + sum_i f_i (X_i) wide "with" wide f_i (x) = EE [Y | X_i = x] - mu $

The coefficient of determination of this additive model equals the sum of the first-order indices, and so quantifies directly the extent to which the one-at-a-time profiles of Supplementary Note S1 describe the joint behavior of the pipeline.

==== Dependence on the swept ranges

Sobol indices are variance ratios with respect to a distribution over the inputs.
They therefore describe sensitivity within the examined region and are not invariant to the choice of that region.
As such, indices were additionally computed over an intentionally over-wide grid extending to $10^6$ bp (Supplementary Table X) in order to compare physical arguments for the previously described parameter boundaries with empirical data.

==== Pareto front

Because precision and recall span very different ranges across the grid, F1 is close to a monotone function of recall alone and obscures the trade-off between the two.
Parameter settings were therefore also summarized by their Pareto front: a setting is dominated if another attains at least equal precision and recall and strictly exceeds it on one, and the front comprises the non-dominated settings.
The front answers a different question from an optimum.
It separates settings that buy precision at a real cost in recall from settings that are simply worse on both, so a dominated setting has no argument in its favor.
Dominance compares measured performance, so it is informative only between settings that score the same events under the same definition of a match.
Three of the four parameters do not satisfy that: the padding and the floor change which intervals enter the comparison, and the classification threshold changes what counts as a match, so for these the front records a change in the comparison rather than an improvement in the pipeline.
Fronts were therefore computed over the whole grid and again with the classification threshold, and then the size floor, held at the values adopted for the pipeline.

=== Results

The profiles of Supplementary Note S1 move one parameter with the other three held fixed.
To ask whether that reading survives elsewhere in the parameter field, all four were varied jointly over the full factorial grid of 22,743 settings and evaluated for each of the three consensus call sets at 30x.
The 2-of-3 consensus is reported here; the other two levels are given in Supplementary Table X.

The four parameters act very nearly independently (Supplementary Figure SensitivityA).
Their first-order indices sum to 0.94 for precision, 0.94 for recall and 0.95 for F1, so an additive model, a grand mean plus one curve per parameter, reproduces the joint field to within 6% of its variance.
The largest single interaction anywhere is between the size floor and the classification threshold, at 2.9% of the variance in F1 and 3.9% in precision.
The one-at-a-time profiles of Supplementary Note S1 are therefore not artifacts of where the remaining parameters were held.

The variance is not shared evenly (Supplementary Figure SensitivityA, B).
The size floor carries 76.6% of the variance in F1 and 76.6% in recall, and the classification threshold carries 81.5% of the variance in precision.
The benchmark padding accounts for no more than 0.1% of the variance in any of the three metrics.
The same ordering holds at the 1-of-3 and 3-of-3 consensus levels, with one shift: the consensus threshold's share of the variance in F1 rises from 0.2% for the 1-of-3 set to 15.1% for the 3-of-3 set, which follows the intuition that for higher stringencies, changes in the threshold more readily remove calls.

These shares per parameter describe the region examined, and are influenced by the ranges we evaluate.
Recomputed over the deliberately over-wide grid, in which the padding and the floor extend to $10^6$ bp, the padding's share of the variance in F1 rises from below 0.1% to 34.2% and the share carried by interactions rises from 4.9% to 25.5% (Supplementary Table X).
The padding's apparent inertness and the additivity of the field are properties of the ranges used, and are physically constrained around the values that we use for the subsequent analyses.

Across the grid, precision exceeds recall at every setting, and F1 is close to a monotone function of recall alone.
The Pareto front comprises 79 of the 22,743 settings, spanning recall from 0.048 to 0.301 and precision from 0.831 to 0.955 (Supplementary Figure SensitivityC).
Every setting on the front uses a classification threshold of 0.05, the loosest value tested in the grid.
The setting attaining the highest F1 of 0.442 uses zero padding, a 2 kb floor, a consensus threshold of 0.05, and a classification threshold of 0.05.
Both metrics decline monotonically in the classification threshold, so every value above the smallest is dominated; what this records is that a looser definition of a match credits more calls, not that the pipeline performs better under it.
The front is read with the crediting rule fixed for that reason.

Held at a classification threshold of 0.5, the front comprises 113 of the 1,197 remaining settings, and the adopted setting is one of them: no setting in the grid reaches both a higher precision and a higher recall (Supplementary Figure SensitivityD).
With the floor also held at the 1 kb the caller resolution fixes, the front comprises 61 of the 171 settings that remain and the adopted setting is again on it.
That front is traced by the consensus threshold, which takes every one of its nineteen values along it; the padding takes five of its nine, and moving from zero padding to 200 bp at the adopted consensus threshold exchanges 0.002 of recall for 0.002 of precision.
The choice of an operating point on this front is therefore a choice of consensus stringency, and the rate at which the two metrics exchange turns at the adopted value: raising the threshold from 0.05 to 0.50 gains 7.1 percentage points of precision for 1.3 of recall, while raising it from 0.50 to 0.95 gains a further 3.8 for 9.2.

== Supplementary Note S3: Breakpoint precedence under union merging <note-s3>

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

== Supplementary Note S4: Independence of the benchmark from the sequencing data <note-s4>

The 1000 Genomes high-coverage SV call set, one of the three benchmark sources, was called from the same 30x alignments as the sequencing call sets in this study.
An artifact that both pipelines draw from those reads could therefore be confirmed by the benchmark, which would score a sequencing call as a true positive while being correctly left undetected by the SNP Array.
This note bounds how much of the comparison could rest on that shared source.

=== Approach

Two measurements were taken at the adopted parameters, separately for deletions and duplications.
First, against the full merged benchmark, the benchmark intervals each call set recovered were divided by whether they are supported by the 1000 Genomes set alone or by at least one of the two long-read sources; recoveries of the first kind are the only ones the shared source could have produced.
Second, every call set was scored again against a benchmark merged from HGSVC3 and ONT Vienna alone, neither of which shares reads with the sequencing data.

=== Results

Deletions depend little on the shared source (Supplementary Table Benchmark Independence).
Of the 16,742 deletion intervals in the full benchmark above the size floor, 8.3% are supported by the 1000 Genomes set alone, and these make up 1.6--1.8% of the deletions recovered by each 2-of-3 consensus call set against 3.0% of those recovered by the SNP array.
Against the long-read-only benchmark of 14,330 deletion intervals, deletion precision fell by 0.030--0.052 for the four consensus call sets and by 0.029 for the array, so every consensus call set retained a deletion precision between 0.838 and 0.896 against 0.570 for the array.
Deletion recall rose slightly for every call set, since the smaller benchmark is the denominator.

Duplications cannot be tested the same way.
Of the 6,451 duplication intervals in the full benchmark above the floor, 89.1% are supported by the 1000 Genomes set alone, and the long-read-only benchmark retains 469 duplication intervals, too few for its duplication metrics to describe the same truth set.
The duplications recovered by each call set rest on the 1000 Genomes set alone in 80.2--85.1% of cases for the consensus call sets and in 87.5% for the array.

=== Conclusion

The advantage of the sequencing call sets over the array on deletions does not depend on the benchmark source that shares their reads, and the array relies on that source slightly more than they do.
The duplication results rest largely on the 1000 Genomes set; the array's equal reliance on it argues against that set merely echoing artifacts of the shared reads, but an independent duplication truth set would be needed to confirm it.

// =============================================================================
#pagebreak()
= Supplementary Figures

#figure(
  image("/results/parameterization/benchmark_padding.png", width: 100%)
)

#cap("Supplementary Figure Benchmark Padding:")[
  Effect of benchmark padding on the truth set and on the comparison, swept from 0 to 100 kb.
  Padding is applied to both ends of every benchmark record before the three benchmark sets are merged, and the 1 kb size floor is applied afterwards.
  (A) Benchmark intervals: all merged intervals, those reaching 1 kb, the manufactured subset of the latter, and those found by the 30x 2-of-3 consensus call set.
  An interval is manufactured when the longest benchmark record within it is itself shorter than 1 kb, so that it clears the floor only through fusion.
  (B) Percentage of benchmark intervals found by the same call set, with native and manufactured intervals scored separately.
  (C) Departures from one-to-one matching, at a permissive 0.1 classification threshold.
  The fixed 0.5 threshold was not used for this comparison because a query call cannot reach half of two benchmark intervals that do not themselves overlap.
  (D) Precision, recall, and F1 for the 30x 2-of-3 consensus at the 0.5 threshold.
  The dotted vertical line marks the 1 kb cap applied to padding in the joint parameter grid.
  Panels C and D are drawn for the 2-of-3 consensus; the same profiles for all six 30x call sets, together with the benchmark size distribution across the sweep, are given in Supplemental Figure Benchmark Padding Profiles.
]

#pagebreak()

#figure(
  image("/results/size_floor/detectable_size_domain_pub.png", width: 100%)
)

#cap("Supplementary Figure Size Floor:")[
  Effect of a size floor applied symmetrically to 30x WGS-derived query call sets and the merged benchmark call set, swept over \[1 bp, 100 kb\] at 80 logarithmically spaced points.
  (A) Number of intervals surviving the floor in each call set, with the merged benchmark as a black dash-dot line.
  (B) Recall, with the maximum attainable recall (query calls divided by truth intervals, capped at one) as a dashed line of the same color.
  (C) Precision.
  (D) F1, with each call set's maximum marked.
  The shaded band marks 1--5 kb, spanning every F1 maximum; the dotted vertical line marks the 1 kb floor adopted for all subsequent analyses.
  In panels C and D, each curve is drawn only while the call set behind it retains at least 100 intervals, since precision estimated from a few dozen calls is not comparable with precision estimated from thousands.
  Duplications are 17.3% of the merged benchmark over all sizes and 27.8% above 1 kb, so the curves describe both CNV types.
]

#pagebreak()

#figure(
  image("/results/parameterization/consensus_overlap.png", width: 100%)
)

#cap("Supplementary Figure Consensus Overlap:")[
  Effect of the query consensus reciprocal-overlap threshold on the three 30x consensus call sets, swept from 0.05 to 0.95.
  The benchmark is held at zero padding, the size floor at 1 kb on both sides, and the classification threshold at 0.5.
  (A) Consensus calls surviving the 1 kb floor.
  (B) Precision.
  (C) Recall, whose denominator is the same 23,193 benchmark intervals at every point.
  (D) F1.
  The dotted vertical line marks the adopted value of 0.5.
]

#pagebreak()

#figure(
  image("/results/parameterization/classification_overlap.png", width: 100%)
)

#cap("Supplementary Figure Classification Overlap:")[
  Effect of the classification reciprocal-overlap threshold, swept from 0 to 0.99.
  The benchmark is held at zero padding, the size floor at 1 kb on both sides, and the consensus threshold at 0.5.
  (A) Precision, recall, and F1 for the 30x 2-of-3 consensus.
  (B) Query calls credited against more than one benchmark interval, and (C) benchmark intervals credited to more than one query call, for all six 30x call sets; both axes are symmetric-log so that zero has a position.
  (D) F1 for the same six call sets.
  Panels B--D share the legend in D.
  The dotted vertical line marks the adopted value of 0.5.
]

#pagebreak()

#figure(
  image("/results/parameterization/sensitivity.png", width: 100%)
)

#cap("Supplementary Figure Sensitivity:")[
  Joint behavior of the four parameters over the full factorial grid, for the 30x 2-of-3 consensus call set.
  (A) Sobol' first-order indices (solid) and total-order indices (pale extension) for precision, recall and F1, with the variance carried only by interactions at the right.
  (B) Marginal F1 of each parameter under the additive model, against position within that parameter's swept range; the dotted line is the grand mean.
  (C) Precision and recall at every setting in the grid (gray), the Pareto front (black), the adopted setting (star), and the setting attaining the highest F1 (diamond).
  (D) The same plane with the classification threshold held at 0.5, colored by size floor.
]

#pagebreak()

#figure(
  image("/results/breakpoint_precedence/breakpoint_precedence.png", width: 100%)
)

#cap("Supplementary Figure Breakpoint Precedence:")[
  Breakpoint precedence among the three callers in the 30x 2-of-3 consensus call set (Supplementary Note S3).
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
  Breakpoint precedence of each caller in the 2-of-3 consensus call set at each coverage, at the adopted 50% reciprocal overlap (Supplementary Note S3).
  Components is the number of consensus components carrying at least two callers in which the caller appears.
  Holds is the percentage of those components in which the caller's call reaches the union start or end, counting boundaries tied with another caller; Holds alone counts only boundaries no other caller reaches.
  Span Ratio is the median ratio of the caller's own extent within a component to the component's union span.
]

#pagebreak()

#block(width: 100%)[
#set text(hyphenate: false, size: 9.5pt)
#show table.cell.where(y: 1): strong
#table(
  columns: (1fr, auto, auto, auto, auto, auto, auto, auto),
  align: (left, right, right, right, right, right, right, right),
  stroke: none,
  table.hline(stroke: 0.6pt),
  table.header(
    table.cell(rowspan: 2, align: left + bottom)[Call Set],
    table.cell(colspan: 3, align: center)[Deletion Precision],
    table.cell(colspan: 2, align: center)[Deletion Recall],
    table.cell(colspan: 2, align: center)[1000G-only Recoveries (%)],
    table.hline(start: 1, end: 4, stroke: 0.3pt),
    table.hline(start: 4, end: 6, stroke: 0.3pt),
    table.hline(start: 6, end: 8, stroke: 0.3pt),
    [Full], [Long-read], [Change], [Full], [Long-read], [DEL], [DUP],
  ),
  table.hline(stroke: 0.6pt),
  [30x -- 2/3 Consensus], [0.921], [0.891], [−0.030], [0.240], [0.271], [1.8], [81.3],
  [6x -- 2/3 Consensus], [0.908], [0.868], [−0.040], [0.073], [0.082], [1.8], [85.1],
  [4x -- 2/3 Consensus], [0.934], [0.896], [−0.038], [0.042], [0.047], [1.6], [81.9],
  [2x -- 2/3 Consensus], [0.890], [0.838], [−0.052], [0.020], [0.022], [1.8], [80.2],
  table.hline(stroke: 0.3pt),
  [SNP Array], [0.599], [0.570], [−0.029], [0.030], [0.033], [3.0], [87.5],
  table.hline(stroke: 0.6pt),
)
]

#cap("Supplementary Table Benchmark Independence:")[
  Dependence of the comparison on the 1000 Genomes high-coverage SV call set, which was called from the same alignments as the sequencing call sets (Supplementary Note S4).
  Full is the merged benchmark of all three sources (16,742 deletion intervals above the 1 kb floor) and Long-read the benchmark merged from HGSVC3 and ONT Vienna alone (14,330); Change is the long-read precision minus the full precision.
  1000G-only Recoveries is the percentage of the benchmark intervals each call set recovered against the full benchmark that are supported by the 1000 Genomes set alone, for deletions and for duplications.
  All values at the adopted parameters.
]
