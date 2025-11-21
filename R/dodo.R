#' dodo: Extinction Models in R
#'
#' dodo implements a selection of published extinction models.
#'
#' @section Model Functions:
#' * [SL88F1()] Springer & Lilje's (1988) "Broken Stick" model.
#' * [SS89F1()] Strauss & Sadler's (1989) "Classical" model.
#' * [SS89B1()] Strauss & Sadler's (1989) "Bayesian" model.
#' * [SO93F1()] Solow's (1993) "Classical" model.
#' * [SO93B1()] Solow's (1993) "Bayesian" model.
#' * [SO93F2()] Solow's (1993) "Declining" model.
#' * [MA94F1()] Marshall's (1994) "Distribution-free" model.
#' * [BU95F1()] Burgman et al.'s (1995) "Discrete-time" model.
#' * [BU95F2()] Burgman et al.'s (1995) "Runs Test" model.
#' * [BU95F3()] Burgman et al.'s (1995) "Empty Cells" model.
#' * [MC98F1()] McCarthy's (1998) "Non-random" model.
#' * [MC99F1()] McFarlane's (1999) "Partial Solow" model.
#' * [RS03F1()] Roberts & Solow's (2003) "Optimal Linear Estimation" model.
#' * [SR03F1()] Solow & Roberts' (2003) "Non-parametric" model.
#' * [MC06F1()] McInerny et al.'s (2006) "Sighting Rate" model.
#' * [RO09B1()] Rout et al.'s (2009) "Declining" model.
#' * [JE10F1()] Jarić & Ebenhard's (2010) "Stationary" model.
#' * [JE10F2()] Jarić & Ebenhard's (2010) "Non-stationary" model.
#' * [BR12F1()] Bradshaw et al.'s (2012) "Inverse-weighted McInerny" model.
#' * [SO12B1()] Solow et al.'s (2012) "Uncertain" model.
#' * [TH13B1()] Thompson et al.'s (2013) "Deterministic" model.
#' * [TH13B2()] Thompson et al.'s (2013) "Quenched" model.
#' * [TH13B3()] Thompson et al.'s (2013) "Annealed" model.
#' * [AL14B1()] Alroy's (2015) "Creeping Shadow of a Doubt" model.
#' * [CB14B1()] Caley & Barry's (2014) "Constant" model.
#' * [CB14B2()] Caley & Barry's (2014) "Non-constant" model.
#' * [JR14F1()] Jarić & Roberts' (2014) "Solow" model.
#' * [LE14B1()] Lee's (2014) "Simple" model.
#' * [SB14B1()] Solow & Beet's (2014) "Model 1" model.
#' * [SB14B2()] Solow & Beet's (2014) "Model 2" model.
#' * [AL15B1()] Alroy's (2015) "Agnostic" model.
#' * [SO16B1()] Solow's (2016) "Retrospective" model.
#' * [SO16B2()] Solow's (2016) "Sequential" model.
#' * [WA16B1()]	Wang et al.'s (2016) "Adaptive" model.
#' * [TH17I1()] Thompson et al.'s (2017) "Iterative" model.
#' * [BR19F1()] Brook et al.'s (2019) "Solow" model.
#' * [BR19F2()] Brook et al.'s (2019) "Solow & Roberts" model.
#' * [BR19F3()] Brook et al.'s (2019) "Roberts & Solow" model.
#' * [BR19F4()] Brook et al.'s (2019) "McInerny" model.
#' * [BR19F5()] Brook et al.'s (2019) "Last Appearance Date" model.
#' * [TH19B1()] Thompson et al.'s (2017) "Bayesian Updating" model.
#' * [KO20B1()] Kodikara et al.'s (2020) "Certain-only" model.
#' * [KO20B2()] Kodikara et al.'s (2020) "Uncertain" model.
#' * [KO21B1()] Kodikara et al.'s (2021) "Change-point" model.
#'
#' @section Other Functions:
#' * [convert_dodo()] Convert basic sightings spreadsheet into all `dodo`
#' formats.
#'
#' @name dodo
"_PACKAGE"
