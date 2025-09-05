#' @title Caribbean Monk Seal sightings
#'
#' @description
#' Caribbean Monk Seal (*Monachus tropicalis*) sighting records from Solow
#' (1993) based on data presented in LeBouef et al. (1986).
#'
#' @format A numeric vector with 5 values.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction from Sighting Data. *Ecology*,
#' 74(3), 962-964. \doi{10.2307/1940821}
#'
#' **Other References**
#'
#' LeBoeuf, B. J., Kenyon, K. W., & Villa‐Ramirez, B. (1986). The Caribbean monk
#' seal is extinct. *Marine Mammal Science*, 2(1), 70-72.
#' \doi{10.1111/j.1748-7692.1986.tb00028.x}
"monk_seal"

#' @title Black-footed Ferret sightings (continuous)
#'
#' @description
#' Black-footed Ferret (*Mustela nigripes*) sighting records from Solow
#' (1993) based on data provided by the US Department of the Interior.
#'
#' @format A numeric vector with 28 values.
#'
#' @note
#' Dates represent months post January 1972 (e.g. June 1972 is 6).
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction in a Declining Population.
#' *Journal of Mathematical Biology*, 32(1), 79-82. \doi{10.1007/Bf00160376}
"ferret1"

#' @title Black-footed Ferret sightings (discrete)
#'
#' @description
#' Black-footed Ferret (*Mustela nigripes*) sighting records from Solow
#' (1993) based on data provided by the US Department of the Interior.
#'
#' @format A `data.frame` with 2 columns and 21 rows.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction in a Declining Population.
#' *Journal of Mathematical Biology*, 32(1), 79-82. \doi{10.1007/Bf00160376}
"ferret2"

#' @title *Metrarabdotos* n. sp. 5 sightings
#'
#' @description
#' Stratigraphic intervals for *Metrarabdotos* n. sp. 5 from Marshall (1994)
#' based on data supplied by A. H. Cheetham.
#'
#' @format A numeric vector with 11 values.
#'
#' @note
#' The oldest sightings are set at time 0.
#'
#' @references
#' **Key Reference**
#'
#' Marshall, C. R. (1994). Confidence intervals on stratigraphic ranges:
#' partial relaxation of the assumption of randomly distributed fossil
#' horizons. *Paleobiology*, 20(4), 459-469. \doi{10.1017/S0094837300012938}
"metrarabdotos"

#' @title Example sightings from Burgman et al.
#'
#' @description
#' Example frequency-based sighting records from Figure 1b in Burgman et al.
#' (1995).
#'
#' @format a `data.frame` with 2 columns and 16 rows.
#'
#' @references
#' **Key Reference**
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring Threat from
#' Scientific Collections. Conservation Biology, 9(4), 923-928.
#' \doi{0.1046/j.1523-1739.1995.09040923.x}
"burgmanf1b"

#' @title Lord Howe Gerygone sightings
#'
#' @description
#' Frequency-based collection-effort sighting record of the Lord Howe Gerygone
#' (*Gerygone insularis*) based on data from the Atlas of Living Australia.
#'
#' @format a `data.frame` with 3 columns and 236 rows.
#'
#' @note
#' Includes all observations of birds within 10km of Lord Howe Island with
#' valid spatial coordinates. The index of collection effort is the total
#' number of bird species recorded in that year.
#'
#' @references
#' **Key Reference**
#' Atlas of Living Australia occurrence download, accessed on 25 August 2025.
#' \doi{10.26197/ala.dad3870b-1b79-4600-ad12-af429feb1ad2}
#'
#' **Other References**
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#' museum data. *Biological Conservation*, 83(1), 9-17.
#' \doi{10.1016/S0006-3207(97)00048-7}
"gerygone"

#' @title Dodo sightings
#'
#' @description
#' Dodo (*Raphus cucullatus*) sighting records from Roberts & Solow (2003).
#'
#' @format A numeric vector with 10 values.
#'
#' @references
#' **Key Reference**
#'
#' Roberts, D. L., & Solow, A. R. (2003). Flightless birds: when did the dodo
#' become extinct? *Nature*, 426(6964), 245. \doi{10.1038/426245a}
"dodos"

#' @title Alaskan Mammoth fossil record
#'
#' @description
#' Woolly Mammoth (*Mammuthus primigenius*) dated fossil records from Solow et
#' al. (2006) based on data presented in Guthrie (2004).
#'
#' @format A numeric vector with 25 values.
#'
#' @note
#' Includes the 25 most recent fossil records, measured in years before present
#' (note that values are all <0). Radiometric dating error estimates are not
#' included.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., Roberts, D. L., & Robbirt, K. M. (2006). On the Pleistocene
#' extinctions of Alaskan mammoths and horses. *Proceedings of the National*
#' *Academy of Sciences of the United States of America*, 103(19), 7351-7353.
#' \doi{10.1073/pnas.0509480103}
#'
#' **Other References**
#'
#' Guthrie, R. D. (2004). Radiocarbon evidence of mid-Holocene mammoths
#' stranded on an Alaskan Bering Sea island. *Nature*, 429(6993), 746-749.
#' \doi{10.1038/nature02612}
#'
#' Bradshaw, C. J. A., Cooper, A., Turney, C. S. M., & Brook, B. W. (2012).
#' Robust estimates of extinction time in the geological record. *Quaternary*
#' *Science Reviews*, 33, 14-19. \doi{10.1016/j.quascirev.2011.11.021}
"mammoth"

#' @title *Anabarella* fossil record
#'
#' @description
#' Siberian and Mongolian records of Cambrian *Anabarella* sp. from Wang et al.
#' (2015) based on data presented in Maloof et al. (2010).
#'
#' @format A numeric vector with 19 values.
#'
#' @note
#' Dates are measured in million years before present (note that values are all
#' <0). Radiometric dating error estimates are not included.
#'
#' @references
#' **Key Reference**
#'
#' Wang, S. C., Everson, P. J., Zhou, H. J., Park, D., & Chudzicki, D. J.
#' (2016). Adaptive credible intervals on stratigraphic ranges when recovery
#' potential is unknown. *Paleobiology*, 42(2), 240-256.
#' \doi{10.1017/pab.2015.37}
#'
#' **Other Reference**
#'
#' Maloof, A. C., Porter, S. M., Moore, J. L., Dudás, F. Ö., Bowring, S. A.,
#' Higgins, J. A., ... & Eddy, M. P. (2010). The earliest Cambrian record of
#' animals and ocean geochemical change. *GSA Bulletin*, 122(11-12), 1731-1774.
#' \doi{10.1130/B30346.1}
"anabarella"

#' @title Eskimo Curlew sightings
#'
#' @description
#' Eskimo Curlew (*Numenius borealis*) sighting records from Jarić & Roberts
#' (2014) based on data from Roberts et al. (2010).
#'
#' @format a `data.frame` with 2 columns and 49 rows.
#'
#' @note
#' Sighting certainty values are the mean values of the reliability intervals
#' (i.e. 0.85, 0.7, 0.25).
#'
#' @references
#' **Key Reference**
#'
#' Jarić, I., & Roberts, D. L. (2014). Accounting for observation reliability
#' when inferring extinction based on sighting records. *Biodiversity and*
#' *Conservation*, 23(11), 2801-2815. \doi{10.1007/s10531-014-0749-8}
#'
#' **Other References**
#'
#' Roberts, D. L., Elphick, C. S., & Reed, J. M. (2010). Identifying anomalous
#' reports of putatively extinct species and why it matters. *Conservation*
#' *Biology*, 24(1), 189-196.\doi{10.1111/j.1523-1739.2009.01292.x}
"eskimo_curlew"

#' @title Ivory-billed Woodpecker sightings (certain)
#'
#' @description
#' Ivory-billed Woodpecker (*Campephilus principalis*) sighting records from
#' Roberts et al. (2010).
#'
#' @format A numeric vector with 22 values.
#'
#' @note
#' Only certain sightings are included.
#'
#' @references
#' **Key Reference**
#'
#' Roberts, D. L., Elphick, C. S., & Reed, J. M. (2010). Identifying anomalous
#' reports of putatively extinct species and why it matters. *Conservation*
#' *Biology*, 24(1), 189-196.\doi{10.1111/j.1523-1739.2009.01292.x}
"woodpecker1"

#' @title Ivory-billed Woodpecker sightings (certain and uncertain)
#'
#' @description
#' Ivory-billed Woodpecker (*Campephilus principalis*) sighting records from
#' Roberts et al. (2010).
#'
#' @format A `data.frame` with 2 columns and 68 rows.
#'
#' @note
#' Sightings are indicated as either certain (`certain = TRUE`) or uncertain
#' (`certain = FALSE`).
#'
#' @references
#' **Key Reference**
#'
#' Roberts, D. L., Elphick, C. S., & Reed, J. M. (2010). Identifying anomalous
#' reports of putatively extinct species and why it matters. *Conservation*
#' *Biology*, 24(1), 189-196.\doi{10.1111/j.1523-1739.2009.01292.x}
"woodpecker2"

#' @title Tasmanian Fox carcass recoveries
#'
#' @description
#' European Red Fox (*Vulpes vulpes*) carcass recovery dates from Caley & Barry
#' (2014).
#'
#' @format A numeric vector with 4 values.
#'
#' @references
#' **Key Reference**
#'
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities from
#' sighting records: inference and uncertainties. *PLoS One*, 9(4), e95857.
#' \doi{10.1371/journal.pone.0095857}
"fox"

#' @title Mauritian orchid sightings
#'
#' @description
#' Mauritian orchid (*Angraecum aff. rutenbergianum*) sighting records from
#' Solow & Roberts (2003).
#'
#' @format A numeric vector with 11 values.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for extinction
#' based on a sighting record. *Ecology*, 84(5), 1329-1332.
#' \doi{10.1890/0012-9658(2003)084[1329:ANTFEB]2.0.CO;2}
"orchid"

