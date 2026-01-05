#' @title Caribbean Monk Seal sightings
#'
#' @description
#' Caribbean Monk Seal (*Monachus tropicalis*) sighting record from Solow
#' (1993) based on data presented in LeBouef et al. (1986). The temporal range
#' (`init.time` to `test.time`) of these data is 1915 to 1992. Format: `ccon`.
#'
#' @format A `numeric` vector with 5 values.
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

#' @title Black-footed Ferret sightings
#'
#' @description
#' Black-footed Ferret (*Mustela nigripes*) sighting record from Solow (1993)
#' based on data provided by the US Department of the Interior, processed into
#' various `dodo` formats. The temporal range (`init.time` to `test.time`) of
#' these data is either 0 to 229 (`ccon`) or 1972 to 1992 (`cdis`). Formats:
#' `ccon`, `cdis`.
#'
#' @format A `list` with 2 elements.
#'
#' @note
#' Times for `ccon` represent months post January 1972 (e.g. June 1972 is 6).
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R. (1993). Inferring Extinction in a Declining Population.
#' *Journal of Mathematical Biology*, 32(1), 79-82. \doi{10.1007/Bf00160376}
"ferret"

#' @title *Metrarabdotos* n. sp. 5 sightings
#'
#' @description
#' Stratigraphic intervals for *Metrarabdotos* n. sp. 5 from Marshall (1994)
#' based on data supplied by A. H. Cheetham. Format: `ccon`.
#'
#' @format A `numeric` vector with 11 values.
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
#' Example frequency-based sighting record from Figure 1b in Burgman et al.
#' (1995). The temporal range (`init.time` to `test.time`) of these data is 0 to
#' 16. Format: `cdis`.
#'
#' @format a `numeric` vector with 16 values.
#'
#' @references
#' **Key Reference**
#'
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring Threat from
#' Scientific Collections. Conservation Biology, 9(4), 923-928.
#' \doi{0.1046/j.1523-1739.1995.09040923.x}
"burgman_figure1b"

#' @title Lord Howe Gerygone sightings
#'
#' @description
#' Lord Howe Gerygone (*Gerygone insularis*) sighting record based on data from
#' the Atlas of Living Australia. The temporal range (`init.time` to
#' `test.time`) of these data is 1788 to 2023. Format: `cdis`.
#'
#' @format a `numeric` vector with 236 values.
#'
#' @references
#' **Key Reference**
#'
#' Atlas of Living Australia occurrence download, accessed on 25 August 2025.
#' \doi{10.26197/ala.dad3870b-1b79-4600-ad12-af429feb1ad2}
#'
#' **Other References**
#'
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#' museum data. *Biological Conservation*, 83(1), 9-17.
#' \doi{10.1016/S0006-3207(97)00048-7}
#'
#' @seealso [gerygone_effort]
"gerygone"

#' @title Lord Howe Gerygone effort
#'
#' @description
#' Lord Howe Gerygone (*Gerygone insularis*) sampling effort proxy data.
#' Following McCarthy et al. (1998), effort is approximated by the number of
#' species recorded in the species' range (within 10km of Lord Howe Island) in
#' each year. Species records come from the Atlas of Living Australia (all
#' records of birds with coordinates). The temporal range (`init.time` to
#' `test.time`) of these data is 1788 to 2023.
#'
#' @format a `numeric` vector with 236 values.
#'
#' @references
#' **Key Reference**
#'
#' Atlas of Living Australia occurrence download, accessed on 25 August 2025.
#' \doi{10.26197/ala.dad3870b-1b79-4600-ad12-af429feb1ad2}
#'
#' **Other References**
#'
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#' museum data. *Biological Conservation*, 83(1), 9-17.
#' \doi{10.1016/S0006-3207(97)00048-7}
#'
#' @seealso [gerygone]
"gerygone_effort"

#' @title Dodo sightings
#'
#' @description
#' Dodo (*Raphus cucullatus*) sighting record from Roberts & Solow (2003). The
#' temporal range (`init.time` to `test.time`) of these data is 1598 to 2002.
#' Format: `ccon.`
#'
#' @format A `numeric` vector with 10 values.
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
#' al. (2006) based on data presented in Guthrie (2004). Format: `ccon`.
#'
#' @format A `numeric` vector with 25 values.
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
#' (2016) based on data presented in Maloof et al. (2010). Format: `ccon`.
#'
#' @format A `numeric` vector with 19 values.
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

#' @title Ivory-billed Woodpecker sightings
#'
#' @description
#' Ivory-billed Woodpecker (*Campephilus principalis*) sighting record from
#' Roberts et al. (2010), processed into various `dodo` formats. The temporal
#' range (`init.time` to `test.time`) of these data is 1897 to 2010. Formats:
#' `ccon`, `cbin`, `cdis`, `ucon`, `ubin`, `umcd`.
#'
#' @format a `list` with 6 elements.
#'
#' @references
#' **Key Reference**
#'
#' Roberts, D. L., Elphick, C. S., & Reed, J. M. (2010). Identifying anomalous
#' reports of putatively extinct species and why it matters.
#' *Conservation Biology*, 24(1), 189-196.
#' \doi{10.1111/j.1523-1739.2009.01292.x}
"woodpecker"

#' @title Tasmanian Fox carcass recoveries
#'
#' @description
#' European Red Fox (*Vulpes vulpes*) carcass recovery data from Caley & Barry
#' (2014). The temporal range (`init.time` to `test.time`) of these data is 2001
#' to 2012. Format: `cbin`.
#'
#' @format A `numeric` vector with 12 values.
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
#' Mauritian orchid (*Angraecum aff. rutenbergianum*) sighting record from
#' Solow & Roberts (2003). The temporal range (`init.time` to `test.time`) of
#' these data is 1933 to 1996. Format: `ccon`.
#'
#' @format A `numeric` vector with 11 values.
#'
#' @references
#' **Key Reference**
#'
#' Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for extinction
#' based on a sighting record. *Ecology*, 84(5), 1329-1332.
#' \doi{10.1890/0012-9658(2003)084[1329:ANTFEB]2.0.CO;2}
"orchid"

#' @title Purple-winged Ground Dove sightings
#'
#' @description
#' Purple-winged Ground Dove (*Paraclaravis geoffroyi*) sighting record from
#' Lees et al. (2021), processed into various `dodo` formats. The temporal
#' range (`init.time` to `test.time`) of these data is 1811 to 2019. Formats:
#' `ccon`, `cbin`, `ubin`.
#'
#' @format A `list` with 3 elements.
#'
#' @note
#' Adapted from Supplementary Table 1 in Lees et al. (2021). Almost all specimen
#' records are considered certain, while sight and literature-based records are
#' all considered uncertain. The controversial specimen from 1991 is classified
#' as uncertain. Records have been deduplicated to ensure they are all
#' independent from one another. Records without a date are not included.
#'
#' @references
#' **Key Reference**
#'
#' Lees, A. C., Devenish, C., Areta, J. I., de Araújo, C. B., Keller, C.,
#' Phalan, B., & Silveira, L. F. (2021). Assessing the Extinction Probability
#' of the Purple-winged Ground Dove, an Enigmatic Bamboo Specialist.
#' *Frontiers in Ecology and Evolution*, 9. \doi{10.3389/fevo.2021.624959}
"ground_dove"

#' @title Example sightings from Thompson et al.
#'
#' @description
#' Example multi-class sighting record from Table 1 in Thompson et al. (2013).
#' The temporal range (`init.time` to `test.time`) of these data is 0 to 7.
#' Format: `umcd`.
#'
#' @format a `data.frame` with 4 columns and 7 rows.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Lee, T. E., Stone, L., McCarthy, M. A., & Burgman, M. A.
#' (2013). Inferring extinction risks from sighting records.
#' *Journal of Theoretical Biology*, 338, 16-22.
#' \doi{10.1016/j.jtbi.2013.08.023}
"thompson_table1"

#' @title Slender-billed Curlew sightings (raw)
#'
#' @description
#' Unprocessed Slender-billed Curlew (*Numenius tenuirostris*) sighting record
#' from Buchanan et al. (2025).
#'
#' @format a `data.frame` with 6 columns and 1174 rows.
#'
#' @references
#' **Key Reference**
#'
#' Buchanan, G. M., Chapple, B., Berryman, A. J., Crockford, N., Jansen, J. J.
#' F. J., & Bond, A. L. (2025). Global extinction of Slender-billed Curlew
#' (*Numenius tenuirostris*). *Ibis*, 167(2), 357-370. \doi{10.1111/ibi.13368}
#'
#' @seealso [curlew]; [convert_dodo()]
"curlew_raw"

#' @title Slender-billed Curlew sightings (processed)
#'
#' @description
#' Slender-billed Curlew (*Numenius tenuirostris*) sighting record from
#' Buchanan et al. (2025), processed into various `dodo` formats. The temporal
#' range (`init.time` to `test.time`) of these data is 1817 to 2022. Formats:
#' `ccon`, `cbin`, `cdis`, `ucon`, `ubin`, `umcd`, `iucn`. An additional
#' "format", `buchanan`, contains the data in the exact `iucn` format as used
#' in Buchanan et al. (2025).
#'
#' @format a `list` with 8 elements.
#'
#' @references
#' **Key Reference**
#'
#' Buchanan, G. M., Chapple, B., Berryman, A. J., Crockford, N., Jansen, J. J.
#' F. J., & Bond, A. L. (2025). Global extinction of Slender-billed Curlew
#' (*Numenius tenuirostris*). *Ibis*, 167(2), 357-370. \doi{10.1111/ibi.13368}
#'
#' @seealso [curlew_raw]; [convert_dodo()]; [curlew_effort]
"curlew"

#' @title Slender-billed Curlew surveys
#'
#' @description
#' Slender-billed Curlew (*Numenius tenuirostris*) survey data from Buchanan et
#' al. (2025). The temporal range (`init.time` to `test.time`) of these data
#' is 1892 to 2022.
#'
#' @format a `data.frame` with 11 columns and 206 rows.
#'
#' @references
#' **Key Reference**
#'
#' Buchanan, G. M., Chapple, B., Berryman, A. J., Crockford, N., Jansen, J. J.
#' F. J., & Bond, A. L. (2025). Global extinction of Slender-billed Curlew
#' (*Numenius tenuirostris*). *Ibis*, 167(2), 357-370. \doi{10.1111/ibi.13368}
"curlew_surveys"

#' @title Slender-billed Curlew effort
#'
#' @description
#' Sampling effort proxy data for the Slender-billed Curlew
#' (*Numenius tenuirostris*). Following McCarthy et al. (1998), effort is
#' approximated by the number of species recorded in the species' range in each
#' year. Species records come from GBIF (all records of Aves with coordinates),
#' and range is defined as all of the countries  with 3+ confirmed (p_ci >= 0.9)
#' records in Buchanan et al. (2025). The temporal range (`init.time` to
#' `test.time`) of these data is 1817 to 2022.
#'
#' @format a `numeric` vector with 206 values.
#'
#' @references
#' **Key Reference**
#'
#' Buchanan, G. M., Chapple, B., Berryman, A. J., Crockford, N., Jansen, J. J.
#' F. J., & Bond, A. L. (2025). Global extinction of Slender-billed Curlew
#' (*Numenius tenuirostris*). *Ibis*, 167(2), 357-370. \doi{10.1111/ibi.13368}
#'
#' **Other References**
#'
#' McCarthy, M. A. (1998). Identifying declining and threatened species with
#' museum data. Biological Conservation, 83(1), 9-17.
#' \doi{10.1016/S0006-3207(97)00048-7}
#'
#' @seealso [curlew]
"curlew_effort"

#' @title Alaotra Grebe sightings
#'
#' @description
#' Alaotra Grebe (*Tachybaptus rufolavatus*) sighting record from Thompson et
#' al. (2017). The temporal range (`init.time` to `test.time`) of these data
#' is 1929 to 2017. Format: `iucn`.
#'
#' @format a `data.frame` with 5 columns and 89 rows.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Koshkina, V., Burgman, M. A., Butchart, S. H. M., & Stone,
#' L. (2017). Inferring extinctions II: A practical, iterative model based on
#' records and surveys. *Biological Conservation*, 214, 328-335.
#' \doi{10.1016/j.biocon.2017.07.029}
"grebe"

#' @title Alaotra Grebe surveys
#'
#' @description
#' Alaotra Grebe (*Tachybaptus rufolavatus*) survey data from Thompson et
#' al. (2017). The temporal range (`init.time` to `test.time`) of these data
#' is 1929 to 2017.
#'
#' @format a `data.frame` with 11 columns and 89 rows.
#'
#' @references
#' **Key Reference**
#'
#' Thompson, C. J., Koshkina, V., Burgman, M. A., Butchart, S. H. M., & Stone,
#' L. (2017). Inferring extinctions II: A practical, iterative model based on
#' records and surveys. *Biological Conservation*, 214, 328-335.
#' \doi{10.1016/j.biocon.2017.07.029}
"grebe_surveys"

#' @title Pohnpei Mountain Starling sightings
#'
#' @description
#' Pohnpei Mountain Starling (*Aplonis pelzelni*) sighting record from Lee
#' (2014). The temporal range (`init.time` to `test.time`) of these data
#' is 1930 to 2013. Format: `iucn`.
#'
#' @format a `data.frame` with 5 columns and 84 rows.
#'
#' @references
#' **Key Reference**
#'
#' Lee, T. E. (2014). A simple numerical tool to infer whether a species is
#' extinct. *Methods in Ecology and Evolution*, 5(8), 791-796.
#' \doi{10.1111/2041-210x.12227}
"starling"

#' @title Pohnpei Mountain Starling surveys
#'
#' @description
#' Pohnpei Mountain Starling (*Aplonis pelzelni*) survey data from Lee (2014).
#' The temporal range (`init.time` to `test.time`) of these data is 1930 to
#' 2013.
#'
#' @format a `data.frame` with 11 columns and 84 rows.
#'
#' @references
#' **Key Reference**
#'
#' Lee, T. E. (2014). A simple numerical tool to infer whether a species is
#' extinct. *Methods in Ecology and Evolution*, 5(8), 791-796.
#' \doi{10.1111/2041-210x.12227}
"starling_surveys"

#' @title Weiss & Marshall surveys
#'
#' @description
#' Palaeontological survey data from Weiss & Marshall (2014).
#'
#' @format a `numeric` vector with 36 values.
#'
#' @references
#' **Key Reference**
#'
#' Weiss, R. E., & Marshall, C. R. (1999). The Uncertainty in the True End
#' Point of a Fossil's Stratigraphic Range When Stratigraphic Sections Are
#' Sampled Discretely. *Mathematical Geology*, 31(4), 435-453.
#' \doi{10.1023/A:1007542725180}
"weissmarshall_surveys"

#' @title Verneuilinoides sp. A records
#'
#' @description
#' Verneuilinoides sp. A (species 25) sighting record from Weiss & Marshall
#' (2014). Format: `cbin`.
#'
#' @format a `numeric` vector with 36 values.
#'
#' @references
#' **Key Reference**
#'
#' Weiss, R. E., & Marshall, C. R. (1999). The Uncertainty in the True End
#' Point of a Fossil's Stratigraphic Range When Stratigraphic Sections Are
#' Sampled Discretely. *Mathematical Geology*, 31(4), 435-453.
#' \doi{10.1023/A:1007542725180}
"verneuilinoides"

#' @title Eggerellina brevis records
#'
#' @description
#' Eggerellina brevis (species 43) sighting record from Weiss & Marshall
#' (2014). Format: `cbin`.
#'
#' @format a `numeric` vector with 36 values.
#'
#' @references
#' **Key Reference**
#'
#' Weiss, R. E., & Marshall, C. R. (1999). The Uncertainty in the True End
#' Point of a Fossil's Stratigraphic Range When Stratigraphic Sections Are
#' Sampled Discretely. *Mathematical Geology*, 31(4), 435-453.
#' \doi{10.1023/A:1007542725180}
"eggerellina_brevis"

#' @title Example sightings from Lee et al.
#'
#' @description
#' Example binary sighting record from Appendix S1 in Lee et al. (2014). The
#' temporal range (`init.time` to `test.time`) of these data is 1 to 9.
#' Format: `ubin`.
#'
#' @format a `data.frame` with 2 columns and 9 rows.
#'
#' @references
#' **Key Reference**
#'
#' Lee, T. E., McCarthy, M. A., Wintle, B. A., Bode, M., Roberts, D. L., &
#' Burgman, M. A. (2014). Inferring extinctions from sighting records of
#' variable reliability. *Journal of Applied Ecology*, 51(1), 251-258.
#' \doi{10.1111/1365-2664.12144}
"lee_s1"

#' @title Bitterweed sightings
#'
#' @description
#' Bitterweed (*Helenium amarum*) sighting record from Rout et al. (2009).
#' Format: `cbin`.
#'
#' @format a `numeric` vector with 178 values.
#'
#' @references
#' **Key Reference**
#'
#' Rout, T. M., Salomon, Y., & McCarthy, M. A. (2009). Using sighting records
#' to declare eradication of an invasive species. *Journal of Applied Ecology*,
#' 46(1), 110-117.\doi{10.1111/j.1365-2664.2008.01586.x}
"bitterweed"
