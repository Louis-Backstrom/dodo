curlew_effort <- c(
  10, 17, 18, 13, 16, 20, 18, 20, 15, 22, 28, 26, 33, 25, 13,
  20, 67, 46, 39, 57, 34, 43, 63, 83, 64, 103, 50, 77, 142,
  106, 134, 94, 96, 136, 112, 124, 107, 143, 155, 148, 135,
  160, 343, 268, 283, 290, 160, 182, 135, 189, 284, 195, 227,
  145, 152, 216, 259, 230, 298, 311, 384, 267, 332, 315, 310,
  356, 347, 358, 343, 369, 348, 323, 302, 333, 261, 276, 304,
  341, 387, 369, 352, 306, 359, 358, 362, 398, 409, 378, 456,
  393, 394, 404, 429, 432, 460, 486, 454, 469, 367, 321, 304,
  282, 288, 333, 373, 334, 360, 391, 401, 471, 466, 431, 409,
  434, 382, 396, 346, 368, 363, 366, 372, 394, 310, 269, 256,
  225, 219, 218, 211, 257, 319, 335, 331, 373, 376, 356, 333,
  423, 378, 376, 419, 437, 435, 390, 414, 368, 434, 397, 404,
  426, 411, 442, 421, 438, 465, 442, 463, 478, 521, 502, 480,
  542, 544, 631, 561, 504, 542, 591, 623, 557, 596, 648, 623,
  596, 676, 682, 727, 700, 736, 670, 666, 677, 749, 732, 770,
  758, 780, 742, 800, 813, 836, 832, 867, 849, 890, 907, 924,
  915, 953, 973, 962, 976, 1024, 959, 1028, 1028
)
usethis::use_data(curlew_effort, overwrite = TRUE)

gerygone_effort <- c(
  7, 0, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 7, 0, 4, 0, 1,
  40, 5, 4, 2, 3, 7, 3, 1, 0, 0, 4, 5, 0, 1, 0, 1, 21, 2, 0,
  4, 17, 9, 8, 14, 6, 1, 20, 15, 9, 1, 2, 1, 1, 3, 14, 16, 4,
  7, 3, 3, 2, 7, 3, 1, 1, 0, 2, 3, 4, 19, 0, 5, 1, 2, 4, 2,
  4, 6, 3, 1, 0, 0, 12, 2, 6, 5, 4, 1, 4, 6, 8, 0, 25, 13, 7,
  17, 26, 10, 5, 2, 23, 14, 14, 9, 38, 20, 11, 21, 37, 4, 36,
  30, 11, 43, 31, 7, 9, 10, 49, 27, 10, 20, 37, 22, 24, 51,
  45, 11, 43, 42, 25, 49, 62, 36, 63, 58, 63, 49, 52, 44, 50,
  40, 48, 48, 47, 63, 62, 70, 56, 59, 65, 66, 70, 44, 51, 70,
  74
)
usethis::use_data(gerygone_effort, overwrite = TRUE)

laughing_owl_effort <- data.frame(year = 1844:2018)

for (taxon_key in c(1, 212, 185019954)) { # Animals, Birds, Morepork
  effort_raw <- data.frame(
    year = integer(),
    species = character()
  )

  for (year in 1844:2018) {
    species_list <- rgbif::occ_search(
      taxonKey = taxon_key,
      country = "NZ",
      year = year,
      facet = "speciesKey",
      facetLimit = 10000,
      limit = 0 # Just get facets, no occurrences
    )$facets$speciesKey$name

    if (!is.null(species_list)) {
      species_list <- data.frame(
        year = year,
        species = species_list
      )
    }

    effort_raw <- rbind(effort_raw, species_list)

    Sys.sleep(0.1)
  }

  effort_full <- tapply(
    effort_raw$species,
    effort_raw$year,
    function(x) length(unique(na.omit(x)))
  )

  effort_full <- effort_full[as.character(1844:2018)]
  effort_full[is.na(effort_full)] <- 0
  effort_full <- as.integer(effort_full)

  laughing_owl_effort <- cbind(laughing_owl_effort, effort_full)
}

colnames(laughing_owl_effort) <- c("year", "animals", "birds", "moreporks")
rownames(laughing_owl_effort) <- laughing_owl_effort$year
laughing_owl_effort$year <- NULL

usethis::use_data(laughing_owl_effort, overwrite = TRUE)
