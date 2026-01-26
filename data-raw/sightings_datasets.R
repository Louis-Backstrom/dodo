monk_seal <- c(1915, 1922, 1932, 1948, 1952)
usethis::use_data(monk_seal, overwrite = TRUE)

ferret <- list(
  ccon = c(
    6, 7, 8, 10, 17, 18, 19, 20, 21, 22, 30, 31, 41, 44, 46, 53, 57, 58,
    66, 90, 117, 118, 122, 123, 127, 139, 151, 153
  ),
  cdis = c(4, 6, 2, 3, 3, 1, 0, 1, 0, 2, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0)
)
usethis::use_data(ferret, overwrite = TRUE)

metrarabdotos <- c(0, 0, 0, 0, 0, 0, 0.05, 0.10, 0.20, 0.30, 1.40)
usethis::use_data(metrarabdotos, overwrite = TRUE)

burgman_figure1b <- c(1, 2, 0, 2, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
usethis::use_data(burgman_figure1b, overwrite = TRUE)

dodos <- c(1598, 1601, 1602, 1607, 1611, 1628, 1628, 1631, 1638, 1662)
usethis::use_data(dodos, overwrite = TRUE)

mammoth <- -c(
  11.50, 11.54, 11.91, 12.12, 12.19, 12.34, 12.43, 12.44, 12.48,
  12.49, 12.51, 12.58, 12.68, 12.88, 13.06, 13.23, 13.29, 13.34,
  13.38, 13.41, 13.44, 13.66, 13.69, 14.02, 14.09
) * 1000
usethis::use_data(mammoth, overwrite = TRUE)

anabarella <- -c(
  533.06579, 531.07032, 530.05212, 530.02947, 529.78322,
  529.47184, 528.11654, 527.84072, 527.72420, 527.68698,
  527.62881, 525.62909, 525.00290, 524.67879, 523.80700,
  523.76623, 523.67817, 522.95226, 522.19971
)
usethis::use_data(anabarella, overwrite = TRUE)

woodpecker <- list(
  ccon = c(
    1897, 1898, 1899, 1900, 1901, 1902, 1904, 1905, 1906, 1907, 1908,
    1909, 1910, 1913, 1914, 1917, 1924, 1925, 1932, 1935, 1938, 1939
  ),
  cbin = c(
    1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ),
  cdis = c(
    1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ),
  ucon = data.frame(
    time = c(
      1897, 1898, 1899, 1900, 1901, 1902, 1904, 1905, 1906, 1907, 1908,
      1909, 1910, 1913, 1914, 1917, 1924, 1925, 1932, 1935, 1938, 1939,
      1911, 1916, 1920, 1921, 1923, 1926, 1929, 1930, 1931, 1933, 1934,
      1936, 1937, 1941, 1942, 1943, 1944, 1946, 1948, 1949, 1950, 1951,
      1952, 1955, 1958, 1959, 1962, 1966, 1967, 1968, 1969, 1971, 1972,
      1973, 1974, 1976, 1981, 1982, 1985, 1986, 1987, 1988, 1999, 2004,
      2005, 2006
    ),
    certainty = c(rep(0.99, 22), rep(0.5, 17), rep(0.01, 29))
  ),
  ubin = data.frame(
    certain = c(
      1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ),
    uncertain = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1,
      1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0,
      0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1,
      0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0
    )
  ),
  umcd = data.frame(
    time = 1897:2010,
    class_1 = c(
      1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ),
    class_2 = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0,
      0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0
    ),
    class_3 = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1,
      0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0,
      1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 1, 1, 1, 0, 0, 0, 0
    )
  )
)
usethis::use_data(woodpecker, overwrite = TRUE)

fox <- c(1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0)
usethis::use_data(fox, overwrite = TRUE)

orchid <- c(1933, 1941, 1952, 1962, 1965, 1966, 1967, 1971, 1974, 1976, 1978)
usethis::use_data(orchid, overwrite = TRUE)

ground_dove <- list(
  ccon = c(
    1820, 1837, 1838, 1848, 1890, 1898, 1899, 1901, 1902, 1921, 1926,
    1929, 1937, 1943, 1946, 1950, 1951, 1953, 1956, 1957, 1959, 1985
  ),
  cbin = c(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0,
    0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0
  ),
  ubin = data.frame(
    certain = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    ),
    uncertain = c(
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1,
      1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1,
      1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0
    )
  )
)
usethis::use_data(ground_dove, overwrite = TRUE)

thompson_table1 <- data.frame(
  time = 1:7,
  class_1 = c(1, 0, 1, 0, 0, 0, 0),
  class_2 = c(1, 0, 0, 1, 0, 1, 0),
  class_3 = c(0, 1, 1, 0, 0, 0, 1)
)
usethis::use_data(thompson_table1, overwrite = TRUE)

curlew_raw <- readxl::read_xlsx("data-raw/ibi13368-sup-0001-tables1.xlsx")
curlew_raw <- as.data.frame(curlew_raw)
usethis::use_data(curlew_raw, overwrite = TRUE)

curlew <- convert_dodo(
  x = curlew_raw, init.time = 1817, test.time = 2022, threshold = 0.9,
  unique = TRUE, aggregate = "max", time = "year", certainty = "p_ci",
  certainty_lower = "p_ci_min", certainty_upper = "p_ci_max"
)

curlew$buchanan <- readxl::read_xlsx("data-raw/Numenius tenuirostris.xlsx",
  sheet = 2, range = "A9:D125"
)
names(curlew$buchanan) <- c(
  "time", "certainty_lower", "certainty",
  "certainty_upper"
)
curlew$buchanan <- merge(curlew$buchanan, data.frame("time" = 1892:2022),
  all = TRUE
)
curlew$buchanan$record <- !is.na(curlew$buchanan$certainty)
curlew$buchanan[is.na(curlew$buchanan)] <- 0
curlew$buchanan <- curlew$buchanan[, c(
  "time", "record", "certainty",
  "certainty_lower", "certainty_upper"
)]

usethis::use_data(curlew, overwrite = TRUE)

curlew_passive <- list(
  epsilon = c(0.7, 0.75, 0.8), p_i = c(0.6, 0.65, 0.7), p_r = c(0.7, 0.8, 0.9)
)

curlew_active <- data.frame(
  time = c(
    1989, 1990, 1992, 1994, 1996, 1997, 1998, 1999,
    2000, 2003, 2004, 2005, 2009, 2010, 2011
  ),
  epsilon = c(
    0.125, 0.125, 0.125, 0.125, 0.01, 0.01, 0.125, 0.01,
    0.0015, 0.5, 0.5, 0.5, 0.125, 0.9, 0.75
  ),
  epsilon_lower = c(
    0.1, 0.1, 0.1, 0.1, 0.001, 0.001, 0.1, 0.001,
    0.001, 0.33, 0.33, 0.33, 0.1, 0.85, 0.7
  ),
  epsilon_upper = c(
    0.15, 0.15, 0.15, 0.15, 0.05, 0.05, 0.25, 0.02,
    0.002, 0.66, 0.66, 0.66, 0.15, 0.95, 0.8
  ),
  p_i = c(
    0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65,
    0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65
  ),
  p_i_lower = c(
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6
  ),
  p_i_upper = c(
    0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
    0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7
  ),
  p_r = c(
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8
  ),
  p_r_lower = c(
    0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
    0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7
  ),
  p_r_upper = c(
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9,
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9
  )
)

curlew_surveys <- data.frame(time = 1817:2022)
curlew_surveys$survey <- curlew_surveys$time %in% curlew_active$time
curlew_surveys <- merge(curlew_surveys, curlew_active,
  by = "time",
  all.x = TRUE
)

curlew_surveys$epsilon_lower[is.na(curlew_surveys$epsilon_lower)] <-
  curlew_passive$epsilon[1]
curlew_surveys$epsilon[is.na(curlew_surveys$epsilon)] <-
  curlew_passive$epsilon[2]
curlew_surveys$epsilon_upper[is.na(curlew_surveys$epsilon_upper)] <-
  curlew_passive$epsilon[3]
curlew_surveys$p_i_lower[is.na(curlew_surveys$p_i_lower)] <-
  curlew_passive$p_i[1]
curlew_surveys$p_i[is.na(curlew_surveys$p_i)] <-
  curlew_passive$p_i[2]
curlew_surveys$p_i_upper[is.na(curlew_surveys$p_i_upper)] <-
  curlew_passive$p_i[3]
curlew_surveys$p_r_lower[is.na(curlew_surveys$p_r_lower)] <-
  curlew_passive$p_r[1]
curlew_surveys$p_r[is.na(curlew_surveys$p_r)] <-
  curlew_passive$p_r[2]
curlew_surveys$p_r_upper[is.na(curlew_surveys$p_r_upper)] <-
  curlew_passive$p_r[3]

rm(curlew_passive, curlew_active)
usethis::use_data(curlew_surveys, overwrite = TRUE)

curlew_effort <- c(
  19, 19, 15, 16, 14, 19, 13, 24, 31, 7, 13, 10, 10, 22, 3, 7,
  20, 20, 20, 13, 14, 19, 18, 36, 16, 40, 15, 18, 111, 28, 36,
  31, 21, 27, 56, 40, 42, 70, 57, 50, 48, 59, 207, 114, 171,
  159, 68, 70, 57, 109, 200, 117, 129, 80, 78, 102, 189, 112,
  209, 192, 278, 178, 207, 203, 186, 172, 151, 144, 136, 227,
  237, 188, 162, 200, 154, 142, 140, 159, 237, 192, 126, 142,
  201, 261, 188, 280, 221, 257, 379, 238, 304, 268, 319, 342,
  316, 361, 333, 356, 288, 274, 225, 225, 230, 245, 318, 275,
  298, 320, 333, 364, 385, 380, 350, 373, 342, 347, 293, 309,
  320, 325, 332, 356, 277, 247, 225, 183, 197, 195, 177, 239,
  290, 293, 279, 329, 326, 309, 295, 397, 351, 349, 382, 411,
  417, 398, 391, 360, 406, 403, 408, 417, 403, 458, 414, 443,
  470, 444, 455, 476, 517, 500, 494, 549, 562, 634, 553, 522,
  566, 614, 641, 584, 630, 682, 645, 626, 715, 719, 753, 713,
  753, 690, 691, 702, 777, 751, 795, 783, 804, 776, 830, 839,
  871, 856, 889, 877, 989, 992, 975, 993, 1010, 1023, 1015,
  1028, 1085, 1020, 1080, 1073
)
usethis::use_data(curlew_effort, overwrite = TRUE)

gerygone <- c(
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
)
usethis::use_data(gerygone, overwrite = TRUE)

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

grebe_records <- data.frame(
  time = c(
    1929, 1960, 1963, 1969, 1970, 1971,
    1972, 1982, 1985, 1986, 1988
  ),
  certainty_lower = c(
    0.99, 0.95, 0.75, 0.75, 0.1, 0.1,
    0.6, 0.6, 0.2, 0.2, 0.2
  ),
  certainty_upper = c(
    1, 0.99, 0.94, 0.94, 0.4, 0.4,
    0.8, 0.8, 0.8, 0.7, 0.5
  )
)

grebe <- data.frame(time = 1929:2017)
grebe$record <- grebe$time %in% grebe_records$time
grebe <- merge(grebe, grebe_records, by = "time", all.x = TRUE)
grebe[is.na(grebe)] <- 0
grebe$certainty <- (grebe$certainty_lower + grebe$certainty_upper) / 2
grebe <- grebe[, c(
  "time", "record", "certainty",
  "certainty_lower", "certainty_upper"
)]
rm(grebe_records)
usethis::use_data(grebe, overwrite = TRUE)

grebe_passive <- list(
  epsilon = c(0, 0.05), p_i = c(0.1, 0.65),
  p_r = c(0.4, 0.6)
)

grebe_active <- data.frame(
  time = c(1989, 1990, 1993, 1994, 1997, 1998, 1999, 2000, 2004, 2009),
  epsilon_lower = c(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7, 0.8, 0.8),
  epsilon_upper = c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.9, 0.95, 0.95),
  p_i_lower = rep(0.9, 10),
  p_i_upper = rep(0.95, 10),
  p_r_lower = rep(0.7, 10),
  p_r_upper = rep(0.9, 10)
)

grebe_surveys <- data.frame(time = 1929:2017)
grebe_surveys$survey <- grebe_surveys$time %in% grebe_active$time
grebe_surveys <- merge(grebe_surveys, grebe_active, by = "time", all.x = TRUE)

grebe_surveys$epsilon_lower[is.na(grebe_surveys$epsilon_lower)] <-
  grebe_passive$epsilon[1]
grebe_surveys$epsilon_upper[is.na(grebe_surveys$epsilon_upper)] <-
  grebe_passive$epsilon[2]
grebe_surveys$epsilon <-
  (grebe_surveys$epsilon_lower + grebe_surveys$epsilon_upper) / 2

grebe_surveys$p_i_lower[is.na(grebe_surveys$p_i_lower)] <-
  grebe_passive$p_i[1]
grebe_surveys$p_i_upper[is.na(grebe_surveys$p_i_upper)] <-
  grebe_passive$p_i[2]
grebe_surveys$p_i <-
  (grebe_surveys$p_i_lower + grebe_surveys$p_i_upper) / 2

grebe_surveys$p_r_lower[is.na(grebe_surveys$p_r_lower)] <-
  grebe_passive$p_r[1]
grebe_surveys$p_r_upper[is.na(grebe_surveys$p_r_upper)] <-
  grebe_passive$p_r[2]
grebe_surveys$p_r <-
  (grebe_surveys$p_r_lower + grebe_surveys$p_r_upper) / 2

grebe_surveys <- grebe_surveys[, c(
  "time", "survey", "epsilon", "epsilon_lower", "epsilon_upper", "p_i",
  "p_i_lower", "p_i_upper", "p_r", "p_r_lower", "p_r_upper"
)]
rm(grebe_passive, grebe_active)
usethis::use_data(grebe_surveys, overwrite = TRUE)

starling <- data.frame(
  time = 1930:2013,
  record = 1930:2013 %in% c(1930, 1995, 2008),
  certainty = ifelse(1930:2013 == 2008, 0.6, 1930:2013 %in%
    c(1930, 1995, 2008)),
  certainty_lower = ifelse(1930:2013 == 2008, 0.4, 1930:2013 %in%
    c(1930, 1995, 2008)),
  certainty_upper = ifelse(1930:2013 == 2008, 0.8, 1930:2013 %in%
    c(1930, 1995, 2008))
)
usethis::use_data(starling, overwrite = TRUE)

starling_surveys <- data.frame(
  time = 1930:2013,
  survey = 1930:2013 %in% c(1983, 2008, 2010)
)
starling_surveys$epsilon_lower <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.8^(1 / 3), 0)
starling_surveys$p_i_lower <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.8^(1 / 3), 0)
starling_surveys$p_r_lower <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.8^(1 / 3), 0)
starling_surveys$epsilon_upper <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.95^(1 / 3), 0)
starling_surveys$p_i_upper <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.95^(1 / 3), 0)
starling_surveys$p_r_upper <- ifelse(starling_surveys$time %in%
  c(1983, 2008, 2010), 0.95^(1 / 3), 0)
starling_surveys$epsilon <- (starling_surveys$epsilon_lower +
  starling_surveys$epsilon_upper) / 2
starling_surveys$p_i <- (starling_surveys$p_i_lower +
  starling_surveys$p_i_upper) / 2
starling_surveys$p_r <- (starling_surveys$p_r_lower +
  starling_surveys$p_r_upper) / 2
starling_surveys <- starling_surveys[, c(
  "time", "survey", "epsilon", "epsilon_lower", "epsilon_upper", "p_i",
  "p_i_lower", "p_i_upper", "p_r", "p_r_lower", "p_r_upper"
)]
usethis::use_data(starling_surveys, overwrite = TRUE)

weissmarshall_surveys <- c(
  0.0, 1.4, 1.6, 2.3, 3.9, 4.1, 4.7, 5.1, 5.4, 5.9, 6.1, 6.6, 7.0, 7.4, 7.9,
  8.6, 9.2, 9.9, 10.2, 10.4, 10.8, 11.2, 11.4, 11.7, 11.9, 12.5, 14.1, 14.9,
  15.4, 16.3, 18.4, 20.7, 22.5, 25.8, 27.6, 31.4
)
usethis::use_data(weissmarshall_surveys, overwrite = TRUE)

verneuilinoides <- c(
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0
)
usethis::use_data(verneuilinoides, overwrite = TRUE)

eggerellina_brevis <- c(
  1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0
)
usethis::use_data(eggerellina_brevis, overwrite = TRUE)

lee_s1 <- data.frame(
  certain = c(1, 0, 0, 0, 0, 1, 0, 0, 0),
  uncertain = c(0, 1, 0, 1, 0, 1, 0, 1, 0)
)
usethis::use_data(lee_s1, overwrite = TRUE)

bitterweed <- c(
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0,
  0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1,
  1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
  1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0,
  1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0
)
usethis::use_data(bitterweed, overwrite = TRUE)
