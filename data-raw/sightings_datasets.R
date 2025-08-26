monk_seal <- c(1915, 1922, 1932, 1948, 1952)
usethis::use_data(monk_seal, overwrite = TRUE)

ferret1 <- c(6, 7, 8, 10, 17, 18, 19, 20, 21, 22, 30, 31, 41, 44, 46, 53, 57, 58,
             66, 90, 117, 118, 122, 123, 127, 139, 151, 153)
usethis::use_data(ferret1, overwrite = TRUE)

ferret2 <- data.frame(
  time = 1972:1992,
  records = c(4, 6, 2, 3, 3, 1, 0, 1, 0, 2, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0)
)
usethis::use_data(ferret2, overwrite = TRUE)

metrarabdotos <- c(0, 0, 0, 0, 0, 0, 0.05, 0.10, 0.20, 0.30, 1.40)
usethis::use_data(metrarabdotos, overwrite = TRUE)

burgmanf1b <- data.frame(
  time = 1:16,
  records = c(1, 2, 0, 2, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0)
)
usethis::use_data(burgmanf1b, overwrite = TRUE)

dodo <- c(1598, 1601, 1602, 1607, 1611, 1628, 1628, 1631, 1638, 1662)
usethis::use_data(dodo, overwrite = TRUE)

mammoth <- -c(11.50, 11.54, 11.91, 12.12, 12.19, 12.34, 12.43, 12.44, 12.48,
              12.49, 12.51, 12.58, 12.68, 12.88, 13.06, 13.23, 13.29, 13.34,
              13.38, 13.41, 13.44, 13.66, 13.69, 14.02, 14.09) * 1000
usethis::use_data(mammoth, overwrite = TRUE)
