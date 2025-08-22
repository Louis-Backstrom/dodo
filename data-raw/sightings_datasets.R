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
