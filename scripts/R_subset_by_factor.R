y <- structure(list(external_gene_name = structure(c(1L, 1L, 1L, 6L, 
      6L, 4L, 3L, 5L, 5L, 2L), .Label = c("FAM87B", "ISG15", "KLHL17", 
      "NOC2L", "PLEKHN1", "SAMD11"), class = "factor"), shortestDistance = c(99L, 
      2L, 1552L, 885L, 1677L, 2160L, 882L, 421L, 497L, 1246L), A = c(8.388, 
      9.913, 22.876, 19.826, 25.163, 122.766, 122.766, 40.414, 16.013, 
      155.554), B = c(3.403, 0.851, 0.851, 33.179, 26.373, 80.821, 
      80.821, 8.507, 17.015, 165.045), C = c(0.541, 1.622, 11.892, 
      2.162, 3.243, 17.298, 17.298, 0.541, 1.081, 33.515)), .Names = c("external_gene_name", 
      "shortestDistance", "A", "B", "C"), row.names = c(5L, 7L, 8L, 
      19L, 20L, 21L, 22L, 23L, 25L, 31L), class = "data.frame")


# external_gene_name shortestDistance       A       B      C
# 5              FAM87B               99   8.388   3.403  0.541
# 7              FAM87B                2   9.913   0.851  1.622
# 8              FAM87B             1552  22.876   0.851 11.892
# 19             SAMD11              885  19.826  33.179  2.162
# 20             SAMD11             1677  25.163  26.373  3.243
# 21              NOC2L             2160 122.766  80.821 17.298
# 22             KLHL17              882 122.766  80.821 17.298
# 23            PLEKHN1              421  40.414   8.507  0.541
# 25            PLEKHN1              497  16.013  17.015  1.081
# 31              ISG15             1246 155.554 165.045 33.515

library(plyr)
ddply(y, .(external_gene_name), summarise, shortestDistance=min(shortestDistance))
# external_gene_name shortestDistance
# 1             FAM87B                2
# 2              ISG15             1246
# 3             KLHL17              882
# 4              NOC2L             2160
# 5            PLEKHN1              421
# 6             SAMD11              885

do.call(rbind, by(y, y$external_gene_name, function(z) z[which.min(z$shortestDistance), ] ))
# external_gene_name shortestDistance       A       B      C
# FAM87B              FAM87B                2   9.913   0.851  1.622
# ISG15                ISG15             1246 155.554 165.045 33.515
# KLHL17              KLHL17              882 122.766  80.821 17.298
# NOC2L                NOC2L             2160 122.766  80.821 17.298
# PLEKHN1            PLEKHN1              421  40.414   8.507  0.541
# SAMD11              SAMD11              885  19.826  33.179  2.162


library(data.table)
setDT(y)
y[, .SD[shortestDistance == min(shortestDistance)], .(external_gene_name)]
# y[, .SD[which.min(shortestDistance)], by=external_gene_name]

#    external_gene_name shortestDistance       A       B      C
# 1:             FAM87B                2   9.913   0.851  1.622
# 2:             SAMD11              885  19.826  33.179  2.162
# 3:              NOC2L             2160 122.766  80.821 17.298
# 4:             KLHL17              882 122.766  80.821 17.298
# 5:            PLEKHN1              421  40.414   8.507  0.541
# 6:              ISG15             1246 155.554 165.045 33.515

setkey(y, external_gene_name, shortestDistance)
y[, head(.SD, 1), .(external_gene_name)]
#    external_gene_name shortestDistance       A       B      C
# 1:             FAM87B                2   9.913   0.851  1.622
# 2:              ISG15             1246 155.554 165.045 33.515
# 3:             KLHL17              882 122.766  80.821 17.298
# 4:              NOC2L             2160 122.766  80.821 17.298
# 5:            PLEKHN1              421  40.414   8.507  0.541
# 6:             SAMD11              885  19.826  33.179  2.162
