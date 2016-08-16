#/usr/bin/env Rscript

library("reshape2")
library("ggplot2")

dates_file <- "/Users/steve/codes/dates2.tsv"

dates_df <- read.delim(dates_file,header = FALSE)

# convert to long format
dates_df_long <- reshape2::melt(dates_df,id.vars='V1', value.name = "userID")

# removed empty cells
dates_df_long <- dates_df_long[!apply(dates_df_long[3] == "", 1, all),]

# remove V2 col
dates_df_long <- dates_df_long[,c(1,3)]

# split the first column into two columns
dates_df_long <- with(dates_df_long, cbind(userID, colsplit(V1, pattern = " ", names = c('month', 'day'))))

# add months as factors
dates_df_long[["month"]] <- factor(dates_df_long[["month"]], levels = month.name)

write.table(dates_df_long,file = "/Users/steve/codes/dates2_long.tsv",quote = FALSE,sep = '\t',row.names = FALSE,col.names = TRUE)


plot_cols <- rainbow(n = length(unique(dates_df_long[["month"]])))

# make a plot
png(filename = "/Users/steve/codes/dates.png",width = 600,height = 600)
ggplot(data = dates_df_long, aes(x = month, fill = month)) + 
    geom_bar() + 
    coord_flip() + 
    ggtitle("Dates") + 
    theme_minimal() #+ scale_fill_manual( values = plot_cols)
dev.off()    
    

