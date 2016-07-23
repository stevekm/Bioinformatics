library("ggplot2")

users_df <- data.frame( Type = c("Active","Semi-Active","Inactive","Missing"),
												Number = c(1011,133,1978,1556))

# regular bar plot
png(filename = "/Users/steve/plot1.png",width = 600,height = 600)
ggplot(data = users_df, aes(x = Type, y = Number,fill = Type)) + 
	geom_bar(stat="identity") + ggtitle("Sample Breakdown")
dev.off()

# single stacked
png(filename = "/Users/steve/plot2.png",width = 600,height = 600)
ggplot(data = users_df, aes(x = 1, y = Number,fill = Type)) + 
	geom_bar(stat="identity") + xlab("") + coord_flip() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) + ggtitle("Sample Breakdown") #+labs(fill="Type")
dev.off()

