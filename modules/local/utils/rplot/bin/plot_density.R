library(ggplot2)

#Input extracted annotation info
df <- read.table('${output_dir}/${filename}.txt', header=FALSE)

#Change colname
colnames(df)[1] <- '${filename}'

#Change to numeric
df[[1]] <- as.numeric(df[[1]])

#Calculate summary stats
summary_stats <- summary(df[[1]])

#Summary string
summary_string <- paste('Summary: Min =', summary_stats[1],
                        '1st Qu. =', summary_stats[2],
                        'Median =', summary_stats[3],
                        'Mean =', round(summary_stats[4], 2),
                        '3rd Qu. =', summary_stats[5],
                        'Max =', summary_stats[6])

#Density plot
p <- ggplot(df, aes_string(x='${filename}')) +
    geom_density(fill='#A7C7E7', colour='black', alpha=0.7) +
    ggtitle(paste('Density Plot for', '${filename}'),
            subtitle = summary_string) +
    labs(x = '${filename}', y = 'Density') +
    theme_minimal(base_size = 14) +
    theme(
        plot.title = element_text(face='bold', hjust = 0.5, size=16),
        axis.title = element_text(face='bold'),
        panel.grid.major = element_line(color='grey85', size=0.3),
        panel.grid.minor = element_line(color='grey90', size=0.2),
        plot.margin = unit(c(1, 1, 1, 1), 'cm')
    )

#Save the plot
ggsave(filename = '${output_dir}/${filename}_density_plot.png', plot = p, width = 15, height = 6, dpi = 300)
