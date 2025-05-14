process R_PLOT_DENSITY {
    label 'process_medium'
    conda "bioconda::r-base=4.3.2"

    input:
        tuple val(meta), val(annotation), file(annotation_txt)

    output:
        tuple val(meta), val(annotation), file("${annotation}_density_plot.png"), file("versions.yml"), emit: plot_png

    script:
    """
    Rscript -e "
    library(ggplot2)
    df <- read.table('${annotation_txt}', header=FALSE)
    colnames(df)[1] <- '${annotation}'
    df[[1]] <- as.numeric(df[[1]])
    summary_stats <- summary(df[[1]])
    summary_string <- paste('Summary: Min =', summary_stats[1],
                            '1st Qu. =', summary_stats[2],
                            'Median =', summary_stats[3],
                            'Mean =', round(summary_stats[4], 2),
                            '3rd Qu. =', summary_stats[5],
                            'Max =', summary_stats[6])
    p <- ggplot(df, aes_string(x='${annotation}')) +
        geom_density(fill='#A7C7E7', colour='black', alpha=0.7) +
        ggtitle(paste('Density Plot for', '${annotation}'),
                subtitle = summary_string) +
        labs(x = '${annotation}', y = 'Density') +
        theme_minimal(base_size = 14) +
        theme(
            plot.title = element_text(face='bold', hjust = 0.5, size=16),
            axis.title = element_text(face='bold'),
            panel.grid.major = element_line(color='grey85', size=0.3),
            panel.grid.minor = element_line(color='grey90', size=0.2),
            plot.margin = unit(c(1, 1, 1, 1), 'cm')
        )
    ggsave(filename = '${annotation}_density_plot.png', plot = p, width = 15, height = 6, dpi = 300)
    "
    cat <<-END_VERS > versions.yml
    R: \$(R --version | head -1)
    END_VERS
    """
}
