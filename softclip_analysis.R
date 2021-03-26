library(tidyverse)
library(filesstrings)
library(TSRexploreR)

options(scipen = 999)

# Generate sample sheet
samples <- data.frame(sample_name = before_last_dot(list.files("BAMs/")),
                                                    file_1 = list.files("BAMs", full.names = TRUE),
                                                    file_2 = NA)

# Store paths to assembly and annotation
assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")

# Generate TSRexploreR object
exp <- tsr_explorer(sample_sheet = samples,
                    genome_assembly = assembly,
                    genome_annotation = annotation)

# Import BAMs
exp <- import_bams(exp, paired = FALSE)

# Soft-clipping analysis

## Generate histograms of soft-clipped base frequency. This approach is taken to enable labeling of columns with 
## numbers of soft-clipped bases.
softclip_histogram(exp, return_table = TRUE) %>%
    group_by(sample) %>%
    mutate(sample_total = sum(total)) %>%
    mutate(n_soft_frac = total/sample_total) %>%
    ggplot(aes(x = n_soft, y = total)) +
    geom_col(aes(fill = sample)) +
    geom_text(aes(label = round(n_soft_frac, 4)), vjust = -1) +
    theme_bw() +
    theme(legend.position = "none") +
    scale_fill_viridis_d() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +        
    facet_wrap(~sample, ncol = 2, scales = "free_y")

softclip_composition(exp, ncol = 1) +
    theme_bw() +    
    theme(axis.text.x = element_text(angle = 90)) +
    scale_fill_viridis_d()
