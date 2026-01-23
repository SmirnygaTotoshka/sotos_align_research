library(dplyr)
library(purrr)
library(tidyr)
tables = Sys.glob("/media/HEAP-EPI/ANTON_ANALYSIS/nextflow/debug/2026/volodin_sotos/QC/statistics/*.csv")
output = "/media/HEAP-EPI/ANTON_ANALYSIS/nextflow/debug/2026/volodin_sotos/"
names = sapply(tables, basename)

total = do.call(rbind.data.frame, map2(tables, names, function(path, name){
    t = read.csv2(path)
    t$sample = name
    t
}, .progress = T))

locuses = unique(total$X3)

for (locus in locuses) {
    total %>%
        filter(X3 == locus) %>% 
        select(sample, mean_depth, coverage_width, uniformity, total_reads, 
               mean_methylation, mean_conversion, meth_by_caller, conversion_by_caller) %>%
        pivot_longer(cols = -sample, names_to = "parameter", values_to = "value") %>%
        pivot_wider(names_from = sample, values_from = value) %>% 
        tibble::column_to_rownames(var = "parameter") %>% 
        write.csv(file.path(output, paste0(locus,".csv")))
}
