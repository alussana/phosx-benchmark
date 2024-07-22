#!/usr/bin/env Rscript

load(file="cptacData.rda")
load(file="cptacGS.rda")

write.table(
    cptacData,
    file = "datasets/cptac/data.tsv",
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
)

cancer_types = c("BRCA", "CCRCC", "GBM", "HNSCC", "LSCC", "LUAD", "UCEC")

for (cancer_type in cancer_types) {
    for (kinase in names(cptacGS[[cancer_type]]$GS_pos_pairs)) {
        for (sample in cptacGS[[cancer_type]]$GS_pos_pairs[[kinase]]) {
            cat(
                sprintf("%s\t%s\t1\n", sample, kinase),
                    file = "datasets/cptac/metadata.tsv",
                    append = TRUE
                )
        }
    }
    for (kinase in names(cptacGS[[cancer_type]]$GS_neg_pairs)) {
        for (sample in cptacGS[[cancer_type]]$GS_neg_pairs[[kinase]]) {
            cat(
                sprintf("%s\t%s\t-1\n", sample, kinase),
                    file = "datasets/cptac/metadata.tsv",
                    append = TRUE
                )
        }
    }
}