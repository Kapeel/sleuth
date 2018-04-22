#!/usr/bin/env Rscript

library("sleuth")
library("getopt")

args = commandArgs(trailingOnly=TRUE)

options<-matrix(c('design', 'd', 1, "character",
                  'fullmodel', 'fm', 1, "character",
                  'datadir', 's',1,"character",
		  'reducedmodel', 'rm', 1, "character",
		  'plotting_covariate', 'pc', 1, "character"),
                   ncol=4, byrow=TRUE)

ret.opts<-getopt(options, args)

base_dir <- ret.opts$datadir
design_matrix <- ret.opts$design
full_model <- ret.opts$fullmodel
reduced_model <- ret.opts$reducedmodel
plot_cov <- ret.opts$plotting_covariate

sample_id <- dir(file.path(base_dir,"kallisto_quant_output"))
print(sample_id)
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_quant_output", id))
print(kal_dirs)

s2c <- read.table(file.path(design_matrix), header = TRUE, stringsAsFactors=FALSE)

print(s2c)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
so <- sleuth_fit(so, as.formula(full_model), 'full')
so <- sleuth_fit(so, as.formula(reduced_model), 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

results_table <- sleuth_results(so, 'reduced:full', 'lrt')
write.table(results_table,"results_trans.tsv",sep="\t")

sleuth_significant <- dplyr::filter(results_table, qval <= 0.05)
write.table(sleuth_significant,"results_trans.sig.tsv",sep="\t")
 
sleuth_save(so, file="sleuth_obj")

jpeg('pca.jpg')
plot_pca(so, color_by = plot_cov, text_labels = TRUE)
dev.off()

jpeg('density.jpg')
plot_group_density(so, use_filtered = TRUE, units = "est_counts",trans = "log", grouping = plot_cov, offset = 1)
dev.off()

jpeg('mean_variance.jpg')
plot_mean_var(so, which_model = "full", point_alpha = 0.4,
  point_size = 2, point_colors = c("black", "dodgerblue"),
  smooth_alpha = 1, smooth_size = 0.75, smooth_color = "red")
dev.off()
