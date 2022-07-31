# Setup the directory structure
#

create_dir <- function(dir) {
  dir.create(dir, FALSE, TRUE)
}

data_dir <- "data"
kallisto <- file.path(data_dir, "kallisto")
hisat <- file.path(data_dir, "hisat2")
figure_dir <- "figures"
results_dir <- "results"
kallisto_res <- file.path(results_dir, "kallisto")
hisat_res <- file.path(results_dir, "hisat2")
kallisto_fig <- file.path(figure_dir, "kallisto")
hisat_fig <- file.path(figure_dir, "hisat2")

dirs <- c(kallisto,
          hisat,
          kallisto_res,
          hisat_res,
          kallisto_fig,
          hisat_fig)

sapply(dirs, create_dir)
