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

dirs <- c(kallisto, hisat, figure_dir, results_dir)

sapply(dirs, create_dir)
