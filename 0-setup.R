# Setup the directory structure
#

create_dir <- function(dir) {
  dir.create(dir, FALSE, TRUE)
}

data_dir <- "data"
figure_dir <- "figures"
results_dir <- "results"

dirs <- c(data_dir, figure_dir, results_dir)

sapply(dirs, create_dir)
