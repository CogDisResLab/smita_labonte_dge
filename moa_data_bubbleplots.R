suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

human_files <- list.files(file.path("raw", "processed_moa"),
  "^Human.*xlsx",
  full.names = TRUE
)

sheets <- human_files |>
  set_names(basename) |>
  map(~ excel_sheets(.x) |> keep(~ str_detect(.x, "[Ff]inal"))) |>
  enframe(name = "file", value = "sheets") |>
  unnest(cols = sheets)

sheet_data <- sheets |>
  mutate(
    data = map2(
      file.path("raw", "processed_moa", file),
      sheets, ~ read_excel(.x, sheet = .y)
    ),
    comparison = str_extract(sheets, "(.*)_[Ff]inal", 1),
    brain_region = str_extract(file, "Human_(.*).xlsx", 1)
  ) |>
  mutate(
    clean_data = map(
      data, ~ .x |>
        rename_with(~ str_replace_all(.x, " ", "_") |>
          str_squish() |>
          str_to_lower()) |>
        slice_max(moa_count, n = 5, na_rm = TRUE) |>
        select(moa, sum_discordant, dis)
    )
  ) |>
  select(-data, -file, -sheets) |>
  unnest(clean_data) |>
  select(brain_region, comparison, moa,
    discordance_score = dis,
    discordant_signature_count = sum_discordant
  ) |>
  filter(str_detect(moa, "[Ee]xperimental", negate = TRUE))

create_bubble_plot <- function(dataset, label) {
  g <- ggplot(
    dataset,
    aes(
      x = brain_region,
      y = moa,
      color = discordance_score,
      size = discordant_signature_count
    )
  )

  p <- g + geom_point(alpha = 0.75) +
    theme_minimal() +
    xlab("Brain Region") + ylab("Mechanism of Action") +
    scale_color_gradient(
      high = "white", low = "darkred",
      name = NULL
    ) +
    scale_size_continuous(range = c(5, 20), name = "# of Signatures") +
    theme(
      panel.background = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
      text = element_text(size = 28)
    ) +
    guides(
      color = guide_colorbar(
        position = "right", reverse = TRUE,
        theme = theme(
          legend.key.width = unit(5, "lines"),
          legend.key.height = unit(50, "lines")
        )
      )
    )


  ggsave(str_glue("{label}_moa_bubble_plot.png"), p,
    bg = "white", width = 6.5 * 3, height = 9 * 3, path = "figures"
  )
  ggsave(str_glue("{label}_moa_bubble_plot.svg"), p,
    bg = "white", width = 6.5 * 3, height = 9 * 3, path = "figures"
  )
}

male_data <- sheet_data |>
  filter(comparison == "Male") |>
  separate_longer_delim(
    cols = moa, delim = "|"
  )

create_bubble_plot(male_data, "male")


female_data <- sheet_data |>
  filter(comparison == "Female") |>
  separate_longer_delim(cols = moa, delim = "|")

create_bubble_plot(female_data, "female")
