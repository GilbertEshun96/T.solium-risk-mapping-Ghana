library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Data: one row per study
studies <- tribble(
  ~StudyID, ~Year, ~Diseases,
  "S1",  2021, "HCC;PCC",
  "S2",  2020, "PCC",
  "S3",  2024, "NCC",
  "S4",  2022, "HTT",
  "S5",  2015, "PCC",
  "S6",  2025, "HTT",
  "S7",  2012, "PCC;HTT",
  "S8",  2018, "HTT",
  "S9",  2025, "HTT",
  "S10", 2013, "HCC-HTT",
  "S11", 1999, "PCC",
  "S12", 2011, "HTT",
  "S13", 2010, "HTT",
  "S14", 2011, "HTT",
  "S15", 2018, "HTT",
  "S16", 2024, "HTT",
  "S17", 2019, "HCC;PCC"
)

# Arrange studies within each year
# Each separate study occupies one full count on the y-axis.
# Studies assessing two diseases are split within the same count.
study_stack <- studies %>%
  arrange(Year, StudyID) %>%
  group_by(Year) %>%
  mutate(
    stack_no = row_number(),
    ymin = stack_no - 1,
    ymax = stack_no
  ) %>%
  ungroup()

# Split bars horizontally where one study assessed two diseases concurrently
plot_data <- study_stack %>%
  separate_rows(Diseases, sep = ";") %>%
  group_by(StudyID, Year) %>%
  mutate(
    disease_no = row_number(),
    disease_total = n(),
    xmin = Year - 0.38 + (disease_no - 1) * (0.76 / disease_total),
    xmax = Year - 0.38 + disease_no * (0.76 / disease_total)
  ) %>%
  ungroup()

disease_cols <- c(
  "HCC" = "#e41a1c",
  "HCC-HTT" = "#4daf4a",
  "HTT" = "#377eb8",
  "NCC" = "#ff7f00",
  "PCC" = "#984ea3"
)

# Main plot
p_main <- ggplot(plot_data) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = Diseases
    ),
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(values = disease_cols) +
  scale_x_continuous(
    breaks = 1999:2025,
    limits = c(1998.5, 2025.5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 0:2,
    limits = c(0, 2.15),
    expand = c(0, 0)
  ) +
  labs(
    x = "Year",
    y = "Count"
  ) +
  theme_bw() +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Manual one-column legend
legend_panel <- ggplot() +

  # Disease legend
  annotate("text", x = 0, y = 9.60, label = "Disease", hjust = 0, size = 4.2) +

  geom_rect(aes(xmin = 0, xmax = 0.35, ymin = 8.95, ymax = 9.30),
            fill = disease_cols["HCC"], colour = "black", linewidth = 0.25) +
  annotate("text", x = 0.48, y = 9.125, label = "HCC", hjust = 0, size = 3.8) +

  geom_rect(aes(xmin = 0, xmax = 0.35, ymin = 8.45, ymax = 8.80),
            fill = disease_cols["HCC-HTT"], colour = "black", linewidth = 0.25) +
  annotate("text", x = 0.48, y = 8.625, label = "HCC-HTT", hjust = 0, size = 3.8) +

  geom_rect(aes(xmin = 0, xmax = 0.35, ymin = 7.95, ymax = 8.30),
            fill = disease_cols["HTT"], colour = "black", linewidth = 0.25) +
  annotate("text", x = 0.48, y = 8.125, label = "HTT", hjust = 0, size = 3.8) +

  geom_rect(aes(xmin = 0, xmax = 0.35, ymin = 7.45, ymax = 7.80),
            fill = disease_cols["NCC"], colour = "black", linewidth = 0.25) +
  annotate("text", x = 0.48, y = 7.625, label = "NCC", hjust = 0, size = 3.8) +

  geom_rect(aes(xmin = 0, xmax = 0.35, ymin = 6.95, ymax = 7.30),
            fill = disease_cols["PCC"], colour = "black", linewidth = 0.25) +
  annotate("text", x = 0.48, y = 7.125, label = "PCC", hjust = 0, size = 3.8) +

  # Study type legend, closer to Disease legend
  annotate("text", x = 0, y = 6.55, label = "Study type", hjust = 0, size = 4.2) +

  # Assessing one disease
  geom_rect(aes(xmin = 0, xmax = 0.45, ymin = 5.70, ymax = 6.30),
            fill = "grey90", colour = "black", linewidth = 0.35) +
  annotate("text", x = 0.65, y = 6.00,
           label = "Assessing\n1 disease",
           hjust = 0, size = 3.6) +

  # Assessing two diseases concurrently
  geom_rect(aes(xmin = 0, xmax = 0.45, ymin = 4.85, ymax = 5.45),
            fill = "grey90", colour = "black", linewidth = 0.35) +
  geom_segment(aes(x = 0.225, xend = 0.225, y = 4.85, yend = 5.45),
               colour = "black", linewidth = 0.35) +
  annotate("text", x = 0.65, y = 5.15,
           label = "Assessing\n2 diseases\nconcurrently",
           hjust = 0, size = 3.6) +

  xlim(0, 2.7) +
  ylim(4.6, 9.9) +
  theme_void()

# Final plot
final_plot <- p_main + legend_panel +
  plot_layout(widths = c(6.3, 0.85))

print(final_plot)

ggsave(
  "publication_year_distribution_final954.png",
  plot = final_plot,
  width = 13,
  height = 7,
  dpi = 300
)
