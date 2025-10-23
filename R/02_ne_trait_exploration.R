# Libraries --------------------------------------------------------------

library(here)
library(tidyverse)

# Read data --------------------------------------------------------------

ne_prefs <- readr::read_tsv(here::here(
  "data",
  "nematode_grazing_traits.tsv"
))

# Nematode preferece distributions ---------------------------------------

df01 <- ne_prefs %>%
  mutate(
    ne_volume_um3_rel_OP50 = ne_volume_um3 / ne_volume_um3_OP50
  ) %>%
  dplyr::select(-ne_volume_um3_OP50)


# Plots ------------------------------------------------------------------

# p01: worm volume
p01 <- df01 %>%
  ggplot() +
  geom_histogram(aes(x = ne_volume_um3, fill = ne_id), bins = 20) +
  scale_x_continuous(guide = guide_axis(angle = 45), trans = "log10") +
  labs(
    x = "Volume (um^3)",
    y = "N at volume bin",
    fill = "Nematode species",
    title = "Distribution of Nematode volumes from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id)

# p02: worm weight
p02 <- df01 %>%
  ggplot() +
  geom_histogram(aes(x = ne_dryweight_ug, fill = ne_id), bins = 20) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(
    x = "Dry weight (ug)",
    y = "N at weight bin",
    fill = "Nematode species",
    title = "Distribution of Nematode dry weight from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id, scales = "free_x")

# p03: worm preference
p03 <- df01 %>%
  ggplot() +
  geom_histogram(aes(x = lfc_spx_over_op50, fill = ne_id), bins = 20) +
  geom_vline(xintercept = 0, lty = 2) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(
    x = "Nematode relative preference (log2 scale)",
    y = "N at preference bin",
    fill = "Nematode species",
    title = "Distribution of Nematode preferences from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id)

# p04: worm preference binned
# We define a log2 fold change 20% around 0 to be the "Mid" range
p04 <- ne_prefs %>%
  mutate(
    pref_bin = cut(
      lfc_spx_over_op50,
      breaks = c(-Inf, log2(4 / 5), log2(5 / 4), Inf),
      labels = c("Low", "Mid", "High")
    )
  ) %>%
  ggplot() +
  geom_histogram(aes(x = pref_bin, fill = ne_id), stat = "count") +
  labs(
    x = "Binned nematode relative preference (log2 scale):\n(Low < log2(4/5) < Mid < log2(5/4) < High)",
    y = "N at preference bin",
    fill = "Nematode species",
    title = "Distribution of Nematode preferences from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id)

# p05: worm clearance
p05 <- df01 %>%
  ggplot() +
  geom_histogram(aes(x = ne_clearance, fill = ne_id), bins = 20) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(
    x = "Bacterial clearance",
    y = "N at clearance bin",
    fill = "Nematode species",
    title = "Distribution of Nematode clearances from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id)

# p06: worm volumes relative to OP50
p06 <- df01 %>%
  ggplot() +
  geom_histogram(
    aes(x = log2(ne_volume_um3_rel_OP50), fill = ne_id),
    bins = 20
  ) +
  scale_x_continuous(guide = guide_axis(angle = 45)) +
  labs(
    x = "Nematode relative size (log2 scale)",
    y = "N at relative size bin",
    fill = "Nematode species",
    title = "Distribution of Nematode relative volumes from feeding on 122 SynCom species"
  ) +
  facet_grid(~ne_id)

# Relationship between nematode preference and size ----------------------

p07 <- df01 %>%
  ggplot(aes(
    x = lfc_spx_over_op50,
    y = log2(ne_volume_um3_rel_OP50),
    color = ne_id
  )) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point() +
  geom_linerange(
    aes(
      xmin = lfc_spx_over_op50 - lfc_spx_over_op50_cl95,
      xmax = lfc_spx_over_op50 + lfc_spx_over_op50_cl95
    ),
    alpha = 0.25
  ) +
  labs(
    x = "Nematode relative preference (log2 scale)",
    y = "Nematode relative size (log2 scale)",
    color = "Nematode"
  ) +
  geom_smooth(method = "lm") +
  facet_grid(~ne_id, scales = "free_x")
