# Libraries --------------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
set.seed(12353)

# Load data --------------------------------------------------------------

ba_carbon <- readr::read_tsv(here::here(
  "data",
  "bacteria_carbon_utilization.tsv"
))

ba_carbon_cl <- readr::read_tsv(here::here(
  "data",
  "bacteria_carbon_utilization_classes.tsv"
))

ba_traits <- readr::read_tsv(here::here(
  "data",
  "bacteria_traits.tsv"
))

ne_prefs <- readr::read_tsv(here::here(
  "data",
  "nematode_grazing_traits.tsv"
))


# PCA of bacterial carbon preferences ------------------------------------

# convert to data.frame
ba_carbon_df <- ba_carbon %>%
  column_to_rownames(var = "ba_id") %>%
  data.frame()

# normalize by maximum in rows. PCA looks less pathological when normalizing
# this way
ba_carbon_df_nm <- ba_carbon_df / apply(ba_carbon_df, 1, max)

# perform the pca
pca_ord <- vegan::rda(ba_carbon_df_nm, scale = TRUE)

# get the loadings
carbon_loadings <- scores(
  pca_ord,
  choices = 1:2,
  display = "species",
  scaling = 0
) %>%
  data.frame() %>%
  rownames_to_column(var = "carbon_source") %>%
  left_join(
    ba_carbon_cl,
    by = join_by(carbon_source)
  )

# get the variance of the loadings
carbon_loadings_sd <- carbon_loadings %>%
  mutate(c = sqrt((PC1)^2 + (PC2)^2)) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    sd = sd(c),
    .by = "carbon_class"
  ) %>%
  # calculate slope of vector from origin
  mutate(m = (PC2 - 0) / (PC1 - 0)) %>%
  # calculate slope of perpendicular line = -m^-1
  mutate(m_perp = -(m^-1)) %>%
  # we want to use the standard deviation of the loading vector magnitudes as a
  # measure of the variance for the mean vector projected onto the PCA:
  # https://math.stackexchange.com/a/5088611 We assume that the sd is the
  # magnitude (hypotenuse) along the triangle starting at PC1, PC2. Thus, we
  # need to calculate dx and dy for the triangle using pythagorean theorem
  mutate(dx = sd / sqrt(1 + m_perp^2)) %>%
  mutate(dy = m_perp * dx) %>%
  mutate(x1 = PC1 + dx, y1 = PC2 + dy, x2 = PC1 - dx, y2 = PC2 - dy)

# make a tibble for plotting the polygon that is the "triangle" of the standard
# deviation around the mean vector
polygon_dat <- carbon_loadings_sd %>%
  dplyr::select(carbon_class, x1:y2) %>%
  pivot_longer(
    cols = c(-carbon_class),
    names_pattern = "(x|y)()",
    names_to = c(".value", "names")
  ) %>%
  bind_rows(
    expand_grid(
      carbon_class = c("amino_acid", "organic_acid", "others", "sugar"),
      x = 0,
      y = 0
    )
  )

# calculate the placement of the species
species_sites <- scores(
  pca_ord,
  choices = 1:2,
  display = "sites",
  scaling = 0
) %>%
  data.frame() %>%
  rownames_to_column(var = "ba_id") %>%
  left_join(ba_traits)

# get percentage of variance explained by each PC
varexplained <- round(pca_ord$CA$eig / pca_ord$tot.chi * 100, 1)

# make the ordination showing species
p01 <- ggplot() +
  geom_point(data = species_sites, aes(x = PC1, y = PC2, color = ba_phylum)) +
  coord_fixed(xlim = c(-0.4, 0.2), ylim = c(-0.25, 0.25)) +
  labs(
    x = paste0("PC1 (", varexplained["PC1"], "%)"),
    y = paste0("PC1 (", varexplained["PC2"], "%)"),
    color = "Phylum"
  )

# make the plot showing loadings
p02 <- ggplot() +
  geom_point(data = species_sites, aes(x = PC1, y = PC2)) +
  geom_segment(
    data = carbon_loadings,
    arrow = grid::arrow(),
    alpha = 0.25,
    show.legend = FALSE,
    aes(x = 0, xend = PC1, y = 0, yend = PC2, color = carbon_class)
  ) +
  geom_polygon(
    data = polygon_dat,
    aes(x = x, y = y, color = NULL, fill = carbon_class, group = carbon_class),
    alpha = 0.25,
    show.legend = FALSE
  ) +
  geom_segment(
    data = carbon_loadings_sd,
    size = 2,
    aes(x = 0, xend = PC1, y = 0, yend = PC2, color = carbon_class)
  ) +
  coord_fixed(xlim = c(-0.4, 0.2), ylim = c(-0.25, 0.25)) +
  labs(
    x = paste0("PC1 (", varexplained["PC1"], "%)"),
    y = paste0("PC1 (", varexplained["PC2"], "%)"),
    color = "Carbon class",
    fill = ""
  )

# pca summary
summary(pca_ord)

# histogram of loadings
p03 <- data.frame(pca_ord$CA$v) %>%
  rownames_to_column(var = "carbon_source") %>%
  left_join(ba_carbon_cl) %>%
  ggplot(aes(x = PC1)) +
  geom_histogram(aes(fill = carbon_class), bins = 10) +
  facet_grid(~carbon_class)

p04 <- ordiplot(pca_ord, type = "n") %>%
  points("sites", pch = 16, col = "red") %>%
  text("species", arrows = T, length = 0.05, col = "blue")

# Sugar amino acid preference --------------------------------------------

# combine the classes and use data and calculate the Sugar Amino Acid Preference
# (SAAP)
saap <- ba_carbon %>%
  pivot_longer(-ba_id, names_to = "carbon_source", values_to = "od600") %>%
  left_join(ba_carbon_cl) %>%
  summarize(
    AA = mean(od600[carbon_class == "amino_acid"]),
    Sug = mean(od600[carbon_class == "sugar"]),
    .by = ba_id
  ) %>%
  # calculate SAAP
  mutate(SAAP = (Sug - AA) / (Sug + AA))

# distribution of SAAP values
p05 <- left_join(saap, ba_traits) %>%
  ggplot(aes(x = SAAP, fill = ba_phylum)) +
  geom_histogram(bins = 20) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  facet_grid(~ba_phylum)

p06 <- left_join(saap, ne_prefs) %>%
  left_join(ba_traits) %>%
  filter(!is.na(ba_phylum)) %>%
  ggplot(aes(
    x = SAAP,
    y = lfc_spx_over_op50,
    color = ba_order,
    shape = ne_id
  )) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_linerange(
    aes(
      ymin = lfc_spx_over_op50 - lfc_spx_over_op50_cl95,
      ymax = lfc_spx_over_op50 + lfc_spx_over_op50_cl95
    ),
    alpha = 0.15
  ) +
  geom_point() +
  labs(
    x = "Sugar-Amino Acid Preference (SAAP)",
    y = "Preference relative to OP50: Log2(SpX/OP50)",
    color = "Taxonomic Order",
    shape = "Nematode"
  ) +
  facet_grid(~ba_phylum)

# join saap to bacteria traits
p07 <- left_join(saap, ba_traits) %>%
  pivot_longer(cols = c(ba_volume_um3:ba_biofilm)) %>%
  ggplot(aes(x = SAAP, y = value)) +
  geom_point(aes(color = ba_phylum)) +
  geom_smooth(method = "lm") +
  scale_y_log10() +
  facet_wrap(name ~ ., scale = "free_y", ncol = 2)
