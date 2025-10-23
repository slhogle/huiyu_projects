# libraries --------------------------------------------------------------

library(here)
library(tidyverse)
library(errors)

# Read raw data ----------------------------------------------------------

# giant single file of all bacteria traits
df01 <- readr::read_csv(here::here(
  "_data_raw",
  "traits",
  "TOTAL_Final_merged.csv"
))
# nematode feeding preference traits
df02 <- readr::read_csv(here::here(
  "_data_raw",
  "traits",
  "nem.pre.size.gra.csv"
))
# nematode feeding preference traits in raw form
df02A <- readr::read_csv(here::here(
  "_data_raw",
  "traits",
  "preference.all.csv"
))
# bacteria taxonomy
df03 <- readr::read_csv(here::here("_data_raw", "traits", "taxonomy.csv"))
# SMILES strings for VOCs
df04 <- readr::read_csv(here::here("_data_raw", "vocs", "vocs_SMILES.csv"))
# bacteria voc production
df05 <- readr::read_csv(here::here("_data_raw", "vocs", "vocs.csv"))
# bacteria use of carbon substrates
df06 <- read_csv(here::here(
  "_data_raw",
  "carbon",
  "carbon.utilization.csv"
))
# carbon subtrate chemical classes
df07 <- read_csv(here::here("_data_raw", "carbon", "carbon_class2.csv"))

# Nematode traits --------------------------------------------------------

# formatting the raw preference data. We calculate the mean and 95% confidence
# interval using function Hmisc::smean.cl.boot. It calculates the mean and then
# does a nonparametric boostrap (n=1000) to estimate the 95% confidence interval

df02A_spx <- summarize(
  df02A,
  ggplot2::mean_cl_boot(spe.ID.num),
  .by = c(nematode, spe.ID)
) %>%
  dplyr::mutate(cl95_spx = (ymax - ymin) / 2) %>%
  dplyr::select(nematode, spe.ID, mn_spx = y, cl95_spx)

df02A_op50 <- summarize(
  df02A,
  ggplot2::mean_cl_boot(spe.control.number),
  nreps = n(),
  .by = c(nematode, spe.ID)
) %>%
  dplyr::mutate(cl95_op50 = (ymax - ymin) / 2) %>%
  dplyr::select(nematode, spe.ID, mn_op50 = y, cl95_op50, nreps)

df02B <- left_join(df02A_spx, df02A_op50, by = join_by(nematode, spe.ID)) %>%
  mutate(
    spx = errors::set_errors(mn_spx, cl95_spx),
    op50 = errors::set_errors(mn_op50, cl95_op50)
  ) %>%
  mutate(
    spx_over_total = spx / (spx + op50),
    lfc_spx_over_op50 = log2(spx / op50)
  ) %>%
  mutate(spx_over_total_cl95 = errors::errors(spx_over_total)) %>%
  mutate(spx_over_total = errors::drop_errors(spx_over_total)) %>%
  mutate(lfc_spx_over_op50_cl95 = errors::errors(lfc_spx_over_op50)) %>%
  mutate(lfc_spx_over_op50 = errors::drop_errors(lfc_spx_over_op50)) %>%
  dplyr::select(
    ne_id = nematode,
    ba_id = spe.ID,
    lfc_spx_over_op50,
    lfc_spx_over_op50_cl95,
    spx_over_total,
    spx_over_total_cl95,
    mn_spx,
    cl95_spx,
    mn_op50,
    cl95_op50,
    nreps
  )

df01 %>%
  dplyr::select(
    ba_id = spe.ID,
    ne_id = nematode,
    ne_volume_um3 = ne.bodysize.um3,
    ne_dryweight_ug = ne.dryweight.ug.,
    ne_preference = ne.preference,
    ne_clearance = ne.grazing
  ) %>%
  left_join(df02, by = join_by(ne_id == nematode, ba_id == spe.ID)) %>%
  dplyr::select(
    ba_id:ne_clearance,
    ne_volume_um3_OP50 = mean_size_OP50
  ) %>%
  left_join(df02B, by = join_by(ba_id, ne_id)) %>%
  relocate(ne_volume_um3_OP50, .after = ne_volume_um3) %>%
  relocate(ne_id) %>%
  arrange(ne_id, ba_id) %>%
  write_tsv(here::here("data", "nematode_grazing_traits.tsv"))

# Bacteria traits --------------------------------------------------------

df01 %>%
  dplyr::select(
    ba_id = spe.ID,
    ba_gram = gram,
    ba_shape = shape,
    ba_volume_um3 = ba.bodysize.um3.,
    ba_dryweight_ug = ba.dryweight.ug.,
    ba_protein = protein,
    ba_sugar = sugar,
    ba_respiration = CO2.activity,
    ba_siderophore = sidero,
    ba_biofilm = biofilm
  ) %>%
  left_join(df03, by = join_by(ba_id == Species)) %>%
  dplyr::rename(
    ba_genus = Genus,
    ba_family = Family,
    ba_order = Order,
    ba_class = Class,
    ba_phylum = Phylum
  ) %>%
  distinct() %>%
  arrange(ba_id) %>%
  write_tsv(here::here("data", "bacteria_traits.tsv"))

# Bacteria bioynthetic gene clusters

df01 %>%
  dplyr::select(
    ba_id = spe.ID,
    acyl_amino_acids:other
  ) %>%
  distinct() %>%
  pivot_longer(-ba_id, names_to = "bgc_cluster", values_to = "gene_count") %>%
  mutate(bgc_cluster = str_replace_all(bgc_cluster, "_", " ")) %>%
  mutate(bgc_cluster = str_replace_all(bgc_cluster, "\\.", " ")) %>%
  arrange(ba_id, bgc_cluster) %>%
  write_tsv(here::here("data", "bacteria_bgc_counts.tsv"))


# Bacteria VOC production ------------------------------------------------

# read smiles file and create a numeric ID that makes later joining much easier
smls <- df04 %>%
  arrange(VOC_name) %>%
  group_by(VOC_name) %>%
  mutate(
    VOC_id = paste0(
      "VOC_",
      str_pad(cur_group_id(), width = 3, side = "left", pad = "0")
    )
  ) %>%
  ungroup()

# save this mapping file for later
readr::write_tsv(
  smls,
  here::here("data", "vocs_SMILES_id.tsv"),
  col_names = FALSE
)

# write a SMI file in way that read.SMIset expects no header and the SMILES
# string first
smls %>%
  dplyr::select(SMILES, VOC_id) %>%
  readr::write_tsv(here::here("data", "vocs_SMILES.smi"), col_names = FALSE)

# also format the VOC concentration data and save for later
vocs <- df05 %>%
  # species PAE.3 has no data so drop it here
  dplyr::filter(spe.ID != "PAE3") %>%
  pivot_longer(
    cols = c(-spe.ID, -Sample_File),
    names_to = "VOC_name",
    values_to = "VOC_counts"
  ) %>%
  arrange(spe.ID, VOC_name) %>%
  group_by(VOC_name) %>%
  mutate(
    VOC_id = paste0(
      "VOC_",
      str_pad(cur_group_id(), width = 3, side = "left", pad = "0")
    )
  ) %>%
  ungroup() %>%
  # there are duplicates of some compounds for some species. Take the mean of
  # those here
  summarize(
    VOC_counts = mean(VOC_counts, na.rm = TRUE),
    .by = c(spe.ID, VOC_name, VOC_id)
  ) %>%
  replace_na(list(VOC_counts = 0)) %>%
  dplyr::select(ba_id = spe.ID, VOC_id, VOC_counts, VOC_name)

# save VOC counts data
readr::write_tsv(vocs, here::here("data", "bacteria_vocs.tsv"))


# Bacteria carbon utlization ---------------------------------------------

# read in the carbon use data and perform some formatting/renaming
cutil <- df06 %>%
  rename(Carbon = sources) %>%
  mutate(Carbon = str_replace_all(Carbon, "-", "_")) %>%
  mutate(Carbon = str_replace_all(Carbon, " ", "."))

ccombo <- left_join(cutil, df07) %>%
  pivot_longer(c(-Carbon, -Class), names_to = "spe.ID", values_to = "od600") %>%
  dplyr::select(
    ba_id = spe.ID,
    carbon_source = Carbon,
    carbon_class = Class,
    od600
  ) %>%
  mutate(
    carbon_source = stringr::str_to_lower(carbon_source),
    carbon_class = stringr::str_to_lower(carbon_class)
  ) %>%
  mutate(
    carbon_source = str_replace_all(carbon_source, "\\.", "_"),
    carbon_class = str_replace_all(carbon_class, "\\.", "_")
  ) %>%
  mutate(
    carbon_source = if_else(
      carbon_source == "2_oxoglutaric",
      "two_oxoglutaric_acid",
      carbon_source
    )
  ) %>%
  arrange(ba_id, carbon_class, carbon_source)

ccombo %>%
  dplyr::select(carbon_source, carbon_class) %>%
  distinct() %>%
  readr::write_tsv(here::here(
    "data",
    "bacteria_carbon_utilization_classes.tsv"
  ))

ccombo %>%
  dplyr::select(-carbon_class) %>%
  pivot_wider(
    id_cols = c(ba_id),
    names_from = carbon_source,
    values_from = od600
  ) %>%
  readr::write_tsv(here::here(
    "data",
    "bacteria_carbon_utilization.tsv"
  ))
