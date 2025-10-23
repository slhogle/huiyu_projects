# libraries
library(here)
library(tidyverse)
library(ChemmineR)

set.seed(45617)
# read the tab-delimited file containing SMILES strings
smiset <- read.SMIset(here::here("data", "vocs_SMILES.smi"))

# convert the SMILES dataset to SDF format
sdfset <- smiles2sdf(smiset)

# make atom pair descriptor database for searching
apset <- sdf2ap(sdfset)

# make fingerprints from descriptor vectors such as atom pairs stored in APset
fpset <- desc2fp(apset)

# cluster the fingerprints using the Tanimoto coefficient. This operates at cutoff levels of 0.5
# 0.7 and 0.9. You can then choose which is best
clusters <- cmp.cluster(
  fpset,
  cutoff = c(0.7, 0.8, 0.9),
  method = "Tanimoto",
  quiet = FALSE
)

clusters_fmt <- clusters %>%
  dplyr::select(VOC_id = ids, CLID_0.7, CLID_0.8, CLID_0.9) %>%
  mutate(
    CLID_0.7 = paste0(
      "cut70_",
      str_pad(CLID_0.7, width = 3, side = "left", pad = "0")
    ),
    CLID_0.8 = paste0(
      "cut80_",
      str_pad(CLID_0.8, width = 3, side = "left", pad = "0")
    ),
    CLID_0.9 = paste0(
      "cut90_",
      str_pad(CLID_0.9, width = 3, side = "left", pad = "0")
    )
  )

# save the formatted cluster file for later
readr::write_tsv(
  clusters_fmt,
  here::here("data", "vocs_fingerprint_clusters.tsv")
)
