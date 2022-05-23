devtools::load_all(here::here("../WMGSPipeline"))

mymre <- filter_samples(mre = mymre, sample_ids = metadata$SampleID)

mymre@metadata@metadata_df <- as.tibble(metadata)
mymre@taxa@metaphlan@phyloseq@sam_data <- metadata %>%
  phyloseq::sample_data(.)
mymre@taxa@metaphlan@phyloseq_sec@sam_data <- metadata %>%
  phyloseq::sample_data(.)




mymre@metadata@categorical_vals <- here::here("Metadata", "CategoricalVariables.txt") %>%
  read.delim(., header = T) %>%
  tibble()

mymre@metadata@numeric_vals <- here::here("Metadata", "NumericalVariables.txt") %>%
  read.delim(., header = T) %>%
  tibble()

mymre@metadata@longitudinal_vals <- here::here("Metadata", "LongitudinalVariables.txt") %>%
  read.delim(., header = T) %>%
  tibble()

