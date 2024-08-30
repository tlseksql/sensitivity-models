ctrpv2 <- read.csv(file = 'pset_data/treatment_response/metadata/manual_annotations/ctrpv2_annotated.csv')
prism <- read.csv(file = 'pset_data/treatment_response/metadata/manual_annotations/prism_annotated.csv')
gdsc1 <- read.csv(file = 'pset_data/treatment_response/metadata/manual_annotations/gdsc1_annotated.csv')
gdsc2 <- read.csv(file = 'pset_data/treatment_response/metadata/manual_annotations/gdsc2_annotated.csv')

library(dplyr)
library(tidyr)

# Function to merge all unique manual annotation(s)
annotation_merger <- function(annotation) {
  unique_annotation <- unique(unlist(strsplit(annotation, ", ")))
  paste(unique_annotation, collapse = ", ")
}

# Bind datasets
binder <- bind_rows(
  ctrpv2 %>% select(treatment, manual_annotation),
  prism %>% select(treatment, manual_annotation),
  gdsc1 %>% select(treatment, manual_annotation),
  gdsc2 %>% select(treatment, manual_annotation)
)

# Group by treatment to merge annotation(s)
metadata <- binder %>%
  group_by(treatment) %>%
  summarise(manual_annotation = annotation_merger(manual_annotation)) %>%
  ungroup()

saveRDS(metadata, file = 'pset_data/treatment_response/metadata/manual_annotations/COMPLETE LIST.rds')
