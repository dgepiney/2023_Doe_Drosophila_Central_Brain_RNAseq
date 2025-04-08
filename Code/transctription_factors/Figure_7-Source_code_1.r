library(Seurat)

my_data <- readRDS("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/r_data/t2_rnaseq.rds")

# List of transcription factor datasets from flybase
gene_lists <- list(
  bdtfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/basic_domain_transcription_factors.csv"),
  hdtfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/hdtf.csv"),
  hmgbtfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/high_mobility_group_box_transcription_factors.csv"),
  hthtfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/helix_turn_helix_transcription_factors.csv"),
  udbdtfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/unclassified_dna_binding_domain_transcription_factors.csv"),
  zftfs = read_csv("/Users/derek/University of Oregon Dropbox/Derek Epiney/lab_stuffs/gene_lists/zinc_finger_transcription_factors.csv")
)

# Iterate over the gene lists, run FindAllMarkers and save each output to a separate CSV file
for (list_name in names(gene_lists)) {
  gene_list <- gene_lists[[list_name]]
  
  markers <- FindAllMarkers(my_data, assay="RNA", only.pos = TRUE, features=intersect(rownames(my_data), gene_list$SYMBOL))
  
  output_file <- paste0("/Users/derek/Desktop/RNAseq/preprocessed_files/", list_name, "_markers_preprocessed.csv")
  
  write.csv(markers, output_file)
  
  print(paste0("Processed and saved: ", output_file))
}

