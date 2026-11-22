args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) {
  stop("Usage: Rscript reproduce.R [output-directory]", call. = FALSE)
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
root <- dirname(script)
out <- if (length(args) == 1L) args[[1]] else file.path(root, "results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Could not create the output directory.", call. = FALSE)

read_input <- function(name) {
  read.csv(file.path(root, "data", name), na.strings = c("", "NA"),
           check.names = FALSE, stringsAsFactors = FALSE)
}
save_table <- function(value, name) {
  write.csv(value, file.path(out, name), row.names = FALSE)
}

# Figure 1D and E: intersect the two saved cNMF rank lists.
wound <- read_input("Wound_W2_top200.csv")
pdac <- read_input("PDAC_F1_top200.csv")
annotations <- read_input("matrisome_annotations.csv")
stopifnot(identical(wound$Rank, seq_len(200L)), identical(pdac$Rank, seq_len(200L)),
          !anyDuplicated(wound$Gene), !anyDuplicated(pdac$Gene), !anyDuplicated(annotations$gene))
sets <- list(Shared74 = intersect(wound$Gene, pdac$Gene),
             Wound126 = setdiff(wound$Gene, pdac$Gene),
             PDAC126 = setdiff(pdac$Gene, wound$Gene))
membership <- stack(sets)
names(membership) <- c("gene", "program")
membership$division <- annotations$division[match(membership$gene, annotations$gene)]
stopifnot(!anyNA(membership$division))
save_table(membership, "shared_program_membership.csv")
composition <- as.data.frame(xtabs(~program + division, membership))
names(composition) <- c("program", "division", "n_genes")
save_table(composition, "matrisome_composition.csv")

# Figure 3A to D: summarize the saved section results and D30 fitted contrasts.
moran <- read_input("pdac_section_moran.csv")
spatial <- read_input("pdac_section_contractile_changes.csv")
stopifnot(nrow(moran) == 12L, nrow(spatial) == 12L,
          setequal(moran$sample_id, spatial$sample_id),
          all(is.finite(moran$Moran_I)), all(is.finite(spatial$adjusted_delta_q75_q25)))
spatial_summary <- data.frame(
  n_sections = nrow(spatial),
  f1_moran_median = median(moran$Moran_I),
  contractile_change_min = min(spatial$adjusted_delta_q75_q25),
  contractile_change_median = median(spatial$adjusted_delta_q75_q25),
  contractile_change_max = max(spatial$adjusted_delta_q75_q25),
  n_positive_changes = sum(spatial$adjusted_delta_q75_q25 > 0)
)
save_table(spatial_summary, "spatial_summary.csv")
d30_summary <- do.call(rbind, lapply(c("Discovery", "CRA001160"), function(cohort) {
  z <- read_input(paste0(cohort, "_D30_contrast_curve.csv"))
  stopifnot(all(is.finite(z$estimate)), all(diff(z$scoreA) > 0))
  data.frame(cohort = cohort, contrast_min = min(z$estimate),
             contrast_max = max(z$estimate), ecm_score_min = min(z$scoreA),
             ecm_score_max = max(z$scoreA))
}))
save_table(d30_summary, "D30_contrast_summary.csv")

# Figure 3E: compare the fixed elastic net and baseline within each external patient.
scores <- read_input("CRA001160_fixed_predictions.csv")
stopifnot(!anyNA(scores), length(unique(scores$patient)) == 16L)
performance <- do.call(rbind, lapply(split(scores, scores$patient), function(z) {
  data.frame(
    patient = z$patient[1],
    n_cells = nrow(z),
    elastic_net_spearman = cor(z$contractile, z$predicted_full, method = "spearman"),
    baseline_spearman = cor(z$contractile, z$predicted_baseline, method = "spearman")
  )
}))
performance <- performance[order(performance$elastic_net_spearman), ]
performance$spearman_improvement <- with(performance, elastic_net_spearman - baseline_spearman)
save_table(performance, "external_performance_by_patient.csv")
prediction_summary <- data.frame(
  n_patients = nrow(performance),
  n_cells = sum(performance$n_cells),
  median_elastic_net = median(performance$elastic_net_spearman),
  median_baseline = median(performance$baseline_spearman),
  n_patients_improved = sum(performance$spearman_improvement > 0)
)
save_table(prediction_summary, "external_prediction_summary.csv")

# Figure 4D and E: use the fixed discovery target sets and saved external contrasts.
targets <- read_input("CRA001160_fixed_TGFB_targets.csv")
patient_cols <- grep("^T[0-9]+$", names(targets), value = TRUE)
stopifnot(length(patient_cols) == 16L, !anyDuplicated(targets[c("ligand", "gene")]))
matrix <- as.matrix(targets[patient_cols])
stopifnot(is.numeric(matrix))
targets$measurable <- rowSums(is.finite(matrix)) > 0L
stopifnot(all(rowSums(is.finite(matrix)) %in% c(0L, length(patient_cols))))
ligand_responses <- lapply(split(targets, targets$ligand), function(z) {
  stopifnot(nrow(z) == 30L, any(z$measurable))
  data.frame(ligand = z$ligand[1], patient = patient_cols,
             n_targets_locked = nrow(z), n_targets_measured = sum(z$measurable),
             mean_log2cpm_difference = colMeans(z[z$measurable, patient_cols, drop = FALSE]))
})
save_table(do.call(rbind, ligand_responses), "external_target_response_by_patient.csv")
target_summary <- do.call(rbind, lapply(ligand_responses, function(z) {
  data.frame(ligand = z$ligand[1], n_targets_locked = z$n_targets_locked[1],
             n_targets_measured = z$n_targets_measured[1], n_patients = nrow(z),
             n_positive = sum(z$mean_log2cpm_difference > 0),
             mean_difference = mean(z$mean_log2cpm_difference))
}))
save_table(target_summary, "external_target_response_summary.csv")
display <- targets[targets$measurable, ]
display <- display[order(-display$regulatory_potential, display$gene, display$ligand), ]
display <- head(display[!duplicated(display$gene), ], 12)
save_table(display[c("gene", "regulatory_potential", patient_cols)], "displayed_12_target_responses.csv")

message("Wrote ", length(list.files(out, pattern = "[.]csv$")), " CSV results to ", out)
