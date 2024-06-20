#' Genelist
#'
#' @return a named list where indexes are gene names used in the original publication, and elements represents a vector of valid gene names for that gene.
#' @param published_only forces return of a single vector containing the originally published rig gene names only (no alts)
#' @export
#'
#' @examples
#' rig_genes()
rig_genes <- function(published_only = FALSE){
 df <- read.csv(system.file("RIG_expression_profile_genes.txt", package = "rig"), sep = "\t")
 ls <- strsplit(df[["ValidGeneNames"]], split = ",")
 names(ls) <- df[["Gene"]]

 if(published_only)
    return(df[["Gene"]])

 return(ls)
}

#' Example Dataset
#'
#' @return example dataset that could be analysed with rig_predict
#' @export
#'
#' @examples
#' rig_example_data()
rig_example_data <- function(){
  # Simulated dataset
  data <- data.frame(
    sample_id = c('sample1', 'sample1', 'sample1', 'sample1', 'sample1', 'sample1', 'sample1', 'sample1', 'sample1', 'sample1',
                  'sample2', 'sample2', 'sample2', 'sample2', 'sample2', 'sample2', 'sample2', 'sample2', 'sample2', 'sample2',
                  'sample3', 'sample3', 'sample3', 'sample3', 'sample3', 'sample3', 'sample3', 'sample3', 'sample3', 'sample3',
                  'sample4', 'sample4', 'sample4', 'sample4', 'sample4', 'sample4', 'sample4', 'sample4', 'sample4', 'sample4',
                  'sample5', 'sample5', 'sample5', 'sample5', 'sample5', 'sample5', 'sample5', 'sample5', 'sample5', 'sample5',
                  'sample6', 'sample6', 'sample6', 'sample6', 'sample6', 'sample6', 'sample6', 'sample6', 'sample6', 'sample6',
                  'sample7', 'sample7', 'sample7', 'sample7', 'sample7', 'sample7', 'sample7', 'sample7', 'sample7', 'sample7',
                  'sample8', 'sample8', 'sample8', 'sample8', 'sample8', 'sample8', 'sample8', 'sample8', 'sample8', 'sample8',
                  'sample9', 'sample9', 'sample9', 'sample9', 'sample9', 'sample9', 'sample9', 'sample9', 'sample9', 'sample9',
                  'sample10', 'sample10', 'sample10', 'sample10', 'sample10', 'sample10', 'sample10', 'sample10', 'sample10', 'sample10'),
    gene = rep(c("MYT1", "PCSK2", "KCNJ10", "PLP1", "OLIG2", "CROC4", "COL11A1", "PMP2", "PCDH15", "CA10"), 10),
    expression = c(16.49, 12.33, 15.23, 18.00, 21, 20, 12.78, 15, 12.88, 19,
                   50, 12, 24, 40, 21, 15.12, 24, 30, 13.56, 12.99,
                   12.93, 12.34, 13.36, 12.98, 11.76, 11.23, 12.69, 12.41, 11.84, 12.39,
                   8.80,  8.39,  9.97, 10.74,  9.23,  9.83,  9.94, 10.09,  9.57,  9.37,
                   11.48, 10.92, 11.19,  9.40, 10.80, 10.53, 11.59, 10.22, 11.42, 11.00,
                   11.67,  9.55, 11.24, 11.66, 11.37,  9.19,  9.10,  9.21,  8.47, 10.10,
                   9.99,  9.40,  9.77,  9.29,  8.90,  9.43,  9.56,  9.89,  8.87,  9.29,
                   10.87, 10.57, 11.32, 11.28, 10.77, 10.92, 11.26, 11.09, 10.69, 11.23,
                   9.15,  9.20,  9.49,  9.57,  9.68,  9.30,  9.21,  9.79,  9.16,  9.29,
                   8.79,  8.49,  9.02,  8.89,  8.61,  8.47,  8.37,  8.67,  8.39,  8.67)
  )

  return(data)

}


#' Identify Radiation Induced Glioma Signatures
#'
#' This function identifies Radiation Induced Glioma (RIG) signatures by looking for overexpression of RIG genes within a given cohort.
#' Predictions will only be accurate if RIGs make up <50% of the samples supplied.
#'
#' @param data A data.frame representing gene expression results for either a single sample or many.
#' It must contain at least 3 columns: one for gene names, one for gene expression measures (ideally TPM), and one for sample identifiers.
#' @param col_samples The name of the column containing sample identifiers.
#' @param col_genes The name of the column containing gene names.
#' @param col_expression The name of the column containing gene expression measures (ideally TPM).
#' @param signature A named list where each element is a character vector of valid gene names for a single gene, and names are the common names for those genes. Defaults to `rig_genes()`.
#' @param classification_threshold The proportion of RIG genes that must be overexpressed before the sample is classified as a RIG
#' @param expression_threshold The robust Zscore above which a rig gene is considered overexpressed.
#' @param require_all_genes_found A boolean indicating whether an error should be thrown if there are genes in the signature gene set that cannot be found in the expression dataset. Defaults to FALSE.
#' @return A data.frame with sample-level summaries including the number and proportion of overexpressed RIG genes, median z-scores, and a prediction of Radiation Induced Glioma status.
#' @export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   sample_id = c('sample1', 'sample1', 'sample2', 'sample2'),
#'   gene = c('geneA', 'geneB', 'geneA', 'geneC'),
#'   expression = c(5, 3, 4, 7)
#' )
#' rig_identification <- rig_predict(data, 'sample_id', 'gene', 'expression')
#' print(rig_identification)
#' }
rig_predict <- function(data, col_samples, col_genes, col_expression,
                        signature = rig_genes(),
                        threshold = 0.7,
                        expression_threshold = 1.3,
                        require_all_genes_found = FALSE
                        ) {

  # Validate input types and contents
  assertions::assert_dataframe(data)
  assertions::assert_string(col_samples)
  assertions::assert_string(col_genes)
  assertions::assert_string(col_expression)
  assertions::assert_names_include(data, col_genes)
  assertions::assert_names_include(data, col_expression)
  assertions::assert_names_include(data, col_samples)
  assertions::assert_list(signature)
  assertions::assert_no_duplicates(unlist(signature), msg = "The signature gene set contains duplicate gene names across different genes, which makes it impossible to automatically find the best gene names to use.")

  # Tailor gene names to match the signature gene set
  sig_genes <- tailor_gene_names_to_signature(
    dataset_genes = data[[col_genes]],
    signature_gene_set = signature,
    require_all_genes_found = require_all_genes_found
  )

  n_signature_genes_found <- length(sig_genes)

  # Filter data to include only the genes found in the signature gene set
  data <- data[data[[col_genes]] %in% sig_genes, , drop = FALSE]

  # Compute robust Z-scores for each gene in the signature
  data[["zscore"]] <- ave(data[[col_expression]], data[[col_genes]], FUN = distrikit::compute_zscore_robust)

  # Compute per-sample summaries

  # Number of RIG genes with Z-score > expression_threshold
  number_of_rig_genes_with_zscore_above_threshold <- tapply(
    data[["zscore"]],
    data[[col_samples]],
    function(z) { sum(z > expression_threshold, na.rm = TRUE) }
  )

  # Proportion of RIG genes with Z-score > expression_threshold
  proportion_of_rig_genes_with_zscore_above_threshold <- number_of_rig_genes_with_zscore_above_threshold / n_signature_genes_found

  # Median Z-score of all RIG genes per sample
  median_zscore_per_sample <- tapply(
    data[["zscore"]],
    data[[col_samples]],
    function(z) { median(z, na.rm = TRUE) }
  )

  # Compile per-sample results into a summary data.frame
  sample_summary <- data.frame(
    sample = names(number_of_rig_genes_with_zscore_above_threshold),
    n_overexpressed_rig_genes = number_of_rig_genes_with_zscore_above_threshold,
    prop_overexpressed_rig_genes = proportion_of_rig_genes_with_zscore_above_threshold,
    median_rig_zscore = median_zscore_per_sample
  )

  # Identify RIG based on the proportion of overexpressed RIG genes
  sample_summary[["predicted_to_be_rig"]] <- sample_summary[["prop_overexpressed_rig_genes"]] >= threshold


  return(sample_summary)
}


#' Tailor Gene Names from a Signature Gene Set to Your Dataset
#'
#' This function takes a list of gene names from an expression dataset and tailors them to match a signature gene set.
#' The function returns the gene names from the dataset that correspond to the signature gene set.
#'
#' @param dataset_genes A character vector of gene names from your expression dataset.
#' @param signature_gene_set A named list where each element is a character vector of valid gene names for a single gene, and names are the common names for those genes.
#' @param require_all_genes_found A boolean indicating whether an error should be thrown if there are genes in the signature gene set that cannot be found in the expression dataset. Defaults to TRUE.
#' @return A character vector of gene names that are present in 'dataset_genes' and correspond to the signature_gene_set.
#' @export
#'
tailor_gene_names_to_signature <- function(dataset_genes, signature_gene_set = rig_genes(), require_all_genes_found = TRUE, verbose = TRUE){

  # Validate input types
  assertions::assert_character(dataset_genes)
  assertions::assert_list(signature_gene_set)
  assertions::assert_no_duplicates(unlist(signature_gene_set), msg = "The signature gene set contains duplicate gene names across different genes, which makes it impossible to automatically find the best gene names to use.")

  # Remove duplicate genes from the input list
  dataset_genes <- unique(dataset_genes)

  # Find the best matching gene names from the expression dataset for each gene in the signature gene set
  matched_genes <- lapply(signature_gene_set, function(signature_genes) {
    dataset_genes[match(signature_genes, dataset_genes)[1]]
  })
  matched_genes = unlist(matched_genes)

  # Identify signature genes that are not present in the expression dataset
  missing_signature_genes <- unique(names(signature_gene_set)[sapply(matched_genes, function(g) all(is.na(g)))])
  n_missing_genes = length(na.omit(missing_signature_genes))
  n_total_sig_genes = length(signature_gene_set)

  # If any signature genes are not found and the requirement is to find all signature genes, throw an error
  if(n_missing_genes > 0 & require_all_genes_found) {
    cli::cli_abort("Failed to find any of the following genes (or any of their alternative names) in your {.strong dataset_genes} vector: [{missing_signature_genes}]")
  }

  if(n_missing_genes > 0 & !require_all_genes_found & verbose) {
    cli::cli_alert_info("Successfully tailored {n_missing_genes}/{n_total_sig_genes} to your dataset")
  }

  # Return the matched gene names
  return(na.omit(matched_genes))
}


