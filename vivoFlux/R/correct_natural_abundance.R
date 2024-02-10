#' Correction of isotopic natural abundance.
#'
#' This is a wrapper function of "natural_abundance_correction" in the accucor package. It functions to correct the 13-carbon natural abundance to unveil the net 13-C enrichment due to experimental 13-C tracer labeling.
#'
#' @param dataset The input may be a dataset after background subtraction, or without subtraction. In the latter case, the blank samples, if any in the dataset input, will be processed just as a normal sample.
#' @param resolution # Mass spectrometry scan resolution. Default value set to 120000
#' @param resolution_defined_at # Molecular mass (Dalton) at which the mass spectrometry resolution is defined. Default value set to 200 Dalton.
#'
#' @return A list of four data frames: the original data set, corrected intensity, corrected and normalized intensity, and the pool.
#' @export
#'
#' @examples



flx.correct_natural_abundance = function(dataset, resolution = 120000, resolution_defined_at = 200){
  natural_abundance_correction(dataset, resolution = resolution, resolution_defined_at = resolution_defined_at)
}



# test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>-

# library(vivoFlux)
# library(accucor)
# d1 = read_excel("/Users/boyuan/Desktop/Harvard/Research/vivoFlux/example_data_parasite_labeling.xlsx", sheet = "i_see_fantacy")
# d1.backgroundSubtracted = flx.subtract_background(d1, blankSamples = c("blank_r1", "blank_r2", "blank_r3"))
#
# list.corrected = d1.backgroundSubtracted %>%
#   flx.correct_natural_abundance()

# d.normalized = list.corrected$Normalized
# d.corrected = list.corrected$Corrected
# d.original = list.corrected$Original


