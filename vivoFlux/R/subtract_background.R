#' Background subtraction.
#'
#' Background level of each isotopomer of each metabolite is averaged from designated blank sample(s), which are then subtracted from the associated signal intensity in biological samples.
#' @param dataset Dataframe. The first three columns are reserved to "isotopeLbel", "compound", and "formula", respectively. The following columns contain the sample signal intensities as numeric values; non-numeric values, such as characters or symbols, must be avoided.
#' @param blankSamples A vector of blank samples' names.It indicates the blank samples whose average signal should be subtracted from real sample signal.
#' @return A data frame, with background subtracted. The blank samples (columns) are removed from the output dataset.
#' @export
#'
#' @examples


flx.subtract_background = function(dataset, # first three columns must be isotopeLbel, compound, and formula, followed by sample peak area in numeric type
                                  blankSamples = NA){ # need to manually specify blank sample names


  # Check the identity of the thirst three columns (order does not matter)
  if (sum(names(dataset)[1:3] %in% c("isotopeLabel", "compound", "formula")) != 3) {
    stop ("The first three columns should be \"isotopeLabel\", \"compound\", and \"formula.\" ") }

  # Check the remaining columns are numeric numbers
  if (class(dataset[, -c(1:3)] %>% sum()) != "numeric") {
    stop ("Data of peak intensity should start from the 4th column, and must be numeric values (characters not allowed).")
  }

  if (is.na(blankSamples[[1]])) { # the argument passed to "is.na" is of length one, thus take the first element of vector blankSamples
    stop("The blank (background) sample names must be specified. Use argument \"blankSamples = ...\" to specify the blank samples. If no blank samples are used, you can skip this step.") }



  # when specified blank sample name in the argument is not in the input dataset, report error
  logic.notIn = blankSamples %in% colnames(dataset)
  n.logic.notIn = sum (logic.notIn == F)
  names.notIn = blankSamples[!logic.notIn] %>% toString()


  if (! is.na(blankSamples[1]) & (n.logic.notIn >=1)   ) {
    stop (paste0("The designated blank (background) sample(s) -- ", paste(names.notIn), ", -- ",
                ifelse(n.logic.notIn > 1, "are", "is"),  " not found in the input dataset.") )
  }


  # subset, with columns containing background (blank) data removed
  d.blankRemoved = dataset %>% select(-contains(blankSamples))
  # subset of background data
  d.blank = dataset %>% select(contains(blankSamples))
  # calculate average background signal
  blk = rowMeans(d.blank)

  # subtract the background
  # the first three columns must be isotopeLbel, compound, and formula! Function assumes all the other following columns are samples data in numeric type
  d.blankRemoved[, -c(1:3)] = apply( d.blankRemoved[, -c(1:3)], 2, function(x) {x - blk} )

  # Force negative values (due to background subtraction) to zero (All numbers plus its absolute value then divided by 2)
  d.blankRemoved[, -c(1:3)] =  ( d.blankRemoved[, -c(1:3)] + abs( d.blankRemoved[, -c(1:3)] ) )/2

  # output background-subtracted dataset
  d.blankRemoved %>% return()

}




# test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>-

# library(readxl)
# library(tidyverse)
# library(accucor)
# library(rebus)
# library(RColorBrewer)
# library(cowplot)
#
#
# d1 = read_excel("/Users/boyuan/Desktop/Harvard/Research/vivoFlux/example_data_parasite_labeling.xlsx", sheet = "i_see_fantacy")
# d1.backgroundSubtracted = flx.subtract_background(d1, blankSamples = c("blank_r1", "blank_r2", "blank_r3"))
#
# # testing error report
# d1 %>% flx.subtract_background(blankSamples =  c("blank_r1", "blank_r2", "blank_r4")) # some of the specified blank names not in the input dataset
# d1 %>% flx.subtract_background(blankSamples =  c("blank_r1", "blank_r2", "blank_r4", "BLANK_r1")) # some of the specified blank names not in the input dataset
# d1 %>% flx.subtract_background(blankSamples =  c(" ")) # some of the specified blank names not in the input dataset
# d1 %>% flx.subtract_background(blankSamples =  c("BLK1", NA)) # testing NA
# d1 %>% flx.subtract_background() # no blank argument

