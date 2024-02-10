#' Plot labeling enrichment
#'
#' This function takes as input the output list from the function of "flx.correct_natural_abundance", and plot the labeling enrichment for any specified compound. When the "sample_group" argument is specified, the plot is faceted into panels based on the specified group category.
#'
#' @param listed.data A list of datasets. It is the output from the function of "flx.correct_natural_abundance". The two subsets "Corrected" and "Normalized" datasets are in particular used for plotting under the hood.
#' @param compound A character vector of length one to specify the compound name.
#' @param enrichment_lower_bound A single number to specify the lower boundary of the y-axis in the bar plot. By default, an enrichment range of 0 - 1 is plotted on the y-axis in the bar plot.
#' @param sample_group A dataset of sample information, optional argument. If specified, the plot is faceted into panels based on sample groups. The input dataset contains two essential columns: the specified column "sample" contains all sample names as contained in the "listed.data" argument input, and the column of "group" denotes the category by which the plot is faceted into group-panel. Extra columns can be stored in this dataset for other analysis, but will not be utilized for this function.
#'
#' @return A ggplot.
#' @export
#'
#' @examples

flx.plot_labeling_enrichment = function(listed.data = NULL,
                                       compound = NULL,
                                       enrichment_lower_bound = 0,
                                       sample_group = NULL,
                                       facet.nrow = 1,
                                       orderedSample = NULL) {

  # Check there is listed data input
  if (is.null(listed.data)) {
    stop("Please specify the dataset: the output from function \"flx.correct_natural_abundance\", which is in the format of a list.\n\n")
  }

  # Extract dataset from the list
  d.corrected = listed.data$Corrected
  d.normalized = listed.data$Normalized


  cmpd.all = d.corrected$Compound %>% unique()

  # Check there is compound input
  if (is.null(compound )) {
    stop("A compound name must be specified in the argument of \"compound = ... \" to make a plot.\n ")
  }

  # Check compound is of length of one
  if (length(compound) > 1) {
    stop("Please input only one compound in the argument of \"compound = ... \" when making a plot.\n\n")
  }
  # Check compound matches compound names in the input listed.data
  if (! compound %in% cmpd.all) {
    stop("Compound", " \"", compound, "\"", " is not found in the input dataset.\n" )
  }



  # Convert C-label into ordered factor for visualization purpose
  ordered.C_label = d.normalized $ C_Label %>% unique() %>% rev()
  d.normalized$C_Label = factor(d.normalized$ C_Label, levels = ordered.C_label, ordered = T)
  d.corrected$C_Label = factor(d.corrected$C_Label, levels = ordered.C_label, ordered = T)


  # color setup for C-label: same color assignment rule for all compounds
  colors = c ("grey", "firebrick", "yellow", "lightgreen",  "skyblue", "steelblue", "black")  # Specify the first seven colros (M+0, ...M + 6)
  colors.more = colorRampPalette(brewer.pal(8, "Dark2"))( (nmax <- d.corrected$C_Label %>% as.numeric() %>% max()) - length(colors) + 1)
  colors = c(colors, colors.more)
  names(colors) = 0:nmax %>% factor(ordered = T)


  # Select specified compound and tidy up
  d.corrected.tidy = d.corrected %>% filter(Compound == compound) %>%
    gather(-c(Compound, C_Label), key = sample, value = intensity)

  d.TIC.tidy = d.corrected.tidy %>% group_by(sample, Compound) %>%
    summarise(TIC = sum(intensity, na.rm = T))

  d.normalized.tidy = d.normalized %>%  filter(Compound == compound) %>%
    gather(-c(Compound, C_Label), key = sample, value = enrichment)



  # If sample info dataset is input as argument, then combine with sample info dataset, and plot with groups in separate panel
  if ( !is.null(sample_group) ) {

    # Check there is "group" column
    if (! "group" %in% colnames(sample_group)) {
      stop ("A", ' "group"', " column needs to be specified in the dataset for the argument of \"sample_group\".\n\n"  )
    }

    # Check the sample name match each other
    n1 = d.normalized.tidy$sample %>% unique()
    n2 = sample_group$sample
    if (sum(n1 %in% n2) == length(n1)) { # all sample names in the listed data has a match in the sample information dataset (no problem if on the other hand the "sample info" dataset has extra sample names not found in the "listed dataset")
      d.normalized.tidy = d.normalized.tidy %>% left_join(sample_group, by = c("sample"))
      d.TIC.tidy = d.TIC.tidy %>% left_join(sample_group, by = "sample")

    } else {
      notFoundInSampleInfo = n1[! n1 %in% n2] # Check listed data samples names that have no match in the "sample info"; but the sample info dataset may have extra (and unused) sample names and this is not problematic!
      stop("The following sample", ifelse(length(notFoundInSampleInfo) == 1, " is", "s are"),
           " not found in the dataset for the argument of \"sample_group\":\n", str_c(notFoundInSampleInfo, "\n "), "\n" )
    }
  }



  # max TIC intens.
  TIC.max = (d.TIC.tidy)$TIC %>% max(na.rm = T)
  # Number of color
  C_label.i = d.normalized.tidy$C_Label %>% unique() %>% as.character()



  # Define plotting function
  theme_set(theme_bw() +
              theme(#axis.text.x = element_text(angle = 45, hjust = 1),
                strip.background = element_blank(),
                strip.text = element_text(size = 14, face = "bold"),
                panel.grid = element_blank(),
                panel.border = element_rect(colour = "black", size = 1),
                axis.text = element_text(colour = "black", size = 9),
                title = element_text(face = "bold", hjust = .5, size = 12),
                plot.title = element_text(hjust = .5, size = 16),
                legend.title = element_blank()) )


  # specify sample order (on the x-axis)
  if( is.null(orderedSample) == F) {
    d.normalized.tidy$sample = factor(d.normalized.tidy$sample, levels = orderedSample, ordered = T)
    d.TIC.tidy$sample = factor(d.TIC.tidy$sample, levels = orderedSample, ordered = T)
  }


  # Label fraction bar plot
  plt.bar = d.normalized.tidy %>%
    filter(Compound == compound) %>%
    ggplot(aes(x = sample)) +
    geom_bar(aes(y = enrichment, fill = C_Label, color = C_Label),
             stat = "identity", alpha = .8) +


    scale_color_manual(values = colors [C_label.i %>% as.character()] ) +
    scale_fill_manual(values = colors [C_label.i %>% as.character()] ) +
    labs(x = " ", y = "Labelling fraction\n", title = compound)  + # \n increase y axis - text gap
    scale_x_discrete(expand = c(0, 0)) +
    coord_cartesian(ylim = c(enrichment_lower_bound, 1)) +
    scale_y_continuous(expand = c(0, 0),
                       # Need this right side axis as place holder
                       sec.axis = sec_axis(~.*TIC.max, name = "Total isotopologues ion counts\n")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, color = "white"),
          axis.text.y.right = element_text(color = NA),  # turn off TIC axis text
          axis.title.y.right = element_text(colour = NA),
          axis.ticks.y.right = element_blank() ) # turn off TIC ticks



  # TIC dot plot
  plt.TIC = d.TIC.tidy %>%
    ggplot(aes(x = sample)) +

    # set bar plot (transparent) so as to keep the legend as a place holder to allow overlay with real bar plot
    geom_bar(aes(y = TIC),
             stat = "identity", alpha = 0, color = NA, position = "fill") +



    # scale_color_manual(values = colors) +
    # scale_fill_manual(values = colors) +
    labs(x = " ", y = "Labelling fraction\n", title = compound)  + # \n increase y axis - text gap
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0),
                       sec.axis = sec_axis(~.*TIC.max, name = "Total isotopologues ion counts\n")) +

    theme(panel.background = element_blank(),

          # turn on bar plot as place holder for legend, but not show legend here
          # so as to show legend of the bar plot
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.text = element_blank(),

          plot.background = element_blank(),
          axis.ticks.y.left = element_blank(), # turn off label fraction axis ticks
          axis.text.x = element_text(angle = 60, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.text.y.left = element_text(colour = NA), # turn off label fraction axis text
          axis.title.y.left = element_text(colour = NA) # turn off label fraction axis title
    )  +

    # TIC counts
    geom_point(aes(y = TIC  / TIC.max), color = "black", alpha = 1, shape = 18, size = 3)
    # geom_line(aes(y = TIC  / TIC.max, group = 1), color = "black", alpha = .5) # not draw line


  # if sample info is input, then plot by faceting groups
  if (!is.null(sample_group)){
    plt.TIC = plt.TIC + facet_wrap( ~ group, scales = "free_x", nrow = facet.nrow) + theme(panel.spacing = unit(1, "lines"))
    plt.bar = plt.bar + facet_wrap( ~ group, scales = "free_x", nrow = facet.nrow) + theme(panel.spacing = unit(1, "lines"))
  }

  # Overlay and align up the two plots
  aligned_plots <- align_plots(plt.bar, plt.TIC, align="hv", axis="tblr")
  ggdraw(aligned_plots[[1]]) + draw_plot(aligned_plots[[2]])

}




# test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>- test -<>-
# library(readxl)
# library(cowplot)
# library(RColorBrewer)
# library(tidyverse)
# library(vivoFlux)
# library(accucor)
#
# # read labeling dataset
# d1 = read_excel("/Users/boyuan/Desktop/Harvard/Research/vivoFlux/example_data_parasite_labeling.xlsx", sheet = "i_see_fantacy")
#
# # read sample information dataset
# sample_group = read_excel("/Users/boyuan/Desktop/Harvard/Research/vivoFlux/example_data_parasite_labeling.xlsx", sheet = "sample_info")
#
# # background subtraction
# d1.backgroundSubtracted = flx.subtract_background(d1, blankSamples = c("blank_r1", "blank_r2", "blank_r3"))
#
# # natural abundance correction
# list.corrected = d1.backgroundSubtracted %>% flx.correct_natural_abundance()
#
# flx.plot_labeling_enrichment(listed.data =  list.corrected, compound = c("Alanine"),
#                             enrichment_lower_bound = 0)
#
# flx.plot_labeling_enrichment(listed.data =  list.corrected, compound = c("Alanine"),
#                             enrichment_lower_bound = 0, sample_group = d.sample)
#
# flx.plot_labeling_enrichment(listed.data =  list.corrected, compound = c("Alanine"),
#                             enrichment_lower_bound = 0, sample_group = d.sample, facet.nrow = 2)


# # Check error report -<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>--<>-
#
# # No dataset input
# flx.plot_labeling_enrichment()
# flx.plot_labeling_enrichment(compound = "Alanine")
#
# # No compound input
# flx.plot_labeling_enrichment(listed.data = list.corrected)
#
# # input compound not found in the dataset
# flx.plot_labeling_enrichment(listed.data =  list.corrected,
#                             compound = c("XXXXX"))
#
# # more than one compound in the compound arguement
# flx.plot_labeling_enrichment(listed.data = list.corrected,
#                             compound = c("Alanine", "df"),
#                             enrichment_lower_bound = 0)
#
# # sample information dataset has no group column
# flx.plot_labeling_enrichment(listed.data = list.corrected,
#                             compound = c("Alanine"),
#                             sample_group = d.sample %>% select(-group))
#
#
#
# # Not all listed dataset sample names are specified with a group category in the sample_group argument
# # 1) one not matched
# flx.plot_labeling_enrichment(listed.data = list.corrected,
#                             compound = c("Alanine"),
#                             sample_group = d.sample[-c(1), ] )
#
# # 2) Multiple not matched (check the plural form in the error message)
# flx.plot_labeling_enrichment(listed.data = list.corrected,
#                             compound = c("Alanine"),
#                             sample_group = d.sample[-c(1:3), ] )
