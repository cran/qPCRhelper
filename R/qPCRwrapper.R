#' @title qPCRhelper
#'
#' @description Takes in user input table of per gene raw Ct values and computes delta Ct, delta delta Ct
#'     and log2(2^(-1*delta delta Ct)). It then tests for significance of normalized Ct values (delta Ct)
#'     across groups per gene using t.test. Plots log2(2^(-1*delta delta Ct)) across groups per gene,
#'     and displays respective p-values.
#' @param data.dir a path, file path of raw qPCR Ct values per gene per sample, expects "Sample",
#'     "Group", and "Gene" columns.
#' @param plot.title a string, title of the plot, optional
#' @param ref.gene a string, expects a character string matching the one of the column names of input
#'     table to be used for delta Ct computation.
#' @param ref.group a string, expects a character string matching the Group column of input table to
#'     be used for delta delta Ct computation.
#' @param plot.ref.group a string, expects a character string matching the Group column of input table
#'     to be used for setting reference in plotting
#' @param plot.nrow numeric, optional, number of rows for plotting
#' @export
#' @return a dataframe of delta delta Ct computation across genes, also returns a plot with test of
#'    significance using t.test
#' @importFrom dplyr "%>%"
qPCRhelper <- function(data.dir = NULL, ref.gene = NULL, ref.group = NULL,
                 plot.ref.group = NULL, plot.nrow = 1,
                 plot.title = NULL) {
  Sample <- p <- NULL

  if (requireNamespace("magrittr", quietly = TRUE)) {
    magrittr::"%>%"
  }

  # Read the user data
  mydata <- read.table(file = data.dir, header = TRUE)

  # Find unique groups
  unique_groups <- unique(mydata[, 2])

  # Add blank rows after each unique group
  new_data <- data.frame(matrix(nrow = 0, ncol = ncol(mydata)))
  for (group in unique_groups) {
    group_indices <- which(mydata[, 2] == group)
    last_index <- group_indices[length(group_indices)]
    group_data <- mydata[group_indices,]
    blank_row <- c(paste0("Average_", group), group, rep(NA, ncol(mydata) - 2))
    blank_row[2] <- group
    new_group_data <- rbind(group_data, blank_row)
    new_data <- rbind(new_data, new_group_data)
  }

  # Compute column averages per 'Group' and add them to blank rows
  for (i in 3:ncol(mydata)) {
    for (group in unique_groups) {
      group_indices <- which(mydata[, 2] == group)
      group_data <- mydata[group_indices, i]
      avg <- mean(group_data, na.rm = TRUE)
      blank_row_index <- which(is.na(new_data[, i]) & new_data[, 2] == group)
      new_data[blank_row_index, i] <- avg
    }
  }

  # Subtract ref.gene value from other columns and store in new columns
  if (!is.null(ref.gene) & ref.gene %in% colnames(mydata)) {
    ref_gene_col <- which(colnames(mydata) == ref.gene)
    for (i in 3:ncol(new_data)) {
      if (colnames(new_data)[i] != ref.gene) {
        new_colname <- paste0("dCT_", colnames(new_data)[i])
        new_data[, new_colname] <- as.numeric(new_data[, i]) - as.numeric(new_data[, ref_gene_col])
      }
    }
  }

  # Compute ddCT based on dCT values and average dCT of ref.group
  if (!is.null(ref.group) & ref.group %in% mydata$Group) {
    ref_group_avg <- subset(new_data, Sample == paste0("Average_", ref.group))
    for (i in 3:ncol(new_data)) {
      if (startsWith(colnames(new_data)[i], "dCT_")) {
        ddCT_colname <- paste0("ddCT_", colnames(new_data)[i])
        new_data[, ddCT_colname] <- as.numeric(new_data[, i]) - as.numeric(ref_group_avg[, colnames(new_data)[i]])
      }
    }
  }

  # Compute 2^(-ddCT) and log2(2^(-ddCT)) based on ddCT values for each gene
  for (i in 3:ncol(new_data)) {
    if (startsWith(colnames(new_data)[i], "ddCT_")) {
      relExp_colname <- paste0("log2RelExp_", colnames(new_data)[i])
      new_data[, relExp_colname] <- log2(2^(-1*(as.numeric(new_data[, i]))))
    }
  }

  remove_avg <- subset(new_data, !startsWith(new_data$Sample, "Average_"))
  remove_avg <- remove_avg[,!names(remove_avg) %in% c(ref.gene)]
  #print(remove_avg)

  # Find unique gene names
  gene_names <- unique(gsub("log2RelExp_ddCT_dCT_", "", grep("log2RelExp_ddCT_dCT_", colnames(remove_avg), value = TRUE), fixed = TRUE))

  # Initialize empty data frame to store the final result
  final_data <- data.frame(Sample = character(), Group = character(), Gene = character(), dCT = numeric(), log2RelExp = numeric(), stringsAsFactors = FALSE)

  # Loop over each unique gene name, extract the dCT and log2RelExp columns, and append them to the final data frame
  for (gene_name in gene_names) {
    dCT_colname <- paste0("dCT_", gene_name)
    log2RelExp_colname <- paste0("log2RelExp_ddCT_dCT_", gene_name)

    # Extract the columns with the current gene name
    current_data <- subset(remove_avg, select = c("Sample", "Group", dCT_colname, log2RelExp_colname))

    # Create a new column called 'Gene' and set it to the current gene name
    current_data$Gene <- rep(gene_name, nrow(current_data))

    # Rename the dCT and log2RelExp columns to 'dCT' and 'log2RelExp'
    colnames(current_data)[3] <- "dCT"
    colnames(current_data)[4] <- "log2RelExp"

    # Append the current data to the final data frame
    final_data <- rbind(final_data, current_data)
  }

  # Order the final data frame by Sample and Group
  final_data <- final_data[order(final_data$Gene, final_data$Group, final_data$Sample), ]


  #Plot
  anno_df <- ggpubr::compare_means(dCT ~ Group, method = "t.test", data = final_data,
                                   group.by = "Gene", ref.group = plot.ref.group) %>%
    rstatix::add_xy_position(x = "Group", data=final_data, formula=dCT ~ Group) %>%
    dplyr::mutate(myformatted.p = paste0("p = ", format.pval(p, digits = 1)))

  #Create the fancy plot with formatted p-values
  y.position <- (max(final_data$log2RelExp) + 0.5)
  ggpubr::ggboxplot(final_data, x = "Group", y = "log2RelExp",
                         palette = "RdBu", add = "jitter", fill = "Group", bxp.errorbar = T,
                         ylab = "log2(Relative expression)", xlab = F, title = plot.title) +
    ggplot2::facet_wrap(~Gene, nrow = plot.nrow) +
    ggpubr::stat_pvalue_manual(anno_df,
                               label = "myformatted.p",
                               step.increase = 0.09,
                               step.group.by = "Gene",
                               y.position = y.position,
                               bracket.shorten = 0.1,
                               bracket.nudge.y = 0.5) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, face = "bold"),
          plot.title = ggplot2::element_text(face = "bold"),
          strip.text.x = ggplot2::element_text(face = "bold.italic", size = 12),
          panel.background = ggplot2::element_rect(fill = "#DCD6D0"))

  # Return the final data frame
  return(final_data)
}

