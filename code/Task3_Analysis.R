
# BVG7003 - ASSIGNMENT 1
# Analysis of REM16 and FT1 Expression in Cannabis sativa

# Run from the repository root using:
# source("code/Task3_Analysis.R")


# 1. LOAD REQUIRED PACKAGES

library(data.table)
library(ggplot2)


# 2. IMPORT DATASET

# Define relative file path (relative to the repository root)
data_path <- "data/2_Data_RNASeq_Cannabis_Sex.csv"

# Locate the repository root.
# If the script is run from inside code/ (e.g., after opening the
# script directly or "Set Working Directory > To Source File Location"),
# move up one level so all relative paths (data/, figures/) work.
if (!file.exists(data_path) && file.exists(file.path("..", data_path))) {
  setwd("..")
  message("Working directory changed to repository root: ", getwd())
}

# Stop with a clear message if the dataset still cannot be found
if (!file.exists(data_path)) {
  stop(
    "Input file not found: ", data_path,
    "\nCurrent working directory: ", getwd(),
    "\nSet the working directory to the repository root ",
    "(the folder containing code/, data/ and README.md), ",
    "then run: source(\"code/Task3_Analysis.R\")",
    call. = FALSE
  )
}

# Import normalized RNA-seq dataset
data <- fread(data_path, check.names = FALSE)

# Display dataset dimensions
print(dim(data))

# Note: the assignment describes >34,000 genes, but the supplied
# file contains 20,500 gene/feature rows (see README, Section 2).
# No expression-based gene filtering is performed in this script.
cat("Gene/feature rows in file:", nrow(data), "\n")
stopifnot(!anyDuplicated(data$V1))


# 3. EXAMINE DATASET STRUCTURE

# Check column names
print(names(data)[1:10])

# Display the first 3 rows and 6 columns
print(head(data[, 1:6, with = FALSE], 3))


# 4. DEFINE TARGET GENES

# REM16: LOC115699937
# FT1:   LOC115696989

genes <- c(
  "LOC115699937",
  "LOC115696989"
)


# 5. DATA QUALITY CHECK - GENE IDS

# Initial attempt using exact Gene ID matching
initial_selected <- data[V1 %in% genes]

# Display initially identified genes
print(initial_selected[, .(V1)])

# Identify unmatched target genes
print(setdiff(genes, initial_selected$V1))


# DATASET ISSUE 1:
# REM16 contains an additional "Cannabis_" prefix.
#
# Expected: LOC115699937
# Actual:   Cannabis_LOC115699937
#
# Exact matching therefore fails to identify REM16.

# Standardize Gene IDs.
# Remove leading/trailing whitespace and the
# "Cannabis_" prefix when present at the beginning.
# Preserve the original V1 column.

data[, Gene_ID := sub("^Cannabis_", "", trimws(V1))]


# 6. EXTRACT REM16 AND FT1

# Extract the two target genes
selected <- data[Gene_ID %in% genes]

# Display original and standardized IDs
print(selected[, .(V1, Gene_ID)])

# Check the number of selected rows
print(nrow(selected))

# Verify that both genes were identified exactly once
stopifnot(identical(
  sort(selected$Gene_ID),
  sort(genes)
))


# 7. CHECK SAMPLE SEX DISTRIBUTION

# Identify sample columns.
# Exclude the original and standardized Gene ID columns.

sample_names <- setdiff(
  names(data),
  c("V1", "Gene_ID")
)

# Extract sex information from sample names
sex <- sub(".*_", "", sample_names)

# Display sex distribution
print(table(sex))

# Check total number of samples
print(length(sample_names))

# Verify sample information
stopifnot(length(sample_names) == 138)
stopifnot(all(sex %in% c("XX", "XY")))
stopifnot(!anyDuplicated(sample_names))

# Count missing values across the entire expression matrix
total_na <- sum(is.na(data[, ..sample_names]))
cat("Missing values in full expression matrix:", total_na, "\n")


# 8. PREPARE EXPRESSION DATA FOR PLOTTING

# Verify that both target genes are available
stopifnot(nrow(selected) == 2)
stopifnot(setequal(selected$Gene_ID, genes))

# Convert Wide Format into Long Format
long_data <- melt(
  selected,
  id.vars = c("V1", "Gene_ID"),
  measure.vars = sample_names,
  variable.name = "Sample_ID",
  value.name = "Expression",
  variable.factor = FALSE
)

# Assign readable gene names
long_data[, Gene := fifelse(
  Gene_ID == "LOC115699937",
  "REM16",
  "FT1"
)]

# Extract sex labels from sample names
long_data[, Sex := sub(
  ".*_",
  "",
  Sample_ID
)]

# Verify that expression values are numeric
stopifnot(is.numeric(long_data$Expression))

# Display first observations
print(head(
  long_data[, .(
    Gene,
    Sample_ID,
    Sex,
    Expression
  )]
))

# Check observations per gene and sex
print(table(long_data$Gene, long_data$Sex))

# Check total number of observations
print(nrow(long_data))

stopifnot(nrow(long_data) == 276)


# 9. DATA QUALITY CHECK - MISSING VALUES

# DATASET ISSUE 2:
# One missing FT1 expression value was identified
# in sample 011_XX.
#
# Missing expression values are excluded from
# plotting rather than replaced with estimated values.

# Identify missing observations
print(long_data[
  is.na(Expression),
  .(Gene, Sample_ID, Sex, Expression)
])

# Count missing observations for each gene
print(long_data[, .(
  Missing = sum(is.na(Expression))
), by = Gene])

# Count total valid expression values
print(sum(!is.na(long_data$Expression)))


# Prepare dataset for plotting
# Keep original long_data unchanged

plot_data <- copy(long_data[
  !is.na(Expression)
])

# Convert sex labels into readable factor levels
plot_data[, Sex := factor(
  Sex,
  levels = c("XX", "XY"),
  labels = c("Female (XX)", "Male (XY)")
)]

# Set consistent gene order
plot_data[, Gene := factor(
  Gene,
  levels = c("REM16", "FT1")
)]

# Verify valid sample numbers
print(table(plot_data$Gene, plot_data$Sex))
stopifnot(nrow(plot_data) == 275)


# 10. DESCRIPTIVE STATISTICS AND STATISTICAL TESTS

# Median, IQR and sample size per gene and sex
summary_stats <- plot_data[, .(
  n      = .N,
  Median = round(median(Expression), 2),
  Q1     = round(quantile(Expression, 0.25), 2),
  Q3     = round(quantile(Expression, 0.75), 2),
  Min    = round(min(Expression), 2),
  Max    = round(max(Expression), 2)
), by = .(Gene, Sex)][order(Gene, Sex)]

print(summary_stats)

# Wilcoxon rank-sum test (non-parametric; no normality assumption)
for (g in levels(plot_data$Gene)) {
  cat("\nWilcoxon rank-sum test:", g, "\n")
  print(wilcox.test(
    Expression ~ Sex,
    data = plot_data[Gene == g]
  ))
}

# REM16 group separation: do female and male values overlap?
rem16_female_min <- plot_data[Gene == "REM16" & Sex == "Female (XX)", min(Expression)]
rem16_male_max   <- plot_data[Gene == "REM16" & Sex == "Male (XY)",   max(Expression)]
cat("\nREM16 lowest female value:", round(rem16_female_min, 2),
    "| highest male value:", round(rem16_male_max, 2),
    "| complete separation:", rem16_female_min > rem16_male_max, "\n")

# FT1 low-expression check: values sharing the most common (floor) value
ft1_values <- round(plot_data[Gene == "FT1", Expression], 3)
ft1_mode   <- as.numeric(names(which.max(table(ft1_values))))
cat("FT1 values equal to", ft1_mode, ":", sum(ft1_values == ft1_mode),
    "of", length(ft1_values), "\n")



# COMMON FIGURE STYLE


# Consistent colors across all three figures

sex_colors <- c(
  "Female (XX)" = "#800020",
  "Male (XY)"   = "#66B3FF"
)

# Common visual theme

theme_bvg <- theme_classic(
  base_size = 14,
  base_family = "sans"
) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )

# Y-axis label used in all figures (see README, Section 2)
y_label <- "Normalized expression (log2 scale, presumed)"

# Jittered points: white fill keeps them visible on dark boxes
point_layer_args <- list(
  shape = 21, fill = "white", colour = "black",
  size = 1.6, stroke = 0.4, alpha = 0.85
)

# Axis labels with sample size, e.g. "Female (XX)\nn = 68"
n_labels <- function(dt) {
  counts <- dt[, .N, by = Sex]
  setNames(
    paste0(counts$Sex, "\nn = ", counts$N),
    counts$Sex
  )
}

# Create figure output directory
dir.create("figures", showWarnings = FALSE)


# 11. FIGURE 1 - REM16 EXPRESSION

# Extract REM16 observations
rem16_data <- copy(plot_data[
  Gene == "REM16"
])

# Verify sample numbers
print(table(rem16_data$Sex))

# Create REM16 box plot

p1 <- ggplot(
  rem16_data,
  aes(
    x = Sex,
    y = Expression,
    fill = Sex
  )
) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 1,
    outlier.shape = NA
  ) +
  
  do.call(geom_point, c(
    list(position = position_jitter(width = 0.12, seed = 7003)),
    point_layer_args
  )) +
  
  scale_fill_manual(
    values = sex_colors
  ) +
  
  scale_x_discrete(labels = n_labels(rem16_data)) +
  
  labs(
    title = "REM16 Expression by Plant Sex",
    x = "Plant Sex",
    y = y_label
  ) +
  
  theme_bvg +
  
  theme(
    legend.position = "none"
  )

# Display Figure 1
print(p1)


# 12. FIGURE 2 - FT1 EXPRESSION

# Extract valid FT1 observations
ft1_data <- copy(plot_data[
  Gene == "FT1"
])

# Verify sample numbers
print(table(ft1_data$Sex))

# Create FT1 box plot

p2 <- ggplot(
  ft1_data,
  aes(
    x = Sex,
    y = Expression,
    fill = Sex
  )
) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 1,
    outlier.shape = NA
  ) +
  
  do.call(geom_point, c(
    list(position = position_jitter(width = 0.12, seed = 7003)),
    point_layer_args
  )) +
  
  scale_fill_manual(
    values = sex_colors
  ) +
  
  scale_x_discrete(labels = n_labels(ft1_data)) +
  
  labs(
    title = "FT1 Expression by Plant Sex",
    x = "Plant Sex",
    y = y_label
  ) +
  
  theme_bvg +
  
  theme(
    legend.position = "none"
  )

# Display Figure 2
print(p2)


# 13. FIGURE 3 - COMBINED EXPRESSION

# Prepare data containing both genes
combined_data <- copy(plot_data[
  Gene %in% c("REM16", "FT1")
])

# Verify sample distribution
print(table(combined_data$Gene, combined_data$Sex))

# Sample size labels placed below each box
n_text <- combined_data[, .(
  N = .N,
  y = min(Expression) - 0.25
), by = .(Gene, Sex)]
n_text[, label := paste0("n = ", N)]

# Create Combined Figure
# Four box plots:
# REM16 Female / REM16 Male
# FT1 Female   / FT1 Male

p3 <- ggplot(
  combined_data,
  aes(
    x = Gene,
    y = Expression,
    fill = Sex
  )
) +
  
  geom_boxplot(
    width = 0.5,
    alpha = 1,
    outlier.shape = NA,
    position = position_dodge(
      width = 0.7
    )
  ) +
  
  # position_jitterdodge() needs a fill/colour aesthetic to dodge by;
  # fill is fixed to white, so colour is mapped to Sex (both black).
  geom_point(
    aes(colour = Sex),
    position = position_jitterdodge(
      jitter.width = 0.12, dodge.width = 0.7, seed = 7003
    ),
    shape = 21, fill = "white",
    size = 1.6, stroke = 0.4, alpha = 0.85
  ) +
  
  scale_colour_manual(
    values = c("Female (XX)" = "black", "Male (XY)" = "black"),
    guide = "none"
  ) +
  
  geom_text(
    data = n_text,
    aes(x = Gene, y = y, label = label, group = Sex),
    position = position_dodge(width = 0.7),
    size = 3.5,
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(
    values = sex_colors
  ) +
  
  labs(
    title = "REM16 and FT1 Expression by Plant Sex",
    x = "Gene",
    y = y_label,
    fill = "Plant Sex"
  ) +
  
  theme_bvg +
  
  theme(
    legend.position = "bottom",
    legend.title = element_text()
  )

# Display Figure 3
print(p3)


# 14. SAVE FINAL FIGURES

# Save Figure 1 - REM16
ggsave(
  filename = "figures/REM16_boxplot.png",
  plot = p1,
  width = 7,
  height = 5,
  dpi = 300
)

# Save Figure 2 - FT1
ggsave(
  filename = "figures/FT1_boxplot.png",
  plot = p2,
  width = 7,
  height = 5,
  dpi = 300
)

# Save Figure 3 - Combined
# Square aspect ratio

ggsave(
  filename = "figures/Combined_boxplot.png",
  plot = p3,
  width = 7,
  height = 7,
  dpi = 300
)


# 15. VERIFY FINAL OUTPUTS

# Confirm that all three plot objects exist
stopifnot(
  inherits(p1, "ggplot"),
  inherits(p2, "ggplot"),
  inherits(p3, "ggplot")
)

# Confirm all PNG figures were saved
figure_files <- c(
  "figures/REM16_boxplot.png",
  "figures/FT1_boxplot.png",
  "figures/Combined_boxplot.png"
)

stopifnot(all(file.exists(figure_files)))

# Display generated figure files
print(list.files("figures"))

# Record software versions for reproducibility
print(sessionInfo())

# Final confirmation
message("Task 3 analysis completed successfully.")