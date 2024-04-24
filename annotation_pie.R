setwd(".")

df <- read.csv("tblastn_log.csv")

annotations <- df[, 3]

not_found <- grepl("not found by blast", annotations)
overlap <- grepl(";blast_id", annotations)
no_overlap <- !overlap

df_categorized <- data.frame(
  Annotation = annotations,
  Category = ifelse(not_found, "Not Found",
                    ifelse(overlap, "Overlap", "No Overlap"))
)

category_counts <- table(df_categorized$Category)
labels_with_counts <- paste(names(category_counts), ": ", category_counts)
colors <- c("#A3DEF8","#DB2D43","#FBF9AF")

png("tblastn_pie.png", width = 520, height = 450)

pie(category_counts, labels = labels_with_counts,
    col = colors, main = "tBLASTn Results")

dev.off()
