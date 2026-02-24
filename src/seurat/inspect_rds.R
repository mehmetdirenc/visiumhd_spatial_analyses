library(Seurat)
setwd("/mnt/storage2/users/ahmungm1/research_projects/xenium_spatial/")
obj <- readRDS("Data_Michigan_02_2026.RDS")

obj


colnames(obj@meta.data)
head(obj@meta.data, 6)

summ_uniq <- function(x, name, n_show = 6){
  u <- sort(unique(as.character(x)))
  cat("\n", name, ":", length(u), "\n", sep="")
  if (length(u) < 10) {
    cat("  values: ", paste(u, collapse = ", "), "\n", sep="")
  } else {
    cat(paste(head(u, n_show), collapse = ", "), "\n", sep="")
  }
}

# candidates for patient/sample IDs
summ_uniq(obj$batch,      "batch (candidate patient/sample ID)")
summ_uniq(obj$orig.ident, "orig.ident")
summ_uniq(obj$temp,       "temp")

# biology groupings
summ_uniq(obj$disease,    "disease")
summ_uniq(obj$lesional,   "lesional")
summ_uniq(obj$site,       "site")

# cell types / clusters
summ_uniq(Idents(obj),         "celltypes (Idents)")
summ_uniq(obj$seurat_clusters, "seurat_clusters")


# how many temp per batch?
table(obj$batch, obj$temp) |> (function(x) rowSums(x>0))() |> sort(decreasing=TRUE)

# or simpler:
tapply(obj$temp, obj$batch, function(v) length(unique(v))) |> sort(decreasing=TRUE)

# how many orig.ident per batch?
tapply(obj$orig.ident, obj$batch, function(v) length(unique(v))) |> sort(decreasing=TRUE)

tapply(obj$temp, obj$disease, function(v) length(unique(v)))
tapply(obj$orig.ident, obj$disease, function(v) length(unique(v)))

sort(unique(as.character(obj$celltype)))


# unique values to TSV (long format: category + value)
get_uniques <- function(vec) sort(unique(as.character(vec)))

df_out <- rbind(
  data.frame(category="celltype",   value=get_uniques(obj$celltype)),
  data.frame(category="disease",    value=get_uniques(obj$disease)),
  data.frame(category="batch",      value=get_uniques(obj$batch)),       # candidate patient/donor ID
  data.frame(category="orig.ident", value=get_uniques(obj$orig.ident)),  # sample/library name
  data.frame(category="temp",       value=get_uniques(obj$temp))         # sample key
)

write.table(df_out,
            file="Data_Michigan_02_2026_unique_names.tsv",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
