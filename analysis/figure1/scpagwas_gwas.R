# scpagwas gwas
source(file.path(Sys.getenv("AM_REPO_DIR", unset = "."), "R", "paths.R"))
suppressPackageStartupMessages(library(data.table))
out_dir <- am_data("vibe/191_scpagwas_hg19"); t0 <- Sys.time()
log_phase <- function(p) {
  z <- sprintf("%s\t%s\t%.1f", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), p, as.numeric(difftime(Sys.time(), t0, units = "secs"))); cat(z, "\n"); cat(z, "\n", file = file.path(out_dir, "run.log"), append = TRUE)
}
S <- am_data("HTAN/scPagwas/scPagwas_lc")
g38 <- fread(file.path(S, "lcgwas_Pagwas.txt"))
g37 <- fread(file.path(S, "hg37_lcgwas_GCST90043864.tsv.gz"), select = c("chromosome", "variant_id", "base_pair_location"))
setnames(g37, c("chrom37", "rsid", "pos37")); g37 <- unique(g37, by = "rsid")
m <- merge(g38, g37, by = "rsid", all.x = TRUE)
audit <- data.table(
  n_input = nrow(g38), n_matched = sum(!is.na(m$pos37)), pct_matched = 100 * mean(!is.na(m$pos37)),
  n_chrom_mismatch = sum(!is.na(m$pos37) & as.character(m$chrom) != as.character(m$chrom37)),
  rs1042522_hg38 = g38[rsid == "rs1042522", pos], rs1042522_hg19 = m[rsid == "rs1042522", pos37]
)
print(audit); fwrite(audit, file.path(out_dir, "data", "gwas_liftover_audit.csv"))
out <- m[!is.na(pos37) & as.character(chrom) == as.character(chrom37)][, pos := pos37][, c("pos37", "chrom37") := NULL]
setcolorder(out, names(g38)); setorder(out, chrom, pos)
fwrite(out, file.path(out_dir, "data", "lcgwas_Pagwas_hg19.txt"), sep = "\t")
log_phase("gwas_hg19")
