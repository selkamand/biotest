# Libraries
library(here)
library(VariantAnnotation)
library(StructuralVariantAnnotation)
library(rtracklayer)

# Paths
path_sv_vcf <- here("mutations/structural_variants/tumor_normal.2sample.purple.sv.vcf")
outdir <- here("mutations/structural_variants/")

vcf <- readVcf(path_sv_vcf)

# Export breakpoints to BEDPE
bpgr <- breakpointRanges(vcf)
write.table(breakpointgr2bedpe(bpgr), file = paste0(outdir, "purple.sv.breakpoints.bedpe"), sep = "\t", quote = FALSE, col.names = FALSE)

# Export single breakends to BED
begr <- breakendRanges(vcf)
begr$score <- begr$QUAL
begr$score[is.na(begr$score)] <- 0 # Set inferred breakends score to 0
export(begr, con = paste0(outdir, "purple.sv.single_breakends.bed"))
