# Run LDSC with custom functional categories (made in make-annot.R + run_LDScore_estimation.R)
# Code adapted from Peter Hickey (https://github.com/hansenlab/BrainEpigenomeNN/blob/master/SLDSR/scripts/


# NOTE: On JHPCE, requires `module load python/2.7.6`

library(parallel)
options("mc.cores" = 10)

args <- commandArgs(TRUE)
i <- as.integer(args[1])
message("i = ", i)

categories <- readRDS("../objects/categories.rds")
gwasss <- list.files("../extdata/munge_sumstats/Phase1", full.names = TRUE,
                     pattern = glob2rx("*.sumstats.gz"))

seqlevels <- 1:22

message("category = ", names(categories)[i])

### ============================================================================
### Adjusting for baseline
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (cn == "CNS") {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../extdata/Phase1/cell_type_groups/CNS.,",
                    "../extdata/Phase1/baseline/baseline. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/CNS.", bn, ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/", cn, ".Phase1.,",
                    "../extdata/Phase1/baseline/baseline. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, ".", bn, ".Phase1 ",
                    "--print-coefficients")
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

### ============================================================================
### Marginal analysis (feature + no baseline features except base)
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (cn == "CNS") {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/base.Phase1.,",
                    "../extdata/Phase1/cell_type_groups/CNS. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, "_no_baseline_except_base.", bn,
                    ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../output/ldsc/base.Phase1.,",
                    "../output/ldsc/", cn, ".Phase1. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn, "_no_baseline_except_base.", bn,
                    ".Phase1 ",
                    "--print-coefficients")
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

### ============================================================================
### Adjusting for baseline + CNS + chromHMM_union + H3K27ac
###

lapply(names(categories)[i], function(cn) {
  mclapply(gwasss, function(x) {
    bn <- sub(".sumstats.gz", "", basename(x))
    if (!cn %in% c("CNS", "chromHMM_union", "H3K27ac")) {
      cmd <- paste0("python ",
                    "/users/phickey/software/ldsc/ldsc.py ",
                    "--h2 ", x, " ",
                    "--w-ld-chr ../extdata/Phase1/weights_hm3_no_hla/weights. ",
                    "--ref-ld-chr ../extdata/Phase1/baseline/baseline.,",
                    "../extdata/Phase1/cell_type_groups/CNS.,",
                    "../output/ldsc/chromHMM_union.Phase1.,",
                    "../output/ldsc/H3K27ac.Phase1.,",
                    "../output/ldsc/", cn, ".Phase1. ",
                    "--overlap-annot ",
                    "--frqfile-chr ../extdata/Phase1/1000G_frq/1000G.mac5eur. ",
                    "--out ../output/ldsc/", cn,
                    ".adjusting_for_baseline_and_CNS_and_chromHMM_union_and_H3K27ac.",
                    bn,
                    ".Phase1 ",
                    "--print-coefficients")
    } else {
      cmd <- paste0("echo Nothing to do for ", cn)
    }
    print(cmd)
    system(cmd)
  }, mc.cores = 10)
})

