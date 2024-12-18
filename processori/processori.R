#Script to obtain the execution time

#Start from the background 

#source("../BG_scripts/load_packages.R")
#source("../R/Build_Background.R")
#source("../R/Helper_functions.R")
#source("../R/Scan_and_post_processing.R")
#source("../R/PSMatrix_class.R")

#txdb <- makeTxDbFromUCSC(genome="hg38", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
#seqlevels(txdb) <- seqlevels(txdb)[1:24]

#saveRDS(txdb, "txdb.rds")

#txdb_path <- system.file("extdata", "txdb.rds", package = "PscanR")
#txdb <- readRDS(txdb_path)

#opts <- list()
#opts[["collection"]] <- "CORE"
#opts[["tax_group"]] <- "vertebrates"

#J2020 <- getMatrixSet(JASPAR2020, opts)

#prom_rng <- promoters(txdb, upstream = 950, downstream = 50, use.names = TRUE)
#prom_rng$tx_name_clean <- sub("\\..*$", "", prom_rng$tx_name)
#bg_hg38_2020 <- ps_build_bg_from_file("BG_scripts/hg38_J2020/J2020_hg38_950u_50d_UCSC.psbg.txt", J2020)

#Get the target sequences

#target <- read.csv("processori/liver.txt", header = F) 

#prom_range <- prom_rng[prom_rng$tx_name_clean %in% target[,1]]
#prom_seq <- getSeq(x = BSgenome.Hsapiens.UCSC.hg38, prom_range)

#Create the execution time function
iterations <- 100
core_numbers <- c(1,2,4,6,8,12,18,24)

execution_time <- function(prom_seq, bg, ncores, output_file){
temp <- numeric(iterations)

for (i in 1:iterations){
  timing <- system.time({
    results <- PscanR::pscan(prom_seq, bg, BPPARAM = BiocParallel::SnowParam(ncores))}) 
  if (i == 1){
    saveRDS(PscanR::ps_results_table(results), output_file) 
  }
  temp[i] <- timing["elapsed"]
}
return(temp)
}

all_cores <- function(core_numbers, prom_seq, bg, iterations) {
  results_list <- list()
  
  for (core in core_numbers) {
    output_file <- paste0("core", core, "_liver_950d_50u.rds")
    results_list[[as.character(core)]] <- execution_time(prom_seq, bg, core, output_file)
  }
  
  return(as.data.frame(results_list))
}

#Calculate the execution times

df_processori <- all_cores(core_numbers, prom_seq, bg_hg38_2020, iterations)
colnames(df_processori) <- core_numbers

saveRDS(df_processori, "df_processori_SnowParam.rds")
library(knitr)
kable(df_processori, format = "markdown", file = "df_processori_SnowParam.md")

#Create the boxplot

#boxplot(df_processori, main = "Execution Times of PscanR with SnowParam",
#        sub = "Liver dataset (950d-50u), JASPAR2020 (746 matrices), 329 promoters (1000bp)", 
#        xlab = "Number of Cores", 
#        ylab = "Execution Time (s)",
#        las = 1,
#        cex.axis = 0.8,
#        at = 1:length(core_numbers),  
#        names = paste0(core_numbers, " core", ifelse(core_numbers > 1, "s", "")))

#grid(nx = NA, ny = NULL, lty = "dotted", col = "lightgray")

#Create the regression curve 

#medians <- apply(df_processori, 2, median)
#core_numbers <- as.numeric(colnames(df_processori))
#model_quadratic <- lm(medians ~ poly(core_numbers, 2))
#predicted <- predict(model_quadratic, newdata = data.frame(core_numbers = core_numbers))
#lines(core_numbers, predicted, col = "blue", lwd = 2)
#legend("topright", legend = c("Regression Line"), col = "blue", lwd = 2)

#compare the different models with anova test

#model_linear <- lm(medians ~ core_numbers)
#model_cubic <- lm(medians ~ poly(core_numbers, 3))
#model_quartic <- lm(medians ~ poly(core_numbers, 4))
#model_quint <- lm(medians ~ poly(core_numbers, 5))

#anova_result <- anova(model_linear, model_quadratic, model_cubic, model_quartic, model_quint)
#print(anova_result) #the quadratic regression seems the best one 



