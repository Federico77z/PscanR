## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
file_path <- system.file("extdata", "df_processori.rds", package = "PscanR")
data<- readRDS(file_path)

data <- t(data)

N <- as.numeric(colnames(data))
means <- colMeans(data, na.rm = TRUE)

# Power regression
power_model <- lm(log(means) ~ log(N))  

smooth_N <- seq(min(N), max(N), length.out = 500)  # 500 punti tra 1 e 24
smooth_pred <- exp(predict(power_model, newdata = data.frame(N = smooth_N)))  # Previsioni

ss_total <- sum((means - mean(means))^2)
ss_residual <- sum((means - exp(predict(power_model)))^2)
r_squared <- 1 - (ss_residual / ss_total) 


## ----boxplot-example, fig.width=6, fig.height=4-------------------------------
boxplot(data,
        main = "PscanR Execution Time",
        xlab = "Cores (N)\nLiver dataset (950d−50u), JASPAR2020 (746 matrices), 329 promoters (1000bp)",
        cex.lab = 0.8,
        ylab = "Time (s)",
        col = "lightblue",
        las = 1, names = N, outline = FALSE)

lines(smooth_N, smooth_pred, col = "blue", lwd = 2)

legend("topright", legend = bquote("Power regression, " ~ R^2 == .(round(r_squared, 3))),
       col = "blue", lwd = 2)

## ----session-info, echo=FALSE-------------------------------------------------
sessioninfo::session_info()

