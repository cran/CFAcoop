## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----loadData-----------------------------------------------------------------
library(CFAcoop)
data("CFAdata")
summary(CFAdata)

## ----show1, fig.width=7, fig.height=5-----------------------------------------
data("CFAdata")
data1 <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "T47D")
data2 <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "BT20")
SF <- vector("list", 2)
SF[[1]] <- analyze_survival(
  RD = data1[, c("Cells seeded","0 Gy","1 Gy","2 Gy","4 Gy","6 Gy","8 Gy")], 
  name = as.character(data1[1,1]), 
  xtreat = c(0, 1, 2, 4, 6, 8), 
  c_range = c(5, 20, 100))
SF[[2]] <- analyze_survival(
  RD = data2[,-1], 
  name = as.character(data2[1,1]), 
  xtreat = c(0, 1, 2, 4, 6, 8))
plot_sf(SF = SF)

## ----export_sf----------------------------------------------------------------
summary_df <- export_sf(SF)

## ----summary_export_sf--------------------------------------------------------
colnames(summary_df)
head(summary_df)
summary(summary_df)

## ----sf_details---------------------------------------------------------------
SF[[1]]$fit[1]
SF[[2]]$SF

## ----PowerReg, fig.width=5, fig.height=4--------------------------------------
data <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "BT20")
data <- aggregate(x = data[, -1], by = list(data$`Cells seeded`), FUN = "mean", na.rm = TRUE)
par_0 <- pwr_reg(seeded = data$`Cells seeded`, counted = data$`0 Gy`)
par_0$coefficients
plot(x = log10(data$`Cells seeded`), y = log10(data$`0 Gy`),xlim = c(2,3.5))
abline(a = log10(exp(1)) * par_0$coefficients[1, 1], b = par_0$coefficients[2, 1])

## ----TestingCooperation-------------------------------------------------------
p_value <-
  (1 - pt(
    q = (par_0$coefficients[2, 1] - 1) / par_0$coefficients[2, 2],
    df = par_0$df[2]
  ))

## ----calculateSF--------------------------------------------------------------
data <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "BT20")
data <- aggregate(x = data[, -1], by = list(data$`Cells seeded`), FUN = "mean", na.rm = TRUE)
par_0 <- pwr_reg(seeded = data$`Cells seeded`, counted = data$`0 Gy`)
par_4 <- pwr_reg(seeded = data$`Cells seeded`, counted = data$`4 Gy`)
calculate_sf(par_ref = par_0, par_treat = par_4, c_range = c(5, 100))

## ----bootstrapUncertainty-----------------------------------------------------
ParSet_0 <- mvtnorm::rmvnorm(n = 10^4,mean = par_0$coefficients[,1],
                             sigma = vcov(object = par_0))
ParSet_4 <- mvtnorm::rmvnorm(n = 10^4,mean = par_4$coefficients[,1],
                             sigma = vcov(object = par_4))
parbootstrap <- calculate_sf(par_ref = ParSet_0, 
                             par_treat = ParSet_4, 
                             c_range = c(5, 100))
cat(c(mean(parbootstrap),sd(parbootstrap)))
cat(quantile(x = parbootstrap,probs = c(0.025,0.5,0.975)))

## ----PEfail1a, fig.width=5, fig.height=4--------------------------------------
data(CFAdata)
data <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "T47D")
data <- aggregate(x = data[, -1], by = list(data$`Cells seeded`), FUN = "mean", na.rm = TRUE)
PE_x <- data$`4 Gy` / data$`Cells seeded`
PE_0 <- data$`0 Gy` / data$`Cells seeded`
plot(x = rep(c(0, 1), each = 18), y = c(PE_0, PE_x), lty = 0, ylim = c(0,0.5),xlim = c(-0.1,1.1),
     xlab = "treatment", ylab = "C / S", axes = FALSE, main = "T47D")
axis(side = 1,at = c(0,1),labels = c("0 Gy","4 Gy"))
axis(side = 2,at = seq(0,0.5,0.1))

## ----PEfail1b, fig.width=5, fig.height=4--------------------------------------
SF_resample <- rep(PE_x, each = length(PE_0)) / rep(PE_0, times = length(PE_x))
hist(SF_resample, breaks = 25,xlim = c(0.12,0.25),xlab = "(C_4/S_4) / (C_0/S_0)",
     main = "valid PE-based SF'-values")

## ----PEfail2------------------------------------------------------------------
range(SF_resample,na.rm = TRUE)
as_nc_0 <- analyze_survival(RD = data[,-1],c_range = c(20))
as_nc_0$uncertainty[[3]]

## ----PEfailCoop, fig.width=5, fig.height=4------------------------------------
data(CFAdata)
data <- subset.data.frame(x = CFAdata, subset = CFAdata$cell.line == "BT20")
data <- aggregate(x = data[, -1], by = list(data$`Cells seeded`), FUN = "mean", na.rm = TRUE)
PE_x <- data$`4 Gy` / data$`Cells seeded`
PE_0 <- data$`0 Gy` / data$`Cells seeded`
plot(x = rep(c(0, 1), each = length(PE_x)), 
     y = c(PE_0, PE_x), lty = 0, ylim = c(0,0.08),xlim = c(-0.1,1.1),
     xlab = "treatment", ylab = "C / S", axes = FALSE, main = "BT20")
axis(side = 1,at = c(0,1),labels = c("0 Gy","4 Gy"))
axis(side = 2,at = seq(0,0.08,0.02),las = 1)
SF_resample <- rep(PE_x, each = length(PE_0)) / rep(PE_0, times = length(PE_x))
hist(SF_resample, breaks = 100,xlim = c(0,10),xlab = "(C_4/S_4) / (C_0/S_0)",
     main = "valid PE-based SF'-values")

## ----PEfailCoop2, fig.width=6, fig.height=5-----------------------------------
range(SF_resample,na.rm = TRUE)
as_c_4 <- analyze_survival(RD = data[,-1],c_range = c(20))
as_c_4$uncertainty[[3]]

## ----PEfailCoop2posta---------------------------------------------------------
sem_SF20 <- (as_c_4$uncertainty[[3]][2]-as_c_4$uncertainty[[3]][1])/4
cat("    SF(20) : ",round(as_c_4$SF[[3]],digits = 3),"\n")
cat("SEM SF(20) : ",round(sem_SF20,digits = 3),"\n")

## ----PEfailCoop2postb---------------------------------------------------------
cat("    SF'    : ",round(mean(SF_resample,na.rm = TRUE),digits = 3),"\n")
cat("SEM SF'    : ",round(sd(SF_resample,na.rm = TRUE)/
                            sqrt(length(SF_resample)),digits = 3),"\n")

