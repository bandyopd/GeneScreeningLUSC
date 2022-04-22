A Prognostic Gene Signature in Lung Squamous Cell Carcinoma
================

Lung cancer is the leading cause of cancer-related death worldwide, and
non-small cell lung cancer (NSCLC) is the most common histological cell
type and often presents in an advanced stage. Lung squamous cell
carcinoma (LUSC) is a NSCLC subtype associated with poor clinical
prognosis and lacks available targeted therapy as it is characterized by
high biological and clinical heterogeneity. A crucial component of
contemporary research with high-throughput data lies in identifying
novel biomarkers that can be predictive for LUSC patient prognosis. In
this document, we demonstrate the implementation of a recently proposed
gene screening approach, ESIS (Ke et al. 2022), for identifying a
prognostic gene signature of LUSC.

## TCGA lung squamous cell carcinoma dataset

Gene expression data and clinical data for patients with LUSC were
acquired from TCGA lung cancer cohort. We only used samples of primary
tumors with unique case ID. For a patient without an event (death), the
overall survival time from first diagnosis was censored by the last
follow-up date. Aside from 20,531 genes, 4 clinical covariates were also
included in the analysis: age at diagnosis, gender, smoking history and
tumor stage. In total, 473 cases with completed data were studied. 1,862
genes were excluded prior to screening due to complete missing or low
expression (with an interquartile range of 0).

``` r
### load the cleaned dataset and ESIS functions
lusc <- read.table(file = "/Users/kec2/GeneScreeningDemo1/LUSC.csv", header = T)
source("/Users/kec2/GeneScreeningDemo1/source.R")
```

## Statistical Analysis

Data were divided into a training set and a testing set in a ratio of
4:1 by stratified randomization based on censoring. The training set was
comprised of 379 samples and the testing set was comprised of 95
samples.

``` r
### training - testing split
grp <- table(lusc$OS_STATUS)
tr.prop <- 0.8 #training proportion
set.seed(321) #reproducibility
#stratified split by censoring
tr.ind <- c(sample(c(rep(T,round(grp[1]*0.8)),rep(F,round(grp[1]*0.2)))),
            sample(c(rep(T,round(grp[2]*0.8)),rep(F,round(grp[2]*0.2)))))
train.x <- lusc[tr.ind,-(1:9)]
train.y <- lusc[tr.ind,c("OS_MONTHS","OS_STATUS")]
colnames(train.y) <- c('time', 'status')
test.x <- lusc[!tr.ind,-(1:9)]
test.y <- lusc[!tr.ind,c("OS_MONTHS","OS_STATUS")]
colnames(test.y) <- c('time', 'status')
```

We first performed ESIS on the training set and pre-selected 379 − 1 =
378 genes. The characteristic kernel as well as the smoothing kernel
were both chosen to be the Gaussian kernel.

``` r
### ESIS
ix0 <- esis(train.x, train.y$time, train.y$status, d = 378)$ix
```

A Penalized Cox model with LASSO regularization (abbreviated as PenCox
henceforth) was then applied to the reduced training data for further
gene selection and prognosis simultaneously. The optimal tuning
parameter was determined through 10-fold cross validation.

``` r
### PenCox
rtrain.x <- as.matrix(scale(train.x[,ix0])) #reduced training data
rtest.x <- as.matrix(scale(test.x[,ix0])) #reduced testing data
library(glmnet)
set.seed(1476) #reproducibility
cvcox <- cv.glmnet(rtrain.x, as.matrix(train.y), 
                   family = "cox", lambda = exp(seq(0,-10,-0.1)))
coef <- coef(cvcox, s = cvcox$lambda.min)
ix <- ix0[which(coef[,1]!=0)] #selected indices
colnames(train.x)[ix] #selected genes
```

    ## [1] "NACC2"     "FAM65A"    "LOC641845" "MON1B"     "IBTK"      "SDHAF3"

Six genes were selected ultimately: NACC2, FAM65A, LOC641845, MON1B,
IBTK and SDHAF3. A patient’s risk score was calculated as the linear
predictor of the fitted PenCox model. Patients were classified as having
a high-risk gene signature or a low-risk gene signature, with the median
risk score of the training group being the cutoff. The same cutoff value
was also applied to assign the test samples into two risk groups. To
evaluate the predictive performance of the PenCox model built upon the
ESIS-selected genes, the Kaplan-Meier curves of the two groups for
overall survival were compared using the log-rank test.

``` r
### prognosis based on the 6-gene signature
# risk score
fitted <- rtrain.x%*%coef
pred <- predict(cvcox, s = cvcox$lambda.min, newx = rtest.x, type = "link")
cutoff <- median(fitted) #cutoff
# assign risk group
risk.train <- as.vector(fitted > cutoff)
risk.test <- pred > cutoff
# log-rank test
logrank.train <- survdiff(Surv(train.y$time, train.y$status) ~ risk.train)
logrank.test <- survdiff(Surv(test.y$time, test.y$status) ~ risk.test)
# Kaplan-Meier curve
par(mfrow = c(1, 2))
plot(survfit(Surv(train.y$time,train.y$status) ~ risk.train), mark.time = T,
     col = c("blue","red"), xlab = "Time in months", ylab = "Overall survival", xlim = c(0, 175),
     main = "ESIS + PenCox, Training Cohort")
legend(125, 1, col = c("red", "blue"), lty = 1, legend = c("high-risk", "low-risk"), bty = "n")
legend("bottomleft", 
       legend = bquote(p-value == .(format(1-pchisq(logrank.train$chisq, 1), digits = 4))),
       bty = "n")
plot(survfit(Surv(test.y$time,test.y$status) ~ risk.test), mark.time = T,
     col = c("blue","red"), xlab = "Time in months", ylab = "Overall survival", xlim = c(0, 175),
     main = "ESIS + PenCox, Testing Cohort")
legend(125, 1, col = c("red", "blue"), lty = 1, legend = c("high-risk", "low-risk"), bty = "n")
legend("bottomleft", 
       legend = bquote(p-value == .(format(1-pchisq(logrank.test$chisq, 1), digits = 4))), 
       bty = "n")
```

<img src="GeneScreeningLUSC_files/figure-gfm/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

The ESIS+PenCox model achieved nice separation of the two risks groups.
Patients with a high-risk gene signature had a shorter median overall
survival than those with a low-risk gene signature (34.7 months vs. 71.3
months; p-value = 0.078). Finally, a Cox model was fitted on the entire
dataset to make inference about independent prognostic factors
associated with survival, and the selected gene signature, age, gender,
tumor stage, and smoking history were used as covariates.

``` r
### multi-variable Cox model
genesig <- drop(as.matrix(lusc[,-(1:9)][, ix0])%*%as.matrix(coef)) #gene signature
# clinical covariates (create dummy variables for factors)
lusc[,2:4] <- lapply(lusc[,2:4], factor)
lusc$AJCC_PATHOLOGIC_TUMOR_STAGE <- relevel(lusc$AJCC_PATHOLOGIC_TUMOR_STAGE, "4")
lusc$TOBACCO_SMOKING_HISTORY_INDICATOR <- relevel(lusc$TOBACCO_SMOKING_HISTORY_INDICATOR, "5")
clinical0 <- model.matrix(~AGE+SEX.x+AJCC_PATHOLOGIC_TUMOR_STAGE+TOBACCO_SMOKING_HISTORY_INDICATOR-1,lusc)
clinical <- clinical0[, c("AGE", "SEX.x0", "AJCC_PATHOLOGIC_TUMOR_STAGE1", "AJCC_PATHOLOGIC_TUMOR_STAGE2",
                          "TOBACCO_SMOKING_HISTORY_INDICATOR1", "TOBACCO_SMOKING_HISTORY_INDICATOR2")]
colnames(clinical) <- c("age", "genderM", "stage1", "stage2", "non-smoker", "current-smoker")
coxdata <- data.frame(lusc[,6:7], gene=genesig, clinical)
coxfit <- coxph(Surv(OS_MONTHS, OS_STATUS)~., coxdata)
coxfit
```

    ## Call:
    ## coxph(formula = Surv(OS_MONTHS, OS_STATUS) ~ ., data = coxdata)
    ## 
    ##                     coef exp(coef)  se(coef)      z        p
    ## gene            2.532868 12.589557  0.571149  4.435 9.22e-06
    ## age             0.023065  1.023333  0.008745  2.637  0.00835
    ## genderM        -0.080908  0.922279  0.167398 -0.483  0.62887
    ## stage1         -0.530645  0.588225  0.177207 -2.994  0.00275
    ## stage2         -0.462956  0.629420  0.195019 -2.374  0.01760
    ## non.smoker      0.662908  1.940426  0.433477  1.529  0.12619
    ## current.smoker  0.429529  1.536533  0.152071  2.825  0.00474
    ## 
    ## Likelihood ratio test=46.77  on 7 df, p=6.186e-08
    ## n= 473, number of events= 207

According to the multivariable Cox regression, the 6-gene signature
selected by ESIS+PenCox was a strong predictor with an hazard ratio of
12.59 (p-value
![\<](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3C "<")
0.001) after controlling for other clinical covariates. There was a 2%
increase in the expected hazard relative to a one year increase in age
(p-value=0.008). Having the first or the second stage cancer reduced the
hazard by 41% (p-value=0.003) and 37% (p-value=0.018), respectively,
compared to later stage cancers. Current smokers were associated with
worse prognosis (p-value=0.005).

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ke2022gene" class="csl-entry">

Ke, Chenlu, Dipankar Bandyopadhyay, Mark Acunzo, and Robert Winn. 2022.
“Gene Screening in High-Throughput Right-Censored Cancer Data.”

</div>

</div>
