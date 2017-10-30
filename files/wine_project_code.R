# ==============================================
# project code
# A multiple regression model of white wine quality
# by Weiqi Weng, 04/13/2016 
# ==============================================

# load data
winedata <- read.csv(file="winequality-white-data.csv", header=TRUE)
summary(winedata)

# ==============================================

# data preprocessing

# outlier elimination
# fixed acidity is fine
# volatile.acidity
attach(winedata)
winedata$v.a[v.a <= 0.1 | v.a >= 1] <- NA
# citric.acid
winedata$c.a[c.a > 0.7] <- NA
# residual.sugar
winedata$r.s[r.s > 60] <- NA
# chlorides is fine
# total.sulfur.dioxide
winedata$t.s.d[t.s.d > 350] <- NA
# density
winedata$de[de > 1] <- NA
# pH and alcohol are fine
# no reference on sulphates
# quality 
winedata$qu[qu > 10 | qu < 0] <- NA
# finally delete outliers
# fix skewness of residual sugar
winedata$r.s <- log(r.s)
detach(winedata)
winedata <- na.omit(winedata)

# add 2 variable (relative index)
attach(winedata)
# free.sulfur.dioxide.proportion
winedata$f.s.d.p <- f.s.d/t.s.d
#citric.acid.proportion
winedata$c.a.p <- c.a/f.a
detach(winedata)

# ==============================================

# model selection

# data splitting
library(caret)
sample_size <- dim(winedata)[1]
set.seed(sample_size)
index <- createDataPartition(winedata$qu, p = 0.5, list = FALSE, times = 1)
build_winedata <- winedata[index, ]
val_winedata <- winedata[-index, ]

# correlation transformation
std_winedata <- as.data.frame(1/sqrt(sample_size - 1)*scale(winedata[-1]))
build_winedata_std <- std_winedata[index, ]
val_winedata_std <- std_winedata[-index, ]

# correlation matrix
cor(std_winedata)
# pairs(std_winedata)

# box plot for explanatory variable candidates
opar <- par(no.readonly = TRUE)
par(font.main = 1, cex = 0.8, mai = c(0.5,0.5,0.4,0.2))
boxplot(std_winedata[-c(1, 12)],
        col = "lightgray",
        main = "data of explanatory variables")
# Q-Q plot for response variable
qqnorm(scale(winedata$qu))
abline(0,1)
par(opar)

# exhaustive subset selection
library(leaps)
library(car)
winefit_full <- regsubsets(qu ~ ., data = build_winedata[-1], force.in = 1)
par(font.main = 1, cex = 0.8, mai = c(0.7,0.5,0.3,0.2))
plot(winefit_full, scale = "bic", main = "model selection with BIC")
par(opar)
subsets(winefit_full, statistic = "bic", main = "BIC plot for All Subsets Regression")
reg_summary_full <- summary(winefit_full)
coef(winefit_full, which.min(reg_summary_full$bic))

# ==============================================

# multicollinearity

# variance inflation factor
winefit1 <- lm(qu~f.a+v.a+r.s+f.s.d+de+pH+al+f.s.d.p,data = build_winedata)
summary(winefit1)
vif(winefit1)

# standardized model
winefit1_std <- lm(qu~0+f.a+v.a+r.s+f.s.d+de+pH+al+f.s.d.p,data = build_winedata_std)
summary(winefit1_std)

# rigid regression
library(MASS)
winefit_rigid <- lm.ridge(qu~., data = build_winedata[-1], lambda = seq(0,100,0.01))
par(font.main = 1, cex = 0.8, mai = c(1,0.9,0.5,0.2))
plot(winefit_rigid)
abline(h = 0, col = "black")
par(opar)

# ==============================================

# Least Square Method and ANOVA
winefit2 <- lm(qu~f.a+v.a+r.s+pH+al+f.s.d.p,data = build_winedata)
summary(winefit2)
anova(winefit2)

# lack of fit test
winefit2_lof <- lm(qu~factor(f.a)+factor(v.a)+factor(r.s)+
                     factor(pH)+factor(al)+factor(f.s.d.p),
                   data = build_winedata)
summary(winefit2_lof)
anova(winefit2_lof)

# ==============================================

# curvature and interaction investigation

# partial residual plot
par(font.main = 1, cex = 0.8)
crPlots(winefit2, main = "")
par(opar)

# introduce all the second-order terms and interactions
winefit_cur <- lm(qu~f.a+v.a+r.s+pH+al+f.s.d.p+I(f.a^2)+I(pH^2)+I(r.s^2)+f.a:pH+f.a:r.s+pH:r.s,
                  data = build_winedata)
summary(winefit_cur)

# add fixed.acidity^2, fixed.acidity*residual.sugar and residual.sugar^2
build_winedata$fars <- build_winedata$f.a*build_winedata$r.s
build_winedata$fa2 <- build_winedata$f.a*build_winedata$f.a
build_winedata$rs2 <- build_winedata$r.s*build_winedata$r.s

val_winedata$fars <- val_winedata$f.a*val_winedata$r.s
val_winedata$fa2 <- val_winedata$f.a*val_winedata$f.a
val_winedata$rs2 <- val_winedata$r.s*val_winedata$r.s

build_winedata_std$fars <- build_winedata_std$f.a*build_winedata_std$r.s
build_winedata_std$fa2 <- build_winedata_std$f.a*build_winedata_std$f.a
build_winedata_std$rs2 <- build_winedata_std$r.s*build_winedata_std$r.s

# fit the model
winefit_cur1 <- lm(qu~f.a+v.a+r.s+pH+al+f.s.d.p+fa2+fars+rs2,
                  data = build_winedata)
summary(winefit_cur1)
anova(winefit_cur1)

# lack of fit test
winefit_cur1_lof <- lm(qu~factor(f.a)+factor(v.a)+factor(r.s)+factor(pH)+
                          factor(al)+factor(f.s.d.p)+factor(fa2)+factor(fars)+factor(rs2),
                          data = build_winedata)
summary(winefit_cur1_lof)
anova(winefit_cur1_lof)

# ==============================================

# Diagnostics

# diagnostic plot
par(mfrow = c(2,2), font.main = 1, cex.main = 0.8, cex = 0.8, mai = c(0.4,0.5,0.4,0.2))
plot(winefit_cur1)
par(opar)

# studentized residual plot
# published by Dr. Robert I. Kabacoff in his excellent                                                       
# book, R in Action.
residplot <- function(fit, nbreaks=10) {
  z <- rstudent(fit)
  hist(z, breaks=nbreaks, freq=FALSE,
       xlab="Studentized Residual",
       main="Distribution of Errors")
  rug(jitter(z))
  curve(dnorm(x, mean=mean(z), sd=sd(z)),
        add=TRUE, lwd=2)
  lines(density(z)$x, density(z)$y, lwd=2, lty=2)
  legend("topleft",
         legend = c( "Normal", "Kernel"),
         lty=1:2, cex=.7)
}

par(font.main = 1, cex = 0.8)
residplot(winefit_cur1,15)
par(opar)

# influence plot
par(font.main = 1, cex = 0.8, mai = c(0.9,0.9,0.4,0.2))
influencePlot(winefit_cur1, 
              id.method = "observation id", 
              main = "Influence Plot", 
              sub = "Circle size proportional to Cook's D")
par(opar)

# Bonferroni outlier test
outlierTest(winefit_cur1)

# Cook's distance
build_data_size <- nrow(build_winedata)
cookD <- as.matrix(cooks.distance(winefit_cur1))
cookD_id <- which(cookD > 4/(build_data_size-10-1))
hist(cookD[cookD_id, ])
cookD_id <- which(cookD > 0.004)

# hat value
hat_val <- as.matrix(hatvalues(winefit_cur1))
hat_val_id <- which(hat_val > 20/build_data_size)
hist(hat_val[hat_val_id, ])
hat_val_id <- which(hat_val > 0.02)

# suspicious data points
sus_id <- c(cookD_id, hat_val_id)

# eliminate suspicious data points and remodel
new_build_winedata <- build_winedata[-sus_id, ]
new_build_winedata_std <- build_winedata_std[-sus_id, ]

new_winefit_cur1 <- lm(qu~f.a+v.a+r.s+pH+al+f.s.d.p+fa2+fars+rs2,
                       data = new_build_winedata)
summary(new_winefit_cur1)
anova(new_winefit_cur1)
par(mfrow = c(1,2), font.main = 1, cex = 0.8, mai = c(0.9,0.9,0.4,0.2))
influencePlot(winefit_cur1, 
              id.method = "observation id", 
              main = "Influence Plot", 
              sub = "Circle size proportional to Cook's D")
residplot(new_winefit_cur1,15)
par(opar)

# standardized model
new_winefit_cur1_std <- lm(qu~0+f.a+v.a+r.s+pH+al+f.s.d.p+fa2+fars+rs2,,
                           data = new_build_winedata_std)
summary(new_winefit_cur1_std)

# ==============================================

# sample for prediction
ori_y <- val_winedata["qu"]
index <- c("f.a","v.a","r.s","pH","al","f.s.d.p","fa2","fars","rs2")
val_winedata_size <- nrow(val_winedata)
set.seed(val_winedata_size)
pred_sample_id <- sample(1:val_winedata_size, 5)
val_set <- as.data.frame(val_winedata[index])
val_set <- val_set[-pred_sample_id,]
pre_sample <- predict(new_winefit_cur1, val_winedata[pred_sample_id, ],interval="predict")

# Validation
pre_y <- predict(new_winefit_cur1, val_set)
diff <- ori_y - pre_y
mspe <- sum(diff^2)/nrow(val_winedata)

# effect plot
library(effects)
plot(Effect("v.a",new_winefit_cur1,list(r.s=seq(1.2,1.6,0.1)),multiline = TRUE))
plot(Effect("f.a",new_winefit_cur1,list(r.s=seq(1.2,1.6,0.1)),multiline = TRUE))
plot(Effect("r.s",new_winefit_cur1,list(f.a=seq(6.5,6.9,0.1)),multiline = TRUE))
plot(Effect("fars",new_winefit_cur1,list(r.s=seq(1.2,1.6,0.1)),multiline = TRUE))

# ==============================================

# relative weight

# The relweights function determines the relative importance of each   
# independent variable to the dependent variable in an OLS regression. 
# The code is adapted from an SPSS program generously provided by      
# Dr. Johnson and published by Dr. Robert I. Kabacoff in his excellent                                                       
# book, R in Action.                                                                     
# See Johnson (2000, Multivariate Behavioral Research, 35, 1-19) for   
# an explanation of how the relative weights are derived.

relweights <- function(fit, ...) {
  R <- cor(fit$model)
  nvar <- ncol(R)
  rxx <- R[2:nvar, 2:nvar]
  rxy <- R[2:nvar, 1]
  svd <- eigen(rxx)
  evec <- svd$vectors
  ev <- svd$values
  delta <- diag(sqrt(ev))
  
  # correlations between original predictors and new orthogonal variables
  lambda <- evec %*% delta %*% t(evec)
  lambdasq <- lambda^2
  
  # regression coefficients of Y on orthogonal variables
  beta <- solve(lambda) %*% rxy
  rsquare <- colSums(beta^2)
  rawwgt <- lambdasq %*% beta^2
  import <- (rawwgt/rsquare) * 100
  lbls <- names(fit$model[2:nvar])
  rownames(import) <- lbls
  colnames(import) <- "Weights"
  
  # plot results
  barplot(t(import), names.arg = lbls, ylab = "% of R-Square", 
          xlab = "Predictor Variables", main = "Relative Importance of Predictor Variables", 
          sub = paste("R-Square = ", round(rsquare, digits = 3)), 
          ...)
  return(import)
}

relweights(new_winefit_cur1)

# ==============================================

# hypothesis testing 1
x1 <- winedata$f.a
x3 <- (-0.011*winedata$r.s-0.378)/0.064
x1bar <- mean(x1)
x3bar <- mean(x3)
s1 <- sd(x1)
s3 <- sd(x3)
up1 <- (s1^2/sample_size+s3^2/sample_size)^2
down1 <- (s1^2/sample_size)^2/(sample_size-1)+(s3^2/sample_size)^2/(sample_size-1)
df1 <- round(up1/down1)
T1 <- (x1bar-x3bar)/sqrt(s1^2/sample_size+s3^2/sample_size)
tq1 <- qt(0.95,df)

# hypothesis testing 2
x2 <- (0.174*winedata$r.s+0.033)/0.011
x2bar <- mean(x2)
s2 <- sd(x2)
up2 <- (s1^2/sample_size+s2^2/sample_size)^2
down2 <- (s1^2/sample_size)^2/(sample_size-1)+(s2^2/sample_size)^2/(sample_size-1)
df2 <- round(up2/down2)
T2 <- (x2bar-x1bar)/sqrt(s1^2/sample_size+s2^2/sample_size)
tq2 <- qt(0.95,df)
