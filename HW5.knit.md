---
title: "HW5"
output: pdf_document
date: "2023-03-31"
---



### QUESTION 1

#### (a) Model Fitting and Interpretation

```r
## fit a poisson log-link model
crab = read.table("HW5-crab.txt", header = TRUE)
M1 = glm(Sa ~ W, family=poisson(link=log), data=crab)
summary(M1)
```

```
## 
## Call:
## glm(formula = Sa ~ W, family = poisson(link = log), data = crab)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.8526  -1.9884  -0.4933   1.0970   4.9221  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)    
## (Intercept) -3.30476    0.54224  -6.095  1.1e-09 ***
## W            0.16405    0.01997   8.216  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for poisson family taken to be 1)
## 
##     Null deviance: 632.79  on 172  degrees of freedom
## Residual deviance: 567.88  on 171  degrees of freedom
## AIC: 927.18
## 
## Number of Fisher Scoring iterations: 6
```

```r
# goodness of fit
res.pearson = residuals(M1, type="pearson")
G.stat = sum(res.pearson^2)
res.dev = residuals(M1, type="deviance")
D.stat = sum(res.dev^2)
1-pchisq(G.stat, df=173-2)
```

```
## [1] 0
```

```r
1-pchisq(D.stat, df=173-2)
```

```
## [1] 0
```
**Goodness of fit:**

The Pearson's Chi-square test statistic is 544.16, deviance test statistic is 567.88. The p-values of both of them are less than 0.05 under 171 df, indicating the fit is not good here.

**Interpretation of the model:**

We fitted a poisson log linear model here where:

log(E(Y)) = β0 + β1·X

with Y=Sa (number of satellites), X=W (carapace width).

Here β0 = -3.3, β1 =  0.16. Indicating that when the carapace width is 0, the expected number of satellites is 0.04 . When there's one unit change in carapace width, the expected change in the number of satellites is 1.18.


#### (b) Model Fitting with More Predictors and Interpretation

```r
# fit again with 2 predictors
M2 = glm(Sa ~ W + Wt, family=poisson(link=log), data=crab)
summary(M2)
```

```
## 
## Call:
## glm(formula = Sa ~ W + Wt, family = poisson(link = log), data = crab)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.9308  -1.9705  -0.5481   0.9700   4.9905  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)   
## (Intercept) -1.29168    0.89929  -1.436  0.15091   
## W            0.04590    0.04677   0.981  0.32640   
## Wt           0.44744    0.15864   2.820  0.00479 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for poisson family taken to be 1)
## 
##     Null deviance: 632.79  on 172  degrees of freedom
## Residual deviance: 559.89  on 170  degrees of freedom
## AIC: 921.18
## 
## Number of Fisher Scoring iterations: 6
```

```r
# goodness of fit
res.pearson2 = residuals(M2, type="pearson")
G.stat2 = sum(res.pearson2^2)
res.dev2 = residuals(M2, type="deviance")
D.stat2 = sum(res.dev2^2)
1-pchisq(G.stat2, df=173-3)
```

```
## [1] 0
```

```r
1-pchisq(D.stat2, df=173-3)
```

```
## [1] 0
```



```r
# compare 2 models
test.stat = M1$deviance - M2$deviance
df = 171-170
pv3 =  1-pchisq(test.stat, df = df)
pv3 # rej, choose larger model
```

```
## [1] 0.004694838
```
The Pearson's Chi-square test statistic is 536.6, deviance test statistic is 559.89. The p-values of both of them are less than 0.05 under 170 df, indicating the fit is not good here. However, the p-value of the deviance test comparing the two models is less than 0.05 under 1 df, indicating that we should choose the larger model over the smaller one.


#### (c) Overdispersion Adjustment


```r
# plot half-normal plot
abs_res = abs(res.pearson2)
plot(qnorm((173+1:173+0.5)/(2*173+1.125)),sort(abs_res))
abline(a=0,b=1)
```

![](HW5_files/figure-latex/unnamed-chunk-4-1.pdf)<!-- --> 

```r
# estimate dispersion parameter
res.p1 = residuals(M2, type='pearson', data=crab)
G1 = sum(res.p1^2)
phi = G1/(173-3)
phi
```

```
## [1] 3.156449
```

```r
# refit with phi
summary(M2, dispersion=phi)
```

```
## 
## Call:
## glm(formula = Sa ~ W + Wt, family = poisson(link = log), data = crab)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -2.9308  -1.9705  -0.5481   0.9700   4.9905  
## 
## Coefficients:
##             Estimate Std. Error z value Pr(>|z|)
## (Intercept) -1.29168    1.59771  -0.808    0.419
## W            0.04590    0.08309   0.552    0.581
## Wt           0.44744    0.28184   1.588    0.112
## 
## (Dispersion parameter for poisson family taken to be 3.156449)
## 
##     Null deviance: 632.79  on 172  degrees of freedom
## Residual deviance: 559.89  on 170  degrees of freedom
## AIC: 921.18
## 
## Number of Fisher Scoring iterations: 6
```

```r
# goodness of fit
pv.od = 1-pchisq(173-3, df=173-3)
```
The dispersion parameter is 3.16 in M2. There is a deviation from the reference line with slope 1, indicating constant over-dispersion. After adjusting for dispersion=3.16, the model parameters are the same as before, but the p-value for the Chi-square test is 0.49 here, greater than 0.05, indicating a good fit.


### QUESTION 2

#### (a) Model Fitting and Interpretation

```r
# fit a poisson log-link model
fish = read.table("HW5-parasite.txt", header = TRUE)
fish = 
  fish %>% 
  mutate(Area = as.factor(Area),
         Year = as.factor(Year)) 
M1.fish = glm(Intensity ~ Area + Year + Length, family=poisson(link=log), data=fish)
summary(M1.fish)
```

```
## 
## Call:
## glm(formula = Intensity ~ Area + Year + Length, family = poisson(link = log), 
##     data = fish)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -9.3632  -2.7158  -2.0142  -0.4731  30.2492  
## 
## Coefficients:
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  2.6431709  0.0542838  48.692  < 2e-16 ***
## Area2       -0.2119557  0.0491691  -4.311 1.63e-05 ***
## Area3       -0.1168602  0.0428296  -2.728  0.00636 ** 
## Area4        1.4049366  0.0356625  39.395  < 2e-16 ***
## Year2000     0.6702801  0.0279823  23.954  < 2e-16 ***
## Year2001    -0.2181393  0.0287535  -7.587 3.29e-14 ***
## Length      -0.0284228  0.0008809 -32.265  < 2e-16 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for poisson family taken to be 1)
## 
##     Null deviance: 25797  on 1190  degrees of freedom
## Residual deviance: 19153  on 1184  degrees of freedom
##   (63 observations deleted due to missingness)
## AIC: 21089
## 
## Number of Fisher Scoring iterations: 7
```
**Interpretation of the model:**

We fitted a poisson log linear model here where:

log(E(Y)) = β0 + β1·X1 + β2·X2 + β3·X3

with Y = Intensity, X1 = Area2, X2 = Area3, X3 = Area3, X4 = Year2000, X5 = Year2001, X5 = Length.

- Here exp(β1) = 0.81, indicating that adjusting for year and length, the expected intensity of parasite in area2 is 0.81 times of that in area1.

- exp(β2) = 0.89, indicating that adjusting for year and length, the expected intensity of parasite in area3 is 0.89 times of that in area1.

- exp(β3) = 4.08, indicating that adjusting for year and length, the expected intensity of parasite in area4 is 4.08 times of that in area1.

- exp(β4) = 1.95, indicating that adjusting for area and length, the expected intensity of parasite in year 2000 is 1.95 times of that in year 1999.

- exp(β5) = 0.8, indicating that adjusting for area and length, the expected intensity of parasite in year 2001 is 0.8 times of that in year 1999.

- exp(β6) = 0.97, indicating that adjusting for area and year, with one unit increase in length, the expected intensity of parasite would decrease by 0.97 times.


#### (b) Goodness of Fit


```r
# goodness of fit
res.pearson.fish = residuals(M1.fish, type="pearson")
G.stat.fish = sum(res.pearson.fish^2)
res.dev.fish = residuals(M1.fish, type="deviance")
D.stat.fish = sum(res.dev.fish^2)
1-pchisq(G.stat.fish, df=1184)
```

```
## [1] 0
```

```r
1-pchisq(D.stat.fish, df=1184)
```

```
## [1] 0
```
**Goodness of fit:**

The Pearson's Chi-square test statistic is \ensuremath{4.216497\times 10^{4}}, deviance test statistic is \ensuremath{1.91528\times 10^{4}}. The p-values of both of them are less than 0.05 under 1184 df, indicating that the fit is not good here.

#### (c) Account for ZIP

```r
# analyze data for ZIP regression
M0 = zeroinfl(Intensity ~ Area + Year + Length, data = fish)
summary(M0)
```

```
## 
## Call:
## zeroinfl(formula = Intensity ~ Area + Year + Length, data = fish)
## 
## Pearson residuals:
##     Min      1Q  Median      3Q     Max 
## -2.1278 -0.8265 -0.5829 -0.1821 25.4837 
## 
## Count model coefficients (poisson with log link):
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  3.8431720  0.0583793  65.831  < 2e-16 ***
## Area2        0.2687838  0.0500467   5.371 7.84e-08 ***
## Area3        0.1463174  0.0439485   3.329 0.000871 ***
## Area4        0.9448070  0.0368342  25.650  < 2e-16 ***
## Year2000     0.3919828  0.0282952  13.853  < 2e-16 ***
## Year2001    -0.0448457  0.0296057  -1.515 0.129831    
## Length      -0.0368067  0.0009747 -37.762  < 2e-16 ***
## 
## Zero-inflation model coefficients (binomial with logit link):
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)  0.552579   0.275762   2.004  0.04509 *  
## Area2        0.718680   0.189552   3.791  0.00015 ***
## Area3        0.657710   0.167402   3.929 8.53e-05 ***
## Area4       -1.022864   0.188201  -5.435 5.48e-08 ***
## Year2000    -0.752121   0.172965  -4.348 1.37e-05 ***
## Year2001     0.456533   0.143962   3.171  0.00152 ** 
## Length      -0.009889   0.004629  -2.136  0.03266 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
## 
## Number of iterations in BFGS optimization: 17 
## Log-likelihood: -6950 on 14 Df
```
**Interpretation of the parameters: **

**1. If the fish are susceptible to parasites:**

- exp(β1) = 1.31, indicating that adjusting for year and length, the expected intensity of parasite in area2 is 1.31 times of that in area1.

- exp(β2) = 1.16, indicating that adjusting for year and length, the expected intensity of parasite in area3 is 1.16 times of that in area1.

- exp(β3) = 2.57, indicating that adjusting for year and length, the expected intensity of parasite in area4 is 2.57 times of that in area1.

- exp(β4) = 1.48, indicating that adjusting for area and length, the expected intensity of parasite in year 2000 is 1.48 times of that in year 1999.

- exp(β5) = 0.96, indicating that adjusting for area and length, the expected intensity of parasite in year 2001 is 0.96 times of that in year 1999.

- exp(β6) = 0.96, indicating that adjusting for area and year, with one unit increase in length, the expected intensity of parasite would decrease by 0.96 times.

**2. Variables of whether fish are susceptible to parasites may depend on:**

- exp(z1) = 2.05, indicating that adjusting for year and length, the odds ratio of having no parasite in area2 is 2.05 times of that in area1.

- exp(z2) = 1.93, indicating that adjusting for year and length, the odds ratio of having no parasite in area3 is 1.93 times of that in area1.

- exp(z3) = 0.36, indicating that adjusting for year and length, the odds ratio of having no parasite in area4 is 0.36 times of that in area1.

- exp(z4) = 0.47, indicating that adjusting for area and length, the odds ratio of having no parasite in 2000 is 0.47 times of that in 1999.

- exp(z5) = 1.58, indicating that adjusting for area and length, the odds ratio of having no parasite in 2001 is 1.58 times of that in 1999.

- exp(z6) = 0.99, indicating that adjusting for year and area, the odds ratio of having no parasite in area2 would decrease by 0.99 times.

