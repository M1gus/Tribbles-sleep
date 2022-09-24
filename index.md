Popovic et al. 2022
================
Yizhou Yu, MRC Toxicology Unit, University of Cambridge, (yzy21 at
mrc-tox.cam.ac.uk)
updated: <i>Sep-24-2022</i>
</h4>

# General information

This pipeline relates to data related to the UK Biobank and human
genetics analyses in our paper. It also contains additional
supplementary materials.

The input files for this analysis pipeline are on the master branch of
this GitHub page (link: <https://github.com/M1gus/Tribbles-sleep/>)

The UK Biobank data is not available in the repository as they require
separate application. Please visit ukbiobank.ac.uk for more information.

# Analyses

## Generate the genetic dataset

Load the snp data downloaded from <https://www.ncbi.nlm.nih.gov/snp>

Note: compiled SNPs from TRIB1, TRIB2, TRIB3

``` r
trib_file_path = "dt/TRIB_snps/"
trib_files = list.files(trib_file_path)

trib_snps = data.frame()
for (i in c(1:length(trib_files))){
  tmp = read.csv(paste0(trib_file_path,trib_files[i]), sep = "\t")
  trib_snps = rbind(trib_snps,tmp)
}
```

``` r
trib_snps$snp = paste0("rs",trib_snps$snp_id)
length(trib_snps$snp)
```

    ## [1] 14676

``` r
length(unique(trib_snps$snp))
```

    ## [1] 14145

Load phenotype data

``` r
library(BGData)
```

    ## Loading required package: BEDMatrix

    ## Loading required package: LinkedMatrix

    ## Loading required package: symDMatrix

    ## 
    ## Attaching package: 'BGData'

    ## The following object is masked from 'package:graphics':
    ## 
    ##     segments

``` r
library(BGLR)
library(data.table) #requires brew update && brew install llvm
library(qqman)
```

    ## 

    ## For example usage please run: vignette('qqman')

    ## 

    ## Citation appreciated but not required:

    ## Turner, (2018). qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. Journal of Open Source Software, 3(25), 731, https://doi.org/10.21105/joss.00731.

    ## 

``` r
#Load each chromosome as a BEDMatrix object and link them by columns to a ColumnLinkedMatrix object:
genoPath = "/home/yizhouyu/raid/UKB_dt_18-4-2020/"
full_ukb_gen_dt <- as.ColumnLinkedMatrix(lapply(c(1:22), function(chrom) {
  BEDMatrix(paste0(genoPath, "ukb_cal_chr", chrom,"_v2.bed"))
}))  
```

    ## Extracting number of samples and rownames from ukb_cal_chr1_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr1_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr2_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr2_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr3_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr3_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr4_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr4_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr5_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr5_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr6_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr6_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr7_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr7_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr8_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr8_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr9_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr9_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr10_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr10_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr11_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr11_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr12_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr12_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr13_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr13_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr14_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr14_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr15_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr15_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr16_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr16_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr17_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr17_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr18_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr18_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr19_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr19_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr20_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr20_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr21_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr21_v2.bim...

    ## Extracting number of samples and rownames from ukb_cal_chr22_v2.fam...

    ## Extracting number of variants and colnames from ukb_cal_chr22_v2.bim...

``` r
rownames(full_ukb_gen_dt) <- sapply(strsplit(rownames(full_ukb_gen_dt), "_"), `[`, 1) 
# convert names from FID_II D to eid
```

``` r
ukbSNP_without_variant = gsub("_.*","",colnames(full_ukb_gen_dt))
trib_snp_intersect = intersect(unique(trib_snps$snp),ukbSNP_without_variant)

#get the index values 
trib_snp_intersect_df = data.frame(row_num = 1:length(colnames(full_ukb_gen_dt)), 
                                          raw_var = colnames(full_ukb_gen_dt), 
                                          SNP_only = ukbSNP_without_variant, 
                                          selected_var = ukbSNP_without_variant %in% trib_snp_intersect)

trib_snp_intersect_df = subset(trib_snp_intersect_df, selected_var == TRUE)

trib_snps = merge(trib_snps[!duplicated(trib_snps$snp),], trib_snp_intersect_df,
                              by.x="snp", by.y = "SNP_only")

write.csv(trib_snps,"dt_out/trib_snp_in_ukbdt.csv", row.names = FALSE)

BGD_trib_snps = full_ukb_gen_dt[,na.omit(trib_snps$raw_var)]
write.csv(BGD_trib_snps, "dt_out/ukb_trib_geno_dt.csv")
```

### Curate the phenotype data

``` r
ukb_pheno = read.csv("dt/ukb_pheno_curated.csv")
colnames(ukb_pheno)
```

    ##  [1] "eid"           "sleep_dur"     "wake"          "chronotype"   
    ##  [5] "nap"           "insomnia"      "daysleep"      "townsend"     
    ##  [9] "ethnicity"     "sex"           "bpdias"        "bpsys"        
    ## [13] "edu_level"     "folate"        "cog_rt"        "cog_pal"      
    ## [17] "cog_symbol"    "cog_numeric"   "cog_tmta"      "cog_tmtb"     
    ## [21] "cog_matrix"    "peri_greyM"    "greyM"         "brainV"       
    ## [25] "left_hippVo"   "right_hippVo"  "left_hippGre"  "right_hippGre"
    ## [29] "ad_diag"       "age"           "whr"           "log_normCRP"  
    ## [33] "log_normVitaD"

``` r
ukb_pheno$ethnicity = ifelse(ukb_pheno$ethnicity == "British", "British", "Not British")
ukb_pheno = subset(ukb_pheno, select = c(
  eid, age, townsend, ethnicity, sex, edu_level, bpdias, bpsys, 
  whr, sleep_dur, nap, daysleep,ad_diag
))
```

### Baseline models

``` r
whr.baseline = lm(whr ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + sleep_dur, data = ukb_pheno)

MASS::stepAIC(whr.baseline, direction = "both")
```

    ## Start:  AIC=-2553205
    ## whr ~ age + townsend + ethnicity + sex + edu_level + bpdias + 
    ##     bpsys + sleep_dur
    ## 
    ##             Df Sum of Sq    RSS      AIC
    ## <none>                   2007.1 -2553205
    ## - sleep_dur  1      0.03 2007.1 -2553200
    ## - bpsys      1      0.06 2007.2 -2553193
    ## - ethnicity  1      1.80 2008.9 -2552786
    ## - townsend   1     15.08 2022.2 -2549702
    ## - edu_level  7     21.36 2028.5 -2548262
    ## - bpdias     1     29.94 2037.1 -2546273
    ## - age        1     49.69 2056.8 -2541755
    ## - sex        1   1434.01 3441.1 -2300758

    ## 
    ## Call:
    ## lm(formula = whr ~ age + townsend + ethnicity + sex + edu_level + 
    ##     bpdias + bpsys + sleep_dur, data = ukb_pheno)
    ## 
    ## Coefficients:
    ##                                                      (Intercept)  
    ##                                                        6.525e-01  
    ##                                                              age  
    ##                                                        1.391e-03  
    ##                                                         townsend  
    ##                                                        1.932e-03  
    ##                                             ethnicityNot British  
    ##                                                        6.344e-03  
    ##                                                          sexMale  
    ##                                                        1.133e-01  
    ##                            edu_levelCollege or University degree  
    ##                                                       -6.255e-03  
    ##                                      edu_levelCSEs or equivalent  
    ##                                                        7.619e-03  
    ##                                       edu_levelNone of the above  
    ##                                                        1.604e-02  
    ##                         edu_levelNVQ or HND or HNC or equivalent  
    ##                                                        7.338e-03  
    ##                            edu_levelO levels/GCSEs or equivalent  
    ##                                                        1.264e-04  
    ## edu_levelOther professional qualifications eg: nursing, teaching  
    ##                                                       -5.811e-05  
    ##                                    edu_levelPrefer not to answer  
    ##                                                        1.294e-02  
    ##                                                           bpdias  
    ##                                                        1.088e-03  
    ##                                                            bpsys  
    ##                                                       -2.813e-05  
    ##                                                        sleep_dur  
    ##                                                        2.244e-04

For WHR, all of the parameters were significant

sleep\_dur, nap, daysleep

``` r
sleep_dur.baseline = lm(sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys+whr, data = ukb_pheno)

MASS::stepAIC(sleep_dur.baseline, direction = "both")
```

    ## Start:  AIC=99927.98
    ## sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpdias + 
    ##     bpsys + whr
    ## 
    ##             Df Sum of Sq    RSS    AIC
    ## - bpdias     1      0.70 579634  99927
    ## <none>                   579633  99928
    ## - bpsys      1      4.92 579638  99930
    ## - whr        1      8.43 579642  99933
    ## - sex        1     63.20 579697  99977
    ## - edu_level  7    131.97 579765 100021
    ## - ethnicity  1    328.55 579962 100191
    ## - townsend   1    946.21 580580 100690
    ## - age        1   1180.42 580814 100879
    ## 
    ## Step:  AIC=99926.55
    ## sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpsys + 
    ##     whr
    ## 
    ##             Df Sum of Sq    RSS    AIC
    ## <none>                   579634  99927
    ## + bpdias     1      0.70 579633  99928
    ## - bpsys      1      5.13 579639  99929
    ## - whr        1      9.16 579643  99932
    ## - sex        1     63.46 579698  99976
    ## - edu_level  7    131.92 579766 100019
    ## - ethnicity  1    328.31 579963 100190
    ## - townsend   1    946.32 580581 100688
    ## - age        1   1269.21 580903 100949

    ## 
    ## Call:
    ## lm(formula = sleep_dur ~ age + townsend + ethnicity + sex + edu_level + 
    ##     bpsys + whr, data = ukb_pheno)
    ## 
    ## Coefficients:
    ##                                                      (Intercept)  
    ##                                                        6.7277021  
    ##                                                              age  
    ##                                                        0.0068063  
    ##                                                         townsend  
    ##                                                       -0.0153509  
    ##                                             ethnicityNot British  
    ##                                                       -0.0856389  
    ##                                                          sexMale  
    ##                                                       -0.0311951  
    ##                            edu_levelCollege or University degree  
    ##                                                        0.0418780  
    ##                                      edu_levelCSEs or equivalent  
    ##                                                       -0.0084888  
    ##                                       edu_levelNone of the above  
    ##                                                        0.0042242  
    ##                         edu_levelNVQ or HND or HNC or equivalent  
    ##                                                       -0.0201886  
    ##                            edu_levelO levels/GCSEs or equivalent  
    ##                                                       -0.0026132  
    ## edu_levelOther professional qualifications eg: nursing, teaching  
    ##                                                       -0.0031313  
    ##                                    edu_levelPrefer not to answer  
    ##                                                       -0.0328688  
    ##                                                            bpsys  
    ##                                                       -0.0001884  
    ##                                                              whr  
    ##                                                        0.0670571

### Prepare models

``` r
BGD_trib_snps_df = as.data.frame(BGD_trib_snps)
snp_list = colnames(BGD_trib_snps_df)
BGD_trib_snps_df$eid = row.names(BGD_trib_snps_df)
ukb_pheno_geno = merge(ukb_pheno, BGD_trib_snps_df, all.x = TRUE)
```

### WHR model

``` r
whr_lm_df = data.frame()

for (i in 1:length(snp_list)){
  formula_whr = paste0("whr ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + ",snp_list[i])
  #print(formula_whr)
  whr_lm = lm(formula = formula_whr, data = ukb_pheno_geno)
  whr_lm.name = snp_list[i]
  
  whr_lm.coef = as.data.frame(summary(whr_lm)$coefficients)$Estimate[length(as.data.frame(summary(whr_lm)$coefficients)$Estimate)]
  
  whr_lm.025 = as.data.frame(confint(whr_lm)[nrow(confint(whr_lm)),1])
  whr_lm.975 = as.data.frame(confint(whr_lm)[nrow(confint(whr_lm)),2])
  
  whr_lm.p = as.data.frame(summary(whr_lm)$coefficients)$`Pr(>|t|)`[length(as.data.frame(summary(whr_lm)$coefficients)$`Pr(>|t|)`)]
  
  whr_lm_df = rbind(whr_lm_df, unname(c(whr_lm.name,whr_lm.coef,whr_lm.025,whr_lm.975,whr_lm.p)))
}

colnames(whr_lm_df) <- c("snp","coef","025","975","p")
whr_lm_df$fdr = p.adjust(whr_lm_df$p, method = "BH")
```

### Night sleep duration

``` r
nightsleep_lm_df = data.frame()

for (i in 1:length(snp_list)){
  formula_nightsleep = paste0("sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + ",snp_list[i])
  nightsleep_lm = lm(formula = formula_nightsleep, data = ukb_pheno_geno)
  
  nightsleep_lm.name = snp_list[i]
  
  nightsleep_lm.coef = as.data.frame(summary(nightsleep_lm)$coefficients)$Estimate[length(as.data.frame(summary(nightsleep_lm)$coefficients)$Estimate)]
    
  nightsleep_lm.025 = as.data.frame(confint(nightsleep_lm)[nrow(confint(nightsleep_lm)),1])
  nightsleep_lm.975 = as.data.frame(confint(nightsleep_lm)[nrow(confint(nightsleep_lm)),2])
  
  nightsleep_lm.p = as.data.frame(summary(nightsleep_lm)$coefficients)$`Pr(>|t|)`[length(as.data.frame(summary(nightsleep_lm)$coefficients)$`Pr(>|t|)`)]
  
  nightsleep_lm_df = rbind(nightsleep_lm_df, unname(c(nightsleep_lm.name,nightsleep_lm.coef,nightsleep_lm.025,nightsleep_lm.975,nightsleep_lm.p)))
}

colnames(nightsleep_lm_df) <- c("snp","coef","025","975","p")

nightsleep_lm_df$fdr = p.adjust(nightsleep_lm_df$p, method = "BH")
```

### Common significant SNPs between night sleep and WHR

``` r
nightsleep_lm_df.sig = subset(nightsleep_lm_df, fdr <= 0.05)
whr_lm_df.sig = subset(whr_lm_df, fdr <= 0.05)

sleep_whr_snp_noFDR = merge(whr_lm_df, nightsleep_lm_df, by = "snp")

sleep_whr_snp_withFDR = merge(whr_lm_df.sig, nightsleep_lm_df.sig, by = "snp")
```

#### output snp data

``` r
nightsleep_lm_df.sig$phenotype = "sleep"
whr_lm_df.sig$phenotype = "obesity"
trib_out = rbind(whr_lm_df.sig,nightsleep_lm_df.sig)
trib_out$rs_id = gsub("_.*","",trib_out$snp)
trib_out = merge(trib_out, trib_snps, by.x = "rs_id", by.y = "snp", all.x = T)
write.csv(trib_out, "dt_out/trib_snp_analysis_output.csv", row.names = F)
```

### MR approach using significant SNPs between night sleep and WHR

``` r
nightsleep_lm_df_mr = nightsleep_lm_df
colnames(nightsleep_lm_df_mr) <- c("snp","coef.sleep","025.sleep","975.sleep","p.sleep","fdr.sleep")

whr_lm_df_mr = whr_lm_df
colnames(whr_lm_df_mr) <- c("snp","coef.whr","025.whr","975.whr","p.whr","fdr.whr")

sleep_whr_mr_merge = merge(nightsleep_lm_df_mr,
                           whr_lm_df_mr, by = "snp")
```

Quick visualisation

``` r
library(ggpubr)
```

    ## Loading required package: ggplot2

``` r
library(ggplot2)
ggplot(data = sleep_whr_mr_merge, aes(x = coef.whr, y = coef.sleep)) + 
  geom_smooth(method = "lm") + 
  geom_point() +
  stat_cor(label.y = 0.2)+ 
  stat_regline_equation(label.y = 0.15)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
ggsave("fig/linear_regression_with_coefficients_sleep_whr_trib.pdf")
```

    ## Saving 7 x 5 in image

    ## `geom_smooth()` using formula 'y ~ x'

p = 0.11 for lm - this is promosing. I will do a further analysis using
the MendelianRandomisation package

Do another visualisation with more significant snps

``` r
sleep_whr_mr_merge_p70 = subset(sleep_whr_mr_merge, p.whr < 0.7)
sleep_whr_mr_merge_p70 = subset(sleep_whr_mr_merge_p70, p.sleep < 0.7)

sleep_whr_mr_merge_p60 = subset(sleep_whr_mr_merge, p.whr < 0.5)
sleep_whr_mr_merge_p60 = subset(sleep_whr_mr_merge_p60, p.sleep < 0.5)

ggplot(data = sleep_whr_mr_merge_p60, aes(x = coef.whr, y = coef.sleep)) + 
  geom_smooth(method = "lm") + 
  geom_point() +
  stat_cor(label.y = 0.2)+ 
  stat_regline_equation(label.y = 0.15)
```

    ## `geom_smooth()` using formula 'y ~ x'

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
sleep_whr_mr_merge_p05 = subset(sleep_whr_mr_merge, p.whr < 0.05)
sleep_whr_mr_merge_p05 = subset(sleep_whr_mr_merge_p05, p.sleep < 0.05)
sleep_whr_mr_merge_p05
```

    ##            snp   coef.sleep   025.sleep    975.sleep    p.sleep  fdr.sleep
    ## 16 rs2980874_G -0.005880482 -0.01057599 -0.001184973 0.01410479 0.05611286
    ##        coef.whr      025.whr      975.whr        p.whr      fdr.whr
    ## 16 0.0006712768 0.0003952805 0.0009472731 1.870213e-06 3.833936e-05

``` r
library(gridExtra)

pdf("fig/sleep_whr_mr_merge_p05.pdf",width = 12, height = 5)
grid.table(sleep_whr_mr_merge_p05)
dev.off()
```

    ## png 
    ##   2

#### Create dataset for an MRInput object

Re-run all the models and save relevant data

``` r
### WHR model

whr_lm_df_MR = data.frame()

for (i in 1:length(snp_list)){
  formula_whr = paste0("whr ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + ",snp_list[i])
  #print(formula_whr)
  whr_lm = lm(formula = formula_whr, data = ukb_pheno_geno)
  whr_lm.name = snp_list[i]
  
  whr_lm.coef = as.data.frame(summary(whr_lm)$coefficients)$Estimate[length(as.data.frame(summary(whr_lm)$coefficients)$Estimate)]
  
  whr_lm.se = as.data.frame(summary(whr_lm)$coefficients)$"Std. Error"[length(as.data.frame(summary(whr_lm)$coefficients)$"Std. Error")]
  
  whr_lm_df_MR = rbind(whr_lm_df_MR, unname(c(whr_lm.name,whr_lm.coef,whr_lm.se)))
}

colnames(whr_lm_df_MR) <- c("snp","exposure.beta","exposure.se")

nightsleep_lm_df_MR = data.frame()

for (i in 1:length(snp_list)){
  formula_nightsleep = paste0("sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + ",snp_list[i])
  nightsleep_lm = lm(formula = formula_nightsleep, data = ukb_pheno_geno)
  
  nightsleep_lm.name = snp_list[i]
  
  nightsleep_lm.coef = as.data.frame(summary(nightsleep_lm)$coefficients)$Estimate[length(as.data.frame(summary(nightsleep_lm)$coefficients)$Estimate)]
    
  nightsleep_lm.se = as.data.frame(summary(nightsleep_lm)$coefficients)$"Std. Error"[length(as.data.frame(summary(nightsleep_lm)$coefficients)$"Std. Error")]
  
  nightsleep_lm_df_MR = rbind(nightsleep_lm_df_MR, unname(c(nightsleep_lm.name,nightsleep_lm.coef,nightsleep_lm.se)))
}

colnames(nightsleep_lm_df_MR) <- c("snp","outcome.beta","outcome.se")
```

Create MRInputObject

Note: MR assumes that the SNPs are not correlated by dfault, so I need
to provide a correlation matrix. Spearman correlation: Spearman
correlation evaluates the monotonic relationship. The Spearman
correlation coefficient is based on the ranked values for each variable
rather than the raw data.

``` r
ukb_pheno_geno_sleep_whr_mr_merge_p60 = subset(ukb_pheno_geno, select = sleep_whr_mr_merge_p60$snp)

trib_snp_cor <- Hmisc::rcorr(as.matrix(ukb_pheno_geno_sleep_whr_mr_merge_p60),
                             type = "spearman")

corrplot::corrplot(trib_snp_cor$r, method="number")
```

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

Take away rs6076473\_G because its highly correlated with rs6076472\_G

``` r
sleep_whr_MR_merge = merge(nightsleep_lm_df_MR,
                           whr_lm_df_MR, by = "snp")
sleep_whr_MR_merge = sleep_whr_MR_merge[sleep_whr_MR_merge$snp %in% whr_lm_df.sig$snp, ]

trib_snp_cor_rho = trib_snp_cor$r

trib_snp_cor_rho[is.na(trib_snp_cor_rho)] <- 0

library(MendelianRandomization)
sleep_whr_mrinput_weak = mr_input(bx = as.numeric(sleep_whr_MR_merge$exposure.beta),
                             bxse = as.numeric(sleep_whr_MR_merge$exposure.se),
                             by = as.numeric(sleep_whr_MR_merge$outcome.beta),
                             byse = as.numeric(sleep_whr_MR_merge$outcome.se),
                             #correlation = trib_snp_cor_rho,
                             exposure = "WHR",
                             outcome = "Sleep duration",
                             snps = sleep_whr_MR_merge$snp)
```

#### Run MR models with MendelianRandomisation

Inverse variance weighted

``` r
sleep_whr_mrinput_weak_mrmedian = MendelianRandomization::mr_median(sleep_whr_mrinput_weak, weighting = "weighted")
sleep_whr_mrinput_weak_mrmedian
```

    ## 
    ##  Weighted median method 
    ## 
    ## Number of Variants : 3 
    ## ------------------------------------------------------------------
    ##                  Method Estimate Std Error   95% CI        p-value
    ##  Weighted median method   -6.644     2.824 -12.179, -1.109   0.019
    ## ------------------------------------------------------------------

``` r
mr_plot(mr_allmethods(sleep_whr_mrinput_weak),error = TRUE, labels = TRUE)
```

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

``` r
ggsave("fig/mrPlot_sleep_whr_trib.pdf")
```

    ## Saving 7 x 5 in image

Calculate wald ratio for rs2980874(G)

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:BGData':
    ## 
    ##     summarize

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
ukb_pheno_geno_scaled = subset(ukb_pheno_geno, select = c(whr, sleep_dur, age, townsend, ethnicity, sex, edu_level, bpdias, bpsys, rs2980874_G))
ukb_pheno_geno_scaled = ukb_pheno_geno_scaled %>%
    mutate_if(is.numeric, scale)

whr_lm_rs2980874 = lm(formula = "whr ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + rs2980874_G", data = ukb_pheno_geno_scaled)
whr_lm_rs2980874.sum = summary(whr_lm_rs2980874)$coefficients

whr_lm_rs2980874.sum_cleaned = whr_lm_rs2980874.sum["rs2980874_G",]

sleep_lm_rs2980874 = lm(formula = "sleep_dur ~ age + townsend + ethnicity + sex + edu_level + bpdias + bpsys + rs2980874_G", data = ukb_pheno_geno_scaled)
sleep_lm_rs2980874.sum = summary(sleep_lm_rs2980874)$coefficients
sleep_lm_rs2980874.sum_cleaned = sleep_lm_rs2980874.sum["rs2980874_G",]
```

``` r
library(TwoSampleMR)
```

    ## TwoSampleMR version 0.5.6 
    ## [>] New: Option to use non-European LD reference panels for clumping etc
    ## [>] Some studies temporarily quarantined to verify effect allele
    ## [>] See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details

    ## 
    ## Attaching package: 'TwoSampleMR'

    ## The following objects are masked from 'package:MendelianRandomization':
    ## 
    ##     mr_ivw, mr_median

``` r
rs2980874_wald = 
mr_wald_ratio(b_exp= as.numeric(whr_lm_rs2980874.sum_cleaned["Estimate"]), 
                           b_out= as.numeric(sleep_lm_rs2980874.sum_cleaned["Estimate"]), 
                           se_exp= as.numeric(whr_lm_rs2980874.sum_cleaned["Std. Error"]), 
                           se_out= as.numeric(sleep_lm_rs2980874.sum_cleaned["Std. Error"]))
```

``` r
generate_odds_ratios(rs2980874_wald)
```

    ## $b
    ## [1] -0.7029324
    ## 
    ## $se
    ## [1] 0.2863743
    ## 
    ## $pval
    ## [1] 0.01410442
    ## 
    ## $nsnp
    ## [1] 1
    ## 
    ## $lo_ci
    ## [1] -1.264226
    ## 
    ## $up_ci
    ## [1] -0.1416387
    ## 
    ## $or
    ## [1] 0.4951313
    ## 
    ## $or_lci95
    ## [1] 0.2824578
    ## 
    ## $or_uci95
    ## [1] 0.8679348

``` r
library(TwoSampleMR)

whr_sleep_ivw = 
mr_wald_ratio(b_exp= as.numeric(whr_lm_rs2980874.sum_cleaned["Estimate"]), 
                           b_out= as.numeric(sleep_lm_rs2980874.sum_cleaned["Estimate"]), 
                           se_exp= as.numeric(whr_lm_rs2980874.sum_cleaned["Std. Error"]), 
                           se_out= as.numeric(sleep_lm_rs2980874.sum_cleaned["Std. Error"]))
```

``` r
generate_odds_ratios(rs2980874_wald)
```

    ## $b
    ## [1] -0.7029324
    ## 
    ## $se
    ## [1] 0.2863743
    ## 
    ## $pval
    ## [1] 0.01410442
    ## 
    ## $nsnp
    ## [1] 1
    ## 
    ## $lo_ci
    ## [1] -1.264226
    ## 
    ## $up_ci
    ## [1] -0.1416387
    ## 
    ## $or
    ## [1] 0.4951313
    ## 
    ## $or_lci95
    ## [1] 0.2824578
    ## 
    ## $or_uci95
    ## [1] 0.8679348

#### output figures

For the figure:

plot a forest plot of the absolute value of the coefficients +
confidence intervals (one for sleep and another one for WHR) –&gt; this
will show the effect size.

##### WHR figure

``` r
library(ggplot2)
ggplot(whr_lm_df.sig, 
       aes(x=reorder(snp,-coef), y=abs(coef))) + 
  geom_errorbar(aes(ymin=abs(`025`), ymax=abs(`975`)),
                width=0,                    # Width of the error bars
                position=position_dodge(.9), color = "grey", size = 1) +
  geom_point(shape=21, size = 2, color = 'red', fill = "red") +
  theme_classic() + 
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip()+ylab("Coefficients") + xlab("")
```

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
ggsave('fig/WHR_TRIB.pdf', width = 4, height = 2)
```

##### sleep figure

``` r
library(ggplot2)
ggplot(nightsleep_lm_df.sig, 
       aes(x=reorder(snp,coef), y=abs(coef))) + 
  geom_errorbar(aes(ymin=abs(`025`), ymax=abs(`975`)),
                width=0,                    # Width of the error bars
                position=position_dodge(.9), color = "grey", size = 1) +
  geom_point(shape=21, size = 2, color = 'red', fill = "red") +
  theme_classic() + 
  geom_hline(yintercept = 0, linetype="dotted") +
  coord_flip()+ylab("Coefficients") + xlab("")
```

![](TRIB_UKB_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
ggsave('fig/sleep_TRIB.pdf', width = 4, height = 2)
```

##### Figure of all SNPS

``` r
library(gridExtra)
trib_table = readxl::read_excel("dt_out/trib_snp_analysis_output_forRP.xlsx", sheet = "figure")
pdf("fig/trib_snp_table.pdf",width = 12, height = 5)
grid.table(trib_table)
dev.off()
```

    ## png 
    ##   2

## Descriptive data

``` r
#library(gtsummary)
#library(dplyr)
#library(gtsummary)
#ukb_pheno_na = merge(ukb_pheno, BGD_trib_snps_df)
#ukb_pheno_na = na.omit(ukb_pheno_na[,2:10])
#ukb_pheno_na %>% tbl_summary() %>% bold_labels()%>%
#  as_gt() %>%
#  gt::gtsave(filename = "fig/descrip_data.pdf")
```

## TRIB gene expression

``` r
library(TwoSampleMR)
library(MRInstruments)
load("dt/gtex_eqtl.RData")
#head(gtex_eqtl)

TRIB1_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    gene_name == "TRIB1" 
  ))
```

    ## Separating the entries into the following phenotypes:
    ## TRIB1 (Skin Not Sun Exposed Suprapubic)
    ## TRIB1 (Thyroid)

    ## Warning in format_data(gtex_eqtl_subset, type = type, phenotype_col = type, : The following columns are not present but are helpful for harmonisation
    ## eaf

    ## Inferring p-values

``` r
TRIB2_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    gene_name == "TRIB2" 
  ))
```

    ## Separating the entries into the following phenotypes:
    ## TRIB2 (Adipose Subcutaneous)
    ## TRIB2 (Adipose Visceral Omentum)
    ## TRIB2 (Artery Aorta)
    ## TRIB2 (Colon Transverse)
    ## TRIB2 (Esophagus Muscularis)
    ## TRIB2 (Heart Left Ventricle)
    ## TRIB2 (Lung)
    ## TRIB2 (Muscle Skeletal)
    ## TRIB2 (Nerve Tibial)
    ## TRIB2 (Skin Not Sun Exposed Suprapubic)
    ## TRIB2 (Thyroid)
    ## TRIB2 (Whole Blood)

    ## Warning in format_data(gtex_eqtl_subset, type = type, phenotype_col = type, : The following columns are not present but are helpful for harmonisation
    ## eaf

    ## Inferring p-values

``` r
TRIB3_exp_dat <-
  format_gtex_eqtl(subset(
    gtex_eqtl,
    gene_name == "TRIB3"
  ))
```

    ## Separating the entries into the following phenotypes:
    ## TRIB3 (Adipose Subcutaneous)
    ## TRIB3 (Adrenal Gland)
    ## TRIB3 (Artery Aorta)
    ## TRIB3 (Artery Coronary)
    ## TRIB3 (Artery Tibial)
    ## TRIB3 (Brain Caudate basal ganglia)
    ## TRIB3 (Brain Cerebellum)
    ## TRIB3 (Brain Cortex)
    ## TRIB3 (Brain Nucleus accumbens basal ganglia)
    ## TRIB3 (Brain Putamen basal ganglia)
    ## TRIB3 (Colon Sigmoid)
    ## TRIB3 (Esophagus Gastroesophageal Junction)
    ## TRIB3 (Esophagus Muscularis)
    ## TRIB3 (Heart Atrial Appendage)
    ## TRIB3 (Heart Left Ventricle)
    ## TRIB3 (Muscle Skeletal)
    ## TRIB3 (Nerve Tibial)
    ## TRIB3 (Ovary)
    ## TRIB3 (Testis)
    ## TRIB3 (Thyroid)
    ## TRIB3 (Whole Blood)

    ## Warning in format_data(gtex_eqtl_subset, type = type, phenotype_col = type, : The following columns are not present but are helpful for harmonisation
    ## eaf

    ## Inferring p-values

``` r
ao <- available_outcomes()
```

    ## API: public: http://gwas-api.mrcieu.ac.uk/

``` r
ao_EUR = subset(ao, population == "European")
# overweight:   ieu-a-93 or     weight: ukb-b-11842
# sleep duration: ukb-b-4424
```

``` r
all_trib_exposure = rbind(TRIB1_exp_dat,TRIB2_exp_dat,TRIB3_exp_dat)

overweight_out_dat <- extract_outcome_data(
    snps = TRIB2_exp_dat$SNP,
    outcomes = 'ieu-b-40'
)
```

    ## Extracting data for 8 SNP(s) from 1 GWAS(s)

    ## Finding proxies for 4 SNPs in outcome ieu-b-40

    ## Extracting data for 4 SNP(s) from 1 GWAS(s)

``` r
overweight_harmonise_dat <- harmonise_data(
    exposure_dat = TRIB2_exp_dat, 
    outcome_dat = overweight_out_dat
)
```

    ## Harmonising TRIB2 (Artery Aorta) (zm6Tma) and body mass index || id:ieu-b-40 (ieu-b-40)

    ## Harmonising TRIB2 (Heart Left Ventricle) (davY4h) and body mass index || id:ieu-b-40 (ieu-b-40)

    ## Removing the following SNPs for being palindromic with intermediate allele frequencies:
    ## rs1544856

    ## Harmonising TRIB2 (Thyroid) (2e0hDe) and body mass index || id:ieu-b-40 (ieu-b-40)

    ## Harmonising TRIB2 (Colon Transverse) (xsVxMk) and body mass index || id:ieu-b-40 (ieu-b-40)

    ## Harmonising TRIB2 (Whole Blood) (iCLR4M) and body mass index || id:ieu-b-40 (ieu-b-40)

``` r
overweight_mr_dat = mr(overweight_harmonise_dat)
```

    ## Analysing '2e0hDe' on 'ieu-b-40'

    ## No SNPs available for MR analysis of 'davY4h' on 'ieu-b-40'

    ## Analysing 'iCLR4M' on 'ieu-b-40'

    ## Analysing 'xsVxMk' on 'ieu-b-40'

    ## Analysing 'zm6Tma' on 'ieu-b-40'

``` r
subset(overweight_mr_dat, pval < 0.06)
```

    ##   id.exposure id.outcome                        outcome        exposure
    ## 1      2e0hDe   ieu-b-40 body mass index || id:ieu-b-40 TRIB2 (Thyroid)
    ##       method nsnp            b          se       pval
    ## 1 Wald ratio    1 -0.009739874 0.004988716 0.05089299

``` r
generate_odds_ratios(subset(overweight_mr_dat, pval < 0.06))
```

    ##   id.exposure id.outcome                        outcome        exposure
    ## 1      2e0hDe   ieu-b-40 body mass index || id:ieu-b-40 TRIB2 (Thyroid)
    ##       method nsnp            b          se       pval       lo_ci        up_ci
    ## 1 Wald ratio    1 -0.009739874 0.004988716 0.05089299 -0.01951776 3.800926e-05
    ##          or  or_lci95 or_uci95
    ## 1 0.9903074 0.9806715 1.000038

``` r
sleep_out_dat <- extract_outcome_data(
    snps = TRIB2_exp_dat$SNP,
    outcomes = 'ebi-a-GCST003839'
)
```

    ## Extracting data for 8 SNP(s) from 1 GWAS(s)

``` r
sleep_harmonise_dat <- harmonise_data(
    exposure_dat = TRIB2_exp_dat, 
    outcome_dat = sleep_out_dat
)
```

    ## Harmonising TRIB2 (Lung) (HstCE2) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Artery Aorta) (zm6Tma) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Heart Left Ventricle) (davY4h) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Removing the following SNPs for being palindromic with intermediate allele frequencies:
    ## rs1544856

    ## Harmonising TRIB2 (Thyroid) (2e0hDe) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Colon Transverse) (xsVxMk) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Whole Blood) (iCLR4M) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Skin Not Sun Exposed Suprapubic) (JiIvjb) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Adipose Subcutaneous) (JNKtbc) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Adipose Visceral Omentum) (y3iePg) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Muscle Skeletal) (nMIeBS) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Esophagus Muscularis) (Ztwyoa) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

    ## Harmonising TRIB2 (Nerve Tibial) (uJy4ue) and Sleep duration || id:ebi-a-GCST003839 (ebi-a-GCST003839)

``` r
sleep_mr_dat = mr(sleep_harmonise_dat)
```

    ## Analysing '2e0hDe' on 'ebi-a-GCST003839'

    ## No SNPs available for MR analysis of 'davY4h' on 'ebi-a-GCST003839'

    ## Analysing 'HstCE2' on 'ebi-a-GCST003839'

    ## Analysing 'iCLR4M' on 'ebi-a-GCST003839'

    ## Analysing 'JiIvjb' on 'ebi-a-GCST003839'

    ## Analysing 'JNKtbc' on 'ebi-a-GCST003839'

    ## Analysing 'nMIeBS' on 'ebi-a-GCST003839'

    ## Analysing 'uJy4ue' on 'ebi-a-GCST003839'

    ## Analysing 'xsVxMk' on 'ebi-a-GCST003839'

    ## Analysing 'y3iePg' on 'ebi-a-GCST003839'

    ## Analysing 'zm6Tma' on 'ebi-a-GCST003839'

    ## Analysing 'Ztwyoa' on 'ebi-a-GCST003839'

``` r
subset(sleep_mr_dat, pval < 0.06)
```

    ##   id.exposure       id.outcome                               outcome
    ## 8      xsVxMk ebi-a-GCST003839 Sleep duration || id:ebi-a-GCST003839
    ##                   exposure     method nsnp          b         se       pval
    ## 8 TRIB2 (Colon Transverse) Wald ratio    1 0.04758451 0.02307697 0.03920859

``` r
generate_odds_ratios(subset(sleep_mr_dat, pval < 0.06))
```

    ##   id.exposure       id.outcome                               outcome
    ## 8      xsVxMk ebi-a-GCST003839 Sleep duration || id:ebi-a-GCST003839
    ##                   exposure     method nsnp          b         se       pval
    ## 8 TRIB2 (Colon Transverse) Wald ratio    1 0.04758451 0.02307697 0.03920859
    ##         lo_ci      up_ci       or or_lci95 or_uci95
    ## 8 0.002353644 0.09281538 1.048735 1.002356 1.097259

Mendelian randomisation analysis variants correlated with increased TRB2
expression, and their effect on weight or sleep duration. The SNP
instrument used for increased TRIB2 in thyroid and weight is rs17390839
(A). The SNP instrument used for increased TRIB2 in transverse colon and
sleep duration is rs72773697 (A).
