##################### load data and packages #################

.libPaths("C:/Users/Missy Dreier/Documents") # set library path

setwd("C:/Users/Missy Dreier/Box/The Hamilton Lab/Ongoing Studies/PLUS2/Data/Final Datasets/EMA/Raw Data Merged")

data <- read.csv("PLUS-2 EMA Long Format.csv") # load main dataset

# baseline data

BL <- read.csv("C:/Users/Missy Dreier/Box/The Hamilton Lab/Ongoing Studies/PLUS2/Data/Final Datasets/Baseline and CSSRS/Raw Data/PLUS-2_Baseline_Survey.csv")

# CSSRS

cssrs <- read.csv("C:/Users/Missy Dreier/Box/The Hamilton Lab/Ongoing Studies/PLUS2/Data/Final Datasets/Baseline and CSSRS/Raw Data/PLUS-2_Baseline_CSSRS.csv")

# packages

library(tidyverse)
library(MplusAutomation)
library(summarytools)
library(brms)
library(parameters)
library(lme4)
library(lmerTest)
library(naniar)
library(effectsize)
library(sjPlot)

##################################################################################################################
######################################### MPLUS Data preparation ##################################################
##################################################################################################################



data$pid <- as.factor(data$pid) # make ID a factor


## combine date and time to create "hours since study start" for tinterval in MPLUS

data$submit_datetime_str <- paste(data$date, data$Survey.Submitted.Time) # bind date and time

data$submit_datetime <- as.POSIXct(data$submit_datetime_str, format = "%Y-%m-%d %H:%M:%S") # ensure correct date format

## create hours since study start variable

data <- data %>% 
  group_by(pid) %>%
  mutate(study_start_datetime = min(submit_datetime, na.rm=T),
         hours_since_study_start = floor(as.numeric(difftime(submit_datetime, study_start_datetime, units = "hours")))) %>%
  ungroup()

## sort data by pid and chronologically for each participant 

data <- data %>%
  group_by(pid) %>%
  arrange(pid, hours_since_study_start)


data <- data %>%
  filter(pid!="2109") # remove participant with no data

#################### additional cleaning for later ###################################################

########## cleaning of emotion reg vars #################

data$sm_emoreg <- as.numeric(data$sm_emoreg)

data$ERStrategy <- ifelse(data$sm_emoreg==11, NA, data$sm_emoreg) # code "did not have a negative SM experience" as NA

data$ERStrategy <- as.factor(data$ERStrategy)

data$pid <- as.factor(data$pid) # change var name of ID to ID (and make it a factor again)


################### ERSM/ERIRL vars #######################

data$irl_neg_exp <- as.numeric(data$irl_neg_exp) #make numeric
data$irl_pos_exp <- as.numeric(data$irl_pos_exp) #make numeric


## score ERSM negative
data <- data %>%
  mutate(across(starts_with("ERSM_"), ~ as.numeric(as.character(.x)))) # make ERSM vars numeric


items <- data[, c("ERSM_1", "ERSM_2", "ERSM_4", "ERSM_6")]

data$ERSM_sum <- rowSums(items, na.rm = TRUE)
data$ERSM_sum[rowSums(!is.na(items)) == 0] <- NA # imputes NA if none were answered

ERSM <- data %>% select("ERSM_1", "ERSM_2", "ERSM_4", "ERSM_6", "ERSM_sum")
#View(ERSM)

### score positive ERSM vars

items <- data[, c("ERSM_3", "ERSM_5", "ERSM_7")]

data$ERSM_pos_sum <- rowSums(items, na.rm = TRUE)
data$ERSM_pos_sum[rowSums(!is.na(items)) == 0] <- NA # imputes NA if none were answered



## create global PA variable ###
data$PA_HAPPY <- as.numeric(data$PA_HAPPY)

data$PA_HOPEFUL <- as.numeric(data$PA_HOPEFUL)

data$PA_SATISFIED <- as.numeric(data$PA_SATISFIED)

items <- data[, c("PA_HAPPY", "PA_HOPEFUL", "PA_SATISFIED")]

data$Global_PA <- rowMeans(items, na.rm = TRUE)
data$Global_PA[rowSums(!is.na(items)) == 0] <- NA # imputes NA if none were answered


####################### add week in study var ########################

## dates in correct format
data <- data %>%
  mutate(
    date = ymd(date),
    aa_start_date = ymd(aa_start_date) # parses YYYY-MM-DD safely
  )




data <- data %>%
  ungroup() %>%
  mutate(
    days_since_start = as.numeric(date - aa_start_date),
    week_in_study = if_else(
      days_since_start >= 0,
      floor(days_since_start / 7) + 1,
      NA_real_
    )
  )


### create max days/weeks for checking later ##################
data <- data %>%
  group_by(pid) %>%
  mutate(max_week = max(week_in_study),
         max_days = max(day_in_study))


####################################### test stationarity assumption ##################################################

df2 <- data[complete.cases(data[, c("Depression","day_in_study","pid")]), ] # subset to necessary vars


dep <- lmer(Depression ~ 1 + (1|pid), data = df2, REML = T) # null model

dep_day <- lmer(Depression ~ day_in_study + (1|pid), data = df2, REML = T) # model with day as predictor

s2_0 <- attr(VarCorr(dep), "sc")^2 # obtain variance
s2_1 <- attr(VarCorr(dep_day), "sc")^2 # obtain variance

100 * (s2_0 - s2_1) / s2_0 # test difference = .13% change (negligible, depressed mood does not systematically change throughout study) - assumption = met!




####################################### test ICCs/fit for MLM ##################################################


### write function to obtain ICCs
icc_lmer <- function(data, var, id = "pid") {
  f <- as.formula(paste0(var, " ~ 1 + (1|", id, ")"))
  m <- lmer(f, data = data, REML = TRUE, na.action = na.omit)
  vc <- as.data.frame(VarCorr(m))
  between <- vc$vcov[vc$grp == id]      # random intercept variance
  within  <- vc$vcov[vc$grp == "Residual"]
  icc <- between / (between + within)
  return(icc)
}



## continuous vars for which to obtain ICCs
vars <- c("Depression", "ERSM_sum", "ERSM_pos_sum", "irl_neg_exp", "irl_pos_exp", "Global_PA")

icc_results <- data.frame(
  variable = vars,
  icc = sapply(vars, \(v) icc_lmer(data, v, id = "pid"))
)

icc_results[order(-icc_results$icc), ]

## ICCs between ~ .54-.67 - good fits for MLM 

#######################################MCAR test missing data##################################################

## full dataset
mcar_df <- data %>%
  ungroup() %>%
  select(-c(1:29, 37:52, 64:100)) # subset to question vars asked at all prompts - e.g., removing vars like pid, date, timestamp which technically always have a response and make it look like missingness not random

mcar_df[] <- lapply(mcar_df, as.numeric)


 

mcar_test(mcar_df) # p = .412, missing completely at random
vis_miss(mcar_df) # visualize missingness



###### missingness for 2nd half of study (ER/IRL questions)
## filter dataset

data$pid <- as.numeric(as.character(data$pid))

data_miss <- subset(data, pid >= 2223)

## 
mcar_df <- data_miss %>%
  ungroup() %>%
  select(-c(1:29, 37:52, 64, 65, 70:100)) # subset to question vars asked at all prompts - once ER prompts added


mcar_df[] <- lapply(mcar_df, as.numeric)



mcar_test(mcar_df) # p = .09, missing completely at random
vis_miss(mcar_df) # visualize missingness

############################# 

############### subset data ###########################
mplus <- data %>% select("pid", "hours_since_study_start", "ERSM_sum", "Depression")



####### write data ###############

setwd("C:/Users/Missy Dreier/OneDrive - Rutgers University/Documents/Dissertation/Dissertation Analyses/MPLUS")

write.table(mplus, "mplus_data.dat", sep = " ", na = "-999", col.names = F, row.names = FALSE, quote = F, fileEncoding="UTF-8")



#########################################################################################################################################
################################################ import MPLUS output for further analyses #################################################
#########################################################################################################################################


################################ visual and analyze person-specific effects#####################

setwd("C:/Users/Missy Dreier/OneDrive - Rutgers University/Documents/Dissertation/Dissertation Analyses/MPLUS")


ERSM=readModels("ersm mood_2.25.26.out",what="all") #import the MPlus file

ERSM1=ERSM[["parameters"]][["wilevel.standardized"]][["stdyx.standardized"]] # extract the within-level stdyx standardized parameters

ERSM2=ERSM1[which(ERSM1$paramHeader=="BETA|DEP.ON"),] #To extract the beta parameter (ERSM & Mood)

x=as.numeric(ERSM2$est) # make the estimate numeric

ERSM2 <- ERSM2 %>% filter(cluster != "2139", cluster != "2252",  cluster != "2172", cluster != "2227") # remove case with < 50 observations

################# frequencies ####################

summary(x)
summarytools::freq(ERSM2$sig, order = "default") #Frequencies for significance, to get an overview of how many of the person-specific effects are significant
summarytools::freq(x, order = "default") #Frequencies for effect sizes, to get an overview of all person-specific effect sizes


# create groups based on effect
ERSM2$effect_group <- ifelse(ERSM2$est < -.05, "Inverse", 
                             ifelse(ERSM2$est > .05, "Expected", "Neutral"))


ERSM2$effect_group <- factor(ERSM2$effect_group, levels = c("Inverse", "Neutral", "Expected"))

# see groups

table(ERSM2$effect_group) # 7 inverse, 47 neutral, 73

############### plots ######################

histERSM= hist(x, breaks = 50, plot=FALSE) #Create the histogram

ccat = cut(histERSM$breaks, c(-Inf, -0.06, 0.04, Inf)) #Divide the histogram in three parts, based on the effect sizes (including positive, negative, and no effects)

#Plot the histogram
plot(histERSM, 
     xlim = c(-.67,.67),
     ylim = c(0,30),
     main="",
     las=1,
     xaxt='n',
     labels = FALSE,
     xlab="Person-Specific Beta Values for Negative ERSM and Mood",
     ylab="Number of Participants",
     col=c("brown1","gray","darkseagreen3")[ccat])
abline(v= 0.1523, col = "black", lwd=2)#To add the sample mean of the person-specific effects (beta=.096))
axis(side=1, at=seq(-0.67, 0.67, 0.05))
legend("topright", cex=0.90, box.lty = 0, pch=15, c("Negative (Inverse) Effect (5.5%)", "No Effect (Neutral) (37%)", "Positive (Expected) Effect (57.5%)"), col=c("brown1","gray","darkseagreen3"))

####################### merge onto data ##################################

### pull just necessary params ###
ERSM_params <- ERSM2 %>% select("est", "cluster", "effect_group") %>%
  rename("ERSM_Mood_effect" = "est", "pid" = "cluster")

ERSM_params$pid <- as.factor(ERSM_params$pid) # make ID factor

data$pid <- as.numeric(as.character(data$pid)) # now make ID numeric

data_ER <- subset(data, pid >= 2223) # for ER data, subset to only those with data

data$pid <- as.factor(data$pid) # now make ID a factor

effects <- right_join(data, ERSM_params, by = "pid") # join full dataset

data_ER$pid <- as.factor(data_ER$pid)

effects_ER <- merge(data_ER, ERSM_params, by = "pid") # now just just the ER dataset - include only people in both datasets

#### datasets without inverse group - for sensitivity analyses ####

effects_sens <- subset(effects, effect_group!="Inverse")
effects_ER_sens <- subset(effects_ER, effect_group!="Inverse")


############################## basic data checks ##########################################

### days in study of those included in individual models###

data_days <- data %>%
  group_by(pid) %>%
  summarise(max_days = mean(max_days))

mean(data_days$max_days) # 52.51 (7.5 weeks)
sd(data_days$max_days) # 10.98

min(data_days$max_days) # 3 days (.43 weeks)
max(data_days$max_days) # 64 days (9.14 weeks)

freq(data_days$max_days) # 93.89% have at least 4 weeks of data, 82.44% have 8 weeks or more

table(data$Response.Type) # 15017 submitted prompts, 20,491 total = 73.3 % compliance overall


### days in study of those included in individual models###

mean(effects$max_days) # 55 (7.86 weeks)
min(effects$max_days) # 16 days - 2.29 weeks
max(effects$max_days) # 64 days - 10.86 weeks


#### check randomized prompts

table(effects_ER$randomized_yn) # total = 8790, 52.5% not randomized, 47.5% randomized

completed <- subset(effects_ER, Response.Type == "Submission")

table(completed$randomized_yn) # total = 6744, 53.1% not randomized, 46.9% randomized


### overall endorsements/completion

effects$ERSM_yn <- ifelse(effects$ERSM_sum > 0, 1, 0)
table(effects$ERSM_yn) # 14991 total prompts, 9352 endorsed something or 62.4% of prompts, endorsed ANY ERSM, 45.6% of possible prompts (including missed)



# aren't as good/attractive/accomplished as others
length(which(effects$ERSM_1 > 0)) # 6516 - 31.8% of prompts, 43.5% of submitted prompts, 69.7% of endorsements

# FOMO
length(which(effects$ERSM_2 > 0)) # 7364, 35.9% all prompts, 49.1% of submitted prompts, 78.7% of endorsements

# nervous to post
length(which(effects$ERSM_4 > 0)) # 5526, 27% of all prompts, 36.9% of submitted prompts, 59.1% of endorsements

# sad or hurt bc of negative interaction
length(which(effects$ERSM_6 > 0)) # 3586, 17.5% of all prompts, 23.9% of submitted prompts, 38.3% of endorsements

######### check ER strategy endorsement rates #####################
nrow(effects_ER) # 8790 total prompts
table(effects_ER$Response.Type) #6744 submitted prompts (76.7% compliance)

length(which(!is.na(effects_ER$ERStrategy))) # 1849,  21% of possible prompts endorsed some coping strategy, 27.4% of submitted prompts
length(which(effects_ER$ERStrategy==10)) # 498 = "nothing", so  true strategies, 15% of possible prompts and 20% of submitted prompts

###################################################################################################################
###################################################################################################################
###################################################################################################################
########################################               Aim 2          ###################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################

################# cleaning for Aim 2 ##############################

############ compare baseline characteristics ##################

###### merge BL & CSSRS #########
cssrs <- cssrs %>% rename("pid" = "ID")
BL <- BL %>% rename("pid" = "ID")

cssrs$pid <- as.factor(cssrs$pid)
BL$pid <- as.factor(BL$pid)

effects_cssrs <- left_join(ERSM_params, cssrs, by = "pid")
effects_bl <- left_join(effects_cssrs, BL, by = "pid")



## create scores for measures
effects_bl <- effects_bl %>%
  mutate(
    MASC_total = rowSums(
      select(., "MASC_1":"MASC_39"),
      na.rm = TRUE),
    MFQ_total = rowSums(
      select(., "MFQ_1":"MFQ_33"),
      na.rm = TRUE
    )
    
  )

########################## create wide dataset - 1 row per person ################

effects_bl_wide <- effects_bl %>%
  group_by(pid) %>%
  slice(1) %>%
  ungroup()

### cssrs scoring BXS #

### create yes/no variable for any SI ##

effects_bl_wide <- effects_bl_wide %>%
  mutate(
    cssrs_SI_any = case_when(
      rowSums(select(., matches("^cssrs_SI.*\\.")) == 1, na.rm = TRUE) > 0 ~ 1,
      rowSums(!is.na(select(., matches("^cssrs_SI.*\\.")))) > 0 ~ 0,
      TRUE ~ NA_real_
    )
  )

#make factor & change levels to be readable

effects_bl_wide$cssrs_SI_any <- as.factor(effects_bl_wide$cssrs_SI_any)
levels(effects_bl_wide$cssrs_SI_any) <-  c("No", "Yes")



### create yes/no variable for attempts ###

effects_bl_wide$attempt <- ifelse(effects_bl_wide$SB1_AA==1 | effects_bl_wide$SB1_IA1==1 | 
                                    effects_bl_wide$SB_ABA1==1, "Yes", "No")

##### change NSSI factor levels #######
effects_bl_wide$NSSI_1A <- as.factor(effects_bl_wide$NSSI_1A)
levels(effects_bl_wide$NSSI_1A) <-  c("No", "Yes")


################ create subset dataset for sensitivity analyses ##########


effects_bl_sens <- subset(effects_bl_wide, effect_group!="Inverse")


################################### Baseline Depression ########################################################


cor.test(effects_bl_sens$MFQ_total, effects_bl_sens$ERSM_Mood_effect)

# r = .48, p < .001


ggplot(effects_bl_sens, aes(x = ERSM_Mood_effect, y = MFQ_total)) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "Beta Value",
    y = "MFQ Score",
    title = "Baseline Depression x Beta Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")


########### plot full sample #####################

### create plot dataset for error bars  
plot_df <- effects_bl_wide %>%
  group_by(effect_group) %>%
  summarise(
    mean_mfq = mean(MFQ_total, na.rm = TRUE),
    se = sd(MFQ_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

### order groups as desired
plot_df$effect_group <- factor(plot_df$effect_group, levels = c("Inverse", "Neutral", "Expected"))

ggplot(plot_df, aes(x = effect_group, y = mean_mfq, fill = effect_group)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_mfq - 1.96 * se,
        ymax = mean_mfq + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    Inverse   = "brown1",
    Neutral = "gray",
    Expected   = "darkseagreen3"
  )) +
  labs(
    x = "Effect Group",
    y = "MFQ Score",
    fill = "Effect Group",
    title = "Baseline Depression by Effect Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


################################### Baseline Anxiety ########################################################


cor.test(effects_bl_sens$MASC_total, effects_bl_sens$ERSM_Mood_effect)

# r = .29, p = .002


ggplot(effects_bl_sens, aes(x = ERSM_Mood_effect, y = MASC_total)) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "Beta Value",
    y = "MASC Score",
    title = "Baseline Anxiety x Beta Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")


########### plot full sample #####################

### create plot dataset for error bars  
plot_df <- effects_bl_wide %>%
  group_by(effect_group) %>%
  summarise(
    mean_masc = mean(MASC_total, na.rm = TRUE),
    se = sd(MASC_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

### order groups as desired
plot_df$effect_group <- factor(plot_df$effect_group, levels = c("Inverse", "Neutral", "Expected"))

ggplot(plot_df, aes(x = effect_group, y = mean_masc, fill = effect_group)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_masc - 1.96 * se,
        ymax = mean_masc + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    Inverse   = "brown1",
    Neutral = "gray",
    Expected   = "darkseagreen3"
  )) +
  labs(
    x = "Effect Group",
    y = "MASC Score",
    fill = "Effect Group",
    title = "Baseline Anxiety by Effect Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


#################################### Baseline SI ###########################################


t.test(ERSM_Mood_effect ~ cssrs_SI_any, data = effects_bl_sens)
# t = -4.47, df = 112.61, p < .001

aggregate(ERSM_Mood_effect ~ cssrs_SI_any, data = effects_bl_sens, sd)


### create plot dataset for error bars  
plot_df <- effects_bl_sens %>%
  filter(!is.na(cssrs_SI_any)) %>%
  group_by(cssrs_SI_any) %>%
  summarise(
    mean_effect = mean(ERSM_Mood_effect, na.rm = TRUE),
    se = sd(ERSM_Mood_effect, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(plot_df, aes(x = cssrs_SI_any, y = mean_effect, fill = cssrs_SI_any)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_effect - 1.96 * se,
        ymax = mean_effect + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    No = "gray",
    Yes   = "#227C9D"
  )) +
  labs(
    x = "History of SI (Yes/No)",
    y = "Beta value",
    fill = "Effect Group",
    title = "Beta Value x SI History"
  ) +
  theme_minimal() +
  theme(legend.position = "none")




#### visualize full sample ###
# create table for test & plot
tab <- table(effects_bl_wide$effect_group, effects_bl_wide$cssrs_SI_any)


mosaicplot(
  tab,
  main = "Effect Group × Suicidal Ideation",
  xlab = "Effect Group",
  ylab = "Any Lifetime SI",
  color = TRUE
)


#################################### Baseline Attempt ###########################################



t.test(ERSM_Mood_effect ~ attempt, data = effects_bl_sens)
# not significant - t = -1.66, df = 55.57, p = .10

aggregate(ERSM_Mood_effect ~ attempt, data = effects_bl_sens, sd)


### create plot dataset for error bars  
plot_df <- effects_bl_sens %>%
  filter(!is.na(attempt)) %>%
  group_by(attempt) %>%
  summarise(
    mean_effect = mean(ERSM_Mood_effect, na.rm = TRUE),
    se = sd(ERSM_Mood_effect, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(plot_df, aes(x = attempt, y = mean_effect, fill = attempt)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_effect - 1.96 * se,
        ymax = mean_effect + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    No = "gray",
    Yes   = "#227C9D"
  )) +
  labs(
    x = "History of Suicide Attempt (Yes/No)",
    y = "Beta value",
    title = "Beta Value x Attempt History"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# plot full sample
tab <- table(effects_bl_wide$effect_group, effects_bl_wide$attempt)


mosaicplot(
  tab,
  main = "Effect Group × Suicide Attempt",
  xlab = "Effect Group",
  ylab = "Any Lifetime Attempts",
  color = TRUE
)


#################################### Baseline NSSI ###########################################


######### NSSI
t.test(ERSM_Mood_effect ~ as.factor(NSSI_1A), data = effects_bl_sens) 
# significant difference in BL NSSI - t -4.04, df = 110.13, p < .001
## visualize

aggregate(ERSM_Mood_effect ~ as.factor(NSSI_1A), data = effects_bl_sens, sd) 



### create plot dataset for error bars  
plot_df <- effects_bl_sens %>%
  filter(!is.na(NSSI_1A)) %>%
  group_by(NSSI_1A) %>%
  summarise(
    mean_effect = mean(ERSM_Mood_effect, na.rm = TRUE),
    se = sd(ERSM_Mood_effect, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(plot_df, aes(x = NSSI_1A, y = mean_effect, fill = NSSI_1A)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_effect - 1.96 * se,
        ymax = mean_effect + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    No = "gray",
    Yes   = "#227C9D"
  )) +
  labs(
    x = "History of NSSI (Yes/No)",
    y = "Beta value",
    title = "Beta Value x NSSI History"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


############# plot full sample ###########

tab <- table(effects_bl_wide$effect_group, effects_bl_wide$NSSI_1A)

mosaicplot(
  tab,
  main = "Effect Group × NSSI History",
  xlab = "Effect Group",
  ylab = "NSSI History (Yes/No)",
  color = TRUE
)












################################################################################################################
#################################################################################################################################################
##########################################  Negative and Positive Experiences  ###################################################
################################################################################################################

#### create binary yes/no variable for experiences ###

library(dplyr)

## full sample
effects <- effects %>%
  ungroup() %>%
  mutate(
    ERSM_yn_score = case_when(
      rowSums(!is.na(select(., ERSM_1, ERSM_2, ERSM_4, ERSM_6))) == 0 ~ NA_real_,
      TRUE ~ rowSums(select(., ERSM_1, ERSM_2, ERSM_4, ERSM_6) >= 1, na.rm = TRUE)
    ),
    ERSM_pos_yn_score = case_when(
      rowSums(!is.na(select(., ERSM_3, ERSM_5, ERSM_7))) == 0 ~ NA_real_,
      TRUE ~ rowSums(select(., ERSM_3, ERSM_5, ERSM_7) >= 1, na.rm = TRUE)
    )
  )


# remove outliers
effects_sens <- effects_sens %>%
  ungroup() %>%
  mutate(
    ERSM_yn_score = case_when(
      rowSums(!is.na(select(., ERSM_1, ERSM_2, ERSM_4, ERSM_6))) == 0 ~ NA_real_,
      TRUE ~ rowSums(select(., ERSM_1, ERSM_2, ERSM_4, ERSM_6) >= 1, na.rm = TRUE)
    ),
    ERSM_pos_yn_score = case_when(
      rowSums(!is.na(select(., ERSM_3, ERSM_5, ERSM_7))) == 0 ~ NA_real_,
      TRUE ~ rowSums(select(., ERSM_3, ERSM_5, ERSM_7) >= 1, na.rm = TRUE)
    )
  )


#### create yn vars for in-person exp ####

# full sample
effects_ER$IRL_yn <- (ifelse(effects_ER$irl_neg_exp > 0, 1, 0)) ## negative IRL

effects_ER$IRL_pos_yn <- (ifelse(effects_ER$irl_pos_exp > 0, 1, 0)) ## pos IRL


# remove outliers
effects_ER_sens$IRL_yn <- (ifelse(effects_ER_sens$irl_neg_exp > 0, 1, 0)) ## negative IRL

effects_ER_sens$IRL_pos_yn <- (ifelse(effects_ER_sens$irl_pos_exp > 0, 1, 0)) ## pos IRL


#################################### Online ERSM ##########################################################

###### ----------------- negative y/n -----------------------------####

ERSM_br <- lmer(ERSM_yn_score ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(ERSM_br, show.std = T) # no relationship


###### ----------------- negative ERSM -----------------------------####

ERSM_sum <- lmer(ERSM_sum ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(ERSM_sum, show.std = T) # no relationship


###### ----------------- positive y/n -----------------------------####

ERSM_pos_br <- lmer(ERSM_pos_yn_score ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(ERSM_pos_br, show.std = T) # no relationship


###### ----------------- positive ERSM -----------------------------####

ERSM_pos_sum <- lmer(ERSM_pos_sum ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(ERSM_pos_sum, show.std = T) # no relationship




#################################### IRL Emotional Reactions##########################################################

###### ----------------- negative y/n -----------------------------####

ERIRL_br <- glmer(IRL_yn ~ ERSM_Mood_effect + (1 | pid), family = binomial, data = effects_ER_sens)
tab_model(ERIRL_br, show.std = T) # no relationship


###### ----------------- negative ERIRL -----------------------------####

ERIRL_sum <- lmer(irl_neg_exp ~ ERSM_Mood_effect + (1 | pid), data = effects_ER_sens)
tab_model(ERIRL_sum, show.std = T) # no relationship


###### ----------------- positive y/n -----------------------------####

ERIRL_pos_br <- lmer(IRL_pos_yn ~ ERSM_Mood_effect + (1 | pid), data = effects_ER_sens)
tab_model(ERIRL_pos_br, show.std = T) # no relationship


###### ----------------- positive ERIRL -----------------------------####

ERIRL_pos_sum <- lmer(irl_pos_exp ~ ERSM_Mood_effect + (1 | pid), data = effects_ER_sens)
tab_model(ERIRL_pos_sum, show.std = T) # no relationship



####################################### Negative ERSM + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_ERSM = mean(ERSM_sum, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_ERSM = mean(daily_ERSM, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_ERSM, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_ERSM),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Negative ERSM",
    title = "Daily Negative ERSM Trajectories and Group Average Over Time"
  ) +
  theme_minimal()



####################################### positive ERSM + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_ERSM = mean(ERSM_pos_sum, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_ERSM = mean(daily_ERSM, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_ERSM, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_ERSM),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Positive ERSM",
    title = "Daily Positive ERSM Trajectories and Group Average Over Time"
  ) +
  theme_minimal()


####################################### negative  ERIRL + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_ER_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_ERIRL = mean(irl_neg_exp, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_ERIRL = mean(daily_ERIRL, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_ERIRL, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_ERIRL),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Negative Emotional Reactions to In-Person Negative Experiences",
    title = "Daily Negative IRL Reactions and Group Average Over Time"
  ) +
  theme_minimal()

####################################### positive  ERIRL + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_ER_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_ERIRL = mean(irl_pos_exp, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_ERIRL = mean(daily_ERIRL, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_ERIRL, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_ERIRL),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Positive Emotional Reactions to In-Person Positive Experiences",
    title = "Daily Positive IRL Reactions and Group Average Over Time"
  ) +
  theme_minimal()


##########################################################################################################################
##########################################################################################################################
##################################################### Mood  #############################################
##########################################################################################################################
##########################################################################################################################
##########################################################################################################################


########### negative mood ##################

neg_mood <- lmer(Depression ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(neg_mood, show.std = T) # associated with neg mood, b = 4.90, p < .001


#################### mood volatility ######################


# 1) RMSSD of Depression per person ----
rmssd_df <- effects_sens %>%
  arrange(pid, hours_since_study_start) %>%
  group_by(pid) %>%
  summarise(
    RMSSD_Depression = {
      x <- Depression
      d <- diff(x)
      sqrt(mean(d^2, na.rm = TRUE))
    },
    ERSM_Mood_effect = mean(ERSM_Mood_effect, na.rm=T),
    .groups = "drop"
  )


cor.test(rmssd_df$RMSSD_Depression, rmssd_df$ERSM_Mood_effect)

# r = .70, p < .001






################### visualize #############################


####################################### Depressed mood + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_mood = mean(Depression, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_mood = mean(daily_mood, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_mood, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_mood),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Depressed Mood",
    title = "Daily Depressed Mood Trajectories and Group Average Over Time"
  ) +
  theme_minimal()



############## positive mood ###########################


pos_mood <- lmer(Global_PA ~ ERSM_Mood_effect + (1 | pid), data = effects_sens)
tab_model(pos_mood, show.std = T) # inversely with pos mood, B = -.15, p = .02


###################### positive mood volatility ##############################


# 1) RMSSD of positive mood per person ----
rmssd_df <- effects_sens %>%
  arrange(pid, hours_since_study_start) %>%
  group_by(pid) %>%
  summarise(
    RMSSD_PA = {
      x <- Global_PA
      d <- diff(x)
      sqrt(mean(d^2, na.rm = TRUE))
    },
    ERSM_Mood_effect = mean(ERSM_Mood_effect, na.rm=T),
    .groups = "drop"
  )

cor.test(rmssd_df$RMSSD_PA, rmssd_df$ERSM_Mood_effect)
### also associated with volatility in PA r = .27, p = .003


####################################### positive mood + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  summarise(daily_mood = mean(Global_PA, na.rm = TRUE), .groups = "drop")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_mood = mean(daily_mood, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = daily_mood, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_mood),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Positive Affect",
    title = "Daily Positive Affect Trajectories and Group Average Over Time"
  ) +
  theme_minimal()






#################### daily SI #############################
### clean daily SI var ###


effects_sens$SI_clean = ifelse(effects_sens$SI_Passive==1 | effects_sens$SI_Active==1 | effects_sens$SI_Active==2 | effects_sens$SI_Method==1 | effects_sens$SI_Method==2, 1, 0)



SI <- glmer(as.factor(SI_clean) ~ ERSM_Mood_effect + (1|pid), family = binomial, data = effect_sens)

tab_model(SI, show.std = T) # B = 2.84, p < .001



####################################### positive ERSM + group ####################

# 1) Person-level daily averages (collapsing 3x/day -> 1 value/day/person)
person_daily <- effects_sens %>%
  group_by(effect_group, pid, day_in_study) %>%
  filter(Source.Name=="Evening")

# 2) Group-level daily mean across people
group_daily <- person_daily %>%
  group_by(effect_group, day_in_study) %>%
  summarise(mean_SI = mean(SI_clean, na.rm = TRUE), .groups = "drop")

# 3) Plot: faint individual trajectories + bold group mean, faceted by group
ggplot() +
  geom_line(
    data = person_daily,
    aes(x = day_in_study, y = SI_clean, group = pid),
    alpha = 0.15,
    linewidth = 0.4
  ) +
  geom_line(
    data = group_daily,
    aes(x = day_in_study, y = mean_SI),
    linewidth = 1.2
  ) +
  facet_wrap(~ effect_group) +
  labs(
    x = "Day in Study",
    y = "Suicidal Thoughts Yes/No",
    title = "Daily SI Trajectories and Group Average Over Time"
  ) +
  theme_minimal()

###################################################################################################################
###################################################################################################################
###################################################################################################################
########################################               Aim 3 (ER)          ###################################################
###################################################################################################################
###################################################################################################################
###################################################################################################################



#################### create dummy variables for individual analyses ############

library(fastDummies)

# Strategy names in order
strategy_names <- c(
  "Acceptance",
  "Problem_Solving",
  "Social_support",
  "Reappraisal",
  "Emotional_expression",
  "Rumination",
  "Cognitive_avoidance",
  "Behavioral_avoidance",
  "Something_else",
  "Nothing"
)

# Create dummy variables
effects_ER <- effects_ER %>%
  mutate(SMER = factor(ERStrategy, levels = 1:10, labels = strategy_names)) %>%
  fastDummies::dummy_cols(select_columns = "SMER", remove_first_dummy = FALSE) %>%
  # Rename columns and convert to yes/no
  rename_with(~ paste0(., "_yn"), starts_with("SMER_")) %>%
  mutate(across(ends_with("_yn")))




###### subset
effects_ER_sens <- subset(effects_ER, effect_group!="Inverse")

effects_ER_sens$randomized_yn <- as.factor(effects_ER_sens$randomized_yn)

### individual strategies ###

##----------------------- acceptance ---------------------------##

accept <- glmer(SMER_Acceptance_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(accept, show.std = T) 


##----------------------- behavioral avoidance ---------------------------##

ba <- glmer(SMER_Behavioral_avoidance_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(ba, show.std = T) 


##----------------------- cognitive avoidance ---------------------------##

ca <- glmer(SMER_Cognitive_avoidance_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(ca, show.std = T) 


##----------------------- emotional expression ---------------------------##

ee <- glmer(SMER_Emotional_expression_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(ee, show.std = T) 

### higher beta value associated with emotional expression (p = .023 - prior to correction)



##----------------------- nothing ---------------------------##

nothing <- glmer(SMER_Nothing_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(nothing, show.std = T) 



##----------------------- problem solving ---------------------------##

ps <- glmer(SMER_Problem_Solving_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(ps, show.std = T) #



##----------------------- reappraisal ---------------------------##

reap <- glmer(SMER_Reappraisal_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(reap, show.std = T) 


##----------------------- rumination ---------------------------##

rum <- glmer(SMER_Rumination_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(rum, show.std = T)



##----------------------- social support ---------------------------##

ss <- glmer(SMER_Social_support_yn ~ ERSM_Mood_effect + randomized_yn + (1 | pid), family = binomial, data = effects_ER_sens)

tab_model(ss, show.std = T) 


####################################### adjust p values #####################################


ps <- c(.482, .139, .346, .024, .583, .558, .966, .056, .882) # list of ps

p.adjust(ps, method = "BH", n = length(ps)) # Benjamini-Hochberg correction

# 0.7495714 0.4170000 0.7495714 0.2160000 0.7495714 0.7495714 0.9660000 0.2520000 0.9660000

################################### entropy ############################################
library(dplyr)
library(tidyr)

### 0) Ensure variables are in the right type
effects_ER_flex <- effects_ER %>%
  mutate(
    pid = as.character(pid),
    randomized_yn = as.integer(randomized_yn),   # or factor; but we use 0/1 below
    ERStrategy = factor(
      ERStrategy,
      levels = c(1,2,3,4,5,6,7,8,9,10),
      labels = c(
        "Acceptance",
        "Problem Solving",
        "Social Support",
        "Reappraisal",
        "Emotional Expression",
        "Rumination",
        "Cognitive Avoidance",
        "Behavioral Avoidance",
        "Something Else",
        "Nothing"
      )
    )
  )

### 1) Count unique strategies per person (raw)
unique_counts <- effects_ER_flex %>%
  group_by(pid) %>%
  summarise(
    n_unique_ER_strategies = n_distinct(ERStrategy[!is.na(ERStrategy)]),
    .groups = "drop"
  )

### 2) Create person x strategy counts (this is your missing `er_counts`)
er_counts <- effects_ER_flex %>%
  filter(!is.na(ERStrategy), !is.na(randomized_yn)) %>%
  group_by(pid, randomized_yn, ERStrategy) %>%
  summarise(n = dplyr::n(), .groups = "drop")

### 3) Get overall strategy proportions within randomized_yn groups
# (overall across people, using counts)
p_yes <- er_counts %>%
  filter(randomized_yn == 1) %>%
  group_by(ERStrategy) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(p_yes = n / sum(n)) %>%
  select(ERStrategy, p_yes)

p_no <- er_counts %>%
  filter(randomized_yn == 0) %>%
  group_by(ERStrategy) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(p_no = n / sum(n)) %>%
  select(ERStrategy, p_no)

### 4) Weight map: how common strategy is in randomized vs nonrandomized
# Use a small epsilon to avoid divide-by-zero if a strategy is absent in one group.
eps <- 1e-8

w_map <- full_join(p_no, p_yes, by = "ERStrategy") %>%
  mutate(
    p_no  = tidyr::replace_na(p_no,  0),
    p_yes = tidyr::replace_na(p_yes, 0),
    w = (p_yes + eps) / (p_no + eps)
  ) %>%
  select(ERStrategy, w)

### 5) Apply weights to NONrandomized counts only
er_counts_adj <- er_counts %>%
  left_join(w_map, by = "ERStrategy") %>%
  mutate(
    w = tidyr::replace_na(w, 1),
    n_adj = if_else(randomized_yn == 0, n * w, as.numeric(n))
  )

### 6) Compute entropy from adjusted counts (person-level)
er_entropy_weighted <- er_counts_adj %>%
  group_by(pid, ERStrategy) %>%
  summarise(n_adj = sum(n_adj), .groups = "drop") %>%
  group_by(pid) %>%
  mutate(total_ER = sum(n_adj),
         p = n_adj / total_ER) %>%
  summarise(
    n_ER_unique = dplyr::n(),                # number of strategies used (after dropping NA)
    total_ER = first(total_ER),
    ERStrategy_entropy = -sum(p * log2(p), na.rm = TRUE),
    .groups = "drop"
  )

### 7) Keep ERSM_Mood_effect and assign effect_group (person-level)
effect_summary <- effects_ER_flex %>%
  group_by(pid) %>%
  summarise(
    ERSM_Mood_effect = mean(ERSM_Mood_effect, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    effect_group = case_when(
      ERSM_Mood_effect < -0.05 ~ "Inverse",
      ERSM_Mood_effect >  0.05 ~ "Expected",
      TRUE ~ "Neutral"
    )
  )


person_summary <- left_join(er_entropy_weighted, effect_summary, by = "pid") ### merge onto beta values/groups


ER_summary <- left_join(er_counts_adj, effect_summary, by = "pid") ## do this with ER strategies

ER_sens <- ER_summary %>% filter(effect_group!="Inverse")
## subset to without inverse group ##

person_sub <- person_summary %>% filter(effect_group!="Inverse")

### ensure factor levels are correctly presented
person_summary$effect_group <- factor(person_summary$effect_group, levels = c("Inverse", "Neutral", "Expected"))
ER_summary$effect_group <- factor(ER_summary$effect_group, levels = c("Inverse", "Neutral", "Expected"))


##################################### ER Strategy diversity ######################


cor.test(person_sub$ERStrategy_entropy, person_sub$ERSM_Mood_effect) # r = .18, p = .21


###########################################################################################
###########################################################################################
################################################# plot   ################################
###########################################################################################
###########################################################################################



######### plot strategy use by group ############

alpha_levels <- sort(unique(as.character(ER_sens$ERStrategy)))

plot_df_adj <- ER_sens %>%
  group_by(effect_group) %>%
  mutate(
    total_n_adj_group = sum(n_adj),
    prop_adj = n_adj / total_n_adj_group
  ) %>%
  ungroup() %>%
  mutate(ERStrategy = factor(as.character(ERStrategy), levels = rev(alpha_levels)))

ggplot(plot_df_adj, aes(x = prop_adj, y = ERStrategy)) +
  geom_col() +
  facet_wrap(~ effect_group) +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    x = "Proportion of strategy use within group",
    y = NULL,
    title = "Proportion of emotion regulation strategies by group"
  ) +
  theme_minimal()




### plot flexibility by group ###############

plot_df <- person_summary %>%
  group_by(effect_group) %>%
  summarise(
    mean_entropy = mean(ERStrategy_entropy, na.rm = TRUE),
    se = sd(ERStrategy_entropy, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(plot_df, aes(x = effect_group, y = mean_entropy, fill = effect_group)) +
  geom_col(width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_entropy - 1.96 * se,
        ymax = mean_entropy + 1.96 * se),
    width = 0.2
  ) +
  scale_fill_manual(values = c(
    Inverse   = "brown1",
    Neutral = "gray",
    Expected   = "darkseagreen3"
  )) +
  labs(
    x = "Effect Group",
    y = "ER Strategy Diversity",
    fill = "Effect Group",
    title = "ER Strategy Diversity by Effect Group"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

########### linear plot flexibilit and beta value ############

ggplot(person_sub, aes(x = ERSM_Mood_effect, y = ERStrategy_entropy)) +
  geom_smooth(method = "lm", color = "black") +
  labs(
    x = "Beta Value",
    y = "ER Strategy Diversity",
    title = "ER Strategy Diversity x Beta Value"
  ) +
  theme_bw() +
  theme(legend.position = "none")




