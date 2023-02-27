#load packages
library(tidyverse)      #for tidying and wrangling
library(lme4)           #for GLMM building
library(lmerTest)       #for better output from GLMMs
library(emmeans)        #for our planned comparisons
library(fitdistrplus)   #for Cullen and Frey plots
library(ggpubr)         #for dataviz.
library(kableExtra)     #for table output

#control scientific notation - makes reading p-values easier
options("scipen"=100000, "digits"=10)

#load in pre-tidied and wrangled dataset
exp1_all <- read_csv("exp1_tidied_ALL-PARTICIPANTS_data.csv")

#check participant no. = 170
exp1_all %>%
  distinct(Participant) %>%
  nrow()

#quick check
head(exp1_all)

#filter NONCE control trials
exp1_all <- exp1_all %>%
  filter(!ProbeType == "NONCE")

#factorising
exp1_all <- exp1_all %>%
  mutate(Participant = factor(Participant),
         ProbeType = factor(ProbeType), 
         ProbePosition = factor(ProbePosition), 
         SubjectType = factor(SubjectType),
         Condition = factor(Condition))

#deviation coding factors
contrasts(exp1_all$ProbeType) <- contr.sum(2) #removed NONCE trials as these are controls
contrasts(exp1_all$ProbePosition) <- contr.sum(2)
contrasts(exp1_all$SubjectType) <- contr.sum(2)

#LOOKING AT QNP AND DNP CONDITIONS

#histogram to check distribution - looks like gamma
exp1_all %>%
  ggplot(aes(x = RT)) + #specify what data to show on x axis
  geom_histogram(binwidth = 200) + #produce histogram with bins that are 200ms wide
  labs(x = "RT (ms)", y ="Frequency") + #x and y axes labels
  ggtitle("Hist. of RTs (ms)") + #add title
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centred
        text = element_text(size = 12)) #adjust rest of chart text - size 12

#Cullen and Frey Graph to check distribution more quantitatively - definitely gamma
descdist(exp1_all$RT,
         discrete = FALSE, #we are not dealing with discrete data - it's continuous RT data
         method = "unbiased", #don't use bootstrapping for sampling - use the exact data without replacement
         graph = TRUE, #produce a graph!
         obs.col = "red", #the colour of the "observation" marker - red in this case
         obs.pch = 15) #the shape of the "observation" marker - square in this case

#summary stats. for experimental conditions
exp1_all %>%
  group_by(Condition) %>% #group by Condition
  summarise("MeanRT" = mean(RT), #mean for each condition
            "StdDevRT" = sd(RT)) %>% #SD for each condition
  arrange(desc(MeanRT)) %>% #arrange descending
  kableExtra::kbl() %>%
  kableExtra::kable_classic(full_width = T, html_font = "Arial")

#summary stats. REL vs. UNREL
exp1_all %>%
  group_by(ProbeType) %>% #group by ProbeType
  summarise("MeanRT" = mean(RT), #mean for each group
            "StdDevRT" = sd(RT)) %>% #SD for each group
  arrange(desc(MeanRT)) %>% #arrange descending
  kableExtra::kbl() %>%
  kableExtra::kable_classic(full_width = T, html_font = "Arial")

#QNP-REL-AFTER-V fastest; QNP-UNREL-AFTER-V slowest.

#viz
set.seed(45422)
exp1_all %>%
  mutate(Condition = fct_reorder(Condition, log(RT), .desc = FALSE)) %>% #reorder log of RT by Condition from lowest to highest from bottom up
  ggplot(aes(x = Condition, y = log(RT), colour = SubjectType)) + #specify what data to show on x and y axes, 
  geom_violin() + #produce a violin plot
  geom_jitter(width = .2, alpha = .08) + #add jitter to minimise overlap of points
  stat_summary(fun.data = mean_se, colour = "black") + 
  guides() + #legend
  labs(x = "ProbeType x ProbePosition x SubjectType", #x axis label
       y = "LogRT (ms)") + #y axis label
  ggtitle("Effect of Probe Type x Probe Position x Subject Type Interaction on RTs") + #add title - \n means 'new line'
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centered
        text = element_text(size = 12)) + #adjust rest of chart text - size 12
  coord_flip() #flip horizontal

#looking at covariates:
#checking for repetition effect/effect of trial number: weak -VE correlation
ggscatter(exp1_all, x = "TrialN", y = "RT", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Trial Number", ylab = "RT (ms)",
          title = "RT-Trial Number Correlation")

#looking at covariates:
#checking for correlation between reading times and RTs: moderate +VE correlation
ggscatter(exp1_all, x = "RT", y = "ReadingTime", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Reading Time (ms)", ylab = "RT (ms)",
          title = "RT-Reading Time Correlation")

#looking at covariates:
#checking for correlation between ProbeFrequency and RTs: weak -VE correlation (NB: not specified in prereg!)
ggscatter(exp1_all, x = "ProbeFrequency", y = "RT", 
          add = "reg.line", #adds regression line
          cor.coef = TRUE, #gives us correlation coefficient
          conf.int = TRUE, #gives us confidence intervals around the line
          add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
          xlab = "Zipf Frequency of Probe Word", ylab = "RT (ms)",
          title = "RT-Frequency of Probe Word Correlation")

#MODEL BUILDING
#gamma model w/ identity link (see Lo & Andrews 2015)
#maximal model: "keep it maximal" (see Bates et al. 2013)
QNPDNP.ALL_PARS.gamma.identity1 <- glmer(RT ~ ProbeType * SubjectType * ProbePosition + #3 factors w/ 3-way interaction
                                  (1+ProbeType * SubjectType * ProbePosition | StimID) + #ranefx
                                  (1+ProbeType * SubjectType * ProbePosition | Participant), #ranefx
                                family = Gamma(link = "identity"), #gamma distribution w/ "identity" link fxn.
                                data = exp1_all)
summary(QNPDNP.ALL_PARS.gamma.identity1) #TO RUN WITH DEFAULT SETTINGS

#reduce ranefx. (see Matuschek et al. 2017; Barr et al. 2013)
QNPDNP.ALL_PARS.gamma.identity2 <- glmer(RT ~ ProbeType * SubjectType * ProbePosition +
                                           (1+ProbeType * ProbePosition | StimID) + #drop SubjectType
                                           (1+ProbeType * SubjectType * ProbePosition | Participant),
                                         family = Gamma(link = "identity"),
                                         data = exp1_all)
summary(QNPDNP.ALL_PARS.gamma.identity2) #singular fit

#reduce ranefx. (see Matuschek et al. 2017; Barr et al. 2013)
QNPDNP.ALL_PARS.gamma.identity3 <- glmer(RT ~ ProbeType * SubjectType * ProbePosition +
                                           (1+ProbeType * ProbePosition | StimID) + 
                                           (1+SubjectType * ProbePosition | Participant), #drop ProbeType
                                         family = Gamma(link = "identity"),
                                         data = exp1_all)

summary(QNPDNP.ALL_PARS.gamma.identity3) #singular fit

#reduce ranefx. (see Matuschek et al. 2017; Barr et al. 2013)
QNPDNP.ALL_PARS.gamma.identity4 <- glmer(RT ~ ProbeType * SubjectType * ProbePosition +
                                           (1+ProbeType | StimID) + #drop ProbePosition
                                           (1+SubjectType * ProbePosition | Participant),
                                         family = Gamma(link = "identity"),
                                         data = exp1_all)

summary(QNPDNP.ALL_PARS.gamma.identity4) #convergence warnings

QNPDNP.ALL_PARS.gamma.identity5 <- glmer(RT ~ ProbeType * SubjectType * ProbePosition +
                                           (1+ProbeType | StimID) +
                                           (1 | Participant), #drop ixn.
                                         family = Gamma(link = "identity"),
                                         data = exp1_all)

summary(QNPDNP.ALL_PARS.gamma.identity5) #converges


#plot it
emmip(QNPDNP.ALL_PARS.gamma.identity5, SubjectType ~ ProbeType | ProbePosition)
emmip(QNPDNP.ALL_PARS.gamma.identity5, ProbePosition ~ ProbeType | SubjectType)

#NEEDS EDITING: significant 3-way ixn. Significant ProbeType:ProbePosition ixn. Significant main effect of ProbeType.
#ixn. plots suggest that in the AFTER-V conditions, RTs are faster in response to REL probes
#when SubjectType is QNP while in the AFTER-VP conditions there's not much difference. There's
#noticeable difference for UNREL probes in the AFTER-VP conditions with UNREL probes responded
#to seemingly slower when the subject is a DNP vs. a QNP. In the AFTER-V conditions there's not much diference.
#Conditioning on SubjectType, a crossover effect is seen for QNP sentences that is not seen for DNP sentences.

#EMMs
emmeans(QNPDNP.ALL_PARS.gamma.identity5, pairwise ~ ProbeType * ProbePosition | SubjectType, adjust = "none")
emmeans(QNPDNP.ALL_PARS.gamma.identity5, pairwise ~ ProbeType, adjust = "none")

#according to predictions, RTs in the QNP-REL-AFTER-VP conditions should be faster than RTs in the
#QNP-REL-AFTER-V conditions, suggesting that the antecedent QNP subject is re-accessed for composition post-VP.
#RTs in the DNP-REL-AFTER-V conditions should be faster than RTs in the DNP-REL-AFTER-VP conditions,
#suggesting that the antecedent DNP subject is re-accessed for composition around the matrix verb.
#EMMs suggest that we find the inverse of our predictions: RTs appear faster in the QNP-REL-AFTER-V conditions
#than in the QNP-REL-AFTER-VP conditions and RTs appear faster in the DNP-REL-AFTER-VP conditions than in the DNP-REL-AFTER-V 
#conditions though not significant. The predictions seem unsupported. We will now build two 2-way models
#as a significant 3-way interaction was observed. Sig. diff for QNP-UNREL-VP vs. V though p = .04 suggesting
#interference effect perhaps...

#LOOKING AT DNP CONDITIONS

#summary stats. for experimental conditions
exp1_all %>%
  filter(SubjectType == "DNP") %>%
  group_by(Condition) %>% #group by Condition
  summarise("MeanRT" = mean(RT), #mean for each condition
            "StdDevRT" = sd(RT)) %>% #SD for each condition
  arrange(desc(MeanRT)) %>% #arrange descending
  kableExtra::kbl(caption = "Mean RTs by Condition") %>%
  kableExtra::kable_classic(full_width = T, html_font = "Arial")

#REL means are the same in AFTER-V and AFTER-VP conditions. Mean RTs faster in UNREL-AFTER-V condition
#than in UNREL-AFTER-VP condition (13ms difference).

#viz
set.seed(45422)
exp1_all %>%
  filter(SubjectType == "DNP") %>%
  mutate(Condition = fct_reorder(Condition, log(RT), .desc = FALSE)) %>% #reorder log of RT by Condition from lowest to highest from bottom up
  ggplot(aes(x = ProbeType:ProbePosition, y = log(RT), colour = ProbeType)) + #specify what data to show on x and y axes, 
  geom_violin() + #produce a violin plot
  geom_jitter(width = .2, alpha = .08) + #add jitter to minimise overlap of points
  stat_summary(fun.data = mean_se, colour = "black") + 
  labs(x = "ProbeType x ProbePosition", #x axis label
       y = "LogRT (ms)") + #y axis label
  ggtitle("Effect of Probe Type x Probe Position Interaction on RTs \n in DNP Conditions") + #add title - \n means 'new line'
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centered
        text = element_text(size = 12)) + #adjust rest of chart text - size 12
  coord_flip() #flip horizontal

#MODEL BUILDING
#gamma model w/ identity link (see Lo & Andrews 2015)
#maximal model: "keep it maximal" (see Bates et al. 2013)
DNP.ALL_PARS.gamma.identity1 <- glmer(RT ~ ProbeType * ProbePosition + #2 factors w/ 2-way interaction
                               (1+ProbeType * ProbePosition | StimID) + #ranefx
                               (1+ProbeType * ProbePosition | Participant), #ranefx
                             family = Gamma(link = "identity"), #gamma distribution w/ "identity" link fxn.
                             data = exp1_all %>% filter(SubjectType == "DNP")) #only DNP condition data
summary(DNP.ALL_PARS.gamma.identity1) #convergence warnings

#maximal model doesn't converge - reduce ranefx. (see Matuschek et al. 2017)
#also increase number of evaluations for good measure to protect against convergence failure
DNP.ALL_PARS.gamma.identity2 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType + ProbePosition | StimID) + #no ixn.
                                        (1+ProbeType * ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity2) #convergence warnings

#simplify ranefx. further
DNP.ALL_PARS.gamma.identity3 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType + ProbePosition | StimID) +
                                        (1+ProbeType + ProbePosition | Participant), #no ixn.
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity3) #convergence warnings

#simplify ranefx. further
DNP.ALL_PARS.gamma.identity4 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType | StimID) + #dropped ProbePosition
                                        (1+ProbeType + ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity4) #convergence warnings

#simplify ranefx. further
DNP.ALL_PARS.gamma.identity5 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) + #dropped ProbeType
                                        (1+ProbeType + ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity5) #convergence warnings

#simplify ranefx. further
DNP.ALL_PARS.gamma.identity6 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) +
                                        (1+ProbeType | Participant), #dropped ProbePosition
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity6) #convergence warnings

#simplify ranefx. further
DNP.ALL_PARS.gamma.identity7 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) +
                                        (1 | Participant), #dropped ProbeType
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "DNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5)))
summary(DNP.ALL_PARS.gamma.identity7) #converges - significant main effect of ProbeType

#plot it
emmip(DNP.ALL_PARS.gamma.identity7, ProbePosition ~ ProbeType)

#EMMs
emmeans(DNP.ALL_PARS.gamma.identity7, pairwise ~ ProbeType * ProbePosition, adjust = "none")

#significant main effect of ProbeType. No main effect of ProbePosition and no ProbeType:ProbePosition ixn.
#According to predictions we should see a significant diff. in RTs in REL-AFTER-V vs. REL-AFTER-VP conditions,
#with the REL-AFTER-V conditions being faster than the REL-AFTER-VP conditions.
#EMMs indicate an 8ms difference in favour of REL probes presented AFTER-VP. This contrast is not significant
#with p = 0.4161. There's therefore statistically no difference in these RTs. 

#EMMs for priming effect
emmeans(DNP.ALL_PARS.gamma.identity7, pairwise ~ ProbeType, adjust = "none")

#A post-hoc test to determine whether there's a priming effect is significant w/ a difference of
#68ms in favour of REL conditions with p <.0001. Effects on RTs therefore appear to be driven by
#ProbeType alone, with REL probes being responded to faster than UNREL probes regardless of ProbePosition.

#LOOKING AT QNP CONDITIONS

#histogram to check distribution - looks like gamma
exp1_all %>% 
  filter(SubjectType == "QNP") %>%
  ggplot(aes(x = RT)) + #specify what data to show on x axis
  geom_histogram(binwidth = 200) + #produce histogram with bins that are 200ms wide
  labs(x = "RT (ms)", y ="Frequency") + #x and y axes labels
  ggtitle("Hist. of RTs (ms)") + #add title
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centred
        text = element_text(size = 12)) #adjust rest of chart text - size 12

#Cullen and Frey Graph to check distribution more quantitatively - definitely gamma
descdist(filter(exp1_all, SubjectType == "QNP")$RT,
         discrete = FALSE, #we are not dealing with discrete data - it's continuous RT data
         method = "unbiased", #don't use bootstrapping for sampling - use the exact data without replacement
         graph = TRUE, #produce a graph!
         obs.col = "red", #the colour of the "observation" marker - red in this case
         obs.pch = 15) #the shape of the "observation" marker - square in this case

#summary stats. for experimental conditions
exp1_all %>%
  filter(SubjectType == "QNP") %>%
  group_by(Condition) %>% #group by Condition
  summarise("MeanRT" = mean(RT), #mean for each condition
            "StdDevRT" = sd(RT)) %>% #SD for each condition
  arrange(desc(MeanRT)) %>% #arrange descending
  kableExtra::kbl(caption = "Mean RTs by Condition") %>%
  kableExtra::kable_classic(full_width = T, html_font = "Arial")


#mean RTs in UNREL-AFTER-V conditions are slowest and RTs in REL-AFTER-V are fastest.

#viz
set.seed(45422)
exp1_all %>%
  filter(SubjectType == "QNP") %>%
  mutate(Condition = fct_reorder(Condition, log(RT), .desc = FALSE)) %>% #reorder log of RT by Condition from lowest to highest from bottom up
  ggplot(aes(x = ProbeType:ProbePosition, y = log(RT), colour = ProbeType)) + #specify what data to show on x and y axes, 
  geom_violin() + #produce a violin plot
  geom_jitter(width = .2, alpha = .08) + #add jitter to minimise overlap of points
  stat_summary(fun.data = mean_se, colour = "black") + 
  guides() + #legend
  labs(x = "ProbeType x ProbePosition", #x axis label
       y = "LogRT (ms)") + #y axis label
  ggtitle("Effect of Probe Type x Probe Position Interaction on RTs \n in QNP Conditions") + #add title - \n means 'new line'
  theme_bw() + #set theme so that background is white not grey and add a border around the chart (aesthetically nicer!)
  theme(plot.title = element_text(size = 14, face = "bold", hjust = .5), #adjust title text - size 14, bold, centered
        text = element_text(size = 12)) + #adjust rest of chart text - size 12
  coord_flip() #flip horizontal

#looking at covariates:
#checking for repetition effect/effect of trial number: weak -VE correlation
exp1_all %>%
  filter(SubjectType == "QNP") %>%
  ggscatter(x = "TrialN", y = "RT", 
            add = "reg.line", #adds regression line
            cor.coef = TRUE, #gives us correlation coefficient
            conf.int = TRUE, #gives us confidence intervals around the line
            add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
            xlab = "Trial Number", ylab = "RT (ms)")

#looking at covariates:
#checking for correlation between reading times and RTs: moderate +VE correlation
exp1_all %>%
  filter(SubjectType == "QNP") %>%
  ggscatter(x = "RT", y = "ReadingTime", 
            add = "reg.line", #adds regression line
            cor.coef = TRUE, #gives us correlation coefficient
            conf.int = TRUE, #gives us confidence intervals around the line
            add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
            xlab = "Reading Time (ms) between Overt Sentential Subject Onset and LDT Probe Onset", ylab = "RT (ms)")

#looking at covariates:
#checking for correlation between ProbeFrequency and RTs: weak -VE correlation (NB: not specified in prereg!)
exp1_all %>%
  filter(SubjectType == "QNP") %>%
  ggscatter(x = "ProbeFrequency", y = "RT", 
            add = "reg.line", #adds regression line
            cor.coef = TRUE, #gives us correlation coefficient
            conf.int = TRUE, #gives us confidence intervals around the line
            add.params = list(color = "blue", fill = "lightgray"), #colour of line + confidence interval
            xlab = "Zipf Frequency of Probe Word", ylab = "RT (ms)")

#MODEL BUILDING
#gamma model w/ identity link (see Lo & Andrews 2015)
#maximal model: "keep it maximal" (see Bates et al. 2013)
QNP.ALL_PARS.gamma.identity1 <- glmer(RT ~ ProbeType * ProbePosition + #2 factors w/ 2-way interaction
                               (1+ProbeType * ProbePosition | StimID) + #ranefx
                               (1+ProbeType * ProbePosition | Participant), #ranefx
                             family = Gamma(link = "identity"), #gamma distribution w/ "identity" link fxn.
                             data = exp1_all %>% filter(SubjectType == "QNP")) #only QNP conditions
summary(QNP.ALL_PARS.gamma.identity1) #convergence warnings

#maximal model doesn't converge - reduce ranefx. (see Matuschek et al. 2017)
#also increase number of evaluations for good measure to protect against convergence failure
QNP.ALL_PARS.gamma.identity2 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType * ProbePosition | StimID) +
                                        (1+ProbeType + ProbePosition | Participant), #no ixn.
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity2) #convergence warnings

#simplify ranefx. further
QNP.ALL_PARS.gamma.identity3 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType + ProbePosition | StimID) + #no ixn.
                                        (1+ProbeType + ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity3) #convergence warnings

#simplify ranefx. further
QNP.ALL_PARS.gamma.identity4 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1+ProbeType | StimID) + #dropped ProbePosition
                                        (1+ProbeType + ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity4) #convergence warnings

#simplify ranefx. further
QNP.ALL_PARS.gamma.identity5 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) + #dropped ProbeType
                                        (1+ProbeType + ProbePosition | Participant),
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity5) #convergence warnings

#simplify ranefx. further
QNP.ALL_PARS.gamma.identity6 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) +
                                        (1+ProbeType | Participant), #dropped ProbePosition
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity6) #convergence warnings

#simplify ranefx. further
QNP.ALL_PARS.gamma.identity7 <- glmer(RT ~ ProbeType * ProbePosition +
                                        (1 | StimID) +
                                        (1 | Participant), #dropped ProbeType
                                      family = Gamma(link = "identity"),
                                      data = exp1_all %>% filter(SubjectType == "QNP"),
                                      control=glmerControl(optCtrl = list(maxfun=2e5))) 
summary(QNP.ALL_PARS.gamma.identity7) #converges - significant ProbeType effect, significant 2-way ixn.

#plot it
emmip(QNP.ALL_PARS.gamma.identity7, ProbePosition ~ ProbeType) #crossover

#EMMs
emmeans(QNP.ALL_PARS.gamma.identity7, pairwise ~ ProbeType * ProbePosition, adjust = "none")

#EMMs for priming effect
emmeans(QNP.ALL_PARS.gamma.identity7, pairwise ~ ProbeType, adjust = "none") #significant priming effect

#significant main effect of ProbeType. Significant ProbeType:ProbePosition interaction.
#According to predictions we should see a significant diff. in RTs in REL-AFTER-V vs. REL-AFTER-VP conditions.
#EMMs indicate a 16ms difference in favour of REL probes presented AFTER-V. This contrast is not significant
#with p = 0.0834.

#there is a priming effect overall with a 67ms difference in favour of REL probes. This is significant
#with p < .0001. The ixn. plot shows a crossover effect: REL probes responded to faster AFTER-V vs. AFTER-VP
#(suggesting facilitated processing) and UNREL probes responded to slower AFTER-V vs. AFTER-VP 
#(suggesting a potential interference effect). A post-hoc test that includes this latter contrast reveals that the
#interference effect is significant with p = 0.0213*2 = 0.0426.

#Overall, it looks like REL probes are responded to faster regardless of ProbePosition in the DNP conditions,
#while REL probes are responded to faster AFTER-V in the QNP conditions. The cross-over effect suggests that
#UNREL probes are responded to slower AFTER-VP in the QNP conditions, though this is the only significant
#contrast once this is added.

#let's look at possible problematic items and do post-hoc analysis w/o
by_items_exp1_all <- exp1_all %>%
  group_by(ProbeType, StimID) %>%
  summarise("meanRT" = mean(RT)) %>%
  pivot_wider(names_from = ProbeType, values_from = c(meanRT)) %>%
  mutate(Difference = UNREL - REL) %>%
  arrange(Difference)

#plot
by_items_exp1_all %>%
  ggplot(aes(x = StimID, y = Difference)) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = StimID), size = 3) +
  geom_hline(yintercept=0)

#filter out problematic item(s) for a post-hoc analysis for just the QNP sentences
#since this is where we found the most interesting effects and cross-over

#QNP MODEL W/O PROBLEMATIC ITEM(S)
QNP.ALL_PARS.mod.items_removed <- update(QNP.ALL_PARS.gamma.identity7, 
                       data = exp1_all %>% 
                         filter(!StimID == "52", 
                                !StimID == "49",
                                SubjectType == "QNP"))
summary(QNP.ALL_PARS.mod.items_removed) #converges - significant ProbeType effect, significant 2-way ixn.

#plot it
emmip(QNP.ALL_PARS.mod.items_removed, ProbePosition ~ ProbeType) #crossover still present

#EMMs
emmeans(QNP.ALL_PARS.mod.items_removed, pairwise ~ ProbeType * ProbePosition, adjust = "none")
emmeans(QNP.ALL_PARS.mod.items_removed, pairwise ~ ProbeType, adjust = "none") #significant priming effect

#predictions unsupported but factoring in a post-hoc test comparing UNREL conditions as well as REL conditions, 
#there is a significant diff. between AFTER-V and AFTER-VP conditions with RTs in response to an 
#UNREL probe AFTER-V being 27ms slower than in response to an UNREL probe AFTER-VP with p = 0.0085*2 = 0.017. 
#No significant difference is seen for the REL conditions with p = 0.1568*2 = 0.3136. There's also a significant 
#priming effect with REL probes responded to 69ms faster with p < .0001. This warrants further investigation 
#into QNP sentences in a follow-up experiment targeting just QNP sentences to determine whether the effect 
#is robust and can be replicated.

#some summary stats. for write-up
exp1_all %>%
  summarise(MeanAge = mean(Age, na.rm = TRUE), StdDevAge = sd(Age, na.rm = TRUE))


