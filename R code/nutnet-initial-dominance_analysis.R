## Load Libraries

## data wrangling
library(tidyverse)
library(broom)

## analysis
library(codyn)
library(nlme)
library(lme4)
library(emmeans)
library(partR2)
library(car)
library(MuMIn)
#library(RRPP)
library(merTools)
library(trajr)

## graphics
library(ggplot2)
library(cowplot)
library(GGally)
library(RColorBrewer)
library(gridExtra)
library(grid)


# create graphics parameters/aesthetics
theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 24), axis.text = element_text(size = 24),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = 'right', legend.justification = c(1, 0.4),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=18,face='bold'),legend.text = element_text(size = 18),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        strip.text =element_text(size=18,face='bold'),
        axis.line = element_line(colour = 'grey60', size = 0.35), axis.ticks = element_line(colour = 'grey60', size = 0.35)) 
theme_set(theme_figs)

# functions to tidy up model output tables
anova.tablefunc <- function(mod) {  
  dat <- data.frame(anova(mod)) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    dplyr::mutate(F.value = round(F.value, digits = 2), 
                  p.value = round(p.value, digits = 3)) %>%
    dplyr::mutate(p.value = ifelse(p.value == 0.000, '< 0.001', p.value)) %>%
    tidyr::unite(DF, numDF, denDF, sep = ',') %>%
    dplyr::rename(`F value` = F.value, `P value` = p.value)
  return(dat)
}

summary.tablefunc <- function(mod) {  
  dat <- data.frame(summary(mod)$tTable) %>%
    tibble::rownames_to_column(var = 'Effect') %>%
    rename_with(stringr::str_replace, 
                pattern = "-", replacement = ".", 
                matches("Length")) %>% 
    dplyr::mutate(Estimate = signif(Value, digits = 3),
                  Std.Error = signif(Std.Error, digits = 2),
                  t.value = signif(t.value, digits = 2),
                  p.value = signif(p.value, digits = 2)) %>%
    dplyr::mutate(p.value = ifelse(p.value <= 0.001, '< 0.001', p.value)) %>% 
    dplyr::select(-Value) %>% 
    relocate(Estimate,.before = Std.Error)
  return(dat)
}


#read in data
setwd("H:/other NutNet papers/peter Wilfahrt/")

dominant_pop <- read_csv('Dominants-through-time_2023-01-08.csv')
dominant_pop # take a look

# create nice labels
dominant_pop$local_lifespan <- gsub('ANNUAL','Annual',dominant_pop$local_lifespan)
dominant_pop$local_lifespan <- gsub('PERENNIAL','Perennial',dominant_pop$local_lifespan)


#### Model of rank decay ####

##### swapping in site richness in year 0 for site richness across all years

dominant_pop$site_richness <- dominant_pop$site_rich_0

### find percentage of species that never dip below 1 (all sites, hen 10+ yr_trt sites)

rank_avg <- dominant_pop %>% 
  group_by(site_code) %>% 
  mutate(max_year = max(year_trt)) %>% ungroup() %>% 
  group_by(site_code, plot, Taxon, max_year) %>% 
  summarize(avg_rank = mean(perc_rank))

nrow(rank_avg %>% filter(avg_rank == 1)) / nrow(rank_avg) ## 22.1%
nrow(rank_avg %>% filter(avg_rank == 1, max_year >= 10)) / nrow(rank_avg %>%filter(max_year >= 10)) ## 9.1%


# alternative code for arcsin transformation
# asinTransform <- function(p) { asin(sqrt(p)) }
# dominant_pop$rank_as <- asinTransform(dominant_pop$perc_rank)

### logit transformation
dominant_pop$rank_adj <- ifelse(dominant_pop$perc_rank == 1, 0.99,
                                ifelse(dominant_pop$perc_rank == 0, 0.01, dominant_pop$perc_rank))
dominant_pop$rank_logit <- logit(dominant_pop$rank_adj)

#### reduce functional groups to two roughly equal (taxon-wise) groups for simplicity
dom0 <- dominant_pop %>% filter(year_trt == 0)
tapply(dom0$functional_group,dom0$functional_group,length)

dominant_pop$func_simple <- ifelse(dominant_pop$functional_group == 'GRAMINOID','Graminoid','Non-graminoid')

dom_complete <- dominant_pop %>% filter(local_provenance != 'UNK')

# why block was not included in your model 
mod_rank <- lme(rank_logit ~ year_trt * initial_rel_cover * NPK * Fence +
                                   local_lifespan * year_trt * NPK * Fence  +
                                   local_provenance * year_trt * NPK * Fence  +
                                   func_simple * year_trt * NPK * Fence +
                                   MAP_v2 * year_trt * NPK * Fence +
                                   MAP_VAR_v2 * year_trt * NPK * Fence +
                                   site_richness * year_trt * NPK * Fence,
                                   random = ~ 1|site_code/plot,
                                 data = dom_complete)


# check normality of residuals
qqnorm(residuals(mod_rank)) #looks fine
plot(mod_rank) 
# QC: check model fit, model does not fit the data well 

 # create summary table
rank_table <- summary.tablefunc(mod_rank)

rank_table <- rank_table %>% 
  mutate(order = if_else(Effect %in% grep('initial',rank_table$Effect, value = TRUE),2,
                   if_else(Effect %in% grep('lifespan',rank_table$Effect, value = TRUE),3,
                     if_else(Effect %in% grep('provenance',rank_table$Effect, value = TRUE),4,
                       if_else(Effect %in% grep('func',rank_table$Effect, value = TRUE),5,
                        if_else(Effect %in% grep('richness',rank_table$Effect, value = TRUE),6,
                          if_else(Effect %in% grep('MAP_v2',rank_table$Effect, value = TRUE),7,
                           if_else(Effect %in% grep('MAP_VAR_v2',rank_table$Effect, value = TRUE),8,1)))))))) %>% 
  arrange(order)

rank_table

#correlation
cor(dominant_pop[!duplicated(dominant_pop$site_code) == FALSE,c('MAP_v2','site_richness','MAP_VAR_v2','MAT_v2')])

r.squaredGLMM(mod_rank)

# write_csv(rank_table,
#           paste0('./Tables/rank-decay-logit-table_',
#                  Sys.Date(),'.csv'))

#* relative importance via AIC ####

AIC(mod_rank) #42439.56
mod.cov <- update(mod_rank, . ~ . - initial_rel_cover - NPK:initial_rel_cover - year_trt:initial_rel_cover
              - Fence:initial_rel_cover -year_trt:NPK:initial_rel_cover - year_trt:Fence:initial_rel_cover
              - year_trt:NPK:Fence:initial_rel_cover) # 42491.69

mod.life <- (update(mod_rank, . ~ . - local_lifespan - NPK:local_lifespan - year_trt:local_lifespan
           - Fence:local_lifespan -year_trt:NPK:local_lifespan - year_trt:Fence:local_lifespan
           - year_trt:NPK:Fence:local_lifespan)) # 42466.72

mod.func <- (update(mod_rank, . ~ . - func_simple - NPK:func_simple - year_trt:func_simple
           - Fence:func_simple -year_trt:NPK:func_simple - year_trt:Fence:func_simple
           - year_trt:NPK:Fence:func_simple)) # 42433.07

mod.prov <- (update(mod_rank, . ~ . - local_provenance - NPK:local_provenance - year_trt:local_provenance
           - Fence:local_provenance -year_trt:NPK:local_provenance - year_trt:Fence:local_provenance
           - year_trt:NPK:Fence:local_provenance)) # 42465.61

mod.rich <- (update(mod_rank, . ~ . - site_richness - NPK:site_richness - year_trt:site_richness
           - Fence:site_richness -year_trt:NPK:site_richness - year_trt:Fence:site_richness
           - year_trt:NPK:Fence:site_richness)) # 42385.85

mod.map <- (update(mod_rank, . ~ . - MAP_v2 - NPK:MAP_v2 - year_trt:MAP_v2
           - Fence:MAP_v2 -year_trt:NPK:MAP_v2 - year_trt:Fence:MAP_v2
           - year_trt:NPK:Fence:MAP_v2)) # 42331.08

mod.var <- (update(mod_rank, . ~ . - MAP_VAR_v2 - NPK:MAP_VAR_v2 - year_trt:MAP_VAR_v2
           - Fence:MAP_VAR_v2 -year_trt:NPK:MAP_VAR_v2 - year_trt:Fence:MAP_VAR_v2
           - year_trt:NPK:Fence:MAP_VAR_v2)) # 42419.36

aic.tab <- data.frame(Predictor = c('Full Model','Initial Cover','Lifespan','Provenance','Functional Group','Site Richness','MAP','Preciptation variability'),
                         AIC = c(AIC(mod_rank),AIC(mod.cov),AIC(mod.life),AIC(mod.prov),AIC(mod.func),AIC(mod.rich),AIC(mod.map),AIC(mod.var))) %>% 
  mutate(deltaAIC = round(AIC - AIC[1],1))
                         
aic.tab
# QC: mod.map has the lowest aic
plot(mod.map)
                         
# to add pseudo-r2
                         # R2 = c(r.squaredGLMM(mod_rank)[1],r.squaredGLMM(mod.cov)[1],r.squaredGLMM(mod.life)[1],r.squaredGLMM(mod.prov)[1],r.squaredGLMM(mod.func)[1],
                         #        r.squaredGLMM(mod.rich)[1],r.squaredGLMM(mod.map)[1],r.squaredGLMM(mod.var)[1]))


#### graphing rank decay ####

## define function for determining graphical display of continuous gradients

# create color palette for graphs
display.brewer.pal(11, "Blues")
brewer.pal(11, "Blues")

pos <- position_dodge(0.5)

ctrl_pal <- c('black',"#D53E4F","#6BAED6","#5E4FA2")

# define estimate of continuous variable to graph
# upper <- function(x) max(x)
# lower <- function(x) min(x)
upper <- function(x) quantile(x, 0.95)
lower <- function(x) quantile(x, 0.05)
# upper <- function(x) mean(x) + (sd(x)*1.5)
# lower <- function(x) mean(x) - (sd(x)*1.5)



#create discrete categories for graphing purporses
dominant_pop <- dominant_pop %>% 
  mutate(
    initial_cover  = if_else(initial_rel_cover < median(initial_rel_cover), 'Low Initial Cover','High Initial Cover'),
    richness  = if_else(site_richness < median(site_richness), 'Low site richness','High site richness'),
    ppt  = if_else(MAP_v2 < median(MAP_v2), 'Low MAP','High MAP'),
    pptvar  = if_else(MAP_VAR_v2 < median(MAP_VAR_v2), 'Low PPT Variation','High PPT Variation')
    )


#* lifespan effect --------

rank_table %>% filter(Effect %in% grep('lifespan',rank_table$Effect, value = TRUE))

### all sig

### get model estimates at specified covariate levels
# QC: according to your model selection process, the model without MAP(?) has the lowest AIC, and therefore should be the final model
# but below you estimate parameters using the original model 
rank_life_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_lifespan * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

# rename treatments
rank_life_emm <- rank_life_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

# code in significance
#rank_life_emm$sig <- ifelse(rank_life_emm$trt == 'NPK+Fence','No','Yes')
rank_life_emm$sig <- 'Yes'

#back transformation for arcsin
# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

#back transformation for logit
rank_life_emm$em_btrans <- invlogit(rank_life_emm$emmean)
rank_life_emm$low_btrans <- invlogit(rank_life_emm$lower.CL)
rank_life_emm$up_btrans <- invlogit(rank_life_emm$upper.CL)

#create graph
gg_rank_year_life_trt <- ggplot(rank_life_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  geom_line(aes(group = trt), size = 2) +
  #geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  #scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  facet_grid(~local_lifespan)
gg_rank_year_life_trt

# code for saving individually, but will stitch these together later on
# ggsave(paste0('Graphs/Presentation/rank-decay-lifespan_',Sys.Date(),'.png'),
#        gg_rank_year_life_trt, width = 12, height = 8, units = "in", dpi = 600)


#* provenance effect --------

rank_table %>% filter(Effect %in% grep('provenance',rank_table$Effect, value = TRUE))
### Fence sig, NPK marginal

rank_prov_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_provenance * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_prov_emm <- rank_prov_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_prov_emm$sig <- ifelse(rank_life_emm$trt %in% c('NPK+Fence'),'No','Yes')

# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

rank_prov_emm$em_btrans <- invlogit(rank_prov_emm$emmean)
rank_prov_emm$low_btrans <- invlogit(rank_prov_emm$lower.CL)
rank_prov_emm$up_btrans <- invlogit(rank_prov_emm$upper.CL)

# rename factors for nicer labels
rank_prov_emm$local_provenance <- ifelse(rank_prov_emm$local_provenance == 'NAT','Native','Non-native')
dom_complete$local_provenance <- ifelse(dom_complete$local_provenance == 'NAT','Native','Non-native')

gg_rank_year_prov_trt <- ggplot(rank_prov_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dom_complete, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  #geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  facet_grid(~local_provenance)
gg_rank_year_prov_trt


#* functional group effect --------


rank_table %>% filter(Effect %in% grep('func',rank_table$Effect, value = TRUE))
### only fence is sig
##QC: from what I see, several terms are significant 

rank_func_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | func_simple * year_trt,
                                            at = list(year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_func_emm <- rank_func_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_func_emm$sig <- ifelse(rank_func_emm$trt %in% c('NPK+Fence','NPK'),'No','Yes')

# rank_life_emm$em_btrans <- sin(rank_life_emm$emmean)^2
# rank_life_emm$low_btrans <- sin(rank_life_emm$lower.CL)^2
# rank_life_emm$up_btrans <- sin(rank_life_emm$upper.CL)^2

rank_func_emm$em_btrans <- invlogit(rank_func_emm$emmean)
rank_func_emm$low_btrans <- invlogit(rank_func_emm$lower.CL)
rank_func_emm$up_btrans <- invlogit(rank_func_emm$upper.CL)


gg_rank_year_func_trt <- ggplot(rank_func_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dom_complete, aes(x=year_trt, y = perc_rank),  width=0.2, alpha = 0.8, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.15,color=NA) +
  #geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  facet_grid(~func_simple)
gg_rank_year_func_trt


#* initial cover effect ----

rank_table %>% filter(Effect %in% grep('initial',rank_table$Effect, value = TRUE))

### all sig

# cov.low <- mean(dominant_pop$initial_rel_cover) - 2 * sd(dominant_pop$initial_rel_cover)
# cov.high <- mean(dominant_pop$initial_rel_cover) + 2 * sd(dominant_pop$initial_rel_cover)

cov.low <- lower(dominant_pop$initial_rel_cover)
cov.high <- upper(dominant_pop$initial_rel_cover)


rank_init_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | initial_rel_cover * year_trt,
                                            at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by =cov.high-cov.low),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_init_emm <- rank_init_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_init_emm$initial_cover <- factor(ifelse(rank_init_emm$initial_rel_cover == cov.low, 'Low Initial Cover','High Initial Cover'), levels = c('Low Initial Cover','High Initial Cover'))

# rank_init_emm$em_btrans <- sin(rank_init_emm$emmean)^2
# rank_init_emm$low_btrans <- sin(rank_init_emm$lower.CL)^2
# rank_init_emm$up_btrans <- sin(rank_init_emm$upper.CL)^2

#### should negative arcsin values be set to 0 on backtransformation (otherwise artificial uptick due to sin)?
rank_init_emm$em_btrans <- invlogit(rank_init_emm$emmean)
rank_init_emm$low_btrans <- invlogit(rank_init_emm$lower.CL)
rank_init_emm$up_btrans <- invlogit(rank_init_emm$upper.CL)



gg_rank_year_init_trt <- ggplot(rank_init_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  geom_line(aes(group = trt), size = 2, color= 'black') +
  geom_line(aes(group = trt), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  #ylim(0,1.05) +
  facet_wrap(~initial_cover)

gg_rank_year_init_trt

# ggsave(paste0('Graphs/Presentation/rank-decay-initcover_',Sys.Date(),'.png'),
#        gg_rank_year_init_trt, width = 12, height = 8, units = "in", dpi = 600)

#* richness effect -----

rank_table %>% filter(Effect %in% grep('richness',rank_table$Effect, value = TRUE))

### only npk sig

rich.low <- lower(dominant_pop$site_richness)
rich.high <- upper(dominant_pop$site_richness)


rank_rich_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | site_richness * year_trt,
                                            at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_rich_emm <- rank_rich_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_rich_emm$sig <- ifelse(rank_rich_emm$trt %in% c('Fence','NPK+Fence'),'No','Yes')


rank_rich_emm$richness <- factor(ifelse(rank_rich_emm$site_richness == rich.low, 'Low site richness','High site richness'), levels = c('Low site richness','High site richness'))

# rank_rich_emm$em_btrans <- sin(rank_rich_emm$emmean)^2
# rank_rich_emm$low_btrans <- sin(rank_rich_emm$lower.CL)^2
# rank_rich_emm$up_btrans <- sin(rank_rich_emm$upper.CL)^2

rank_rich_emm <-  mutate(rank_rich_emm,
                          em_btrans = invlogit(emmean), 
                          low_btrans = invlogit(lower.CL),
                          up_btrans = invlogit(upper.CL))


gg_rank_year_rich_trt <- ggplot(rank_rich_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~richness)

gg_rank_year_rich_trt

# ggsave(paste0('Graphs/Presentation/rank-decay-initcover_',Sys.Date(),'.png'),
#        gg_rank_year_init_trt, width = 12, height = 8, units = "in", dpi = 600)


#* MAP effect ----

rank_table %>% filter(Effect %in% grep('MAP_v2',rank_table$Effect, value = TRUE))

#fence sig

map.low <- lower(dominant_pop$MAP_v2)
map.high <- upper(dominant_pop$MAP_v2)


rank_map_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_v2 * year_trt,
                                            at = list(MAP_v2 = seq(from = map.low, to = map.high, by = map.high-map.low),
                                                      year_trt = seq(from = 0, to = 15, by = 1)),
                                            cov.reduce = TRUE,
                                            type = "response")$emmeans))

rank_map_emm <- rank_map_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))


rank_map_emm$sig <- ifelse(rank_map_emm$trt %in% c('NPK','NPK+Fence'),'No','Yes')


rank_map_emm$ppt <- factor(ifelse(rank_map_emm$MAP_v2 == map.low, 'Low MAP','High MAP'), levels = c('Low MAP','High MAP'))

# rank_map_emm$em_btrans <- sin(rank_map_emm$emmean)^2
# rank_map_emm$low_btrans <- sin(rank_map_emm$lower.CL)^2
# rank_map_emm$up_btrans <- sin(rank_map_emm$upper.CL)^2

rank_map_emm <-  mutate(rank_map_emm,
                         em_btrans = invlogit(emmean), 
                         low_btrans = invlogit(lower.CL),
                         up_btrans = invlogit(upper.CL))

gg_rank_year_ppt_trt <- ggplot(rank_map_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.25,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~ppt)

gg_rank_year_ppt_trt

#* MAP VAR effect ----

rank_table %>% filter(Effect %in% grep('MAP_VAR_v2',rank_table$Effect, value = TRUE))
# fence is marginally significant, year is only significant effect otherwise

mvar.low <- lower(dominant_pop$MAP_VAR_v2)
mvar.high <- upper(dominant_pop$MAP_VAR_v2)


rank_pptvar_emm <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_VAR_v2 * year_trt,
                                           at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = mvar.high-mvar.low),
                                                     year_trt = seq(from = 0, to = 15, by = 1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))

rank_pptvar_emm <- rank_pptvar_emm %>%
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL)))))

rank_pptvar_emm$sig <- ifelse(rank_pptvar_emm$trt %in% c('Fence','NPK','NPK+Fence'),'No','Yes')


rank_pptvar_emm$pptvar <- factor(ifelse(rank_pptvar_emm$MAP_VAR_v2 == mvar.low, 'Low PPT Variation','High PPT Variation'), levels = c('Low PPT Variation','High PPT Variation'))

# rank_pptvar_emm$em_btrans <- sin(rank_pptvar_emm$emmean)^2
# rank_pptvar_emm$low_btrans <- sin(rank_pptvar_emm$lower.CL)^2
# rank_pptvar_emm$up_btrans <- sin(rank_pptvar_emm$upper.CL)^2

rank_pptvar_emm <-  mutate(rank_pptvar_emm,
                        em_btrans = invlogit(emmean), 
                        low_btrans = invlogit(lower.CL),
                        up_btrans = invlogit(upper.CL))


gg_rank_year_pptvar_trt <- ggplot(rank_pptvar_emm, aes(x = year_trt, y = em_btrans, col = trt)) +
  geom_jitter(data = dominant_pop, aes(x=year_trt, y = perc_rank),  width=0.11, alpha = 0.7, col = 'grey90') +
  geom_ribbon(aes(ymin = low_btrans,ymax= up_btrans,fill=trt),alpha=0.35,color=NA) +
  #geom_line(aes(group = trt), size = 3.2, color= 'black') +
  geom_line(aes(group = trt, linetype = sig), size = 1.8) +
  guides(color = guide_legend(title = "Treatment"), fill = guide_legend(title = "Treatment")) +
  ylab('Rank percentage') + xlab('Year after treatment') +
  #scale_color_brewer(palette = "Set1") +
  scale_color_manual(values = ctrl_pal) +
  scale_fill_manual(values = ctrl_pal) +
  scale_linetype_manual(values = c('dashed','solid'), guide = 'none') +
  ylim(0,1.05) +
  facet_wrap(~pptvar)

gg_rank_year_pptvar_trt

#### Fig 2 ####
#### stitch them together

rank_legend <- get_legend(gg_rank_year_func_trt + 
                            theme(legend.direction = "horizontal", legend.box = "horizontal"))

gg_rank_init <- gg_rank_year_init_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
                                             
gg_rank_life <- gg_rank_year_life_trt + theme(legend.position="none", axis.title.y=element_blank())
                                              #axis.title.x=element_blank(),axis.text.x=element_blank())

# Fig 2 with initial rank and lifespan
# Fig2 <- ggdraw() +
#   draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
#   draw_plot(gg_rank_init, x = 0.03, y = 0.52, width = .97, height = .43) +
#   draw_plot(gg_rank_life, x = 0.03, y = 0, width = .97, height = .5) +
#   draw_plot_label(label = c("A", "B"), size = 20,
#                     x = c(0.03, 0.03), y = c(0.97, 0.51)) +
#   draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25)

#Fig 2 just initial cover

Fig2 <- gg_rank_year_init_trt

Fig2


# Fig 3 - species chars

gg_rank_life <- gg_rank_year_life_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_func <- gg_rank_year_func_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_prov <- gg_rank_year_prov_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())

#axis.title.x=element_blank(),axis.text.x=element_blank())

Fig3 <- ggdraw() +
  draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
  draw_plot(gg_rank_life, x = 0.03, y = 0.65, width = .97, height = .3) +
  draw_plot(gg_rank_prov, x = 0.03, y = 0.35, width = .97, height = .3) +
  draw_plot(gg_rank_func, x = 0.03, y = 0.05, width = .97, height = .3) +
  draw_plot_label(label = c("A", "B", "C"), size = 20,
                    x = c(0.03, 0.03, 0.03), y = c(0.95, 0.65, 0.35)) +
  draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25) +
  draw_text('Year after treatment', x = 0.35, y= 0.02, size = 25)

Fig3

# Fig S1 - environmental characteristics

gg_rank_rich <- gg_rank_year_rich_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_ppt <- gg_rank_year_ppt_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())
gg_rank_pptvar <- gg_rank_year_pptvar_trt + theme(legend.position="none", axis.title.y=element_blank(),axis.title.x=element_blank())

#axis.title.x=element_blank(),axis.text.x=element_blank())

Fig4 <- ggdraw() +
  draw_plot(rank_legend, x = -0.3, y = 0.95, width = 1, height = .05) +
  draw_plot(gg_rank_rich, x = 0.03, y = 0.65, width = .97, height = .3) +
  draw_plot(gg_rank_ppt, x = 0.03, y = 0.35, width = .97, height = .3) +
  draw_plot(gg_rank_pptvar, x = 0.03, y = 0.05, width = .97, height = .3) +
  draw_plot_label(label = c("A", "B", "C"), size = 20,
                  x = c(0.03, 0.03, 0.03), y = c(0.95, 0.65, 0.35)) +
  draw_text('Rank Percentile',angle = 90, x = 0.015, size = 25) +
  draw_text('Year after treatment', x = 0.35, y= 0.02, size = 25)

Fig4

Fig2

# save_plot(paste0('Graphs/Fig2-rank-initial-cover_',Sys.Date(),'.png'),
#           Fig2,
#           base_height = 8, base_width = 10)
# # 
# save_plot(paste0('Graphs/Fig3-rank-traits_',Sys.Date(),'.png'),
#           Fig3,
#           base_height = 14, base_width = 10)
# # 
# save_plot(paste0('Graphs/Fig4-rank-env_',Sys.Date(),'.png'),
#           Fig4,
#           base_height = 14, base_width = 10)
# 



# Defining inertia #####

### This creates an 'effect size' table of when model estimates drops below mean rank value for different covariate levels
### output is estimated year (in tenths) at which given conditions are expected to drop below the mean

#define mean rank percentile
cf = mean(dom_complete$perc_rank[dom_complete$year_trt != 0])

inertia_npk <-  data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence |  year_trt,
                                                          at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                                          cov.reduce = TRUE,
                                                          type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL),Predictor = 'Average') %>% 
  group_by(NPK,Fence) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = Predictor) 

inertia_npk

inertia_init <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | initial_rel_cover * year_trt,
#                                           at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by = (cov.high-cov.low)/2),
                                            at = list(initial_rel_cover = seq(from = cov.low, to = cov.high, by = (cov.high-cov.low)),
                                                                         year_trt = seq(from = 0, to = 40, by = 0.1)),
                                                               cov.reduce = TRUE,
                                                               type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, initial_rel_cover) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(initial_cover = factor(if_else(initial_rel_cover == cov.low, 'Low Initial Cover',
                                        if_else(initial_rel_cover == cov.high,'High Initial Cover', 'Average Initial Cover')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = initial_cover) %>% 
  rename(Predictor = 1) 

inertia_init
# qc: error 

map.low <- lower(dom_complete$MAP_v2)
map.high <- upper(dom_complete$MAP_v2)

inertia_map <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_v2 * year_trt,
                                           # at = list(MAP_v2 = seq(from = map.low, to = map.high, by = (map.high-map.low)/2),
                                                     at = list(MAP_v2 = seq(from = map.low, to = map.high, by = (map.high-map.low)),
                                                     year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, MAP_v2) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(MAP = factor(if_else(MAP_v2 == map.low, 'Low MAP',
                                        if_else(MAP_v2 == map.high,'High MAP', 'Average MAP')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = MAP) %>% 
  rename(Predictor = 1)

inertia_map

mvar.low <- lower(dom_complete$MAP_VAR_v2)
mvar.high <- upper(dom_complete$MAP_VAR_v2)

inertia_mapvar <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | MAP_VAR_v2 * year_trt,
                                          # at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = (mvar.high-mvar.low)/2),
                                                    at = list(MAP_VAR_v2 = seq(from = mvar.low, to = mvar.high, by = (mvar.high-mvar.low)),
                                                    year_trt = seq(from = 0, to = 15, by = 0.1)),
                                          cov.reduce = TRUE,
                                          type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, MAP_VAR_v2) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(ppt_variation = factor(if_else(MAP_VAR_v2 == mvar.low, 'Low PPT Variation',
                              if_else(MAP_VAR_v2 == mvar.high,'High PPT Variation', 'Average PPT Variation')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = ppt_variation)%>% 
  rename(Predictor = 1)

inertia_mapvar

rich.low <- lower(dom_complete$site_richness)
rich.high <- upper(dom_complete$site_richness)

inertia_rich <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | site_richness * year_trt,
                                             # at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)/2),
                                                       at = list(site_richness = seq(from = rich.low, to = rich.high, by = (rich.high-rich.low)),
                                                       year_trt = seq(from = 0, to = 15, by = 0.1)),
                                             cov.reduce = TRUE,
                                             type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, site_richness) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  mutate(Site_Richness = factor(if_else(site_richness == rich.low, 'Low Site Richness',
                                        if_else(site_richness == rich.high,'High Site Richness', 'Average Site Richness')))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = Site_Richness)%>% 
  rename(Predictor = 1)

inertia_rich

inertia_life <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_lifespan * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, local_lifespan) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = local_lifespan)%>% 
  rename(Predictor = 1)

inertia_life  

inertia_func <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | func_simple * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, func_simple) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = func_simple)%>% 
  rename(Predictor = 1)

inertia_func
  
inertia_prov <- data.frame(summary(emmeans(mod_rank, pairwise ~ NPK * Fence | local_provenance * year_trt,
                                           at = list(year_trt = seq(from = 0, to = 15, by = 0.1)),
                                           cov.reduce = TRUE,
                                           type = "response")$emmeans))  %>% 
  mutate(trt = if_else(NPK == 0 & Fence == 0, 'Control',
                       if_else(NPK == 0 & Fence == 1, 'Fence',
                               if_else(NPK == 1 & Fence == 0, 'NPK',
                                       if_else(NPK == 1 & Fence == 1, 'NPK+Fence',NULL))))) %>% 
  mutate(cutoff = cf, mean_btrans = invlogit(emmean), upper_btrans = invlogit(upper.CL)) %>% 
  group_by(NPK,Fence, local_provenance) %>% 
  arrange(year_trt) %>% 
  # slice(if(any(upper_btrans < cutoff)) which.max(upper_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  slice(if(any(mean_btrans < cutoff)) which.max(mean_btrans < cutoff) else which.max(year_trt == max(year_trt))) %>% 
  pivot_wider(names_from = trt, values_from = year_trt, id_cols = local_provenance) %>% 
  rename(Predictor = 1)

inertia_prov

#put the tables together
inertia_tab <- bind_rows(inertia_npk, inertia_init,inertia_life,inertia_func,inertia_prov,inertia_map,inertia_mapvar,inertia_rich)

#calculate percentage change
inertia_tab <-inertia_tab %>% mutate(fence_per = -100*(1-Fence/Control),
                       npk_per = -100*(1-NPK/Control),
                       npkf_per = -100*(1-`NPK+Fence`/Control)) %>% 
  mutate(Fence = paste0(Fence, ' (',signif(fence_per,3),'%)'),
         NPK = paste0(NPK, ' (',signif(npk_per,3),'%)'),
         `NPK+Fence` = paste0(`NPK+Fence`, ' (',signif(npkf_per,3),'%)')) %>% 
  dplyr::select(-c(fence_per,npk_per,npkf_per))
  
inertia_tab

#write_csv(inertia_tab,paste0('./Tables/Inertia-effect-size_',Sys.Date(),'.csv'))
