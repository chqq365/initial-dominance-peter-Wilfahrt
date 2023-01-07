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

## graphics
library(ggplot2)
library(cowplot)
library(GGally)

# remove MASS if loaded from analysis script (interferes with dplyr select)
# detach(package:merTools)
# detach(package:arm)
# detach(package:MASS)

theme_figs <- theme_classic() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 18),
        plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = 'right', legend.justification = c(1, 0.4),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size=18,face='bold'),legend.text = element_text(size = 18),
        plot.margin = unit(c(0.35, 0, 0, 0), 'cm'),
        strip.text =element_text(size=18,face='bold'),
        axis.line = element_line(colour = 'grey60', size = 0.35), axis.ticks = element_line(colour = 'grey60', size = 0.35)) 
theme_set(theme_figs)

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

## Read in data
setwd("H:/other NutNet papers/peter Wilfahrt/")

original_coverdat <- read.csv('full-cover-by-date-2022-10-10.csv')
coverdat <- original_coverdat # backup
dim(coverdat)
# 305412   x  21 - sorted by date

original_comb <- read_csv("comb-by-plot-clim-soil-diversity-2022-11-15.csv")
comb<-original_comb

# QC: I heard from Jon Bakker he updates his function for each version of the cover data release and he suggests that we should use the corresponding taxa adjutment file
# there is a new version of taxonomic-adjustments in the dropbox 
# https://www.dropbox.com/home/qingqing%20chen/NutNet%20data?preview=site-taxonomy-2022-07-13.csv
cover_fix <- read_csv('https://github.com/jon-bakker/NutNet_TaxonomicAdjustments/raw/main/taxonomic-adjustments-2022-05-02.csv')

# choose all experimental sites, with at least an N treatment, 1 year post-treatment:
sites <-unique(comb$site_code[comb$trt %in% c("NPK", "Fence") # at least an NPK treatment or fence
                               & (!comb$experiment_type %in% c("Observational")) # not observational site
                               &  comb$year_trt >=1 # has at least one year post-treatment data
                               ]) 


length(sites)

# just the sites we want
working_coverdat <- coverdat[coverdat$site_code %in% sites,]

### lagoas has off-season pre-treatment year, but data looks sufficient to use after checking

working_coverdat$year_trt <- ifelse(working_coverdat$site_code == 'lagoas.br' & working_coverdat$year_trt == -1,
                                    0, working_coverdat$year_trt)


## need to make sure they have pre-treatment data as well
site.miss <- sites[!sites %in%  unique(working_coverdat$site_code[working_coverdat$year_trt == 0])]

sites <- sort(unique(working_coverdat$site_code[working_coverdat$year_trt == 0]))

working_coverdat <- coverdat[coverdat$site_code %in% sites,]


# reduce to just vascular plants
# maybe also change species with cover > 100 to 100 according to Jon Bakker, you also mentioned in your script that cover has a maximum of 100 
working_coverdat <- droplevels(working_coverdat) %>% 
  filter(live == 1, !functional_group %in% c('BRYOPHYTE','LICHEN','LIVERWORT'), max_cover != 0) %>% 
  group_by(site_code) %>% 
  mutate(max_year = max(year_trt))  ##  add in max years
range(working_coverdat$max_cover)

### reduce functional groups
working_coverdat$functional_group <- ifelse(working_coverdat$functional_group %in% c('GRASS','GRAMINOID'), 'GRAMINOID',working_coverdat$functional_group)

### doane.us only sampled block 1 after pre-treatment year

working_coverdat <- working_coverdat %>% 
  filter(!(site_code == 'doane.us' & block != 1))

#### Jena was only sampled in september in pre-treatment, so only using sept data thereafter

unique(working_coverdat[working_coverdat$site_code == 'jena.de',]$date)
dim(working_coverdat[working_coverdat$site_code == 'jena.de',])

working_coverdat <- working_coverdat %>% 
  filter(!(site_code == 'jena.de' & substr(date,6,7) %in% c('05','06')))


### bnch.us: per Seabloom combine all Carex into single taxa
working_coverdat %>%
  filter(site_code == 'bnch.us', substr(Taxon,1,5) == 'CAREX') %>% 
    count(year,Taxon) 

working_coverdat[working_coverdat$site_code == 'bnch.us' & substr(working_coverdat$Taxon,1,5) == 'CAREX',]$Taxon <- 'CAREX SP.'

working_coverdat %>%
  filter(site_code == 'hopl.us', substr(Taxon,1,5) == 'AVENA') %>% 
  count(year,Taxon) 

# hopl.us	all	Avena barbata	I would check to see whether the zeroes are due to misidentifying as A. fatua or being labeled as Avena species
# plot 2

working_coverdat %>%
  filter(site_code == 'hopl.us', substr(Taxon,1,5) == 'AVENA',plot == 2) %>% 
  select(plot, year_trt, Taxon, max_cover) %>% 
  print(n=60)

### sier.us	all	Bromus diandrus	We had problems early on separating B. diandrus and B. sterilis. You might need to combine these species, due uncertainty early in experiment

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS')) %>% 
  count(year,Taxon) %>% 
  print(n=60)

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

### for now we will combine both into BROMUS SP. and later take the maximum value of duplicates (needed to get around Jena fix - see below)
working_coverdat[working_coverdat$site_code == 'sier.us' & working_coverdat$Taxon %in%  c('BROMUS DIANDRUS','BROMUS STERILIS'),]$Taxon <- 'BROMUS SP.'


# doane.us	8,10, 19	Sporabolus compositus	This was from the pre-treatment year from summer of 2012 when Nebraska was facing a drought, so I wonder if I misidentified this species as Sporabolus, but it was Andropogon gerardii, which is more common there.	

working_coverdat %>%
  filter(site_code == 'doane.us', Taxon %in%  c('SPOROBOLUS COMPOSITUS','ANDROPOGON GERARDII')) %>% 
  count(year,Taxon) %>% 
  print(n=60)
# sporobolus completely disappears in Y2

### for now I'll turn SC into AG and later take the maximum value of duplicates (needed to get around Jena fix - see below)
working_coverdat[working_coverdat$site_code == 'doane.us' & working_coverdat$Taxon %in%  c('SPOROBOLUS COMPOSITUS','ANDROPOGON GERARDII'),]$Taxon <- 'ANDROPOGON GERARDII'


### cbgb.us	61	Schizachyrium scoparium	It should be there, but I checked the original field dataasheet for 2012 (year_trt = 3) and it is not listed, probably overlooked at that time

working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt %in% c(2,3,4)) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover)

working_coverdat <- rbind(working_coverdat,
                          working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt ==2) %>% 
  mutate(year_trt = 3, max_cover = 60))

working_coverdat %>%
  filter(site_code == 'cbgb.us', substr(Taxon,1,5) == 'SCHIZ', plot == 61, year_trt %in% c(2,3,4)) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover, local_lifespan)

#### take only maximum cover value when multiple are present in a year (for some taxa fixes and sites with multiple sampling dates)
# QC: not sure about this, I think it makes more sense that cover from all species within plots should be measured at similar time
# QC: it does not make sense to "construct" a community using cover of species A in September and cover of species B in Summer
  group_by(site_code,year_trt,plot,Taxon) %>% 
  slice_max(order_by = max_cover,with_ties = FALSE)

working_coverdat %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS SP.'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

cover_max %>%
  filter(site_code == 'sier.us', Taxon %in%  c('BROMUS SP.'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)
### looks good in this example


### known site with two sampling dates per year
working_coverdat %>%
  filter(site_code == 'bayr.de', Taxon %in%  c('FESTUCA RUBRA'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)


cover_max %>%
  filter(site_code == 'bayr.de', Taxon %in%  c('FESTUCA RUBRA'),plot == 5) %>% 
  dplyr::select(plot, year_trt, Taxon, max_cover,local_lifespan) %>% 
  print(n=60)

## now apply Bakker fixes - previous script ID'ed that 8 initial dominants were in this fix (very low number)

cover_max_fix <- left_join(cover_max,cover_fix %>% dplyr::select(site_code,Taxon,NewTaxon),
                              by = c('site_code','Taxon')) %>% 
  mutate(NewTaxon = ifelse(is.na(NewTaxon), Taxon, NewTaxon)) %>%
  ungroup() %>% 
  dplyr::select(-Taxon) %>% 
  rename(Taxon = NewTaxon) %>% 
  group_by(site_code, block, plot, subplot, trt, year, year_trt, max_year, Taxon,live) %>%
  summarize(max_cover = sum(max_cover), .groups = "keep")

cover_max <- left_join(cover_max_fix,cover_max %>% 
                         ungroup() %>% 
                         dplyr::select(Taxon,local_provenance,functional_group,local_lifespan) %>% 
                         distinct(Taxon,.keep_all = TRUE),
                       by='Taxon')

### use relative cover to determine rank percentile

working_coverdat <- cover_max %>% 
  group_by(site_code, year_trt, plot) %>% 
  mutate(rel_cover = max_cover / sum(max_cover)) %>%  # calc rel cover
  mutate(cover = rel_cover) %>%  #use relative cover
  #mutate(cover = max_cover) %>%  # use absolute cover
  mutate(rank = min_rank(desc(rel_cover)), rich = n()) %>% 
  mutate(perc_rank = (rich + 1 - rank)/rich) %>%  #calculate percentile rank
  filter(year_trt >= 0) %>% 
  group_by(site_code) %>% 
  mutate(nplots = n_distinct(plot)) %>% # number of plots in a given site
  ungroup() %>%   
  group_by(site_code, Taxon,year_trt) %>% 
  mutate(nplot = n_distinct(plot)) %>%  # number of plots taxon is found in
  mutate(plotfreq_yr = nplot/nplots)  # calculate plot frequency

#229472 species records remaining

### calculate the frequency, average cover, average rank of a species in all pre-treatment plots
cover_pre <- working_coverdat[working_coverdat$year_trt == 0 & working_coverdat$live == 1,] %>%
  filter(year_trt == 0) %>%
  group_by(site_code, Taxon) %>%
  mutate(avgCover = mean(cover)) %>%
  mutate(avgRankPerc = mean(perc_rank)) %>%
  mutate(plotfreq = plotfreq_yr) %>%
  ungroup


### create initial dominance data.frame

### this code selects initial dominants based on most dominant species
initial_dominants <- working_coverdat %>%
  filter(year_trt == 0) %>%
  group_by(site_code,plot)  %>%
  slice_max(order_by = max_cover,with_ties = TRUE) %>%   # TRUE allows ties (multiple dominants per plot)
  mutate(initial_rel_cover = cover, initial_abs_cover = max_cover) %>% 
  dplyr::select(site_code,plot,Taxon,initial_rel_cover,initial_abs_cover) 


## merge to multi-year data and expand grid to include zero for years where spp is missing
dominant_year <- droplevels(initial_dominants[,c('site_code','plot','Taxon')]) %>% 
  left_join(working_coverdat[,c('site_code','plot','Taxon','cover','max_cover','rank','perc_rank','year_trt','max_year')],by = c('site_code','plot','Taxon')) %>% 
  filter(!is.na(site_code))

## some sites skipped a sampling year, need to weed those years out (will be all zeros now)
zero_year <- dominant_year %>% 
  group_by(site_code,year_trt) %>% 
  summarize(total_cov = sum(cover)) %>% 
  filter(total_cov == 0) %>% 
  mutate(id = paste(site_code,year_trt,sep='_'))

dominant_year <- dominant_year %>% 
  mutate(id = paste(site_code,year_trt,sep='_')) %>% 
  filter(!id %in% zero_year$id) %>% 
  left_join(
    (working_coverdat[,c('site_code','Taxon','year_trt','plotfreq_yr')] %>% 
       mutate(temp_id = paste(site_code,Taxon, year_trt,sep='_')) %>% 
       filter(!duplicated(temp_id))), by = c('site_code','year_trt','Taxon')) %>% 
  mutate(plotfreq_yr = ifelse(is.na(plotfreq_yr),0,plotfreq_yr))

scale_col <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

### put it all together
dominant_pop <- initial_dominants %>% 
  #mutate(id = paste(site_code,plot,Taxon,sep='.')) %>% 
  left_join(.,dominant_year[,c('site_code','plot','Taxon','max_cover','cover','perc_rank','year_trt')],
            by=c('site_code'='site_code','plot'='plot','Taxon'='Taxon')) %>% 
  left_join(.,cover_pre[, c('site_code','plot','trt','Taxon',"local_provenance", "local_lifespan", "functional_group",'avgCover')],
            by=c('site_code'='site_code','plot'='plot','Taxon'='Taxon')) %>% 
  left_join(.,comb[!duplicated(comb$site_code),
                   c('site_code','site_richness',"MAP_VAR_v2","MAP_v2","MAT_v2","TEMP_VAR_v2")],
            by = 'site_code') %>% 
  left_join(.,comb[comb$year_trt == 0,c('site_code','plot','vascular_live_mass','litter_mass')], by = c('site_code','plot')) %>% 
    rename(initial_live_mass = vascular_live_mass, initial_litter_mass = litter_mass) %>% 
  # left_join(.,comb_soil[,c('site_code','soil_var')],by='site_code') %>% 
  filter(trt %in% c('Control','NPK','Fence','NPK+Fence')) %>% 
  mutate(NPK = if_else(trt %in% c('NPK','NPK+Fence'),1,0), Fence = if_else(trt %in% c('Fence','NPK+Fence'),1,0))


# Can we fill in missing lifespan/provenance/func_group?
dominant_pop[dominant_pop$local_lifespan == 'NULL',] # 305

unique(dominant_pop$Taxon[dominant_pop$local_lifespan == 'NULL']) # 5

dominant_pop[dominant_pop$Taxon %in% c('NARDUS STRICTA','LOLIUM MULTIFLORUM',"ELYMUS SPICATUS"),]$local_lifespan <- 'PERENNIAL'
dominant_pop[dominant_pop$Taxon %in% c('DAUCUS CAROTA'),]$local_lifespan <- 'BIENNIAL'

#for now. poa sp. at koffler is likely perennial (poa pratensis or trivialis)
dominant_pop$local_lifespan <- ifelse(dominant_pop$local_lifespan %in% c('PERENNIAL','INDETERMINATE','NULL'), 'PERENNIAL','ANNUAL')
table(dominant_pop$local_lifespan)
dominant_pop$functional_group <- ifelse(dominant_pop$functional_group %in% c('GRASS','GRAMINOID'), 'GRAMINOID',dominant_pop$functional_group)

#provenance

unique(dominant_pop$local_provenance)
dominant_pop[dominant_pop$local_provenance == 'NULL',] # 526
unique(dominant_pop$Taxon[dominant_pop$local_provenance == 'NULL']) # 9

dominant_pop %>% filter(Taxon == 'UNKNOWN ')
working_coverdat %>% filter(site_code == 'kibber.in' & plot == 2) %>% dplyr::select(year_trt,Taxon, functional_group,rel_cover,cover) %>% 
  print(n = 100)

#going to delete this unknown as it's tied with another species for initial dominance and there's no certainty if it actually disappears in year 1 on
dominant_pop <- filter(dominant_pop, !(Taxon == 'UNKNOWN '))

dominant_pop %>% filter(Taxon == 'CAREX SP.') %>% distinct(site_code)
# no evidence of non-native carex in Oregon, especially common ones
dominant_pop[dominant_pop$Taxon == 'CAREX SP.' & dominant_pop$site_code == 'bnch.us',]$local_provenance <- 'NAT'

# no idea on provenance of Carex in India
dominant_pop[dominant_pop$Taxon == 'CAREX SP.' & dominant_pop$site_code == 'kibber.in',]$local_provenance <- 'UNK'

# daucus carota
dominant_pop %>% filter(Taxon == 'DAUCUS CAROTA') %>% distinct(site_code)

dominant_pop[dominant_pop$Taxon == 'DAUCUS CAROTA' & dominant_pop$site_code == 'mcdan.us',]$local_provenance <- 'INT'

# elymus spicata

dominant_pop %>% filter(Taxon == "ELYMUS SPICATUS") %>% distinct(site_code) %>% print(n=100)
dominant_pop[dominant_pop$Taxon == "ELYMUS SPICATUS",]$local_provenance <- 'NAT'

# "HETEROTHECA VILLOSA"
dominant_pop %>% filter(Taxon == "HETEROTHECA VILLOSA") %>% distinct(site_code) %>% print(n=100)
dominant_pop[dominant_pop$Taxon == "HETEROTHECA VILLOSA",]$local_provenance <- 'NAT'

# "MELICA NUTANS"  

dominant_pop %>% filter(Taxon == "MELICA NUTANS") %>% distinct(site_code) %>% print(n=100)
dominant_pop[dominant_pop$Taxon == "MELICA NUTANS",]$local_provenance <- 'NAT'


# "SETARIA SP."   

dominant_pop %>% filter(Taxon == "SETARIA SP.") %>% distinct(site_code) %>% print(n=100)
# likely setaria parviflora based on my time at this site (native)
dominant_pop[dominant_pop$Taxon == "SETARIA SP.",]$local_provenance <- 'NAT'



#    "NARDUS STRICTA" 
dominant_pop %>% filter(Taxon == "NARDUS STRICTA") %>% distinct(site_code) %>% print(n=100)
dominant_pop[dominant_pop$Taxon == "NARDUS STRICTA",]$local_provenance <- 'NAT'


# "AGROSTIS SP."

dominant_pop %>% filter(Taxon == "AGROSTIS SP.") %>% distinct(site_code) %>% print(n=100)
#common agrostis is native to Europe, so seems reasonable assumption
dominant_pop[dominant_pop$Taxon == "AGROSTIS SP.",]$local_provenance <- 'NAT'

dominant_pop[dominant_pop$local_provenance == 'NULL',] # 0
unique(dominant_pop$Taxon[dominant_pop$local_provenance == 'UNK']) # 4

# "AMBROSIA ARTEMISIIFOLIA"

dominant_pop %>% filter(Taxon == "AMBROSIA ARTEMISIIFOLIA") %>% distinct(site_code) %>% print(n=100)
# definitely native in MN
dominant_pop[dominant_pop$Taxon == "AMBROSIA ARTEMISIIFOLIA",]$local_provenance <- 'NAT'


# "HIERACIUM SP."           
dominant_pop %>% filter(Taxon == "HIERACIUM SP.") %>% distinct(site_code) %>% print(n=100)
sort(unique(working_coverdat$Taxon[working_coverdat$site_code == 'kbs.us']))
# not solvable


# "DANTHONIA SPICATA"    
dominant_pop %>% filter(Taxon == "DANTHONIA SPICATA") %>% distinct(site_code) %>% print(n=100)
dominant_pop[dominant_pop$Taxon == "DANTHONIA SPICATA",]$local_provenance <- 'NAT'

dominant_pop[dominant_pop$local_provenance == 'UNK',] # 23

unique(dominant_pop$functional_group)
dominant_pop[dominant_pop$functional_group == 'NULL',]
unique(dominant_pop$Taxon[dominant_pop$functional_group == 'NULL']) # 1

dominant_pop[dominant_pop$Taxon == "CAREX SP.",]$functional_group <- 'GRAMINOID'

unique(dominant_pop$Taxon[dominant_pop$functional_group == 'WOODY']) # 23
unique(dominant_pop$Taxon[dominant_pop$functional_group == 'LEGUME']) # 12
unique(dominant_pop$Taxon[dominant_pop$functional_group == 'GRAMINOID']) #127
unique(dominant_pop$Taxon[dominant_pop$functional_group == 'FORB']) #89


### add year 0 site richness
dominant_pop <- left_join(dominant_pop, 
                              working_coverdat %>% 
                                filter(year_trt == 0, live == 1, !Taxon %in% c('LICHEN','BRYOPHYTE')) %>% 
                                group_by(site_code) %>% 
                                distinct(Taxon) %>% 
                                summarize(site_rich_0 = n()),
                              by = 'site_code')
# write.csv(dominant_pop,
#                    paste0('Dominants-through-time_',
#                    Sys.Date(),'.csv'),
#                   row.names = F)

### create site charactersistics table

site_tab <- comb %>% filter(site_code %in% dominant_pop$site_code) %>% 
  group_by(site_code) %>%
  mutate(max_year = max(year_trt)) %>% filter(year_trt == max_year) %>%
  distinct(site_code, .keep_all = TRUE) %>% 
  dplyr::select(site_name, site_code, continent, country, experiment_type, latitude, longitude, MAP_v2, MAP_VAR_v2) %>% 
  left_join(.,dominant_pop %>% group_by(site_code) %>% distinct(site_code,.keep_all = TRUE) %>% dplyr::select(site_code,site_rich_0),
            by='site_code')

# write_csv(site_tab,'site_table.csv')

