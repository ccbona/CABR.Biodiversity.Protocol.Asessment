library(ggplot2)
library(vegan)
library(tidyr)
library(tidyverse)
library(NatParksPalettes)
library(patchwork)
library(readxl)
library(janitor)
library(viridis)

#create ggplot theme
clean <- theme(panel.background = element_rect(fill = NA)) + 
  theme(axis.line = element_line(linetype = "solid"))+ 
  theme(text = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))

# Read in and organize data -----------------------------------------------

#set working directory
setwd("~/SIP Cabrillo/Biodiversity Monitoring/SIP Summer 23 - Modified Biodiversity/data/mod_biodiversity_nps/Data")

#read csv
ucsc <- read_csv('ucsc.csv')
colnames(ucsc)[1] ="site" #Rename column
ucsc <- ucsc %>% group_by(species_lump) %>% filter(sum(percent_cover) > 0) %>% ungroup() #Remove zeros

tides <- read_csv('tide height.csv')
tides$'Tidal Height' <- as.numeric(tides$'Tidal Height')


###Widen and clean UCSC dataset (1990)

ucsc_wide1990 <- ucsc %>% 
  select(site, year, name_1990, percent_cover) %>%
  pivot_wider(names_from = name_1990, values_from = percent_cover, values_fn=sum) %>%
  select(where(~ any(. != 0)))


#Create wide "target" dataframe
target_90[is.na(target_90)] <- 0
target_wide <- target_90 %>%
  select(site_name, type, plot_code, survey_year, scientific_name, pct_cover) %>%
  pivot_wider(names_from = scientific_name, values_from = pct_cover, values_fn=sum) 
colnames(target_wide)[1] = "site" #Rename site
colnames(target_wide)[4] = "year" #Rename year

target_wide_comb <- target_90 %>%
  select(site_name, survey_year, scientific_name, pct_cover) %>%
  pivot_wider(names_from = scientific_name, values_from = pct_cover, values_fn=sum)#Create wide "target" dataframe
colnames(target_wide_comb)[1] = "site" #Rename site
colnames(target_wide_comb)[2] = "year" #Rename year
target_wide_comb['program'] = 'LTM' #Add program column to target_wide
ucsc_wide1990['program'] = 'CBS' #Add program column to ucsc_wide1990
ucsc_wide1990['Unscorable'] = 0 #Add Unscorable to ucsc_wide1990
comb <- rbind(target_wide_comb, ucsc_wide1990)

comb_limited <- comb %>% filter(site=="Cabrillo I" | site=="Cabrillo III") %>% filter(year == 2002 | year==2004 & site=="Cabrillo I" | year==2020 & site=="Cabrillo I" | year==2012 & site=="Cabrillo III" | year==2017 & site=="Cabrillo III")

# Create taxa table and figure --------------------------------------------

crosswalk <- read_csv('crosswalk_full.csv') #Read Crosswalk CSV
crosswalk <- crosswalk %>% rename("LTM 1990" = "name_1990", "2000" = "name_2000", "LTM current" = "species_code", "CBS" = "ucsc") #Rename columns
wide <- crosswalk %>% 
  pivot_longer(
    cols = c(10:17), 
    names_to = "program", 
    values_to = "name") #Create wide version of crosswalk

totals <- wide %>% 
  select(mobility, program, name, type) %>% 
  drop_na() %>%  
  distinct(program, name, .keep_all = TRUE) %>% 
  count(program,type) %>% 
  filter(program == "LTM 1990" |program == "CBS" | program == "Zedler" | program == "LTM current") %>% 
  filter(type == "sessile invertebrates" | type == "mobile invertebrates" | type == "algae")%>%
  mutate(type = fct_relevel(type, c("mobile invertebrates", "sessile invertebrates", "algae"))) #Filter for programs of interest and relevel taxa types

taxa_figure <- ggplot(totals, mapping=aes(x= reorder(program, -n), y = n, fill=type)) + 
  geom_col(width=1) + 
  clean+labs(x = "Program", y = "Number of Taxa")+ 
  scale_fill_manual(values=c("#212E52","#FCC893","#D8511D")) #Create figure
taxa_figure


# Create stacked bar figures ----------------------------------------------

ltm <- ggplot(data = target_90, mapping = aes(x = survey_year, y = pct_cover)) + geom_col(position="fill", aes(fill = scientific_name), width=1) + clean + scale_fill_manual(values=natparks.pals("Acadia", 9))+labs(x = "Year", y = "Percent Cover")+labs(title = "Long-Term Monitoring Program")+ scale_y_continuous(labels = scales::percent)+xlim(1990,2021)

cbs <- ggplot((data = ucsc), mapping = aes(x = year, y = percent_cover, fill = name_1990)) + geom_col(position = "fill", width=1) + clean + scale_fill_manual(values=natparks.pals("Acadia", 9))+labs(x = "Year", y = "Percent Cover")+labs(title = "Coastal Biodiversity Survey")+xlim(1990,2021)+ scale_y_continuous(labels = scales::percent)

grid.arrange(ltm,cbs, ncol=1)

# NMDS --------------------------------------------------------------------

comb_limited[is.na(comb_limited)] <- 0
comb_vegan <- comb_limited[,c(3:11)]

NMDScomb=metaMDS(comb_vegan)
stressplot(NMDScomb)
plot(NMDScomb)
scores(NMDScomb)

data.scores <- as.data.frame(scores(NMDScomb, "sites")) #extract site scores
data.scores$site <- rownames(data.scores) #create a column of sites, from rownames of data.scores
head(data.scores) #view data

species.scores <- as.data.frame(scores(NMDScomb, "species"))  #extract species scores
species.scores$species <- rownames(species.scores)  #create a column of sites, from the rownames of species.scores
head(species.scores)  #view data

comb_limited$NMDS1 <- data.scores$NMDS1
comb_limited$NMDS2 <- data.scores$NMDS2

ggplot(comb_limited)+geom_point(mapping=aes(NMDS1, NMDS2, color=program, shape=site)) +geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),size=3,alpha=0.5)+ clean+stat_ellipse(mapping=aes(NMDS1, NMDS2, color=program),type = "norm")+ 
  scale_color_manual(values=c("#212E52","#D8511D"))

adonis2(comb_vegan ~ comb_limited$year * comb_limited$site * comb_limited$program)

