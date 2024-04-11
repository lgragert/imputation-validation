# Manhattan plots of the TRS by AA position per population group by subject

if (!require("tidyverse"))
  install.packages("tidyverse")
  install.packages("ggtext")
if (!require('ggplot2'))
  install.packages('ggplot2')
if (!require('ggrepel'))
  install.packages('ggrepel')

# Start by dividing out the population groups then take the AA positions
setwd("./")

trs_AFA<-read.csv("HLA_AA_TRS_Average_AFA.csv")  # Have separate files for each population AFA, ASN, CAU, HIS, MLT, NAM
trs_ASN<-read.csv("HLA_AA_TRS_Average_ASN.csv")
trs_CAU<-read.csv("HLA_AA_TRS_Average_CAU.csv")
trs_HIS<-read.csv("HLA_AA_TRS_Average_HIS.csv")
trs_MLT<-read.csv("HLA_AA_TRS_Average_MLT.csv")
trs_NAM<-read.csv("HLA_AA_TRS_Average_NAM.csv")

# Subset data based on population and subject type
afa_subjs = trs_AFA %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'AFA'))
afa_recip = trs_AFA %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'AFA'))
afa_donor = trs_AFA %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'AFA'))

asn_subjs = trs_ASN %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'ASN'))
asn_recip = trs_ASN %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'ASN'))
asn_donor = trs_ASN %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'ASN'))

cau_subjs = trs_CAU %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'CAU'))
cau_recip = trs_CAU %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'CAU'))
cau_donor = trs_CAU %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'CAU'))

his_subjs = trs_HIS %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'HIS'))
his_recip = trs_HIS %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'HIS'))
his_donor = trs_HIS %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'HIS'))

mlt_subjs = trs_MLT %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'MLT'))
mlt_recip = trs_MLT %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'MLT'))
mlt_donor = trs_MLT %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'MLT'))

nam_subjs = trs_NAM %>% filter((Subject_Type =='SUBJECT') & (Ethnicity == 'NAM'))
nam_recip = trs_NAM %>% filter((Subject_Type =='RECIP') & (Ethnicity == 'NAM'))
nam_donor = trs_NAM %>% filter((Subject_Type =='DONOR') & (Ethnicity == 'NAM'))


# GGplot function found on https://r-graph-gallery.com/101_Manhattan_plot.html
manhattan_plot <- function(sub_typ, ethnicity, subj, tolerance=0.04) {
  # Make Loci names L1, L2, L3,... so that you can have it in genomic order
  sub_typ$lgroup <- 'NA'
  sub_typ[is.element(sub_typ$Locus,'A'),"lgroup"] <- 'L1'
  sub_typ[is.element(sub_typ$Locus,'C'),"lgroup"] <- 'L2'
  sub_typ[is.element(sub_typ$Locus,'B'),"lgroup"] <- 'L3'
  sub_typ[is.element(sub_typ$Locus,'DRB345'),"lgroup"] <- 'L4'
  sub_typ[is.element(sub_typ$Locus,'DRB1'),"lgroup"] <- 'L5'
  sub_typ[is.element(sub_typ$Locus,'DQA1'),"lgroup"] <- 'L6'
  sub_typ[is.element(sub_typ$Locus,'DQB1'),"lgroup"] <- 'L7'
  sub_typ[is.element(sub_typ$Locus,'DPA1'),"lgroup"] <- 'L8'
  sub_typ[is.element(sub_typ$Locus,'DPB1'),"lgroup"] <- 'L9'
  # Prepare the data
  don <- sub_typ %>%
    
    # Compute locus size by AA position
    group_by(lgroup) %>%
    summarise(chr_len=max(AA_Position)) %>%
    
    # Calculate cumulative position of each loci
    mutate(tot=cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(sub_typ, ., by=c('lgroup'='lgroup')) %>%
    
    # Add a cumulative position of each AA position
    arrange(lgroup, AA_Position) %>%
    mutate(Loci=AA_Position+tot) %>%
    
    # Add highlight to big variance in TRS Average
    mutate(is_annotate=ifelse((1- TRS_Average)>tolerance, "yes", "no"))
  
  # Make the x-axes
  axisdf = don %>%
    group_by(Locus) %>%
    summarize(center=(max(Loci) + min(Loci)) / 2)
  
  max_y_val <- 0.33
  
  # Plot the data
  ggplot(don, aes(x=Loci, y=(1-TRS_Average))) +
    
    # Show all points
    geom_point(aes(color=as.factor(lgroup)), alpha=0.8, size=1.3) +
    scale_color_manual(values=rep(c("blue", "red"), 22 )) +
    
    # custom X axis:
    scale_x_continuous(label=axisdf$Locus, breaks=axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(0, max_y_val + 0.01) ) +     # remove space between plot area and x axis
    
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=AA_Position), size=2, label.size=NA, max.overlaps=nrow(subset(don, is_annotate=="yes"))) +
    
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()) +
    ggtitle(paste("Average Variance in TRS by AA Position per Loci\nfor", ethnicity, subj)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# AFA population
manhattan_plot(afa_subjs, "AFA", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_AFAsubj.jpg", width=18.00, height=4.00)
manhattan_plot(afa_donor, 'AFA', "Donors")
ggsave("AA_TRS_Manhattan_AFAdonors.jpg",width=18.00, height=4.00)
manhattan_plot(afa_recip, 'AFA', 'Recipients')
ggsave("AA_TRS_Manhattan_AFArecips.jpg",width=18.00, height=4.00)

# ASN population
manhattan_plot(asn_subjs, "ASN", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_ASNsubj.jpg", width=18.00, height=4.00)
manhattan_plot(asn_donor, 'ASN', "Donors")
ggsave("AA_TRS_Manhattan_ASNdonors.jpg",width=18.00, height=4.00)
manhattan_plot(asn_recip, 'ASN', 'Recipients')
ggsave("AA_TRS_Manhattan_ASNrecips.jpg",width=18.00, height=4.00)

# CAU population
manhattan_plot(cau_subjs, "CAU", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_CAUsubj.jpg", width=18.00, height=4.00)
manhattan_plot(cau_donor, 'CAU', "Donors")
ggsave("AA_TRS_Manhattan_CAUdonors.jpg",width=18.00, height=4.00)
manhattan_plot(cau_recip, 'CAU', 'Recipients')
ggsave("AA_TRS_Manhattan_CAUrecips.jpg",width=18.00, height=4.00)

# HIS population
manhattan_plot(his_subjs, "HIS", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_HISsubj.jpg", width=18.00, height=4.00)
manhattan_plot(his_donor, 'HIS', "Donors")
ggsave("AA_TRS_Manhattan_HISdonors.jpg",width=18.00, height=4.00)
manhattan_plot(his_recip, 'HIS', 'Recipients')
ggsave("AA_TRS_Manhattan_HISrecips.jpg",width=18.00, height=4.00)

# MLT population
manhattan_plot(mlt_subjs, "MLT", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_MLTsubj.jpg", width=18.00, height=4.00)
manhattan_plot(mlt_donor, 'MLT', "Donors")
ggsave("AA_TRS_Manhattan_MLTdonors.jpg",width=18.00, height=4.00)
manhattan_plot(mlt_recip, 'MLT', 'Recipients')
ggsave("AA_TRS_Manhattan_MLTrecips.jpg",width=18.00, height=4.00)

# NAM population
manhattan_plot(nam_subjs, "NAM", "Donors and Recipients")
ggsave("AA_TRS_Manhattan_NAMsubj.jpg", width=18.00, height=4.00)
manhattan_plot(nam_donor, 'NAM', "Donors")
ggsave("AA_TRS_Manhattan_NAMdonors.jpg",width=18.00, height=4.00)
manhattan_plot(nam_recip, 'NAM', 'Recipients')
ggsave("AA_TRS_Manhattan_NAMrecips.jpg",width=18.00, height=4.00)
