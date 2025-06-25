#Basic Cohort Imaging
require(readr)
require(dplyr)
require(tidyverse)
require(survminer)
require(survival)
require(KMsurv)
require(ggplot2)
require(RColorBrewer)

part <- read_tsv("participants.ji0tedp.tsv")
samples <-read_tsv("samples.jjhlq.tsv")
biospec <-read_tsv("biospecimens.87egnq.tsv")
treat <-read_tsv("treatments.8q5vqe.tsv")


part%>%
  mutate(death=case_when(
    vital_status=="dead" ~ 1,
    vital_status=="unknown"~ 0,
    vital_status=="not reported"~ 0,
    vital_status=="alive"~ 0))->part


part%>%
  mutate(decade=case_when(
    year_of_diagnosis >=1970 & year_of_diagnosis<1980 ~ "70s",
    year_of_diagnosis >=1980 & year_of_diagnosis<1990 ~ "80s",
    year_of_diagnosis >=1990 & year_of_diagnosis<2000 ~ "90s",
    year_of_diagnosis >=2000 & year_of_diagnosis<2010 ~ "00s",
    year_of_diagnosis >=2010 & year_of_diagnosis<2020 ~ "10s",
    year_of_diagnosis >=2020 & year_of_diagnosis<2030 ~ "20s",
  ))->part

part$decade<- factor(part$decade, levels=c("70s","80s","90s","00s","10s","20s"))
part%>%
  filter(tumor_primary_site=="BREAST (C50)")->breast
part%>%
  filter(tumor_primary_site=="GALLBLADDER (C23)")->gallbladder
part%>%
  filter(tumor_primary_site=="COLON (C18)")->colon
part%>%
  filter(tumor_primary_site=="SKIN (C44)")->skin

part%>%
  mutate(fund = if_else(grepl("non-IBM funded", notes),"No", "Yes"))->part
part%>%
  group_by(tumor_primary_site)%>%
  mutate(total=n())%>%
  arrange(desc(total))->tester
km1<-survfit(Surv(death_date_dfd, death)~tumor_primary_site, data=part)
ggsurvplot(km1)
ggsurvplot_facet(km1, data=part, facet.by='gender')

km3<-survfit(Surv(death_date_dfd, death)~tumor_primary_site, data=part)
ggsurvplot_facet(km3, data=part, facet.by='fund')


km2<-survfit(Surv(death_date_dfd, death)~gender, data=part)
ggsurvplot(km2)

km3<-survfit(Surv(death_date_dfd, death)~tumor_morphology, data=breast)
ggsurvplot(km3)
ggplot(data=part)+
  geom_bar(aes(x=tumor_primary_site, fill=tumor_primary_site))+
  scale_fill_manual("Tumor Location", values=c("lightpink","deeppink","indianred","darkred", "darkorange", "orange",
                                               "yellow","lemonchiffon", "greenyellow", "green","palegreen2","turquoise",
                                               "skyblue","dodgerblue", "blueviolet", "darkorchid", "mediumorchid"))+
  labs(x="Tumor Site", y="Number of Patients", title="Number of Patients by Tumor Location")+
  theme(axis.text.x=element_text(angle= 50, hjust = 1))

ggplot(data=tester)+
  stat_summary(aes(x=year_of_diagnosis, y=total, color=tumor_primary_site), geom="line", fun="sum")+
  scale_color_manual("Tumor Location", values=c("lightpink","deeppink","indianred","darkred", "darkorange", "orange",
                                               "yellow","lemonchiffon", "greenyellow", "green","palegreen2","turquoise",
                                               "skyblue","dodgerblue", "blueviolet", "darkorchid", "mediumorchid"))+
  labs(x="Tumor Site", y="Number of Patients", title="Number of Patients by Tumor Location")+
  theme(axis.text.x=element_text(angle= 50, hjust = 1))

ggplot(data=part)+
  geom_bar(aes(x=decade, fill=tumor_primary_site), position="stack")+
  scale_fill_manual("Tumor Location",  values=c("lightpink","deeppink","indianred","darkred", "darkorange", "orange",
                                                 "yellow","lemonchiffon", "greenyellow", "green","palegreen2","turquoise",
                                                 "skyblue","dodgerblue", "blueviolet", "darkorchid", "mediumorchid"))

ggplot(data=part)+
  geom_bar(aes(x=decade, fill=tumor_primary_site), position="fill")+
  scale_fill_manual("Tumor Location",  values=c("lightpink","deeppink","indianred","darkred", "darkorange", "orange",
                                                "yellow","lemonchiffon", "greenyellow", "green","palegreen2","turquoise",
                                                "skyblue","dodgerblue", "blueviolet", "darkorchid", "mediumorchid"))
ggplot(data=breast)+
  geom_bar(aes(x=decade, fill=tumor_morphology), position="fill")

ggplot(data=treat)+
  geom_bar(aes(x=stop_reason, fill=categories))

treat%>%
  filter("Progression and relapse" == stop_reason)->relapse
ggplot(data=treat)+
  geom_boxplot(aes(x=stop_reason))+
  theme(axis.text.x=element_text(angle= 50, hjust = 1))

#carnoma
part%>%
  filter("C22 (8160/3 CHOLANGIOCARCINOMA)" == tumor_morphology)->carcinoma
part%>%
  mutate(carc= case_when(tumor_morphology == "C22 (8160/3 CHOLANGIOCARCINOMA)"~ 1,
                         tumor_morphology != "C22 (8160/3 CHOLANGIOCARCINOMA)"~0))->part
         
carcsurv<-survfit(Surv(death_date_dfd, death)~tumor_primary_site, data=part)
ggsurvplot(carcsurv)
ggsurvplot_facet(carcsurv, data=part, facet.by='carc')

ggplot(data=treat)+
  stat_summary(aes(x))
  
autopsy<-read.table("/Users/ksweetna/Juric_Rapid_Autopsy_MASTER-participants.txt", sep="\t", fill=TRUE, header=TRUE)

part %>%
  mutate(death_date_dfd = replace_na(death_date_dfd, 30000))->part

km1<-survfit(Surv(death_date_dfd, death)~tumor_primary_site, data=part)
ggsurvplot(km1)
ggsurvplot_facet(km1, data=part, facet.by='gender')
