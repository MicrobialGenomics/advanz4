library(reshape2)
library(gtools)
# library(plyr)
# library(gdata)
# library(xlsx)
# install.packages("gdata")
# read table


#Create dummy table<-

Plot_FollowUp<-function(data, patientColumn, TimePoints ,Dates){
  ## May need a bit of tweaking
  
  
  require(lubridate)
  

  nlevs<-data %>% pull(!!TimePoints) %>% unique() %>% length()
  
  myPlot<-data %>%
    select(Patients = patientColumn, 
           Time = TimePoints, 
           Date = Dates) %>%
    group_by(Patients) %>% 
    dplyr::mutate(Complete=ifelse(n() == nlevs, "yes","no")) %>%
    reshape2::melt(id.vars = c("Patients", "Time", "Complete"), measure.vars = "Date", color="Complete") %>%

    ggplot(.,aes(x=as.Date(value, origin= value-365), y=as.character(Patients), group = Patients, color=Complete)) +
    geom_text(aes(label=as.character(Time), size="size"), show.legend = F) + 
    geom_line() +
    theme_classic()+
    scale_color_brewer(palette = "Set1")+
    ggtitle("Samples by patient and date")+
    xlab("Extraction date")+
    ylab("patients") +
    # scale_shape_manual(values=c(4,19,20))+
    scale_linetype_manual(values=c("dotted","solid"))+
    theme(
      panel.grid.major.y = element_line(colour = "grey", size=0.5),
      # axis.line.y= element_line(color="black",),
      panel.border = element_blank(),
      axis.text.y = element_text(size=6, color="black")
    )
  
  return(myPlot)
}







# 
# ggpie <- function (data) 
# {
#   # prepare name
#   name<-deparse( substitute(data) )
#   countDF<-data.frame(Complete = data) %>%
#     group_by(Complete) %>%
#     dplyr::summarise(count = n()) %>%
#     # arrange(desc(count)) %>%
#     mutate(Perc = paste(round(100*count/sum(count), 2), "%", sep=" "),
#            breaks = sapply(sort(count, decreasing = T), function(x){
#              cumsum(x) - x/2
#            }))
#   
#   
#   # prepare percents for legend
#   count<-table(as.factor(data)) 
#   prop.table( tmp.count1 ) * 100 -> tmp.percent1 ;
#   tmp.percent1 <-sapply(tmp.percent1, round,2)
#   paste( tmp.percent1, " %", sep = "" ) -> tmp.percent2 ;
# 
#   as.vector(tmp.count1) -> tmp.count1 ;
#   
#   # find breaks for legend
#   tmp.count2 <- sort(tmp.count1, decreasing = T);
#   rev( cumsum( tmp.count2 ) - (tmp.count2 / 2) ) -> tmp.breaks1 ;
#   
#   # prepare data
#   data.frame( vector1 = tmp.count1, names1 = names(tmp.percent1) ) -> tmp.df1 ;
#   
#   
#   # plot data
#   tmp.graph1 <- ggplot(tmp.df1, aes(x = 1, y = vector1, fill = names1 ) ) +
#     geom_bar(stat = "identity", color = "black" ) +
#     guides( fill = guide_legend(override.aes = list( colour = NA ) ) ) +
#     # coord_polar( theta = "y" ) +
#     theme(axis.ticks = element_blank(),
#           axis.text.y = element_blank(),
#           axis.text.x = element_text( colour = "black"),
#           axis.title = element_blank(),
#           plot.title = element_text( hjust = 0.5, vjust = 0.5) ) +
#     scale_y_continuous( breaks = tmp.breaks1, labels = tmp.percent2 ) +   
#     ggtitle( name ) #+ 
#     # scale_fill_grey( name = "") ;
#   
#   return( tmp.graph1 )
# 
# }
# 
# 
# ggpie(followup)
# ggplot(proportions, aes(x="",y=Freq, fill=Var1)) + geom_bar(width=2, stat="identity") + 
#   theme_void()+
#   ggtitle("proportion of completed patients")+
#   geom_text(aes(y= c(240,260), x=c(0,-2),label = perc))+
#   scale_fill_brewer(type = "qual",name="Completed patients",labels=c("1 visit left", "Complete", "Incomplete"))+
#   theme(axis.text.y = element_blank(),
#         panel.border = element_blank(),
#         axis.title.x = element_blank(),
#         plot.title = element_text(hjust=0.5)) +
#    coord_polar("y", start=0)
# 
# 
# write.xlsx(metadata_f, "Metadata/metadata_Run2.xlsx")
# 
# metadata<- read.xlsx("Metadata/metadata_Run2.xlsx",1)
# 
# metadata_f<-metadata[!is.na(metadata$SampleID),]
# metadata_f$SampleID<-gsub(" ","",metadata_f$SampleID)
# 
# 
# metadata1$SampleID<-gsub(" ","",metadata1$SampleID)
# 
# metadata_f$Run<-ifelse(metadata_f$SampleID %in% metadata1$SampleID, "Run1","Run2")
# metadata_r2<-metadata_f[metadata_f$Run == "Run2",]
# length(metadata_r2$SampleID)
# sort(metadata_r2$SampleID)
# 
# write.csv(metadata_f, file = "Metadata/metadata_Run2.csv",row.names = F, quote=T)
# metatest<-read.csv("metadata (10).csv")
