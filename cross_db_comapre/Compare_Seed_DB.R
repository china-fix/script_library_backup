library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# plot the base figures (output plot_base_stack, plot_base_fill)
Interpro_all <- read.delim("./Interpro_all.Rtab", header = FALSE)
names(Interpro_all) <- c("Tag_name", "DB", "ex_confirmed", "reference")
#plot_base_stack <- ggplot(Interpro_all, aes(DB, fill=ex_confirmed))+geom_bar(position = "stack")
#plot_base_fill <- ggplot(Interpro_all, aes(DB, fill=ex_confirmed))+geom_bar(position = "fill")

# plot annotation modify
#plot_base_stack <- plot_base_stack+theme(axis.text.x = element_text(angle = 45))
#plot_base_fill <- plot_base_fill+theme(axis.text.x = element_text(angle = 45))

# analysis the raw data
#Interpro_all_summarized <- summarise(group_by(Interpro_all, DB, ex_confirmed), number= n(), percent=n()/nrow(DB))
Interpro_all_summarized_1 <- Interpro_all %>% group_by(DB, ex_confirmed) %>% summarise( number =n())
Interpro_all_summarized_stat <- summarise(group_by(Interpro_all_summarized_1, ex_confirmed), median_num=median(number))
median_num_Yes = as.numeric(filter(Interpro_all_summarized_stat, ex_confirmed == "Yes")[1,2])
median_num_No = as.numeric(filter(Interpro_all_summarized_stat, ex_confirmed == "No")[1,2])
Interpro_all_summarized_2 <- spread(Interpro_all_summarized_1, key = ex_confirmed, value = number)
Interpro_all_summarized_2 <- replace_na(Interpro_all_summarized_2, list(No=0, Yes=0))
Interpro_all_summarized_2 <- mutate(Interpro_all_summarized_2, ex_confirmed_rate = Yes/(Yes+No))
median_num_Yes_rate = median(Interpro_all_summarized_2$ex_confirmed_rate)
# plot the base figures using analyzed data
plot_base_stack <- ggplot(Interpro_all_summarized_1, aes(x= DB, y= number, fill= ex_confirmed, label = number))+geom_bar(stat = "identity")+
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45))

plot_base_fill <- ggplot(Interpro_all_summarized_2, aes(x=DB, y= ex_confirmed_rate, fill=DB)) +
  geom_bar(stat = "identity", show.legend = FALSE)+
  #geom_text(aes(label = paste("number",as.character(Yes))), size = 3, position = position_stack(vjust = 0.5))+
  geom_text(aes(label = paste(label_percent()(ex_confirmed_rate),"(",as.character(Yes),")")), size = 3, position = position_stack(vjust = 0.8))+
  theme(axis.text.x = element_text(angle = 45))

# add annotation lines
plot_base_stack <- plot_base_stack+geom_hline(yintercept = median_num_Yes, linetype="dashed", color = "blue", size=0.4)
plot_base_fill <- plot_base_fill+geom_hline(yintercept = median_num_Yes_rate, linetype="dashed", color = "blue", size=0.4)
plot_base_fill <- plot_base_fill + geom_hline(yintercept = 167/4548, linetype="dashed", color = "red", size=0.4)

###manually
#plot_base_stack + coord_cartesian(ylim = c(0,100))
# plot_base_stack + scale_x_discrete(limits = c("GENE3D", "SUPERFAMILY", "PFAM", "SMART", "TIGRFAM", "PIRSF", "SFLD", "HAMAP", "PROSITE_PROFILES", "CDD", "PRINTS", "PROSITE_PATTERNS", "MOBIDB_LITE", "COILS"))