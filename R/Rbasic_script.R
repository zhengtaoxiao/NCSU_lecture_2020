# Read the dataset into R
covid19_data <- read.csv("RawData_plasma metabolites.csv", header=TRUE, row.names=1)

# log2 transformation
covid19_log2_data <- log2(covid19_data+1)

# define a function to perform t-test 
myttest <- function(idx){
           res <- t.test(x = covid19_log2_data[1:26, idx],
                         y = covid19_log2_data[27:76, idx])
           return (res$p.value)
}
           
# apply the function to our data                               
pvals <- sapply(X=1:ncol(covid19_log2_data),FUN=myttest)

# define a function to calculate log2 fold change
cal_log2fc <- function(idx){
           control_mean <- mean(covid19_data[1:26,idx])
           infected_mean <- mean(covid19_data[27:76,idx])
           return (log2(infected_mean / control_mean))
}

# apply the function to our data
log2_fcs <- sapply(X=1:ncol(covid19_data), FUN=cal_log2fc)

# valcano plot
plot(x = log2_fcs , y = -log10(pvals))

#output the dataset 
write.csv(cbind(log2_fcs, pvals), "ttest_result.csv")

#ggplot2
#install.packages("ggplot2")

#load the package
library(ggplot2)

#read the group data
group_data <- read.csv("group_info.csv", header=T,row.names=1)

#generate the data for mapping to ggplot2 aesthetics
plot_data <- data.frame(value=covid19_data[,1],  group=group_data[,1])

#boxplot, barplot, violin, point, jitter
ggplot(plot_data, aes(x=group, y=value)) + 
  geom_boxplot() 
#barplot and add the errorbar
ggplot(plot_data, aes(x=group, y=value,fill=group)) +
  stat_summary(geom = "bar", fun = mean)+
  stat_summary(geom = "errorbar", fun.data = mean_se, width=0.5)
#plot the metabolites in TCA pathway
plot_data2 <- covid19_data[,c("Citric.acid","Isocitric.acid",
                               "Fumaric.acid", "Succinic.acid")]
head(plot_data2)
#load the package for melting the dataframe
library(reshape2)
plot_data2 <- melt(plot_data2) #Convert an object 
head(plot_data2)
plot_data2$group <- rep(group_data[,1],4)
#plot the TCA metabolits
ggplot(plot_data2, aes(x=group, y=value, fill=variable)) +
  stat_summary(geom = "bar", fun = mean, position= "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position="dodge")
#calculating the relative intensities
for(m in unique(plot_data2$variable)){
  control_mean <- mean(plot_data2[(plot_data2$variable == m) & (plot_data2$group=="CONTROL"),"value"])
  plot_data2[plot_data2$variable == m, "relative"] <- 
    plot_data2[plot_data2$variable == m, "value"] / control_mean
}
ggplot(plot_data2, aes(x=variable, y= relative, fill=group)) +
  stat_summary(geom = "bar", fun = mean, position= "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_se,       position="dodge") + 
  scale_fill_manual(values =   c("#ffbfbf","#ff4040","#bf0000","#6f0000"))
