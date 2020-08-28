#pound sign means the text fowlling it is commented, R willl not run it.

#assign a value to an object
#it will create a object
James_weight <- 250 
James_weight = 250 
James_weight

book_name <- "Harry Potter"
book_name

#Creat vectors 
b <- c(2, 4, 5)   # c means concatenat 
b[1]
b[b>1]
b[b==2]
a <- 34
a[1]
a[2]
fruits <- c("apple","orange","strawberry")
Fruits[1]      # object name is case sensitive
# What if I create a vector contains different types of values?
fruits <- c("apple","orange",1,2,TRUE)
fruits

#creat a list
myList <- list(l1=2, l2=c(1,2,4), l3=fruits)
myList
myList$l1
myList[["l2"]]

#create a matrix
m <- matrix(c(1,2,3,4,5,6), nrow=3, ncol=3)
m
m[2,3]
m <- cbind(c(1,2,3), c(4,5,6))
m
m <- rbind(c(1,2,3), c(4,5,6))
m
m <- matrix(c(1,2,3,4,5,6),
            nrow=3, ncol=3,
            dimnames = list(c("r1","r2","r3"), c("c1","c2","c3")))
m[,"c1"]
m["r1",]
m["r1","c1"]
rownames(m)
rownames(m) <- c("m1","m2","m3")
colnames(m) <- c("n1","n2","n3")


#create the dataframe 
d <- data.frame(c1=c(1,2,3), c2=c("a","d","c"))
d
d$c1 
d$c2
d[["c1"]]
d[,1]
d[1,2]
d <- data.frame(c1=c(1,2,3), c2=c("a","d","c"), row.names = c("r1","r2","r3"))
d
d["r1","c1"]
m <- as.matrix(d)
d <- as.data.frame(m)
d[d$c1==1,]


#create a R function
?c
help(list)
help(help) # or ?help

cal_square  <- function(x) {
  y <- x * x 
  return (y)  ## can be omitted.
}
cal_square
cal_square(2)

################# practice part #################
# Read the dataset into R
#check the working directory
getwd()
setwd("/Users/xiao/Desktop/NCSU_R/R")
covid19_data <- read.csv("RawData_plasma metabolites.csv", header=TRUE, row.names=1)
#or using the full path of the file
covid19_data <- read.csv("/Users/xiao/Desktop/NCSU_R/R/RawData_plasma metabolites.csv", header = TRUE, row.names=1)
group_data <- read.csv("group_info.csv", header=T,row.names=1)
#or
rawdata <- read.csv("RawData_plasma metabolites and group.csv",header=T,row.names=1)
covid19_data <- rawdata[,-1]
group_data <- rawdata[,1]
group_data <- rawdata$GroupID
group_data
group_data <- rawdata[,1,drop=FALSE]
group_data


View(group_data)
head(group_data)
tail(group_data, n=10)

covid19_data <- data.matrix(covid19_data)

# log2 transformation
covid19_log2_data <- log2(covid19_data+1)

#summarize the data
summary(covid19_data)
summary(group_data)

#distribution of metabolites to look the range and shape of the data
hist(covid19_data[1,],50)
hist(covid19_log2_data[1,],50)

#correlation 
cor(covid19_log2_data[1,], covid19_log2_data[2,])
cor(covid19_log2_data)

#t-test using R 
t.test(covid19_log2_data[1:26,1], covid19_log2_data[27:76,1])
t.test(covid19_log2_data[group_data$GroupID == "CONTROL", 1],
       covid19_log2_data[group_data$GroupID != "CONTROL", 1])
?t.test

# define a function to perform t-test 
myttest <- function(idx){
           res <- t.test(x = covid19_log2_data[group_data$GroupID == "CONTROL", idx],
                         y = covid19_log2_data[group_data$GroupID != "CONTROL", idx])
           return (res$p.value)
}
# using the "for" loop to perform t-test
pvals <- c()
for(i in 1:ncol(covid19_log2_data)){
  p <- myttest(i)
  pvals <- c(pvals,p)
}
pvals

# apply the function to our data,                             
pvals <- sapply(X=1:ncol(covid19_log2_data),FUN=myttest)
pvals

# define a function to calculate log2 fold change, the magnitude of changes
cal_log2fc <- function(idx){
           control_mean <- mean(covid19_data[group_data$GroupID == "CONTROL",idx])
           infected_mean <- mean(covid19_data[group_data$GroupID != "CONTROL",idx])
           return (log2(infected_mean / control_mean))
}

# apply the function to our data
log2_fcs <- sapply(X=1:ncol(covid19_data), FUN=cal_log2fc)

# valcano plot
plot(x = log2_fcs , y = -log10(pvals))

###the valcanoplot is not pretty? try use the package developped by others
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
plot_data <- data.frame(log2FoldChange = log2_fcs,
                        pvalue=pvals,
                        row.names = colnames(covid19_log2_data))
EnhancedVolcano(plot_data,
                lab = rownames(plot_data),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 1e-2,
                FCcutoff =2,
                colAlpha = 1,
                )
?EnhancedVolcano
#output the dataset 
output_data <- cbind(log2_fcs, pvals)
rownames(output_data) <- colnames(covid19_data)
write.csv(output_data, "ttest_result.csv")

#select the significant molecules
select <- (abs(log2_fcs) >= 1) & (pvals < 0.01)
write.csv(output_data[select,], "ttest_result.csv")

#one-way anova
?oneway.test
anova_test <- function(i){
  anova_data <- data.frame(value = covid19_data[,i],
                           group = group_data$GroupID)
  anova_res <- oneway.test(value~group, anova_data)
  return(anova_res$p.value)
}
anova_pvals <- sapply(X=1:ncol(covid19_log2_data),FUN=anova_test)
#mean value in each group
group_means <- apply(covid19_log2_data,2,function(x) by(x, group_data$GroupID, mean))
#sd in each group
group_sds <- apply(covid19_log2_data,2,function(x) by(x, group_data$GroupID, sd))

#ggplot2
#install.packages("ggplot2")
#more packages on CRAN:https://cran.r-project.org/

#load the package
library(ggplot2)

#generate the data for mapping to ggplot2 aesthetics
plot_data <- data.frame(value=covid19_data[,1],  group=group_data[,1])

#boxplot, barplot, violin, point, jitter
ggplot(plot_data, aes(x=group, y=value)) +
  geom_boxplot() 
#  geom_col()
#  geom_violin() 
#  geom_point()
#  geom_jitter()

#barplot and add the errorbar
ggplot(plot_data, aes(x=group, y=value,fill=group)) +
  stat_summary(geom = "col", fun = mean)+
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
ggplot(plot_data2, aes(x=group, y=value, fill=Var2)) +
  stat_summary(geom = "col", fun = mean, position= "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position="dodge")
#calculating the relative intensities
for(m in unique(plot_data2$Var2)){
  control_mean <- mean(plot_data2[(plot_data2$Var2 == m) & (plot_data2$group=="CONTROL"),"value"])
  plot_data2[plot_data2$Var2 == m, "relative"] <- 
    plot_data2[plot_data2$Var2 == m, "value"] / control_mean
}
ggplot(plot_data2, aes(x=Var2, y= relative, fill=group)) +
  stat_summary(geom = "col", fun = mean, position= "dodge")+
  stat_summary(geom = "errorbar", fun.data = mean_se, position="dodge") + 
  scale_fill_manual(values =   c("#ffbfbf","#ff4040","#bf0000","#6f0000")) 
#  theme_bw()



#creat heatmap using R
## you can try to refine this graph by change the settings.
library(pheatmap)
my_breaks <- seq(-3, 3, length.out = 100)
pheatmap(covid19_log2_data,
         show_colnames = F,show_rownames = T,
         cluster_rows = F,
         annotation_row = group_data,
         scale="column",
         breaks = my_breaks
         )
