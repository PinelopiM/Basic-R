rm(list=ls())

#Load necessary libraries
library(dplyr)        # For data manipulation
library(car)          # For Levene's test
library(ggplot2)      # For plotting
library(pheatmap)     # For heatmap visualization
library(scales)       # For color gradient in ggplot
library(FSA)          # For compact letter display and Dunn's test

#Set working directory
setwd("your_working_directory")

#Filter the data and get some descriptive statistics
my_sum2<- dat %>%
  filter(Genotype %in% c("WT",..))%>%
  filter(Analyte %in% c("4MTB",..))%>%
  filter(Experiment %in% c("1","2",..))%>%
  group_by(Genotype, Analyte, Experiment) %>%
  summarise(
    n=n(),
    mean=mean(CorrectedConcetration,na.rm = T),
    sd=sd(CorrectedConcetration,na.rm = T)
  )


#Save csv with descriptive stats
write.csv(my_sum2,"Absolute_path.csv", row.names = TRUE)


# Statistics --------------------------------------------------------------

#First, we will evaluate whether the data meet the assumptions necessary for parametric statistical tests.

## Examine normality
#Shapiro test:H0=Data are normally distributed.
shapiro <- shapiro.test(df$Concentration)
shapiro

# Creating a QQ plot as a visual assessment of normality in the data
hist(df$Concentration)
qqnorm(df$Concentration, ylab="Sample Quantiles")
qqline(df$Concentration, col="red")

##Examine homogeneity of variances with Levene test. H0=Variances of the groups being compared are equal
levene <- leveneTest(Concentration ~ Genotype*Condition, data = df )
levene

##If the assumptions are satisfied (p-value>0.05), we will proceed with the parametric tests. Otherwise, consider using non-parametric tests.

##Parametric-tests

# Two-sample t-test
ttest<-t.test(Concentration ~ Genotype, data = df,
              var.equal = TRUE)
ttest

# ONE-way ANOVa
res.aov <- aov(Concentration ~ Genotype, data = df)
AnovaTable <- summary(res.aov)
AnovaTable

#Two-way ANOVA
res.aov <- aov(Concentration ~ Genotype*Condition, data = df)
AnovaTable <- summary(res.aov)
AnovaTable

#Post-hoc Tukey test: Performs pairwise comparisons after ANOVA.
tukey.test <- TukeyHSD(res.aov)
tukey.test

tukey <- as.data.frame(tukey.test$`Genotype*Condition`)
write.csv(tukey,"Absolute_path.csv", row.names = TRUE)

#Compact letter display: Provides compact letter display for significant differences.
cldList(p.adj ~ Comparison,
        data      = tukeytest,
        threshold = 0.05)

##Non-parametric tests

#Two-sample Wilcoxon test
wtest <- wilcox.test(Concentration ~ Genotype, data = df,
                     exact = FALSE)
wtest

#Kruskal-Wallis test: non-parametric equivalent of the one-way ANOVA test
kruskal<- kruskal.test(Concentration ~ Genotype, data = df)
kruskal

#Dunn's test: Performs pairwise comparisons after Kruskal-Wallis.
DT = dunnTest(Concentration ~ Genotype,
              data=analyte,
              method="bh")
DT = DT$res
DT

# Visualizations ----------------------------------------------------------

##Boxplot with customized legend
ggplot(df) +
  # Specifying aesthetics for x-axis, y-axis, and fill color
  aes(x = Genotype, y = Concentration , fill=Condition) +
  # Adding geom_boxplot for the boxplot
  geom_boxplot() +
  # Adding labels for x and y axes
  xlab("Genotype") + ylab("Relative Concentration") +
  # Customizing plot theme
  theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    axis.text.x = element_text(face="bold", color="#993333", size=14)
  ) +
  # Setting manual fill colors for the legend
  scale_fill_manual(
    values=c("#CC6666", "#9999CC"),  # Custom fill colors
    name="Condition",  # Legend title
    breaks=c("+Pi", "-Pi"),  # Specifying breaks
    labels=c("+Pi", "-Pi")  # Specifying labels
  ) 


##Grouped barplot for visualizing gene expression data from qPCR
ggplot(df, aes(fill= Condition, y= mean, x=Genotype), title='Pi-responsive induction of GLS biosynthesis genes') + 
  # Specifying geom_bar for bars with dodge position
  geom_bar(position=position_dodge(), 
           stat="identity", width = 0.6) +
  # Adding error bars
  geom_errorbar( aes(x=Genotype, ymin=mean-sd, ymax=mean+sd),width=0.2, 
                 position=position_dodge(.6), colour="black", alpha=0.6, size=0.7) +
  # Faceting by the variable 'Gene'
  facet_wrap("Gene") +
  # Adding labels for x and y axes
  xlab("Glucosinolate Biosynthesis Genes") + ylab("Relative fold change to UBC9" )+
  # Customizing plot theme
  theme(
    plot.title = element_text(color="black", size=14, face="bold"),
    axis.title.x = element_text(color="black", size=14),
    axis.title.y = element_text(color="black", size=14),
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14),
    axis.text.x = element_text(face="bold", color="#993333", 
                               size=11))+
  # Setting fill colors using Brewer palette
  scale_fill_brewer(palette="Paired")

##Gene Expression Heatmap with Clustering for Visualizing Expression Patterns

# Extract just the numeric data into a matrix with named rows by gene
rownames(geneExp_matrix) <- geneExp_matrix$Gene_ID  # Naming rows by gene ID
geneExp_matrix <- as.matrix(geneExp_matrix[2:5])  # Extracting numeric data
head(geneExp_matrix)  # Display

# Create a green and red color palette
color <- colorRampPalette((c("light blue","black","yellow")))(75)

# Draw the heatmap
out = pheatmap(geneExp_matrix, 
               kmeans_k = NA, 
               cutree_rows = 5,
               color = color,
               scale = 'row',
               clustering_distance_rows = "euclidean",
               cluster_rows = T,
               cluster_cols = T,
               display_numbers = F,
               cellwidth = 25, cellheight = 1,
               show_rownames = F, show_colnames= T
)

# Extract clusters
data.frame(cluster = cutree(out$tree_row, k = 5))

##Dot Plot for Visualizing Functional Enrichment Analysis Result
# Reordering the 'Term' factor so that the different GO terms are ordered based on their significance
dat$Term <- factor(dat$Term, levels = dat$Term[order(dat$FDR)])

# Create the a scatter plot
ggplot(dat, aes(x=FDR, y=Term)) +  
  geom_point(aes(size=Enrichment, color=Count), show.legend=TRUE) +
  # Adding aesthetics for point size and color
  scale_color_gradient(low="red3", high="royalblue3") +  # Setting color gradient
  # Customizing axis labels and legend position
  xlab('FDR') +  # X-axis label
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "bottom"  # Legend position
  )
