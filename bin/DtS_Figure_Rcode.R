
# Required Packages

#############################################################################
library(dplyr)
library(ggplot2)
library(ggpubr)
library(eulerr)




# Figure 2. AD SRM Assay

#############################################################################
comp <- read.csv("comparison.csv", header = TRUE)
comp <- as.data.frame(comp)
ADprot <- comp %>% distinct()
ADprot <- ADprot$Protein

### a) Assay overlap
cpep <- comp %>% filter(Group=="C") %>% select(Peptide)
dpep <- comp %>% filter(Group=="D") %>% select(Peptide)

intpep <- merge(cpep,dpep)
uc <- cpep %>% filter(!Peptide %in% dpep$Peptide)
ud <- dpep %>% filter(!Peptide %in% cpep$Peptide)

fit <- euler(c("C"=nrow(uc), "D"=nrow(ud), "C&D"=nrow(intpep)))

F2A <- plot(fit, quantities = TRUE, 
            fills = list(fill = c("#FFCB39", "#1756E2"), alpha = 0.5),
            labels = list(col = "white", font = 4))

### b) Assay %CV
medD <- comp %>% filter(Group=="D") %>% summarise("medCV"=median(SRM_CV))
medC <- comp %>% filter(Group=="C") %>% summarise("medCV"=median(SRM_CV))
colpal <- c("D" = "#13379B", "C"= "#E8920A")

F2B <- ggplot(comp, aes(x=SRM_CV, fill=Group, color=Group)) + 
  geom_histogram(binwidth=1, alpha=0.2, position = "identity", 
                 mapping = aes(y=..density..)) +
  geom_density(alpha=0.4, linewidth=0.5, mapping = aes(y=..density..)) + 
  scale_fill_manual(values=colpal) +
  geom_vline(data=medD, aes(xintercept=medCV), color="#13379B", 
             linetype="dashed", linewidth=0.5) +
  geom_text(aes(x=(medD$medCV+1), label=round(medD$medCV, digits=2), y=0.195), 
            colour="#13379B", hjust = 0,
            text=element_text(size=8)) +
  geom_vline(data=medC, aes(xintercept=medCV), color="#E8920A",
             linetype="dashed", linewidth=0.5) +
  scale_colour_manual(values = colpal) +
  geom_text(aes(x=(medC$medCV+1), label=round(medC$medCV, digits=2), y=0.215), 
            colour="#E8920A", hjust = 0,
            text=element_text(size=8)) +
  labs(x="% coefficient of variation", y="density")+
  theme_bw()

F2 <- ggarrange(F2A, F2B, ncol=2, widths=c(0.5, 1.5))




# Figure 3. Pain SRM Assay

#############################################################################
pain <- read.csv("pain_CV.csv", header = TRUE)

pain <- as.data.frame(pain)

avpC = mean(pain$CV)
medpC = median(pain$CV)

F3 <- ggplot(pain, aes(x=pain$CV)) + 
  geom_histogram(binwidth=1, alpha=0.3, fill="#1756E2", color = "#1756E2",
                 mapping = aes(y=..density..)) +
  geom_density(alpha = 0.3, fill = "#1756E2", color= "#1756E2", 
               linewidth = 0.5) +
  geom_text(aes(x=(medpC+0.5), label=round(medpC, digits=2), y=0.195), 
            colour="grey20", hjust = 0,
            text=element_text(size=8)) +
  labs(x="% coefficient of variation", y="density") +
  geom_vline(data=pain, aes(xintercept=medpC),linetype="dashed", linewidth=0.5, 
             colour="black") +
  theme_bw()





# Supplemental Figure 2. DIA CV distribution

#############################################################################
dia <- read.csv("DtS_DIA_report.csv", header = TRUE)

pep <- dia$Peptide
diaTIC <- dia %>% dplyr:: select(contains("Total."))
diaTIC <- cbind("Peptide"=pep, diaTIC[,2:ncol(diaTIC)])
diaTIC[diaTIC == "#N/A"] <- NA

# calc medians TIC for normalization
medTIC <- diaTIC %>% dplyr::select(contains("Total.Ion.Current.Area")) %>%  
  rowwise() %>%
  dplyr::mutate(med = median(c_across(where(is.numeric)), na.rm = TRUE))

#TIC normalize pepide area  function
tic.normArea <- function(df){
  df <- diaTIC %>% dplyr::select(contains(df))
  df[is.na(df)] <- 0
  
  #Assigning the highest value to be the representative value
  df$Area <- ifelse(as.numeric(df[,1]) > as.numeric(df[,3]),
                    as.numeric(df[,1]), as.numeric(df[,3]))
  
  #Keeping the TIC for the highest value
  df$TIC <- ifelse(as.numeric(df[,1]) > as.numeric(df[,3]), 
                   df[,2], df[,4])
  
  df.Area <- data.frame("norm_Area"=with(df, Area/TIC*medTIC$med))
  return(df.Area)
}

#Normalize for each GPF-Library replicate
A.Area <- tic.normArea("LibA")
colnames(A.Area) <- "A.Area"

B.Area <- tic.normArea("LibB")
colnames(B.Area) <- "B.Area"

C.Area <- tic.normArea("LibC")
colnames(C.Area) <- "C.Area"

#New dataframe for ease of plotting
## Peptide meta: peptide, protein, precursor, charge, length, mean Area, 
## sd, %CV, ntrans count, rank within protein, 

peptide_df <- dia[,1:12]
peptide_df$Fragment.Ion <- NULL
peptide_df$Fragment.Ion.Type <- NULL
peptide_df$Library.Rank <- NULL

## TIC normalized areas
peptide_df <- data.frame(cbind(peptide_df, A.Area, B.Area, C.Area ))
peptide_df <- peptide_df %>% distinct()

## Calc mean, sd and %CV of TIC normalized peptide areas
peptide_df$mean.Area <- apply(peptide_df[,10:12], 1, mean)
peptide_df$sd.Area <- apply(peptide_df[,10:12], 1, sd)
peptide_df$CV.Area <- (peptide_df$sd.Area/peptide_df$mean.Area)*100

medCV.Area <- median(peptide_df$CV.Area, na.rm=TRUE)
SF2 <- ggplot(peptide_df, aes(x=CV.Area)) +
  geom_histogram(colour='#415168', fill='#8794a0', binwidth=2) +
  geom_vline(xintercept=20, color='orange', 
             linetype='longdash', linewidth=0.6) +
  geom_vline(xintercept=medCV.Area, color='#415168', 
             linetype='dashed',linewidth=0.6) +
  geom_text(aes(x=medCV.Area+2, label=round(medCV.Area, digits=2), y=810), 
            colour="gray30", hjust = 0,
            text=element_text(size=8)) +
  geom_text(aes(x=22, label="20% CV", y=910), 
            colour="gray30", hjust = 0,
            text=element_text(size=8)) +
  scale_x_continuous(limits=c(0,125)) +
  theme_bw() +
  labs(title="Peptide reproducibility from DIA", 
       x="% coefficient of variation", y="count")



# Supplemental Figure 3. GPF-DIA number of transitions

#############################################################################

### a) Distribution of transitions/peptide
shist <- ggplot(peptide_df, aes(x=factor(Transition.Count))) +
  geom_bar(color='black', fill='#8794a0') +
  theme_bw() +
  theme(axis.title.x=element_blank()) +
  labs(title="distribution of peptide %CV")

### b) %CV by transitions/peptide
sbox <- ggplot(peptide_df, aes(x=as.factor(Transition.Count), y=CV.Area))  +
  stat_boxplot(geom = "errorbar", width = 0.2, color="black") +
  geom_boxplot(color='black', fill='#8794a0', outlier.shape=NA) +
  scale_y_continuous(limits=c(0,50))+
  theme_bw() +
  geom_hline(yintercept=20, linetype='dashed', color='#E8920A', linewidth=1) +
  labs(x="interference-free transitions", y="% coefficient of variation")

SF3 <- ggarrange(shist, sbox, nrow = 2, align="hv", heights = c(1,2))




# Supplemental Figure 4. GPF-DIA Peptide length

#############################################################################
# a) Distribution of peptide length
lhist <- ggplot(peptide_df, aes(x=factor(Peptide.Sequence.Length))) +
  geom_bar(color='black', fill='#8794a0') +
  theme_bw() +
  theme(axis.title.x=element_blank())

# b) %CV by peptide length
lbox <- ggplot(peptide_df, aes(x=as.factor(Peptide.Sequence.Length), 
                               y=CV.Area))  +
  stat_boxplot(geom = "errorbar", width = 0.2, color="black") +
  geom_boxplot(color='black', fill='#8794a0', outlier.shape=NA) +
  geom_hline(yintercept=20, linetype='dashed', color='#E8920A', linewidth=1) +
  labs(x="Peptide sequence length", y="% coefficient of variation")+
  scale_y_continuous(limits=c(0,100))+
  theme_bw()

SF4 <- ggarrange(lhist, lbox, nrow = 2, align="hv", heights = c(1,2))




#Session info:

#############################################################################
#R version 4.3.1 (2023-06-16 ucrt)
#Running under: Windows 11 x64 (build 22631)

#attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] ggpubr_0.6.0  ggplot2_3.4.2 eulerr_7.0.2  dplyr_1.1.2  
