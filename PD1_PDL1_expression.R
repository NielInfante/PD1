library(FirebrowseR)
cohorts = Metadata.Cohorts(format = "csv") # Download all available cohorts



immuno.Genes = c("PDCD1", "CD274") #The official gene symbol of PD-1 and PD-L1




slide.size = 100 #number of tcga barcodes for Gene expression extract. FirebrowseR always return error if extract expression data for the whole cohort. Therefore I break it into several slides of small number of samples. 


extract_Exp <- function(tmp.Pats, diff.Exp.Genes, cohort_name){ #Given sample barcode, gene id, and cohort, return the gene expression matrix
  all.Found = F
  page.Counter = 1
  mRNA.Exp = list()
  page.Size = 2000 
  while(all.Found == F){
    mRNA.Exp[[page.Counter]] = Samples.mRNASeq(format = "csv",
                                               gene = diff.Exp.Genes,
                                               cohort = cohort_name,
                                               tcga_participant_barcode =
                                                 tmp.Pats$tcga_participant_barcode,
                                               page_size = page.Size,
                                               page = page.Counter)
    if(nrow(mRNA.Exp[[page.Counter]]) < page.Size)
      all.Found = T
    else
      page.Counter = page.Counter + 1
  }
  return(mRNA.Exp)
}

cohorts = cohorts[-10,] #remove FPPP. - no samples available
cohorts = cohorts[-17,] #remove LAML. - no samples available
nCohorts = nrow(cohorts)
for (i in 1:nCohorts){  #loop all the cohorts
  cancer.Type=cohorts[[1]][i]
  #find the tcga barcode for the samples in the specific cancer cohorts
  all.Received = F
  page.Counter = 1
  page.size = 150
  cancer.Pats = list()
  while(all.Received == F){
    cancer.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                                   cohort = cancer.Type,
                                                   page_size = page.size,
                                                   page = page.Counter)
    if(page.Counter > 1)
      colnames(cancer.Pats[[page.Counter]]) = colnames(cancer.Pats[[page.Counter-1]])
    
    if(nrow(cancer.Pats[[page.Counter]]) < page.size){
      all.Received = T
    } else{
      page.Counter = page.Counter + 1
    }
  }
  cancer.Pats = do.call(rbind, cancer.Pats)
  
  nslide = floor(nrow(cancer.Pats) / slide.size) #find the number of slides based on sample size for each slide
  if(exists(deparse(substitute(immuno.Genes.Exp)))){ #remove the immuno.Genes.Exp from last round
    rm(immuno.Genes.Exp)
  }
  if(nslide > 0){
    for (j in 1:nslide){
      tmp.Pats = cancer.Pats[((j-1)*slide.size+1):(j*slide.size),]
      mRNA.Exp = extract_Exp(tmp.Pats, immuno.Genes, cancer.Type)
      if (j==1){
        immuno.Genes.Exp = mRNA.Exp
      }else{
        immuno.Genes.Exp = append(immuno.Genes.Exp,mRNA.Exp)
      }
    }
  }
  if(nrow(cancer.Pats)>nslide*slide.size){
    tmp.Pats=cancer.Pats[(nslide*slide.size+1):nrow(cancer.Pats),]
    mRNA.Exp = extract_Exp(tmp.Pats, immuno.Genes, cancer.Type)
    if(exists(deparse(substitute(immuno.Genes.Exp)))){
      immuno.Genes.Exp = append(immuno.Genes.Exp,mRNA.Exp)
    }else{
      immuno.Genes.Exp = mRNA.Exp
    }
  }
  immuno.Genes.Exp = do.call(rbind,immuno.Genes.Exp)
  assign(cancer.Type,immuno.Genes.Exp)
}
#Calculate the median exp value for PD-1 and PD-L1 for each cohort
PD1_exp = matrix(0,nrow=nCohorts,ncol=1)
PDL1_exp = matrix(0,nrow=nCohorts,ncol=1)
for(i in 1:nCohorts){
  mRNA.exp = get(cohorts[[1]][i])
  PD1_exp[i] = median(as.numeric(mRNA.exp$expression_log2[which(mRNA.exp$gene == "PDCD1" & mRNA.exp$sample_type=="TP")]),na.rm=T)
  PDL1_exp[i] = median(as.numeric(mRNA.exp$expression_log2[which(mRNA.exp$gene == "CD274" & mRNA.exp$sample_type=="TP")]),na.rm=T)
}


plot(PD1_exp,PDL1_exp,main = "Median expression of PD-1 and PD-L1 from 36 TCGA cohorts")

lines(PDL1_exp,PDL1_exp,type = 'l',col="red")
text(PD1_exp,PDL1_exp,labels = cohorts[[1]],pos=1, offset=0.2,cex=0.5)


library(tidyverse)
library(plotly)


mean_expression <- tibble(Cohort = cohorts[[1]], PD1 = PD1_exp[,1], PDL1 = PDL1_exp[,1])


p <- mean_expression %>% ggplot(aes(PD1, PDL1, tooltip = Cohort)) + geom_point() + 
  theme_minimal()

ggplotly(p)



mean_expression <- mean_expression %>% mutate(Combined = PD1 + PDL1) %>% arrange(desc(Combined))

mean_expression %>% print(n=100)
