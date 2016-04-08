#VISUALIZE CORRELATIONS AMONG TREE AND ALIGNMENT PROPERTIES
#Taken from M. Borowiec GitHub (https://github.com/marekborowiec/good_genes)
#Modified for HybPipe
#--------------------------------------------------------------------------

# disable scientific notation
options(scipen=999)

# install required libraries
# uncomment if you don't have these installed
#install.packages("data.table")

# load needed libraries
library("data.table")

#Read table with combined characteristics
gene_stats <- read.table("combined.txt", header=TRUE)

#names(gene_stats)

#select variables to plot
to_plot <-subset(gene_stats, select=c(Aln_length,Missing_perc,Prop_pars_inf,
                                      GC_content,Aln_entropy,trimAl_sct,Bootstrap,Branch_length,
                                      Clocklikeness,P_distance,Satur_slope,Satur_R_sq,
                                      LBscore_SD))

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.99/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * abs(r)) 
  text(.8, .8, Signif, cex=cex, col=2)
}


png(file="genes_corrs.png",width=20*300,height=16*300,res=300)
par(pch=20)
pairs(to_plot, lower.panel=panel.smooth, upper.panel=panel.cor,cex.labels=2)
dev.off()

pdf(file="genes_corrs.pdf",width=24,height=20)
par(pch=20)
pairs(to_plot, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()