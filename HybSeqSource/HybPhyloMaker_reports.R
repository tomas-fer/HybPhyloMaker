#HybPhyloMaker summary tables as formatted XLSX
#Called from HybPhyloMaker_reports.sh
#Tomas Fer, 2018

#Load library
library(openxlsx)
#Read arguments from command line
args <- commandArgs()
nrsamples <- as.numeric(args[5])
nrgenes <- as.numeric(args[6])
nrselected <- as.numeric(args[7])
missingpercent <- as.numeric(args[8])

# nrsamples=96
# nrgenes=977
# nrselected=125
# missingpercent=50

print(paste("Nr. samples is", nrsamples))
print(paste("Nr. genes is", nrgenes))
print(paste("Nr. selected genes is", nrselected))
print(paste("Missing percent level is", missingpercent))

#Read tables
print("Reading tables")
a <- read.csv("1-Reads_summary.csv", check.names=FALSE, stringsAsFactors=FALSE)
b <- read.csv("2-Mapping_summary.csv", check.names=FALSE, stringsAsFactors=FALSE)
c <- read.csv("3-Missing_data.csv", check.names=FALSE, stringsAsFactors=FALSE)
d <- read.csv(paste("4-Missing_data_", missingpercent, ".csv", sep=""), check.names=FALSE, stringsAsFactors=FALSE)
e <- read.csv("5-All_genes.csv", check.names=FALSE, stringsAsFactors=FALSE)
f <- read.csv("6-Selected_genes.csv", check.names=FALSE, stringsAsFactors=FALSE)

sheet4 <- paste("Missing_data_", missingpercent, sep="")
#make list of tables
print("Making list")
list <- list("Reads_summary" = a, "Mapping_summary" = b, "Missing_data" = c, sheet4 = d, "All_genes" = e, "Selected_genes" = f)
#define styles
print("Defining styles")
hs <- createStyle(fontColour = "black", halign = "left", valign = "center", textDecoration = "Bold", wrapText=TRUE)
tablestyle1 <- createStyle(halign = "right", valign = "center", numFmt = "0.00")
tablestyle2 <- createStyle(halign = "right", valign = "center", numFmt = "0.000")
tablestyle3 <- createStyle(halign = "right", valign = "center", numFmt = "0.0000")
tablestyle4 <- createStyle(halign = "right", valign = "center", numFmt = "# ### ##0")
italics <- createStyle(fontColour = "black", halign = "left", valign = "center", textDecoration = "Italic")
#save xlsx to wb
print("Saving workbook")
wb <- write.xlsx(list, file="summary.xlsx", sheetName="HybPhyloMaker SUMMARY", colWidths="auto", wrapText=TRUE, keepNA = TRUE, firstActiveRow=2, firstActiveCol=2, headerStyle = hs)
names(wb)[[4]] <- sheet4 #change name of the sheet 4
#delete empty cells from sheet 4
print("Deleting empty cells")
deleteData(wb, 4, rows = (nrsamples+2):(nrsamples+3), cols = (nrgenes+4):(nrgenes+5), gridExpand = TRUE)
#add styles
print("Adding styles")
print("...sheet1")
addStyle(wb, 1, style = tablestyle4, rows = 2:(nrsamples+4), cols = c(4,5,6,8,9,10,11,13), gridExpand = TRUE)
addStyle(wb, 1, style = tablestyle1, rows = 2:(nrsamples+4), cols = c(7,12,14), gridExpand = TRUE)
addStyle(wb, 1, style = tablestyle4, rows = 2, cols = c(4,5,6,8,9,10,11,13), gridExpand = TRUE)
addStyle(wb, 1, style = tablestyle1, rows = 2, cols = c(7,12,14), gridExpand = TRUE)
addStyle(wb, 1, style = italics, rows = 2:(nrsamples+1), cols = 2:3, gridExpand = TRUE)
addStyle(wb, 1, style = italics, rows = 2, cols = 2:3, gridExpand = TRUE)
print("...sheet2")
addStyle(wb, 2, style = tablestyle4, rows = 2:(nrsamples+4), cols = 4:8, gridExpand = TRUE)
addStyle(wb, 2, style = tablestyle1, rows = 2:(nrsamples+4), cols = 9, gridExpand = TRUE)
addStyle(wb, 2, style = italics, rows = 2:(nrsamples+1), cols = 2:3, gridExpand = TRUE)
addStyle(wb, 2, style = italics, rows = 2, cols = 2:3, gridExpand = TRUE)
print("...sheet3")
addStyle(wb, 3, style = tablestyle1, rows = nrsamples+2, cols = 2:(nrgenes+3), gridExpand = TRUE)
addStyle(wb, 3, style = italics, rows = 2:(nrsamples+1), cols = 2:3, gridExpand = TRUE)
addStyle(wb, 3, style = italics, rows = 2, cols = 2:3, gridExpand = TRUE)
print("...sheet4")
addStyle(wb, 4, style = tablestyle1, rows = nrsamples+2, cols = 2:(nrgenes+3), gridExpand = TRUE)
addStyle(wb, 4, style = tablestyle1, rows = nrsamples+3, cols = 2:(nrgenes+3), gridExpand = TRUE)
addStyle(wb, 4, style = tablestyle1, rows = 2:(nrsamples+1), cols = (nrgenes+4), gridExpand = TRUE)
addStyle(wb, 4, style = tablestyle4, rows = 2:(nrsamples+1), cols = (nrgenes+5), gridExpand = TRUE)
addStyle(wb, 4, style = italics, rows = 2:(nrsamples+1), cols = 2:3, gridExpand = TRUE)
addStyle(wb, 4, style = italics, rows = 2, cols = 2:3, gridExpand = TRUE)
print("...sheet5")
addStyle(wb, 5, style = tablestyle4, rows = 2:(nrgenes+4), cols = c(2:7,9,13:31), gridExpand = TRUE)
addStyle(wb, 5, style = tablestyle2, rows = 2:(nrgenes+4), cols = c(8,10:12), gridExpand = TRUE)
addStyle(wb, 5, style = tablestyle1, rows = 2:(nrgenes+4), cols = c(32:33), gridExpand = TRUE)
print("...sheet6")
addStyle(wb, 6, style = tablestyle1, rows = 2:(nrselected+4), cols = c(32:33,35,37:41), gridExpand = TRUE)
addStyle(wb, 6, style = tablestyle4, rows = 2:(nrselected+4), cols = c(2:7,9,13:31), gridExpand = TRUE)
addStyle(wb, 6, style = tablestyle2, rows = 2:(nrselected+4), cols = c(8,10:12), gridExpand = TRUE)
addStyle(wb, 6, style = tablestyle3, rows = 2:(nrselected+4), cols = c(36), gridExpand = TRUE)
#set height of first row
setRowHeights(wb, 1, rows = 1, heights = 30)
#modify column widths
print("Setting column width")
setColWidths(wb, sheet = 1, cols = c(4,5,6,7,8,9,10,11,12,13,14), widths = c(10,10,11.8,6.3,10,12,12,15.3,14,16,10))
setColWidths(wb, sheet = 2, cols = c(4,5,6,7,8,9), widths = c(10,10,14,14,11,13))
#setColWidths(wb, sheet = 3, cols = c(), widths = c())
#setColWidths(wb, sheet = 4, cols = c(), widths = c())
setColWidths(wb, sheet = 5, cols = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,28,30,31,32,33), widths = c(5.5,10,11,13.5,7.5,10.5,12,15,11,7.5,7.5,8.2,8.2,8.2,8.2,8.2,8.2,8.2,7.5,6))
setColWidths(wb, sheet = 6, cols = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,28,30,31,32,33,36,38,39,40,41), widths = c(5.5,10,11,13.5,7.5,10.5,12,15,11,7.5,7.5,8.2,8.2,8.2,8.2,8.2,8.2,8.2,7.5,6,7,8,6,6,7))
#freeze header rows/columns
print("Freezing header")
freezePane(wb, 1, firstActiveRow = 2, firstActiveCol = 4)
freezePane(wb, 2, firstActiveRow = 2, firstActiveCol = 4)
freezePane(wb, 3, firstActiveRow = 2, firstActiveCol = 4)
freezePane(wb, 4, firstActiveRow = 2, firstActiveCol = 4)
freezePane(wb, 5, firstActiveRow = 2, firstActiveCol = 2)
freezePane(wb, 6, firstActiveRow = 2, firstActiveCol = 2)

#Add sum/average formulae
print("Adding formulas")
#sheet1
clmn=4
for (cell in c("D","E","F","G","H","I", "J","K","L","M","N")){
	v <- c(paste("SUM(",cell,"2:",cell,nrsamples+1,")",sep=""), paste("AVERAGE(",cell,"2:",cell,nrsamples+1,")",sep=""))
	writeFormula(wb, sheet = 1, x = v, startCol = clmn, startRow = nrsamples+3)
	clmn <- clmn + 1
}
writeData(wb, 1, x = c("TOTAL","AVERAGE"), startCol = 1, startRow = nrsamples+3)
deleteData(wb, 1, rows = nrsamples+3, cols = 7)
deleteData(wb, 1, rows = nrsamples+3, cols = 12)
deleteData(wb, 1, rows = nrsamples+3, cols = 14)
#sheet2
clmn=4
for (cell in c("D","E","F","G","H","I")){
	v <- c(paste("SUM(",cell,"2:",cell,nrsamples+1,")",sep=""), paste("AVERAGE(",cell,"2:",cell,nrsamples+1,")",sep=""))
	writeFormula(wb, sheet = 2, x = v, startCol = clmn, startRow = nrsamples+3)
	clmn <- clmn + 1
}
writeData(wb, 2, x = c("TOTAL","AVERAGE"), startCol = 1, startRow = nrsamples+3)
deleteData(wb, 2, rows = nrsamples+3, cols = 9)

#sheet5
clmn=2
for (cell in c("B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG")){
	v <- c(paste("SUM(",cell,"2:",cell,nrgenes+1,")",sep=""), paste("AVERAGE(",cell,"2:",cell,nrgenes+1,")",sep=""))
	writeFormula(wb, sheet = 5, x = v, startCol = clmn, startRow = nrgenes+3)
	clmn <- clmn + 1
}
v <- c(paste('AVERAGEIF(AF2:AF',nrgenes+1,';">=0")',sep=""))
#writeFormula(wb, sheet = 5, x = paste('AVERAGEIF(AF2:AF',nrgenes+1,';">=0")',sep=""), startCol = 32, startRow = nrgenes+4) #does not work!!!
writeData(wb, 5, x = c("TOTAL","AVERAGE"), startCol = 1, startRow = nrgenes+3)
for (delete in c(2,6,8,10,11,12,32,33)){
	deleteData(wb, 5, rows = nrgenes+3, cols = delete)
}

#sheet6
clmn=2
for (cell in c("B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK","AL","AM","AN","AO")){
	v <- c(paste("SUM(",cell,"2:",cell,nrselected+1,")",sep=""), paste("AVERAGE(",cell,"2:",cell,nrselected+1,")",sep=""))
	writeFormula(wb, sheet = 6, x = v, startCol = clmn, startRow = nrselected+3)
	clmn <- clmn + 1
}
writeData(wb, 6, x = c("TOTAL","AVERAGE"), startCol = 1, startRow = nrselected+3)
for (delete in c(2,6,8,10,11,12,32,33,34,35,36,37,38,39,40,41)){
	deleteData(wb, 6, rows = nrselected+3, cols = delete)
}
deleteData(wb, 6, rows = nrselected+4, cols = 34)

#write XLSX to file
print("Writing final file")
saveWorkbook(wb, "summary.xlsx", overwrite = TRUE)
