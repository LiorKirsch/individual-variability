pp <- scan("mouseGeneNames.txt",what=character(0))

mydata = read.csv("human_genes.csv")


humanEntrez = mydata[3]

library(biomaRt)
# define biomart objects
human=useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse=useMart("ensembl", dataset="mmusculus_gene_ensembl")

#entrezgene

out=getLDS(attributes=c("hgnc_symbol","entrezgene"), filters = "entrezgene", values = humanEntrez,
           mart=human,
           attributesL=c("mgi_symbol","entrezgene"), martL=mouse)
# create output file
write.table(out,"Human2Mouse.txt")

# this function uses ensembl gene ids,
# but you can choose any other ID or symbol.
# to see the full list of names and many other attributes run - 
# att<-listAttributes(mouse)
