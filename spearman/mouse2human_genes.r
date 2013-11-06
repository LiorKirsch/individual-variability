pp <- scan("allen_mouse_entrez.txt",what=character(0))


library(biomaRt)
# define biomart objects
human=useMart("ensembl", dataset="hsapiens_gene_ensembl")
mouse=useMart("ensembl", dataset="mmusculus_gene_ensembl")

#entrezgene

out=getLDS(attributes=c("mgi_symbol","wikigene_name","entrezgene"), filters = "entrezgene", values = pp,
           mart=mouse,
           attributesL=c("hgnc_symbol","entrezgene"), martL=human)
# create output file
write.table(out,"Mouse2Human.txt")

# this function uses ensembl gene ids,
# but you can choose any other ID or symbol.
# to see the full list of names and many other attributes run - 
# att<-listAttributes(mouse)


