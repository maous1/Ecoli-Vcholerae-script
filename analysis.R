library(syngulon2)
library(reutils)
library(ape)
library(seqinr)
library(Biostrings)
library(dplyr)


Sys.setenv(PATH = paste(Sys.getenv("PAH"),"/home/jerome/anaconda3/envs/bacterio/bin/",sep=':'))
system2(command='conda',args = 'list')


phylum <- extract.phylum()
dir.create('01-selected.species/')
write.csv(phylum,'01-selected.species/phylum.table.csv')

bacteria.table <- extract.bacteria.table(phylum = phylum)
write.csv(bacteria.table,'01-selected.species/bacteria.table.csv',row.names = F)

##########  download of accessions 


species <- c('Escherichia_coli','Vibrio_cholerae')
dir.create('02-accessions')
download.accession(species = species,bacteria.table = bacteria.table,outDir = '02-accessions/')

########" download annotation

dir.create('03-annotation')
species <- c('Escherichia_coli','Vibrio_cholerae')
download.annotation(species = species,maxOrganism = 100,indextostart =1,accessionDir = '02-accessions/',outDir = '03-annotation/' )

amelioration.annotation(species = species,collicin = collicin,annotationDir = "03-annotation/")

######## download genomes 

dir.create('04-genomes')
species <- c('Escherichia_coli','Vibrio_cholerae')
download.genome(species=species, maxOrganism =100,indextostart = 1,accessionDir = '02-accessions/',outDir = '04-genomes/' )


########## analyze annotation

bacteria.table <- read.csv('01-selected.species/bacteria.table.csv')
collicin <- read.csv('00-collicin/collicin.csv')
result.annotation <- analyse.annotation(bacteria.table = bacteria.table,collicin =collicin,annotationDir = '03-annotation/'  )

summary.species <- result.annotation$`summmary at species level`
presence <- result.annotation$`presence at species level`
presence.ecoli <- presence%>%filter(.,species =="Escherichia_coli")
presence.cholerae <- presence%>%filter(.,species =="Vibrio_cholerae")


summary.species <- summary.species%>%mutate(species_n=paste0(species,' n = ',n)) %>% select(-c(species,n)) %>% arrange(SubGroup)


#### heatmap au niveau species des gammaproteobacteria

summary.species.Gammaproteobacteria <- summary.species%>%filter(SubGroup=='Gammaproteobacteria') %>% select(-c(SubGroup,n))
matrix.to.plot <- as.matrix(summary.species.Gammaproteobacteria[,-grep('species',colnames(summary.species.Gammaproteobacteria))])

pdf('99-results/heatmap.species.Gammaproteobacteria.pdf',width=18,height = 16)
heatmap.2(as.matrix(summary.species.Gammaproteobacteria[,-grep('species',colnames(summary.species.Gammaproteobacteria))]),dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =summary.species.Gammaproteobacteria$species,margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5 )
dev.off()

#### heatmap au niveau species des ecoli
matrix.to.plot <- as.matrix(presence.ecoli[,-grep('species',colnames(presence.ecoli))])

pdf('99-results/heatmap.e.coli.pdf',width=18,height = 30)
heatmap.2(as.matrix(presence.ecoli[,-grep('species',colnames(presence.ecoli))]),dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =row.names(presence.ecoli),margins = c(7,20) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5 )
dev.off()

#### heatmap au niveau species des cholerae

matrix.to.plot <- as.matrix(presence.cholerae[,-grep('species',colnames(presence.cholerae))])

pdf('99-results/heatmap.cholerae.pdf',width=18,height = 30)
heatmap.2(as.matrix(presence.cholerae[,-grep('species',colnames(presence.cholerae))]),dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =row.names(presence.cholerae),margins = c(7,30) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5 )
dev.off()

#################### Extract genes annotation based

species <- c('Escherichia_coli','Vibrio_cholerae')
dir.create('05-genes/')
collicin <- read.csv('collicin.csv')
extract.all.genes.annotationbased(species =species,collicin = collicin,annotationDir = '03-annotation/',genomeDir = '04-genomes/',outDir = '05-genes/' )


###################compute otus
species <- c('Escherichia_coli','Vibrio_cholerae')
collicin <- read.csv('collicin.csv')

otus(species = species,collicin = collicin,geneDir = '05-genes/',outDir = '08-otus/')


################### Plot phylogenetic of otus
bacteria.table <- read.csv('01-selected.species/bacteria.table.csv')
phylo.from.otus(otusDir = "08-otus/",bacteria.table = bacteria.table)

  
  
  #######################   screen blast
  
dir.create('09-screenblast/')
genelist <- list.files('08-otus/',full.names = T,recursive = F)
genelist <- genelist[grep('.fasta',genelist)]
ngenes <- length(genelist)
genomelist <- list.files('04-genomes/',full.names = T,recursive = T)
ngenomes <- length(genomelist)
result.blast <- matrix(nrow=ngenomes,ncol=ngenes)
rownames(result.blast) <- gsub(genomelist,pattern = '04-genomes/',replacement = '')
colnames(result.blast) <- gsub(basename(genelist),pattern = 'fasta',replacement = '')

sequences <- readDNAStringSet(genelist)
writeXStringSet(sequences,'allcollicingene.fasta')  
for(i in 1:ngenomes)
{
  result.blast[i,] <- screenBlast(reference =  'allcollicingene.fasta',querry = genomelist[i],min.pc.ident = 95 ,min.pc.length=50,dir.out = '09-screenblast/')
  print(i)
}
write.csv(result.blast,'99-results/result.blast.csv')



############### heatmap of blast results
bacteria.table <- read.csv('01-selected.species/bacteria.table.csv')
blastresult <- read.csv('99-results/result.blast.csv',row.names = 1)

blastresult <- analyse.blast.table(bacteria.table,blastresult)

summary.species <- blastresult$`summmary at species level`
presence.species <- blastresult$`presence at species level`
presence.blast.ecoli <- presence.species%>%filter(.,species =="Escherichia_coli")
presence.blast.cholerae <- presence.species%>%filter(.,species =="Vibrio_cholerae")

summary.species <- summary.species%>%mutate(species_n=paste0(species,' n = ',n)) %>% select(-c(species,n)) %>% arrange(SubGroup)



######### heatmap au niveau species des ecoli et v cholera
my_palette <- colorRampPalette(c("white", "gray", "black"))(n = 299)

presence.species.ecoli <- presence.species%>%filter(species=='Escherichia_coli')
matrix.to.plot <- as.matrix(presence.species.ecoli[,-grep('species',colnames(presence.species.ecoli))])
pdf('99-results/heatmap.blast.ecoli.pdf',width=18,height = 30)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =rownames(presence.species.ecoli),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5,breaks=seq(0,1,length.out = 300) )
dev.off()


presence.species.cholerae <- presence.species%>%filter(species=='Vibrio_cholerae')
matrix.to.plot <- as.matrix(presence.species.cholerae[,-grep('species',colnames(presence.species.cholerae))])
pdf('99-results/heatmap.blast.cholerae.pdf',width=18,height = 30)
heatmap.2(matrix.to.plot,dendrogram='none', Rowv=FALSE, Colv=FALSE,
          trace='none',col=my_palette,sepcolor=gray(0.3),colsep=seq(0,50),rowsep=seq(0,150),sepwidth=c(0.015,0.015),
          labRow =rownames(presence.species.cholerae),margins = c(7,35) ,cellnote =matrix.to.plot,keysize = 0.2,key=F,cexRow = 1.5,cexCol = 1.5 )
dev.off()


