#Repository for PL PATH 875 Phylogenetics Practicum (Spring 2026)



# 03 Data and Software Setup

- Clone repository
% git clone https://github.com/solislemuslab/phylo-practicum.

- Download data and move into data subfolder and unzip files

--> Installation

(1) RAxML Next Generation

- Download software, place into designated local folder (folder labeled software)

- Make executable runnable (in Terminal)

% cd ~/software/raxml-ng_v1.2.2_macos
% chmod + x raxml-ng
% raxml -ng -v (to verify installation)

-Add to PATH to run raxml-ng from any where 

% nano ~/.zshrc

export PATH="$HOME/software/raxml-ng_v1.2.2_macos:$PATH" (save, exit, and reload to confirm changes)

(2) SuperTriplets 

-Check if Java is installed

% java -h 

-Download SuperTriplets java file and put into software folder inside a subfolder named supertriplets


(3) ASTRAL

-Install with conda

% conda config --add channels bioconda
% conda install aster

-To check installation location and existence

% conda list
% wastral

(4) Julia and packages (in Terminal)

- Install Julia 

% curl -fsSL https://install.julialang.org | sh

-Run this command after Terminal restart

% juliaup add release

-Add other packages

julia then ] to go into package mode

add PhyloNetworks
add SNaQ
add PhyloPlots
add CSV
add DataFrames

Hit backspace to leave package mode and leave Julia by typing exit()


(5) HyDe

-Move to software folder and download 

% git clone https://github.com/pblischak/HyDe.git

-Install HyDe

% cd HyDe
% python3 -m pip install -r requirements.txt
% python3 -m pip install .

-Check installation

% run_hyde.py



# 04 Gene tree estimation

--> Running RAxML on 10 genes 

--Move to data/Wheat_Relative_History_Data_Glemin_et_al folder and create OneCopyGene-trimmed folder for the first 10 genes

% mkdir OneCopyGenes-trimmed
% cd OneCopyGenes
% ls | head -n10 | xargs -I {} cp "{}" ../OneCopyGenes-trimmed

--Move to the code folder to run 04-raxml.sh script

% ./04-raxml.sh

~Analyze inferred trees in R

library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

getwd()
setwd(/Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/results/RAxML) ##set directory to results/RAxML folder

tree_files <-list.files(pattern="\\.raxml.bestTree$") #List all .bestTree files. $ ensures the end of the name

gene_trees<- list() # list with all the trees
class(gene_trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  gene_trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

#Find the appropriate root taxon
root_taxa <- c("H_vulgare_HVens23", "Er_bonaepartis_TB1", "S_vavilovii_Tr279", "Ta_caputMedusae_TB2")

#Extract the first matching species for each tree
gene_tree_outgroup<- rep(NA,length(gene_trees))
for(i in 1:length(gene_trees)){
  gene_tree <- gene_trees[[i]]
  found_taxa <- root_taxa[root_taxa %in% gene_tree$tip.label]
  
  #Return the first one if it exists
  if (length(found_taxa) > 0) {
    gene_tree_outgroup[i]<-found_taxa[1]
  }
}

#check if any gene trees don't have an outgroup taxon
any(is.na(gene_tree_outgroup)) ## should return FALSE

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(gene_trees)){
  gene_trees[[i]]<- root(gene_trees[[i]],
                         outgroup = gene_tree_outgroup[i],
                         resolve.root=TRUE)
  gene_trees[[i]]<-chronos(gene_trees[[i]]) ## make ultrametric for nicer densitree
}

~Visualize the trees 

plot(gene_trees[[2]]) # plots the second gene tree

#see which genes have which taxa

all_labels <- unlist(lapply(gene_trees, function(x) x$tip.label)) ##count how often each individual appears
df <- as.data.frame(table(all_labels) / length(gene_trees)) #organize as a nice data frame for plotting

ggplot(df, aes(x = all_labels, y = Freq)) +
  geom_col() +
  labs(x = "Individual", y = "Proportion")+
  theme(axis.text.x = element_text(angle = 90))

--SuperTree

st<-superTree(gene_trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

--Densitree 

densiTree(gene_trees,consensus=st,scaleX=T,type='cladogram', alpha=0.05)

common_tips <- Reduce(
  intersect,
  lapply(gene_trees, function(tr) tr$tip.label)
)

length(common_tips) ## 6

trees_pruned <- lapply(
  gene_trees[1:10],
  function(tr) drop.tip(tr, setdiff(tr$tip.label, common_tips))
)

densityTree(trees_pruned,use.edge.length=FALSE,type="cladogram",nodes="centered")

--Comparison of ML trees for each of the 20 RAxML runs (visualization of the best tree across each run)

trees = read.tree(file="Ae_bicornis_Tr406_BIS2_Contig10132_simExt_macseNT_noFS_clean.aln.raxml.mlTrees")

rtrees <- lapply(trees, function(tr) {
  root(tr, outgroup = "H_vulgare_HVens23", resolve.root = TRUE)
})

--Plot densiTree before and after rooting

densityTree(trees,type="cladogram",nodes="intermediate")
densityTree(rtrees,type="cladogram",nodes="intermediate")

densityTree(trees,use.edge.length=FALSE,type="cladogram",nodes="centered")
densityTree(rtrees,use.edge.length=FALSE,type="cladogram",nodes="centered")

--> Running RAxML on all on all the genes

-Modify 04-raxml.sh script to run all genes 

DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes"

-Inside code folder run the script

% cd code
% ./04-raxml.sh

~ Visualize the trees 

getwd() #we want to working in the results/RAxML folder 
setwd(/Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/results/RAxML) 

ibrary(ape)
library(phangorn)


tree_files <-list.files(pattern="\\.raxml.bestTree$") #List all .bestTree files. $ ensures the end of the name

gene_trees<- list() # list with all the trees
class(gene_trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  gene_trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

write.tree(gene_trees, file="../04-all_gene_trees.tre")



# 05 Species tree supermatrix (full) 

--> Full concatenation

#Keep in mind incomplete lineages (not as robust).
#10 Mb windows -- snapshot into chromosomes --telling the same story across.If not, consider ILS, lateral gene transfer, recombination. 
#Multiple introgression -- single vs multiple hybridization events (repeated vs single gene flow event)
#specific regions tell certain stories 

--Move to data/Wheat_Relative_History_Data_Glemin_et_al to run RAxML

% cd /Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/data/Wheat_Relative_History_Data_Glemin_et_al 
% raxml-ng --msa triticeae_allindividuals_OneCopyGenes.fasta --model GTR+G4

--Move output files to the results folder

% cd ../../results/RAxML/
% mkdir full-concatenation
% cd ../../data/Wheat_Relative_History_Data_Glemin_et_al
% ls ## to check files are there
% mv *.raxml* ../../results/RAxML/full-concatenation

~ Visualization

--Open R and move to results/RAxML/full-conatenation folder

library(ape)
library(phangorn)
library(phytools)

tre = read.tree(file="triticeae_allindividuals_OneCopyGenes.fasta.raxml.bestTree")
plot(tre)
nodelabels()

--Root at the outgroup clade

rtre = root(tre,node = 51, resolve.root=TRUE)
plot(ladderize(rtre))

--Compare to the one provided by the authors: MLtree_OneCopyGenes. tree

tre2 = read.tree(file="../../../data/Wheat_Relative_History_Data_Glemin_et_al/MLtree_OneCopyGenes.tree")
plot(tre2)
rtre2 = root(tre2,node = 51, resolve.root=TRUE)
plot(ladderize(rtre2))

--Calculate the RF distance (topological dissimilarity) between our full cocatenation tree and theirs

library(phangorn)
RF.dist(tre,tre2) ## 0



# 06 Species tree supermatrix (10 Mb)

#Creating smaller concatentated files (10 Mb windows) to assess tree discordance 

--Unzip 10 mB Concatenated Sequence files in the data/Wheat.. folder

% cd Concatenation10Mb_OneCopyGenes
% ls | wc -l 

--Within code folder, run 

% ./06-raxml-concatenation.sh

~Visualization

--Move to results/RAxML/10Mb-concatenation folder (in R)

library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

getwd()
setwd("/Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/results/RAxML/10Mb-concatenation")

tree_files <-list.files(pattern="\\.raxml.bestTree$") #List all .bestTree files. $ ensures the end of the name

trees<- list() # list with all the trees
class(trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

--Root all trees by "H_vulgare_HVEns23"

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees)){
  trees[[i]]<- root(trees[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees[[i]]<-chronos(trees[[i]]) ## make ultrametric for nicer densitree
}

--Create a consensus parsimony supertree

st<-superTree(trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

-- Plot the same density tree

densiTree(trees,consensus=st,scaleX=T,type='cladogram', alpha=0.1)

--Compare densitree with one made by the authors (Densitree_OneCopyGenes.nex)

--Make some minor modifications (Densitree_OneCopyGenes_modified.nex):

1. Added ; after "begin trees"
2. Had to manually change taxon 0 as 47 

trees2 <- read.nexus("/Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/data/Wheat_Relative_History_Data_Glemin_et_al/Densitree_OneCopyGenes-modified.nex")

#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees2)){
  trees2[[i]]<- root(trees2[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees2[[i]]<-chronos(trees2[[i]]) ## make ultrametric for nicer densitree
}

st2<-superTree(trees2)
st2<-root(st2,"H_vulgare_HVens23",resolve.root = T)
plot(st2)

densiTree(trees2,consensus=st2,scaleX=T,type='cladogram', alpha=0.1)



# 07 Species tree supertree/coalescent

--> Species tree via supertree

--Move to results folder and run this code

% java -jar ~ /Users/thalia/software/supertriplets/SuperTriplets_v1.1.jar 04-all_gene_trees.tre 07-supertree.tre

--Visulaize in R

library(ape)
tre = read.tree(file="07-supertree.tre")
plot(tre)

rtre = root(tre,outgroup="H_vulgare_HVens23", resolve.root=TRUE)
plot(rtre)

--Compare the tree with the more simple parsimony-based supertree that we constructed in 04-gene-trees

library(phangorn)
gene_trees <- read.tree("04-all_gene_trees.tre")
st_parsimony<-superTree(gene_trees)
st_parsimony<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st_parsimony)

--> Species tree via coalescent models

~At the individual level

--Move to results path and run 

% astral4 -i 04-all_gene_trees.tre -o 07-individual-species-tree-astral4.tre

--Visualize in R (in results)

library(ape)
tre = read.tree(file="07-individual-species-tree-astral4.tre")
plot(tre)
rtre = root(tre,outgroup="H_vulgare_HVens23", resolve.root=TRUE)
plot(rtre)

~At the species level

--Create this mapping in R (inside code)

#First we get all the individual names
#Original code gave me issues, so I made some adjustments with help from Gemini.

genes_dir <- "/Users/thalia/Documents/Path875/phylo-practicum/glemin-wheat/data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes"
gene_files <- list.files(genes_dir, pattern = "\\.aln$", full.names = TRUE)

for(f in gene_files){
  #read.dna comes from the 'ape' package
  headers <- rownames(read.dna(f, format = "fasta"))
  all_individuals <- unique(c(all_individuals, headers))
}
all_individuals <- sort(all_individuals)
all_individuals

cleaned_individuals <- sub("_[^_]+$", "", all_individuals)

#map all individuals to species
mapping <- paste(all_individuals,cleaned_individuals)
writeLines(mapping, "../results/07-species_mapping.txt") ## write to file

--Move to results path and run in Terminal: 

% astral4 -i 04-all_gene_trees.tre -a 07-species_mapping.txt -o 07-species-tree-astral4.tre

--Visulaization in R

library(ape)
tre = read.tree(file="07-species-tree-astral4.tre")
plot(tre)
rtre = root(tre,outgroup="H_vulgare", resolve.root=TRUE)
plot(rtre)



# 08 Species tree visualizations

--4 Estimated species trees (in results)


1) full concatenation (RAxML/full-concatenation/triticeae_allindividuals_OneCopyGenes.fasta.raxml.bestTree)
2) consensus of 10Mb window concatenation trees (folder RAxML/10Mb-concatenation)
3) supertree (07-supertree.tre)
4) ASTRAL4 tree (07-individual-species-tree-astral4.tre)

--> Reproducing Figure 1A

--In results (R): 

library(ape)

tree1 = read.tree(file="RAxML/full-concatenation/triticeae_allindividuals_OneCopyGenes.fasta.raxml.bestTree")
tree2 = read.tree(file="07-supertree.tre")

tree1 = root(tree1,outgroup="H_vulgare_HVens23", resolve.root=TRUE)
tree2 = root(tree2,outgroup="H_vulgare_HVens23", resolve.root=TRUE)

#Suppose your trees are called tree1 and tree2
#First, ladderize both trees
tree1 <- ladderize(tree1)
tree2 <- ladderize(tree2[[2]])

#Make sure they have the same tip labels
common_tips <- intersect(tree1$tip.label, tree2$tip.label)

length(common_tips)
length(tree1$tip.label)
length(tree2$tip.label)

#Reorder the second tree to match the tip order of the first
tree2 <- reorder.phylo(tree2, order = "cladewise")
tree2 <- rotateConstr(tree2, tree1$tip.label) # rotate to match tree1 tip order

#Plot side by side
par(mfrow = c(1, 2))
plot(tree1, main = "Full concatenation", cex = 0.8)
plot(tree2, main = "Supertree", cex = 0.8)

--Calculate RF distance

library(phangorn)
RF.dist(tree1, tree2) ## not zero!

--Reproduce the edge colors on the full concatentation tree (tree1):

library(ape)
library(ggtree)
library(dplyr)

species_colors <- c(
  "Ae_umbellulata" = "yellow",
  "Ae_caudata" = "orange",
  "Ae_comosa" = "darkorange3",
  "Ae_uniaristata" = "sienna4",
  "Ae_bicornis" = "purple",
  "Ae_longissima" = "pink",
  "Ae_sharonensis" = "mediumpurple1",
  "Ae_searsii" = "maroon2",
  "Ae_tauschii"    = "red",
  "T_boeoticum" = "darkgreen",
  "T_urartu" = "green",
  "Ae_speltoides" = "blue4",
  "Ae_mutica" = "steelblue1"
)

tip_species <- sapply(tree1$tip.label, function(x) {
  parts <- strsplit(x, "_")[[1]]
  paste(parts[1:2], collapse = "_")  # combine first two parts
})
tip_species <- as.factor(tip_species)


p <- ggtree(tree1)

#Loop over species to color clades
for(sp in names(species_colors)) {
  #Get tips belonging to this species
  tips <- tree1$tip.label[tip_species == sp]
  #Get MRCA node
  node <- getMRCA(tree1, tips)
  if(!is.null(node)) {
    p <- p + geom_hilight(node = node, fill = species_colors[sp], alpha = 0.3)
  }
}

#Add tip labels
p <- p + geom_tiplab()
p

--Manually rotate some clades

p + geom_text2(aes(subset = !isTip, label = node), hjust = -0.3)

--Rotate the following clades

p <- ggtree(tree1)

p2 <- rotate(p, node = 55)
p2 <- rotate(p2, node = 62)
p2 <- rotate(p2, node = 50)
p2 <- rotate(p2, node = 49)
p2 <- rotate(p2, node = 69)
p2 <- rotate(p2, node = 83)
p2 <- rotate(p2, node = 88)
p2 <- rotate(p2, node = 71)
p2 <- rotate(p2, node = 74)


for(sp in names(species_colors)) {
  tips <- tree1$tip.label[tip_species == sp]
  node <- getMRCA(tree1, tips)  # this works on phylo object
  if(!is.null(node)) {
    p2 <- p2 + geom_hilight(node = node, fill = species_colors[sp], alpha = 0.3)
  }
}

p2 <- p2 + geom_tiplab()
p2

-->Reproducing Figure 1B

--Need to be in results/RAxML/10 Mb-concatenation folder 

library(ape)
library(phangorn)
library(phytools)
library(ggplot2)

tree_files <-list.files(pattern="\\.raxml.bestTree$") #List all .bestTree files. $ ensures the end of the name

trees<- list() # list with all the trees
class(trees)<- "multiPhylo" #make it a multiphylo object for ease of use with other 

i<-1
for(tree_file in tree_files){ ##go thru each file and read the tree
  trees[[i]]<- read.tree(tree_file)
  i<-i+1
}

--Root all trees by outgroup
#re-reroot all our gene trees by the respective outgroup
for(i in 1:length(trees)){
  trees[[i]]<- root(trees[[i]],
                         outgroup = "H_vulgare_HVens23",
                         resolve.root=TRUE)
  trees[[i]]<-chronos(trees[[i]]) ## make ultrametric for nicer densitree
}

--Create a consensus parsimony spupertree

st<-superTree(trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)

--Plot the same density tree as Figure 1B

tree1ultra = chronos(tree1)

png(filename="../../../figures/figure1b.png", width = 1800, height = 900, units = "px")
par(mfrow=c(1,2), mar = c(0.1, 0.1, 0.1, 0.1))
plot(tree1ultra, show.tip.label = FALSE)
densiTree(trees,consensus=tree1, direction='leftwards', scaleX=T,type='cladogram', alpha=0.1)

-->New coalescent-based species tree (ASTRAL4)

--Plot the individual level species along with the full concatenation tree (tree1. In results, run this code:

tree3 = read.tree(file="07-individual-species-tree-astral4.tre")
tree3 = root(tree3,outgroup="H_vulgare_HVens23", resolve.root=TRUE)

par(mfrow=c(1,2), mar = c(0.1, 0.1, 0.1, 0.1))
plot(tree1)
plot(tree3)


--Calculate the RF distance

RF.dist(tree1,tree3)

--Compoare the species tree 

tree4 = read.tree(file="07-species-tree-astral4.tre")
tree4 = root(tree4,outgroup="H_vulgare", resolve.root=TRUE)

par(mfrow=c(1,2), mar = c(0.1, 0.1, 0.1, 0.1))
plot(tree1)
plot(tree4)



# 09 Species network

-->Creating input files

--In Julia, in the code folder

using PhyloNetworks
using SNaQ
using CSV, DataFrames

mappingfile = CSV.read("../results/07-species_mapping.txt", DataFrame; header=false, delim=' ')

rename!(mappingfile, :Column1 => :individual)
rename!(mappingfile, :Column2 => :species)

select!(mappingfile,[:species, :individual])

##We want to remove 3 of the 4 outgroups to simplify the analysis:
##We keep `H_vulgare_HVens23`, 
##We remove `Ta_caputMedusae_TB2`, `Er_bonaepartis_TB1`, `S_vavilovii_Tr279`
filter(row -> row.species in ["Ta_caputMedusae", "Er_bonaepartis", "S_vavilovii"], mappingfile)
size(mappingfile) ## (47,2)

mappingfile = filter(row -> !(row.species in ["Ta_caputMedusae", "Er_bonaepartis", "S_vavilovii"]), mappingfile)
size(mappingfile) ## (44, 2)

CSV.write("../results/09-species_mapping.csv", mappingfile)

##mappingfile = CSV.read("../results/09-species_mapping.csv", DataFrame)
taxonmap = Dict(r[:individual] => r[:species] for r in eachrow(mappingfile)) # as dictionary

#Read the gene trees and compute CF table

--> Running SNaQ to infer phylogenetic network

--Create a subfolder snaq in results folder

--Open Terminal (in the code folder) to start a Julia session with multithreads 

julia -t 2 

--Inside Julia: 

using Distributed
addprocs(4)

@everywhere using PhyloNetworks
@everywhere using SNaQ

##read table of CF
d_sp = readtableCF("../results/09-tableCF_species.csv"); # "DataCF" object for use in snaq!
#read in the species tree from ASTRAL as a starting point
T_sp = readnewick("../results/07-species-tree-astral4.tre")

net = snaq!(T_sp, d_sp, runs=100, Nfail=200, filename= "../results/snaq/09-snaq-h1",seed=8485);

--Plot

using PhyloPlots
#net = readnewick("(Ae_sharonensis,Ae_longissima,(Ae_bicornis,(Ae_searsii,((Ae_tauschii,(((Ae_uniaristata,Ae_comosa)1:0.4918206502664954,(Ae_caudata,Ae_umbellulata)1:0.13449338165653227)1:0.00911821493436927,((T_boeoticum,T_urartu)1:1.6460105085783057,(H_vulgare,((Ae_speltoides,Ae_mutica)1:0.07124470208266999)#H26:0.14159810198824307::0.7563186990617421)1:0.1869824746969994)1:0.46325640725730144)0.99651:0.06048455134529327):0.16788568790821257,#H26:0.0::0.2436813009382579):0.45932787904385297)1:0.9296082436533977)1:0.5926597507029276)1;")
plot(net, showedgenumber=true)

--Root on the outgroup

rootonedge!(net, 16)
rotate!(net,22)
rotate!(net,23)
rotate!(net,-6)
rotate!(net,12)
rotate!(net,11)
plot(net, showgamma=true)

--Replicate the same colors seen in Figure 5 

using DataFrames

tipnodes = [n.number for n in net.node if n.leaf]
tipnames = [n.name for n in net.node if n.leaf]

tipcolors = Dict(
    "T_urartu" => "darkolivegreen",
    "T_boeoticum" => "darkolivegreen",
    "Ae_comosa" => "chocolate",
    "Ae_uniaristata" => "chocolate",
    "Ae_caudata" => "khaki",
    "Ae_umbellulata" => "gold",
    "Ae_tauschii" => "red",
    "Ae_longissima" => "mediumorchid",
    "Ae_sharonensis" => "mediumorchid",
    "Ae_bicornis" => "mediumorchid",
    "Ae_searsii" => "mediumorchid",
    "Ae_mutica" => "dodgerblue",
    "Ae_speltoides" => "navy"
)

colors = [get(tipcolors, name, "black") for name in tipnames]

nodelabel = DataFrame(
    number = tipnodes,
    label = tipnames,
    nodelabelcolor = colors
)

plot(net, nodelabel = nodelabel)

#Result doesn't look like Figure 5 



# 10 Hybrid detection with HyDe

--Move to code folder and open Python with python3

-->Data Preparation for running HyDe 

#We need to convert our sequences (in fasta format) to Phylip format by using the BioPython module.

from Bio import AlignIO
from Bio.AlignIO.PhylipIO import SequentialPhylipWriter

data_path = "../data/Wheat_Relative_History_Data_Glemin_et_al/"
concat_file = "triticeae_allindividuals_OneCopyGenes.fasta"

fasta_file = data_path+concat_file
phylip_file = "../results/10-triticeae_allindividuals_OneCopyGenes.phylip"

#For the following codes. A helpful tip is copy line; hit enter; hit tab and copy next line; repeat until entire block is copied then hit return twice (enter x 2).

# Load the alignment
with open(fasta_file, "r") as f_in:
    alignment = AlignIO.read(f_in, "fasta") --> where I'm stuck 

# Find the length of the longest sequence ID to prevent truncation
max_id_len = max(len(record.id) for record in alignment)

# Write out the alignment (here)
with open(phylip_file, "w") as f_out:
    writer = SequentialPhylipWriter(f_out)
    writer.write_alignment(alignment, id_width=max_id_len + 3) # 3 additional padding

-->Running HyDe on the full concatenation 

--Move to the code folder

run_hyde.py -i ../results/10-triticeae_allindividuals_OneCopyGenes.phylip -m ../results/07-species_mapping.txt -o H_vulgare -n 47 -t 17 -s 11354214 --prefix 10-hyde

--Move results to the results folder

mv 10-hyde-out-filtered.txt ../results
mv 10-hyde-out.txt ../results

-->Running HyDe on 10mB windows
















