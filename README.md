#Repository for PL PATH 875 Phylogenetics Practicum (Spring 2026)



# 03 Data and Software Setup

- Clone repository: git clone https://github.com/solislemuslab/phylo-practicum.
- Download data and move into data subfolder and unzip files

--> Installation

(1) RAxML Next Generation

- Download software, place into designated local folder (folder labeled software)

- Make executable runnable (in Terminal)

cd ~/software/raxml-ng_v1.2.2_macos
chmod + x raxml-ng
raxml -ng -v (to very installation)

-Add to PATH to run raxml-ng from any where 

nano ~/.zshrc
export PATH="$HOME/software/raxml-ng_v1.2.2_macos:$PATH" (save, exit, and reload to confirm changes)

(2) SuperTriplets 

-Check if Java is installed

java -h 

-Download SuperTriplets java file and put into software folder inside a subfolder named supertriplets


(3) ASTRAL

-Install with conda

conda config --add channels bioconda
conda install aster

-To check installation location and existence

conda list
wastral

(4) Julia and packages (in Terminal)

- Install Julia 

curl -fsSL https://install.julialang.org | sh

-Run this command after Terminal restart

juliaup add release

-Add other packages

julia then ] to go into package mode

add PhyloNetworks
add SNaQ
add PhyloPlots
add CSV
add DataFrames

Hit backspace to leave package mode and leave Julia by typing exit()


(4) HyDe

-Move to software folder and download 

git clone https://github.com/pblischak/HyDe.git

-Install HyDe

cd HyDe
python3 -m pip install -r requirements.txt
python3 -m pip install .

-Check installation

run_hyde.py



# 04 Gene tree estimation

--> Running RAxML on 10 genes 

-Move to data/Wheat_Relative_History_Data_Glemin_et_al folder and create OneCopyGene-trimmed folder for the first 10 genes

mkdir OneCopyGenes-trimmed
cd OneCopyGenes
ls | head -n10 | xargs -I {} cp "{}" ../OneCopyGenes-trimmed

-Move to the cold foler to run 04-raxml.sh script

./04-raxml.sh

-Analyze inferred trees in R

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

# Extract the first matching species for each tree
gene_tree_outgroup<- rep(NA,length(gene_trees))
for(i in 1:length(gene_trees)){
  gene_tree <- gene_trees[[i]]
  found_taxa <- root_taxa[root_taxa %in% gene_tree$tip.label]
  
  # Return the first one if it exists
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

-Visualize the trees 

plot(gene_trees[[2]]) # plots the second gene tree

#see which genes have which taxa

all_labels <- unlist(lapply(gene_trees, function(x) x$tip.label)) ##count how often each individual appears
df <- as.data.frame(table(all_labels) / length(gene_trees)) #organize as a nice data frame for plotting

ggplot(df, aes(x = all_labels, y = Freq)) +
  geom_col() +
  labs(x = "Individual", y = "Proportion")+
  theme(axis.text.x = element_text(angle = 90))

-SuperTree

st<-superTree(gene_trees)
st<-root(st,"H_vulgare_HVens23",resolve.root = T)
plot(st)

-Densitree 

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

## to remove NULLs:
## trees_pruned <- trees_pruned[!sapply(trees_pruned, is.null)]

densityTree(trees_pruned,use.edge.length=FALSE,type="cladogram",nodes="centered")

-Comparison of ML trees for each of the 20 RAxML runs (visualization of the best tree across each run)

trees = read.tree(file="Ae_bicornis_Tr406_BIS2_Contig10132_simExt_macseNT_noFS_clean.aln.raxml.mlTrees")

rtrees <- lapply(trees, function(tr) {
  root(tr, outgroup = "H_vulgare_HVens23", resolve.root = TRUE)
})

Plot densiTree before and after rooting

densityTree(trees,type="cladogram",nodes="intermediate")
densityTree(rtrees,type="cladogram",nodes="intermediate")

densityTree(trees,use.edge.length=FALSE,type="cladogram",nodes="centered")
densityTree(rtrees,use.edge.length=FALSE,type="cladogram",nodes="centered")

--> Running RAxML on all on all the genes

-Modify 04-raxml.sh script to run all genes 

DATADIR="../data/Wheat_Relative_History_Data_Glemin_et_al/OneCopyGenes"

-Inside code folder run the script

cd code
./04-raxml.sh

-Visualize the trees 

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



#05 Species tree supermatrix (full) 






