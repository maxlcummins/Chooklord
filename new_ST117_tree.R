library(dplyr)
library(ggtree)
library(phytools)
library(readr)
library(magrittr)
#Read
APEC_Origins_BigChook <- read_csv("~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/ST117_Origins.csv")

hit_table2 <- read_csv("~/Dropbox/Doctorate/Manuscripts/BigChook/Figures/Working Figures/BigChook/Supplementary Figures/Supplementary_Table_3_Genotype_Phylogeny.csv")

hit_table2<- hit_table2[,5:ncol(hit_table2)]

gsub(pattern = "-.*", replacement = "", x = hit_table2$phyloname) -> hit_table2$sample_name

genotype_metadata <- left_join(hit_table2, APEC_Origins_BigChook)

ST117 <- filter(genotype_metadata, ST %in% c(117,4045))

rownames(ST117) <- ST117$phyloname

tree <- read.tree(file = "./SNP_ST117_A160_unicycler_minusref.tree")

midpoint.root(tree = tree) -> tree

state <- c("QLD", "NSW", "VIC", "WA", 0)
state_cols <- c("#d7191c",
  "#fdae61",
  "#6B8E23",
  "#2c7bb6",
  "black")

empty_df <- cbind(gsub(pattern = "UTS(A[0-9]+).out.*", replacement = "\\1", x = tree$tip.label), rep(1,times=length(tree$tip.label)))

colnames(empty_df) <- c("sample_name", "x")

genotype <- ST117[,c(1,7:120)]

#splits the gene hit types to separate dfs to be processed
i <- genotype %>% select(starts_with("i_"))
r <- genotype %>% select(starts_with("r_"))
v <- genotype %>% select(starts_with("v_"))
p <- genotype %>% select(starts_with("p_"))

#remove cols with no hits
i <- i[,colSums(i) > 0]
r <- r[,colSums(r) > 0]
v <- v[,colSums(v) > 0]
p <- p[,colSums(p) > 0]

#change gene hits to a different number based on gene type
i[i > 0] <- "I"
i[i == 0] <- "NA"
r[r > 0] <- "R"
r[r == 0] <- "NA"
v[v > 0] <- "V"
v[v == 0] <- "NA"
p[p > 0] <- "P"
p[p == 0] <- "NA"

#combines processed sub-dfs
genotype <- cbind(i,r,v,p)

#orders the columns based on their names
genotype <- genotype[ , order(names(genotype))]

# #trim off the r_ etc from start of each gene hit
colnames(genotype) <- gsub(pattern = "^[r,i,v,p]_","",colnames(genotype), perl = TRUE)

#Define colors for gene-type dependent coloring of gene hits
colorgenotype <- c("NA" = "white", "I" = "#8dd3c7", "R" = "#bebada", "V" = "#fb8072", "P" = "#80b1d3")

df_for_tiplabs <- left_join(as.data.frame(empty_df), ST117)

tree$tip.label <- df_for_tiplabs$phyloname

rownames(ST117) -> save

#ST117 <- ST117[ , colSums(is.na(ST117)) != nrow(ST117)]

ST117_meta <- ST117[,c(1,123:ncol(ST117))]

ST117_geno <- ST117[,7:114]

rownames(ST117) <- save

rownames(ST117_meta) <- save

rownames(genotype) <- save

p <- ggtree(tree) %<+%
  ST117 +
  geom_tiplab(size = 2, hjust = -.125, align = TRUE, aes(color = Origin)) +
scale_color_manual(breaks = c(state), values = c(state_cols), na.value = "black")

a <- gheatmap(p = p, data = genotype,
              colnames_offset_y = -0.4,
              font.size = 1.5,
              hjust = 0,
              colnames_position = "top",
              colnames_angle = 45,
              width = 11,
              color = "black",
              offset = 1) + 
  scale_fill_manual(values = colorgenotype) +
  theme(legend.position = "none")

a