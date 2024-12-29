**GSE ID**: _GSE58831_
**Title:** _Gene expression data bone marrow CD34+ cell of patients with myelodysplastic syndromes (MDS) and healthy controls._
**Organism:** _Homo sapiens_
**Download:** _Series Matrix File(s)_

-------------------------------------------------------------------------------------------------------------------------------------------------------------

Set working directory
> setwd("~/MDS") ## (“~/path/to/your/directory")

Install necessary packages and Load necessary libraries for data analysis
> install.packages(“GEOquery”) # # For accessing Gene Expression Omnibus (GEO) data
> library(“GEOquery”)
> install.packages(“dplyr”) ## For data manipulation and transformation
> library(“dplyr”)
> install.packages(“tidyverse”) ## For data wrangling and visualization
> library(“tidyverse”)
> install.packages(“data.table”) ## For efficient data manipulation with data tables
> library(“data.table”)

Setting an environment variable
Set the environment variable "VROOM_CONNECTION_SIZE" to a specific value
This adjusts the memory allocation for the vroom package, which is used for fast reading and writing of rectangular data files

> Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)

Data Download
Downloading the GEO dataset with the ID "GSE58831" using the getGEO() function and enabling
GSEMatrix retrieval
> geo_id <- "GSE58831"
> gse <- getGEO(geo_id, GSEMatrix = TRUE)

Extract all three types of data
Phenodata: Information about characteristics and attributes of biological samples.
> phenoData <- pData(phenoData(gse[[1]]))
> phenoData_selected <- phenoData[, c(2, 12)]
> dim(phenoData_selected)

Feature data: Numeric or categorical values describing specific attributes or variables.
> featureData <- fData(gse[[1]])
> featureData <- featureData[, c(11, 12)]
> dim(featureData)

Expression data: Measurements indicating gene or protein expression levels in samples.
> expr_data <- exprs(gse[[1]])
> dim(expr_data)

Transpose dataframe
Transpose dataframe
> expr_data.t <- t(expr_data)
> dim(expr_data.t)

Merge datasets
Merge phenoData_selected and expr_data.t
> expr_plus_pheno <- cbind(phenoData_selected, expr_data.t)
> expr_plus_pheno$characteristics_ch1.2 <- gsub("disease status:", "", expr_plus_pheno$characteristics_ch1.2)
> expr_plus_pheno$characteristics_ch1.2 <- gsub("control", "", expr_plus_pheno$characteristics_ch1.2)
> dim(expr_plus_pheno)

Re-transpose it
> expr_plus_pheno.ret <- t(expr_plus_pheno)
> dim(expr_plus_pheno.ret)

Merge featureData and expr_plus_pheno.ret
> expr_plus_pheno_feature <- cbind(featureData, expr_plus_pheno.ret)
> dim(expr_plus_pheno_feature)

Perform the collapse rows function from WGCNA
Prepare dataset for collapse row function
> expr_plus_pheno_feature[expr_plus_pheno_feature$ENTREZ_GENE_ID == '', ] <- NA ### substitute empty cells with NA
> expr_plus_pheno_feature <- na.omit(expr_plus_pheno_feature)
> New <- strsplit(as.character(expr_plus_pheno_feature$ENTREZ_GENE_ID), '///')
> New <- as.data.frame(sapply(New, "[[", 1))
> expr_plus_pheno_feature_new <- cbind(New, expr_plus_pheno_feature)
> expr_plus_pheno_feature_new <- expr_plus_pheno_feature_new[, -(2:3)]
> colnames(expr_plus_pheno_feature_new)[1] <- "entrez_id"
> expr_plus_pheno_feature_new <- data.frame(lapply(expr_plus_pheno_feature_new, function(x)
as.numeric(as.character(x))), check.names = F, row.names =
rownames(expr_plus_pheno_feature_new))
> dim(expr_plus_pheno_feature_new)

Collapse rows in WGCNA
> install.package("WGCNA") ## Co-expression network analysis
> library("WGCNA")
> expr_plus_pheno_feature_new.collapse <- collapseRows(datET = expr_plus_pheno_feature_new, >
rowGroup = (expr_plus_pheno_feature_new[, 1]), rowID =
rownames(expr_plus_pheno_feature_new))
> expr_plus_pheno_feature_new.collapsed 1 <-
expr_plus_pheno_feature_new.collapsed$datacollapsed
> expr_plus_pheno_feature_new.collapsed 1 <-
as.data.frame(expr_plus_pheno_feature_new.collapsed 1)
> dim(expr_plus_pheno_feature_new.collapsed 1)

Generate that box plot
Generate that box plot
> names(expr_plus_pheno_feature_new.collapsed 1)
> boxplot(expr_plus_pheno_feature_new.collapsed 1[, 1:50])

Normalization of dataset
Assuming your gene expression data is stored in a variable named 'expression_data'
Install and load the preprocessCore package
> install.packages("preprocessCore") ## Preprocessing and normalization of data
> library(preprocessCore)

Perform quantile normalization
> df_norm <-as.data.frame(normalize.quantiles(as.matrix(expr_plus_pheno_feature_new.collapsed 1)))

Add column and row names
> colnames(df_norm) <- colnames(expr_plus_pheno_feature_new.collapsed 1)
> rownames(df_norm) <- rownames(expr_plus_pheno_feature_new.collapsed 1)

Box plot after normalization
> boxplot(df_norm[, 1:50])

DEGs analysis
> install.packages("limma") ## Differential expression analysis for high-dimensional genomics data
> library(limma)
> df_norm.t <- t(df_norm)
> labels <- as.data.frame(rownames(df_norm.t))
> labels[, 1] <- gsub("X.", "", labels[, 1])
> labels[, 1] <- gsub("MDS.*", "MDS", labels[, 1])
> labels[, 1] <- gsub("Healthy.*", "Healthy", labels[, 1])
> group <- labels[, 1]
> design <- model.matrix(~0 + group)
> contrast <- makeContrasts(groupMDS - groupHealthy, levels = design)
> fit <- lmFit(df_norm, design)
> fit2 <- contrasts.fit(fit, contrast)
> fit2 <- eBayes(fit2)
> DEGs <- topTable(fit2, adjust = "BH", number = Inf)

Filtering based on Log2FC and p.value or p.adj
> filtered_DEGs <- DEGs[(DEGs$logFC >=1 | DEGs$logFC <=-1) & DEGs$adj.P.Val <=0.05, ]

Filtered out DEGs based on some p value and Fold change criteria. Add a gene symbol to this object.
Feature data
> featureData <- fData(gse[[1]])
> featureData <- featureData[,c(12,11)]
> featureData[featureData$ENTREZ_GENE_ID=='',] <- NA ## substitute empty cells with NA
> featureData <- na.omit(featureData)
> t1 <- strsplit(as.character(featureData$ENTREZ_GENE_ID),'///')
> t1 <- (sapply(t1, "[[", 1))
> t2 <- strsplit(as.character(featureData$`Gene Symbol`),'///')
> t2 <- (sapply(t2, "[[", 1))
> New_merged <- as.data.frame(cbind(t1, t2)) ## This is a list where we have gene names as per entrez id

To get unique gene list
> New_merged1 <- New_merged[!duplicated(New_merged), ]
> rownames(New_merged1) <- New_merged1$t1

Merge it with DEGs object
> BiocManager::install('EnhancedVolcano') ## For advanced volcano plot visualization
> library("EnhancedVolcano")
> names <- rownames(DEGs)
> newww <- as.data.frame(cbind(names, DEGs))
> m1 <- merge(New_merged1, newww, by.x = "t1", by.y = "names")
> rownames(m1) <- m1[,2]

Generate Volcano plot
Generate volcano plot
> EnhancedVolcano(m1,
lab = rownames(m1),
x = 'logFC',
y = 'P.Value')

for more specific cut off
> EnhancedVolcano(m1,
lab = rownames(m1),
x = 'logFC',
y = 'P.Value',
pCutoff = 0.05,
FCcutoff = 1.0,
pointSize = 3.0,
labSize = 3.0)

Heatmap
Fetch DEGs ids gene expression values from our normalized object
> columns2 <- c("8788","5746","3047","1824","1287","10008","3434","3039","3429","63895","10643",
"81857","81831","29968","7164","55084","4199","360","51280","26050","10379","474344","8645","10007","
57110","55384","132430","10964","10924","151176","400916","56000","146849","23158","29125","57211","
140883","79961","6565","3339","57628","4778","8328","55340","4118","10810","8519","11019","8175","556
15","85439","389831","114880","81035","4938","2968","3437","10388","121599","2535","5784","4493","28
5533","4502","5580","64393","83876","10025","55691","29126","4495","375387","1378","8942","7052","10
498","28560","388685","94121","11282","254531","89781","2549","29128","114548","253430","7262","7071
","54478","3036","84186","115350","55013","116966","64651","2272","22822","55278","10892","195828","5
209","2069","54880","3399","931","79805","8404","4086","414332","79884","57217","29933","9134","2846
1","5522","6583","10734","374393","3007","221061","81575","3394","440423","930","467","1393","51237",
"170482","79733","1901","79161","89796","23169","2355","26776","10123","2004","100128590","768211","
4680","6347","199675","8835","1511","57405","55086","971","84288","1400","9262","93474","54206","554
22","2824","140597","150094","1956","30014","5789","619279","3543","3981","58494","646576","144481",
"1289","6696","83737","2669","49854","2012","1326","7852","84419","6003","8553","54436","51655","973"
,"1870","55790","83891","768206","23462","4885","3569","84620","5055","140733","5791","56937","1178",
"3164","144455","1390","4616","7837","640","7466","2308","9953","54463","5329","928","288","54541","3
280","100996579","152189","23764","5996","6328","1960","1490","114757","340152","54210","78986","151
473","6498","387914","5101","753","284207","152687","3669","80243","1880","932","26051","145474","63
64","387694","131583","3662","1991","9783","55824","3126","101929623","90102","257194","57824","4929
","374","3575","94241","60468","4638","4050","8870","5553","285097","5142","9214","22795","100302650"
,"64399","11245","79661","94120","29760","641518","199786","4852","26659","11067","4092","974","8018
3","8974","51363","51176","100133941","26018","8013","28387","5897","5450","1791","4311","100507254"
,"5079","7441","29802","9934","9590","10777","1879","5896")
> selected_df_new <- df_norm.t[, columns2]
> selected_df_new.t <- t(selected_df_new)

Here again to add gene names we have to merge it with the feature data
> t3 <- as.data.frame(rownames(selected_df_new.t))
> selected_df_new.t2 <- cbind(t3, selected_df_new.t

New_merged1 this is a non duplicated feature data contain gene name and entrez id and here we are merging
based on entrez id
> m2 <- merge(New_merged1, selected_df_new.t2, by.x = "t1", by.y = "rownames(selected_df_new.t)")
> rownames(m2) <- m2[,2]
> m2 <- m2[,-(1:2)]
> typeof(m2)
> install.package(“gplots”) ## For advanced plotting functionalities in R.
> install.package("RColorBrewer") ## For advanced color palettes in R.
> library(gplots)
> library("RColorBrewer")

pheatmap(m2,cluster_cols = T,scale = 'row',breaks = seq(2.25, 15.87, length.out = 101))
> m2 <- as.matrix(m2)
> m2.t <- t(m2)
> hmcol = colorRampPalette(brewer.pal(9, "YlOrRd"))(100)
> library("gplots")
> heatmap.2(m2.t, col = hmcol, trace="none", symbreaks = TRUE, dendrogram= "both", symkey =
FALSE, margin=c(7, 5), cexCol=0.6, cexRow=0.45, density.info="none",breaks = seq(2, 16, length.out
= 101))
> dev.off()
