##### Tutorial Sparse Network And Component (SNAC) analysis #####

#Install SNACR package
library(devtools)
devtools::install_github("PiaTio/SNACR")
library(SNACR)

#Initialize random generator
set.seed(12) 
NIMH_stand[,1:26] # bio-data
NIMH_stand[,27:31] # psycho-data
groups <- c(rep("Biology", 26), rep("Psycholgy", 5))
labels <- c("Chloride", "Creatine", "Glucose", "Potassium", "Sodium",
            "Total carbon dioxide (CO2)", "Urea-nitrogen (BUN)",
            "Alanine transferase (ALT)", "Alkine phosphatase (ALP)",
            "Aspertate transferase (AST)", "Bilirubin", "Basophil", "Eosinophil",
            "Hematocrit", "Hemoglobin (Hb)", "Lymphocytes", "Mean corpuscular hemoglobin (MCH)",
            "Mean corpuscular hemoglobin concentration (MCHC)", "Mean corpuscular volume (MCV)",
            "Monocytes", "Mean platelet volume (MPV)", "Neurotrophils", "Red blood count (RBC)",
            "Creatine kinaase (CK)", "C-reactive protein (CRP)", "Thyroid stimulating hormone (TSH)",
            "Alcohol use", "Verbal IQ", "Nonverbal IQ", "Mental effort - fatigue",
            "Mental effort - performance")

### Estimate EBICglasso network
library(bootnet)
vsize <- 4.5
net_ebic <- estimateNetwork(NIMH_stand, "EBICglasso")
ebic = plot(net_ebic, cut = 0, legend.mode = "style2", nodeNames = labels,
            groups = groups, legend.cex = 1, maximum = 0.8, vsize = vsize)

### SNAC’s component step
R <- 3
Target <- matrix(c(1, 0, 0, 1, 1, 1), ncol = R)
com <- comSNAC(data = NIMH_weight, b1 = 26, b2 = 5, R = R, Target = Target)
com$lambda_com
com$Structure
com$Structure == "common"

### SNAC’s network step
net <- netSNAC(Pcommon = com$Pmatrix[,3], Tcommon = com$Tmatrix[,3])
net$lambda_net
net_block1 <- netSNAC(Pcommon = com$Pmatrix[,1], Tcommon = com$Tmatrix[,1])
net_block2 <- netSNAC(Pcommon = com$Pmatrix[,2], Tcommon = com$Tmatrix[,2])

### Plot common and domain-specific networks
library(qgraph)
net_com <- (qgraph(net$graph, nodeNames = labels, legend = FALSE, layout = ebic$layout,
                   palette = "colorblind", theme = "colorblind", vsize = vsize, cut = 0,
                   maximum = 0.8, groups = groups))
net_bio <- (qgraph(net_block1$graph, nodeNames = labels, legend = FALSE,
                   layout = ebic$layout, palette = "colorblind", theme = "colorblind",
                   vsize = vsize, cut = 0, maximum = 0.8, groups = groups))
net_psy <- (qgraph(net_block2$graph, nodeNames = labels, legend = FALSE,
                   layout = ebic$layout, palette = "colorblind", theme = "colorblind",
                   vsize = vsize, cut = 0, maximum = 0.8, groups = groups))

### Boostrap analysis
#Use lambda values for lasso and graphical lasso from the main analysis
Network <- estimateNetwork(data = NIMH_weight, fun = bootSNAC, b1 = 26,
                           b2 = 5, R = 3, Target = matrix(c(1, 0, 0, 1, 1, 1), ncol = 3), comp = 3,
                           lambda_com = com$lambda_com, lambda_net = net$lambda_net)

boot_check <- bootnet(Network, memorysaver = FALSE)
plot(boot_check, order = "sample", labels = FALSE)

bootnet_res <- bootSNAC_sum(boot_check)
bootnet_res$lambdaCom
bootnet_res$lambdaNet
bootnet_res$Structure
bootnet_res$cocorec

#Use lambda value for graphical lasso from the main analysis
Network_net <- estimateNetwork(data = NIMH_weight, fun = bootSNAC, b1 = 26,
                           b2 = 5, R = 3, Target = matrix(c(1, 0, 0, 1, 1, 1), ncol = 3), comp = 3,
                           lambda_com = com$lambda_com)
boot_check_net <- bootnet(Network_net, memorysaver = FALSE)
bootnet_res_net <- bootSNAC_sum(boot_check_net)

#Use lambda value for lasso from main analysis
Network_com <- estimateNetwork(data = NIMH_weight, fun = bootSNAC, b1 = 26,
                           b2 = 5, R = 3, Target = matrix(c(1, 0, 0, 1, 1, 1), ncol = 3), comp = 3,
                           lambda_net = net$lambda_net)
boot_check_com <- bootnet(Network_com, memorysaver = FALSE, nBoots = 200)
bootnet_res_com <- bootSNAC_sum(boot_check_com)

#Do not use any lambda values from main analysis; both lambda's are estimated using 10-fold cross-validation
Network_comnet <- estimateNetwork(data = NIMH_weight, fun = bootSNAC, b1 = 26,
                                  b2 = 5, R = 3, Target = matrix(c(1, 0, 0, 1, 1, 1), ncol = 3), comp = 3)
boot_check_comnet <- bootnet(Network_comnet, memorysaver = FALSE, nBoots = 200)
bootnet_res_comnet <- bootSNAC_sum(boot_check_comnet)