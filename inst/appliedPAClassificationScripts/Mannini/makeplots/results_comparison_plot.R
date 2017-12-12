rm(list=ls())
library("ggplot2")
library("grid")


proportions_correct <- data.frame(row.names = c(paste("Subject", as.character(1:33)), "Overall"))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_linear_MLR.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$linear_MLR <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_quadratic_MLR.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$quadratic_MLR <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_SVM.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$SVM <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_RF.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$RF <- c(temp, mean(temp))


proportions_correct_long <- reshape(proportions_correct, varying = list(names(proportions_correct)), times = names(proportions_correct), ids = c(as.character(1:33), "Overall"), direction = "long")
rownames(proportions_correct_long) <- NULL
colnames(proportions_correct_long) <- c("Classifier", "Proportion_Correct", "Subject")

p_static <- ggplot(proportions_correct_long[proportions_correct_long$Subject != "Overall", ]) + 
  geom_line(aes(x = Subject, y = Proportion_Correct, group = Classifier, colour = Classifier), size = 1) +
  geom_hline(aes(yintercept = Proportion_Correct, colour = Classifier), size = 1.5, linetype = 2, data = proportions_correct_long[proportions_correct_long$Subject == "Overall", ]) +
  xlim(as.character(1:33)) +
  ylim(c(0, 1)) +
  theme_bw()
p_static



proportions_correct <- data.frame(row.names = c(paste("Subject", as.character(1:33)), "Overall"))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_rfCRF_globalSeq_fullTrans.Rdata")
proportions_correct$rfCRF <- c(combined_prop_correct, mean(combined_prop_correct, na.rm = TRUE))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_boostCRF_globalSeq_randomBeta_fullTrans.Rdata")
combined_prop_correct[25] <- 0.7052239
proportions_correct$boost_quadratic_CRF <- c(combined_prop_correct, mean(combined_prop_correct, na.rm = TRUE))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_baggedGradientTreeBoostCRF_globalSeq_fullTrans.Rdata")
proportions_correct$bagged_gradient_tree_boost_CRF <- c(combined_prop_correct, mean(combined_prop_correct, na.rm = TRUE))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_gradientTreeBoostCRF_globalSeq_fullTrans.Rdata")
proportions_correct$gradient_tree_boost_CRF <- c(combined_prop_correct, mean(combined_prop_correct, na.rm = TRUE))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_McShane_RF_HMM_FullTrans.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$McShane_RF_HMM <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_McShane_linear_MLR_HMM_FullTrans.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$McShane_linear_MLR_HMM <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_McShane_quadratic_MLR_HMM_FullTrans.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$McShane_quadratic_MLR_HMM <- c(temp, mean(temp))

load("C:/Stat/HMM/HMMEnsembles/HMMApplication/Mannini/combinedresults/combined_res_SVM.Rdata")
temp <- sapply(combined_results, function(comp) comp$prop_correct)
proportions_correct$SVM <- c(temp, mean(temp))


proportions_correct_long <- reshape(proportions_correct, varying = list(names(proportions_correct)), times = names(proportions_correct), ids = c(as.character(1:33), "Overall"), direction = "long")
rownames(proportions_correct_long) <- NULL
colnames(proportions_correct_long) <- c("Classifier", "Proportion_Correct", "Subject")

p_dynamic <- ggplot(proportions_correct_long[proportions_correct_long$Subject != "Overall", ]) + 
  geom_line(aes(x = Subject, y = Proportion_Correct, group = Classifier, colour = Classifier), size = 1) +
  geom_hline(aes(yintercept = Proportion_Correct, colour = Classifier), size = 1.5, linetype = 2, data = proportions_correct_long[proportions_correct_long$Subject == "Overall", ]) +
  xlim(as.character(1:33)) +
  ylim(c(0, 1)) +
  theme_bw()

p_dynamic



