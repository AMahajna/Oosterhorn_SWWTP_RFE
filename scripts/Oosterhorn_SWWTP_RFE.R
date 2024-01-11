################################################################################
################################################################################
##Exploratory Data Analysis 

################################################################################
##install and load packages 

source(file = "scripts/install_load_packages.r")

################################################################################
##Create folder for project organization

if(!dir.exists("input_data")){dir.create("input_data")}
if(!dir.exists("output_data")){dir.create("output_data")}
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("scripts")){dir.create("scripts")}

################################################################################
##Load data
glycerol_methanol = read_excel("input_data/ZAWZI_INF_glycerol_methanol.xlsx")
additional_data = read_excel("input_data/ZAWZI_additional_data.xlsx")
inf_lab = read_excel("input_data/ZAWZI_INF_lab.xlsx")
eff_data = read_excel("input_data/ZAWZI_INF_lab.xlsx")
online_data = read_excel("input_data/ZAWZI_Online_data.xlsx") 
SVI_data = read_excel("input_data/ZAWZI_AT_SVI_DS.xlsx") 

##Load weekly averaged data and rarity from microbiome analysis 
weekly_data = read_excel("input_data/WWTP_processdata_perweek_20180125_NGS_T.xlsx") 
rarity = read_csv("input_data/Rarity.csv")

################################################################################
##EDA

#Sheet: ZAWZI_INF_glycerol_methanol
glycerol_methanol[ , c("Glycerol_kg","Methanol_kg")] %>%
  create_report(
    output_file = "EDA_ZAWZI_INF_glycerol_methanol_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_INF_glycerol_methanol",
    config = configure_report(add_plot_correlation = FALSE))

#Sheet: ZAWZI_additional_data
additional_data[4:8]%>%
  create_report(
    output_file = "EDA_ZAWZI_additional_data_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_additional_data", 
    config = configure_report(add_plot_correlation = FALSE))

#Sheet: ZAWZI_INF_lab
inf_lab[3:20] %>%
  create_report(
    output_file = "EDA_ZAWZI_INF_lab_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_INF_lab", 
    config = configure_report(add_plot_correlation = FALSE))

#Sheet: ZAWZI_Online_data
online_data[3:14] %>%
  create_report(
    output_file = "EDA_ZAWZI_Online_data_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_Online_data",
    config = configure_report(add_plot_correlation = FALSE))

#Sheet: ZAWZI_AT_SVI_DS.csv
SVI_data[3:9] %>%
  create_report(
    output_file = "EDA_ZAWZI_AT_SVI_DS_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_AT_SVI_DS",
    config = configure_report(add_plot_correlation = FALSE))

#Sheet: ZAWZI_EFF_lab
eff_data[3:20] %>%
  create_report(
    output_file = "EDA_ZAWZI_EFF_lab_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - ZAWZI_EFF_lab", 
    config = configure_report(add_plot_correlation = FALSE))

################################################################################
################################################################################
##Pre-processing: correlation analysis

#Use rough process data to conduct correlational analysis for each group
#of parameters in order to utilize the biggest amount of informational http://127.0.0.1:40819/graphics/plot_zoom_png?width=904&height=778
#without information loss resulting from data averaging or data selection 

#scatter plot can help identify correlation which is skewed due to outliers 

################################################################################
##Calculating correlation matrix 
##Plotting correlation heatmaps and exporting the figures 
##Plotting pairplots and exporting the figures 

################################################################################
##inf_lab
inf_corr <- round(cor(inf_lab[3:20],use = 'pairwise.complete.obs'), 2)
inf_testRes= cor.mtest(inf_lab[3:20],conf.level = 0.95)

png("figures/inf_cor.png", units="in", width=5, height=5, res=1000)
corr_inf_plot  =corrplot(inf_corr, order ='hclust', addrec =11 , rect.col = "green",
                         rect.lwd = 3 ,method = 'circle', addCoef.col = 'white',
                         number.digits = 1,number.cex = 0.5,tl.pos ='lt', tl.srt=45, tl.cex =0.6)
dev.off()
#Alternatively, we can see with the p-test that correlation is significant 
#corr_inf_plot = ggcorrplot(inf_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
#                           outline.col = "white",
#                           ggtheme = ggplot2::theme_gray,
#                           colors = c("#6D9EC1", "white", "#E46726"), p.mat= inf_testRes$p)

#Nitrogen speciations
png(filename="figures/pairplot_N.png" ,units = 'in',width=10, height=8, res=1000)
ggpairs(inf_lab[c(8,10:13)],lower = list(continuous = "smooth"))
dev.off()
#Missing data from influent 
png(filename="figures/missing_inf.png" ,units = 'in',width=9, height=6, res=1000)
plot_missing(inf_lab, title = "missing data profile for influent lab analysis")
dev.off()
#selected parameters
inf_lab_selected=  colnames(inf_lab[c(4,5,7,11,14,16,18,20)])
cat("The selected parameters from the influent lab dataset are: ",inf_lab_selected )

selected_parameters = inf_lab_selected

################################################################################
##glycerol_methanol
png(filename="figures/glycerol_pairpanels.png", units = 'in',height = 5, width = 7.5,res =1000)
ggpairs(glycerol_methanol[3:4],lower = list(continuous = "smooth"))
dev.off()
#Alternatively, we can see with the p-test that correlation is significant 
#glycerol_corr <- round(cor(glycerol_methanol[3:4]), 2)
#glycerol_testRes= cor.mtest(glycerol_methanol[3:4],conf.level = 0.95)
#corr_glycerol_plot = ggcorrplot(glycerol_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
#                                outline.col = "white",
#                                ggtheme = ggplot2::theme_gray,
#                                colors = c("#6D9EC1", "white", "#E46726"), p.mat= glycerol_testRes$p)

plot_missing(glycerol_methanol, title = "missing data profile for glycerol and methanol dataset")

glycerol_methanol_selected = colnames(glycerol_methanol[3])
cat("The selected parameter from the glycerol and methanol dataset is: ", glycerol_methanol_selected)

selected_parameters = c(selected_parameters, glycerol_methanol_selected)
################################################################################
##additional_data
png(filename="figures/additional_data_pairspanels.png",units = 'in', width = 10, height = 5, res = 1000)
pairs.panels(additional_data[4:6], cex.cor=0.3 ,cex.labels=1, cex.axis = 1)
dev.off()
#Alternatively
#additional_corr <- round(cor(additional_data[4:6],use = 'pairwise.complete.obs'), 2)
#additional_testRes= cor.mtest(additional_data[4:6],conf.level = 0.95)
#corr_additional_plot = ggcorrplot(additional_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
#                                  outline.col = "white",
#                                  ggtheme = ggplot2::theme_gray,
#                                  colors = c("#6D9EC1", "white", "#E46726"), p.mat= additional_testRes$p)

additional_selected=  colnames(additional_data[4:6])
cat("The selected parameters from the influent lab dataset are: ",additional_selected )

selected_parameters = c(selected_parameters, additional_selected)

################################################################################
##online_data
online_corr <- round(cor(online_data[c(3:6,10)],use = 'pairwise.complete.obs'), 2)
online_testRes= cor.mtest(online_data[c(3:6,10)],conf.level = 0.95)
png(filename="figures/online_heatmap.png", units = "in",height=5, width=5, res = 1000)
ggcorrplot(online_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
           outline.col = "white",
           ggtheme = ggplot2::theme_gray,
           colors = c("#6D9EC1", "white", "#E46726"))
dev.off()

plot_missing(online_data, title = "missing data profile for online dataset")

online_selected=  colnames(online_data[c(3,5,6,10)])
cat("The selected parameters from the online dataset are: ",online_selected )

selected_parameters = c(selected_parameters, online_selected)

################################################################################
##SVI_data
png(filename="figures/SVI_pairspanels.png",units = 'in', height=10, width=10, res = 1000)
ggpairs(SVI_data[3:9], lower = list(continuous = "smooth"))
dev.off()
#Alternatively, 
#SVI_corr <- round(cor(SVI_data[3:9],use = 'pairwise.complete.obs'), 2)
#SVI_testRes= cor.mtest(SVI_data[3:9],conf.level = 0.95)
#corr_SVI_plot = ggcorrplot(SVI_corr, hc.order = TRUE, type= "upper", lab = TRUE, 
#                           outline.col = "white",
#                           ggtheme = ggplot2::theme_gray,
#                           colors = c("#6D9EC1", "white", "#E46726"), p.mat= SVI_testRes$p)

plot_missing(SVI_data, title = "missing data profile for SVI dataset")

SVI_selected=  colnames(SVI_data[5])
cat("The selected parameters from the SVI dataset is: ",SVI_selected )

selected_parameters = c(selected_parameters, SVI_selected)

################################################################################
##Effluent data is not used in our analysis approach
#eff_corr <- round(cor(eff_data[3:20],use = 'pairwise.complete.obs'), 2)
#eff_testRes= cor.mtest(eff_data[3:20],conf.level = 0.95)

###############################################################################
##subsetting weekly data 
#weekly data is averaged values of process parameters selected aligned with NGS samples
#In case week value is missing and doesn't align with NGS, we take the value before
#Because the process parameters impact the microbiome and vice versa 

#After we have selected parameters using collinearity analysis, we subset 
#weekly averaged data using the selected parameters 
selected_parameters = c(selected_parameters, "T_avg_C")
weekly_data = weekly_data[selected_parameters]

#check 
#dim(weekly_data)
#should be 18

#check 
#colSums(is.na(weekly_data))

plot_missing(weekly_data, title = "missing data profile for selected parameters in weekly data")
weekly_data =  subset(weekly_data, select = -c(INF_K_mg_per_l,INF_pH_pH))

#check
plot_missing(weekly_data, title = "missing data profile for selected parameters in weekly data")

###############################################################################
## Multi-collinearity analysis 

#Check for multi-collinearity problem
weekly_data_cor = cor(weekly_data, use = "pairwise.complete.obs")
weekly_data_testRes= cor.mtest(weekly_data,conf.level = 0.95)

#Condition number ratio of max to min Eigen values of the correlation matrix  
kappa(weekly_data_cor,exact=TRUE)
# A condition number between 10 and 30 indicates the presence of multicollinearity
#and when a value is larger than 30, the multicollinearity is regarded as strong.

#rarity as a target value to check VIF
weekly_data$rarity = rarity$log_modulo_skewness

model1 = lm(rarity~., data = weekly_data)
vif(model1)
mean(vif(model1))
#if bigger than 4 there is concerning multi-collinearity in the data 


################################################################################
#remove parameters with missing data points
weekly_data =  subset(weekly_data, select = -c(Sludge_load_calc_WWTP_kg_COD_per_kg_DW_per_d,Inf_P_tot_kg_per_d,T_avg_C))

#weekly_data = weekly_data[-c(8,9,10)]

#Check for multi-collinearity problem
weekly_data_cor = cor(weekly_data, use = "pairwise.complete.obs")
weekly_data_testRes= cor.mtest(weekly_data,conf.level = 0.95)

#Condition number ratio of max to min Eigen values of the correlation matrix  
kappa(weekly_data_cor,exact=TRUE)

model1 = lm(rarity~., data = weekly_data)
vif(model1)
mean(vif(model1))
plot_missing(weekly_data, title = "missing data profile for selected parameters in weekly data")
dim(weekly_data)
write.csv(weekly_data, "output_data/relevant_process_data.csv", row.names=FALSE)

#EDA for selected parameters in weekly data
weekly_data %>%
  create_report(
    output_file = "weekly_data_final",
    output_dir = "output_data/" ,
    report_title = "EDA Report - weekly_data_final", 
    config = configure_report(add_plot_correlation = FALSE)) 

###############################################################################
###############################################################################
##Recursive Feature Elimination  

#rfe1 is fitted to Random Forest model with with 10 k-fold cross correlation  
#ref2 is fitted to bagged trees model with with leave one out cross correlation
Weekly_Data <- weekly_data

#12 input and rarity column 13 is output 
#normalization of the variable values and splitting of input and target values  
x <-Weekly_Data[,1:13]
normalization <- preProcess(x)
x <- predict(normalization, x)
x <- as.data.frame(x)
y<- Weekly_Data$rarity

#training scheme: setting up the controls for each recursive feature elimination 
control_RF_CV = rfeControl(functions=rfFuncs, method="cv", repeats = 5, number = 10, returnResamp = 'all')
control_TB_LOOCV = rfeControl(functions=treebagFuncs, method="LOOCV", returnResamp = 'all')

#reproducible data 
set.seed(11227430)

#split data- 80% for training and 20% for testing 
inTrain <- createDataPartition(Weekly_Data$rarity, p= .80, list = FALSE)[,1]

x_train <- x[ inTrain, ]
x_test <- x[-inTrain, ]

y_train <- y[ inTrain]
y_test <- y[ -inTrain]

#run RFE
results_rfe1 <- rfe(x =x_train , y= y_train , sizes=c(1:13),rfeControl=control_RF_CV)
sprintf("The optimal number of variables is: %s", results_rfe1$bestSubset)
sprintf("Optimal Variable: %s", results_rfe1$optVariables)

results_rfe2 <- rfe(x =x_train , y= y_train , sizes=c(1:13),rfeControl=control_TB_LOOCV)
sprintf("The optimal number of variables is: %s", results_rfe2$bestSubset)
sprintf("Optimal Variable: %s", results_rfe2$optVariables)
#plot <- ggplot(data = results_rfe1, metric = "RMSE")

png(filename="figures/Variables_RMSE.png", units ='in', height=8, width=8, res = 1000)
trellis.par.set(caretTheme())
plot1 <- plot(results_rfe1, type = c("g", "o"))
plot2 <- xyplot(results_rfe1, 
                type = c("g", "p"), 
                ylab = "RMSE CV Estimates")
print(plot1, split=c(1,1,1,2), more=TRUE)
print(plot2, split=c(1,2,1,2))
dev.off()

#plot for RMSE of rfe2
plot(results_rfe2, type = c("g", "o"))

#plot variable importance for random forest with cross validation 
png(filename="figures/Var_Imp.png", units ='in', height=5, width=5, res = 1000)
varimp_data <- data.frame(feature = row.names(varImp(results_rfe1))[1:2],
                          importance = varImp(results_rfe1)[1:2, 1])
ggplot(data = varimp_data, 
       aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.25, hjust=1.25,color="white", size=4) + 
  theme_bw() + theme(legend.position = "none") +
  coord_flip()
dev.off()

#plot variable importance for bagged trees with "leave-one-out cross-validation" 
png(filename="figures/Var_Imp_LOOCV.png", units ='in', height=4, width=4, res = 1000)
varimp_data <- data.frame(feature = row.names(varImp(results_rfe2))[1],
                          importance = varImp(results_rfe2)[1, 1])
ggplot(data = varimp_data, 
       aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.25, hjust=1.25,color="white", size=4) + 
  theme_bw() + theme(legend.position = "none") +
  coord_flip()
dev.off()

#plot the density distribution of the RMSE for number of variables included in model
png(filename="figures/Density.png", units ='in', height=5, width=8, res = 1000)
#plot5 = stripplot(results_rfe1)
#plot6 = histogram(results_rfe1)
plot7 = densityplot(results_rfe1)
print(plot7)
dev.off()

################################################################################ 
#Alternative approach

Weekly_Data<-as.data.frame(weekly_data)

rfe_fit<- rfeTerminator(Weekly_Data, x_cols= 1:12, y_cols=13, alter_df = TRUE, eval_funcs = rfFuncs)

#Explore the optimal model results
print(rfe_fit$rfe_model_fit_results)
# Explore the optimal variables selected
print(rfe_fit$rfe_model_fit_results$optVariables)
# Explore the original data passed to the frame
print(head(rfe_fit$rfe_original_data))
# Explore the data adapted with the less important features removed
print(head(rfe_fit$rfe_reduced_data))



