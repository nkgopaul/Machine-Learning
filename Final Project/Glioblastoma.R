install.packages('R.matlab')
install.packages('stringr')
install.packages('cvTools')
install.packages('randomForest')
install.packages('ggplot2')
install.packages('miscTools')
library(R.matlab)
library(stringr)
library(cvTools)
library(randomForest)
library(ggplot2)
library(miscTools)

write_gbm_csv <- function(filename){
	
	#first read the raw matlab file in
	rawdata = readMat(filename)

	#acquire just the raw data, the numbers and put it into a data frame
	df = data.frame(rawdata$RawData)

	#assign the appropriate rows and columns to the data frame
	rownames(df) = unlist(rawdata[[2]])
	colnames(df) = unlist(rawdata[[3]])
	df = t(df)
	#write the file to a csv to work with
	write.csv(data.frame(df), file = paste(str_sub(filename, 1, -5), ".csv", sep = ''))
}

prepare_dataset <- function(feature_csv_file){
	
	#Function that preprocesses the data such that the patients are matched up and the separate outcomes and features are in one csv 
	

	#extracts features from feature file
	features = data.frame(read.csv(feature_csv_file, header = TRUE, row.names = 1))
	features$X <- NULL

	#extracts outcomes from outcome file
	outcomes = data.frame(read.csv("TCGA_GBM_PathologyImageData_n244_v2_ForR.csv", header = TRUE, row.names = 1))
	outcomes$X <- NULL

	#outcome names and new colnames produced
	outcome_names = c(colnames(outcomes)[1], colnames(outcomes)[5], colnames(outcomes)[9], colnames(outcomes)[41], colnames(outcomes)[45], colnames(outcomes)[53])
	predictors = colnames(features)
	column_names = c(predictors, outcome_names)

	#initialize empty dataframe
	df = as.data.frame(setNames(replicate(length(colnames(features)) + 6, numeric(0), simplify = FALSE), column_names))
	
	#acquire patient ids and which overlap
	feature_patient_ids = rownames(features)
	outcome_patient_ids = rownames(outcomes)
	common_patients = intersect(feature_patient_ids, outcome_patient_ids)
	common_row_names = character(0)
	common_patient_count = 0
	
	#loop through and make a dataset of just the patients that both the predictors and outcomes are available
	for (i in 1:length(feature_patient_ids)){
		feature_pid = feature_patient_ids[i]
		if(is.element(feature_pid, outcome_patient_ids) == TRUE){
			common_row_names = c(common_row_names, feature_pid)
			common_patient_count = common_patient_count + 1
			observation = c(features[feature_pid,], outcomes[feature_pid,outcome_names])
			df[common_patient_count,] <- observation
		}
	}

	#write it to a file
	rownames(df) = common_row_names
	write.csv(df, file = paste("GBM_", toString(length(colnames(features)) + 6), ".csv", sep = ''))
}

cross_validation <- function(processed_csv, outcome_name, corr_threshold, num_fold){
	#Edge Length Median, Threshold = 0.0001
	#Cellularity Median, Threshold = 0.0001
	#Cell_Voronoi_Area Median, Threshold = 0.00001
	#Nucleus Area, Threshold = 0.01


	#sets seed for reproducibility
	set.seed(92)
	
	#reads data in
	GBM_data <- data.frame(read.csv(processed_csv, header = TRUE, row.names = 1))
	GBM_data$X <- NULL
	
	#calculates basic bins
	num_obs <- length(rownames(GBM_data))
	bin_size <- num_obs/num_fold
	num_predictors <- length(colnames(GBM_data)) - 6

	#makes the folds
	folds <- cvFolds(num_obs, K = num_fold, type = "random")$subset
	cvss <- 0
	#10 fold cross validation
	rf_mse <- 0
	cart_mse <- 0

	for (i in 1:num_fold){
		#test and training set splits
		test_set_indices <- folds[ ((i-1) * bin_size + 1):(i*bin_size)]
		
		test_set <- data.frame(GBM_data[ test_set_indices, ])
		#print(dim(test_set))
		colnames(test_set) <- colnames(GBM_data)
		rownames(test_set) <- rownames(GBM_data)[test_set_indices]

		training_set_indices <- -1*test_set_indices

		training_set <- data.frame(GBM_data[training_set_indices, ])
		
		colnames(training_set) <- colnames(GBM_data)
		rownames(training_set) <- rownames(GBM_data)[training_set_indices]
		#print(length(rownames(training_set)))
		#for(j in (num_predictors+1):(num_predictors+6)){
		#outcome_name <- colnames(GBM_data)[j]
		outcome_index <- match(outcome_name,colnames(GBM_data))
		predictor_list <- univariate_correlation_screening(training_set, outcome_index, corr_threshold)
		print(length(predictor_list))
		lin_formula <- build_linear_formula(predictor_list, outcome_name)
		#print(lin_formula)
		linear_model <- lm(formula = lin_formula, data = training_set)
		#print(linear_model)
		#print(summary(linear_model))
		r_sq <- summary(linear_model)$r.squared
		print(paste("fold", toString(i), " R^2: ", toString(r_sq), sep = ''))	
			#sig_predictors <- extract_significant_predictors(linear_model)
			#sig_formula <- build_linear_formula(sig_predictors, outcome_name)
			#sig_lm <- lm(formula = sig_formula, data = GBM_data)
			
			#predictions <- predict(sig_lm, test_set)
		predictions <- predict(linear_model, test_set)
		#predictions <- vector()
		#for(j in 1:length(rownames(test_set))){
			#a <- summary(linear_model)$coefficients[,1]
			#b <- data.frame(c(1.0, test_set[j,1:num_predictors-1]))
			#print(str(a))
			#print(str(b))
			#predictions <- c(predictions, sum(a*b))
			#}
		#print(length(predictions)) 
		errors <- test_set[,outcome_index] - predictions
		#print(errors)
		cvss <- cvss + sum(errors^2)
		
		#rf<- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], xtest = test_set[,1:num_predictors], 
										#ytest = test_set[,outcome_index], importance = TRUE, ntree = 1000)
		rf<- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 1000)
		rf_pred = predict(rf, test_set)
		r2 <- rSquared(test_set[,outcome_index], abs(test_set[,outcome_index] - rf_pred))
		print(paste("R^2: ", toString(r2), sep=''))
		mse <- (test_set[,outcome_index] - rf_pred)^2
		rf_mse <- rf_mse + sum(mse)
		#print(paste("MSE: ", toString(mse), sep = ''))
		
		cart <- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 1)
		cart_pred <- predict(cart, test_set)
		c_mse <- (test_set[,outcome_index] - cart_pred)^2
		cart_mse <- cart_mse + sum(c_mse)
		#rfcv(trainx = training_set[,1:num_predictors], trainy = training_set[,outcome_index], cv.fold = num_fold, xtest = test_set[,1:num_predictors], ytest = test_set[,outcome_index])		
	}
	print(paste(outcome_name, ' Lin Reg MSE', ': ', toString(cvss/num_obs), sep = ''))
	print(paste(outcome_name, ' Random Forest MSE', ': ', toString(rf_mse/num_obs), sep = ''))
	print(paste(outcome_name, ' Cart MSE', ': ', toString(cart_mse/num_obs), sep = ''))
}

univariate_correlation_screening <- function(training_set, outcome_index, corr_threshold){
	#training_set is a data frame
	num_obs = length(rownames(training_set))
	num_predictors = length(colnames(training_set))-6
	predictor_list = character(0)
	
	for (i in 1:num_predictors){
		spearman_corr_pval <- cor.test(training_set[,i], training_set[,outcome_index], method = "spearman")$p.value
		
		if(spearman_corr_pval < corr_threshold){
			predictor_list = c(predictor_list, colnames(training_set)[[i]])
		}
	}

	return(predictor_list)
}

#function that takes in a predictor list and an outcome and produces the multivariate linear form used in the linear regression models
build_linear_formula <- function(predictor_list, outcome_var){
	str = predictor_list[1]
	for(i in 2:length(predictor_list)){
		str = paste(str, "+", predictor_list[i])
	}
	linear_formula = as.formula(paste(outcome_var, "~", str))
	return (linear_formula)
}

#helper function that extracts the significant predictors from a model at alpha < 0.10
extract_significant_predictors <- function(model){
	finalized_coef = (summary(model))$coefficients
	
	numb_rows = length(finalized_coef)/4
	matrix_coef = matrix(finalized_coef, ncol = 4, nrow = numb_rows)
	significant_predictors = NULL
	for (i in 2:numb_rows) {
		if(matrix_coef[i,4] < 0.1){
			predictor = rownames(finalized_coef)[i]
			significant_predictors = c(significant_predictors, substr(predictor,1,nchar(predictor)))
		}
	}
	significant_predictors = subset(significant_predictors, duplicated(significant_predictors) != TRUE)
	return (significant_predictors)

}