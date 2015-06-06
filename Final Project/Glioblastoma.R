install.packages('R.matlab')
install.packages('stringr')
install.packages('cvTools')
install.packages('randomForest')
install.packages('ggplot2')
install.packages('miscTools')
install.packages('klaR')
install.packages('nnet')
install.packages('SDMTools')
library(R.matlab)
library(stringr)
library(cvTools)
library(randomForest)
library(ggplot2)
library(miscTools)
library(klaR)
library(nnet)
library(SDMTools)

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

full_experiment <- function(){
	csvs <- c("GBM_106.csv", "GBM_206.csv")
	outcomes <- c("cellularity.Median", "cell_voronoi_area.Median", "cytoplasm_background_intensity_mean.Median", "edge_length.Median", "nucleus_area.Median", "nucleus_background_intensity_mean.Median")
	corr_thresholds <- c(0.0001, 0.00001, 0.01, 0.0001, 0.01, 0.01)
	for (k in 1:length(csvs)){
		print(csvs[k])
		print('---------------------------------------------------')
		for (i in 1:length(outcomes)){
			cross_validation(csvs[k], outcomes[i], corr_thresholds[i], 10)
			print('----------------------------------------------------')
		}
		print('--------------------------------------------------')
		print('--------------------------------------------------')
	}
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

	nb_predictions <- vector()
	nn_predictions <- vector()
	class_rf_pred <- vector()
	class_cart_pred <- vector()
	test_outcome <- vector()

	for (i in 1:num_fold){
		#Preprocessing and Crossvalidation Steps
		#
		#
		#
		#test and training set splits
		test_set_indices <- folds[ ((i-1) * bin_size + 1):(i*bin_size)]		
		test_set <- data.frame(GBM_data[ test_set_indices, ])

		colnames(test_set) <- colnames(GBM_data)
		rownames(test_set) <- rownames(GBM_data)[test_set_indices]

		training_set_indices <- -1*test_set_indices
		training_set <- data.frame(GBM_data[training_set_indices, ])
		
		colnames(training_set) <- colnames(GBM_data)
		rownames(training_set) <- rownames(GBM_data)[training_set_indices]

		outcome_index <- match(outcome_name,colnames(GBM_data))
		

		#Regression Steps
		#Model 1: Linear Regression with Correlation Screening
		#Model 2: Random Forest of Regression Trees
		#Model 3: Cart Tree (Regression)


		predictor_list <- univariate_correlation_screening(training_set, outcome_index, corr_threshold)
		#print(length(predictor_list))
		lin_formula <- build_linear_formula(predictor_list, outcome_name)

		linear_model <- lm(formula = lin_formula, data = training_set)
		r_sq <- summary(linear_model)$r.squared
		#print(paste("fold", toString(i), " R^2: ", toString(r_sq), sep = ''))	
			

		predictions <- predict(linear_model, test_set)
		errors <- test_set[,outcome_index] - predictions
		cvss <- cvss + sum(errors^2)
		
		#rf<- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], xtest = test_set[,1:num_predictors], 
										#ytest = test_set[,outcome_index], importance = TRUE, ntree = 1000)
		rf <- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 2000)
		rf_pred <- predict(rf, test_set)
		r2 <- rSquared(test_set[,outcome_index], abs(test_set[,outcome_index] - rf_pred))
		#print(paste("R^2: ", toString(r2), sep=''))
		mse <- (test_set[,outcome_index] - rf_pred)^2
		rf_mse <- rf_mse + sum(mse)
		#print(paste("MSE: ", toString(mse), sep = ''))
		
		cart <- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 1)
		cart_pred <- predict(cart, test_set)
		c_mse <- (test_set[,outcome_index] - cart_pred)^2
		cart_mse <- cart_mse + sum(c_mse)
		#rfcv(trainx = training_set[,1:num_predictors], trainy = training_set[,outcome_index], cv.fold = num_fold, xtest = test_set[,1:num_predictors], ytest = test_set[,outcome_index])		

		#Classification Portion
		#Step 1: Discretize the Outcome
		#Model 1: Naive Bayes
		#Model 2: Neural Network
		#Model 3: Random Forest (Classification)
		#Model 4: Cart Tree (Classification)

		#calculates the thresholds
		thresholds <- quantile(training_set[,outcome_index], .75)
		train_disc_outcome <- vector()
		test_disc_outcome <- vector()
		
		#makes a new vector for train set with the discrete values
		for (i in 1:length(rownames(training_set))){
			
			attr_val <- training_set[i,outcome_index]
			
			if(attr_val > thresholds[1]){
				
				train_disc_outcome <- c(train_disc_outcome, 1)
				
			} else{
				train_disc_outcome <- c(train_disc_outcome, 0)
			}

		}

		#makes a new vector for test set with discrete values
		for (j in 1:length(rownames(test_set))){
			attr_val <- test_set[j,outcome_index]
			
			if(attr_val > thresholds[1]){
				
				test_disc_outcome <- c(test_disc_outcome, 1)
				
			} else{
				test_disc_outcome <- c(test_disc_outcome, 0)
			}
		}

		#put it into the data frame
		training_set[,outcome_index] <- as.factor(train_disc_outcome)
		test_set[,outcome_index] <- as.factor(test_disc_outcome)
		test_outcome <- c(test_outcome, as.numeric(test_set[,outcome_index])-1)
		#Naive Bayes
		nb <- NaiveBayes(x = training_set[,1:num_predictors], grouping = training_set[,outcome_index])
		nb_predictions <- c(nb_predictions,as.numeric(predict(nb, test_set, type = "class")$class)-1)
		#print(nb_predictions)
		#print(nb_probs)
		#nb_predictions <- vector()
		
		#for (i in 1:length(row(nb_probs))){
			#nb_predictions <- c(nb_predictions, which.max(nb_probs[i,])-1)
		#}

		#Actual values are columns
		#Predicted are rows
		#nb_auc <- auc(test_outcome, vector(nb_predictions))

		#Logistic Regression
		#predictor_list <- univariate_correlation_screening(training_set, outcome_index, corr_threshold)
		#print(length(predictor_list))
		#logit_formula <- build_linear_formula(predictor_list, outcome_name)
		#FIX THIS MODEL
		#logistic_model <- glm(formula = logit_formula, family = binomial, data = training_set)
		#predictions <- predict(logistic_model, test_set)


		#Neural Network
		nn_form <- as.formula(training_set[,outcome_index]~., training_set[,1:num_predictors])
		nn <- nnet(nn_form, training_set, size = 4, trace = FALSE)
		nn_predictions <- c(nn_predictions,as.numeric(predict(nn, test_set, type="class")))
		#print(nn_predictions)
		#nn_auc <- auc(test_outcome, vector(nn_predictions))

		#Random Forest
		class_rf<- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 2000)
		class_rf_pred <- c(class_rf_pred,as.numeric(predict(class_rf, test_set, type="class"))-1)
		#print(class_rf_pred)
		#class_rf_auc <- auc(test_outcome, as.vector(class_rf_pred))
		
		
		#Cart
		class_cart <- randomForest(x = training_set[,1:num_predictors], y = training_set[,outcome_index], importance = TRUE, ntree = 1)
		class_cart_pred <- c(class_cart_pred,as.numeric(predict(class_cart, test_set, type="class"))-1)
		



	}
	print(paste(outcome_name, ' Lin Reg MSE', ': ', toString(cvss/num_obs), sep = ''))
	print(paste(outcome_name, ' Random Forest MSE', ': ', toString(rf_mse/num_obs), sep = ''))
	print(paste(outcome_name, ' Cart MSE', ': ', toString(cart_mse/num_obs), sep = ''))
	print(paste(outcome_name, ' Naive Bayes AUC: ', toString(auc(test_outcome, nb_predictions))))
	print(paste(outcome_name, 'Neural Network AUC: ', toString(auc(test_outcome, nn_predictions))))
	print(paste(outcome_name, 'Random Forest AUC: ', toString(auc(test_outcome, class_rf_pred))))
	print(paste(outcome_name, 'Cart AUC: ', toString(auc(test_outcome, class_cart_pred))))

}

calc_auc <- function(predictions, test_outcomes){
	for (i in 1:length(test_outcomes)){
		predicted_val <- predictions[i]
		true_val <- test_set[i]
		two_by_two <- matrix(c(0,0,0,0),2,2)	
		if(true_val == 0 & predicted_val == 0){
			two_by_two[1,1] = two_by_two[1,1] + 1
		} else if (true_val == 0 & predicted_val == 1){
			two_by_two[2,1] = two_by_two[2,1] + 1
		} else if (true_val == 1 & predicted_val == 0){
			two_by_two[1,2] = two_by_two[1,2] + 1
		} else if (true_val == 1 & predicted_val == 1){
			two_by_two[2,2] = two_by_two[2,2] + 1
		} 
	}
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