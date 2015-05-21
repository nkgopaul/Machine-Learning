import csv, numpy, math, copy, random

def read_csv_data(s):
	'''
	Function that takes a string of a csv file name (with its .csv extension) and imports it into a list of lists
	where each list in the list of lists is a row in the csv file

	Input: CSV File Name
	Output: List of Lists Training Data Object
	'''
	ftrain = open(s)
	training_set = csv.reader(ftrain)
	training_data = []
	for row in training_set:
		training_data.append(row)
	return training_data

def impute_median(data): #pass in the entire dataset
	'''
	Function that takes in one list of lists object named data where the first list, data[0] is a header list that
	indicates the name of each of the attributes. This function imputes the medians of each of the predictors
	wherever a missing value exists. It returns both the header list and the data list of lists with all missing
	values resolved.

	Input: Raw Data with Header of type List of Lists that has missing values
	Output: Header List, Imputed Data without Header of type List of Lists without missing values
	'''
	header = data[0]
	col_data = numpy.transpose(data)
	types = type_of_attr(data)
	medians=[]
	
	#col_data_copy = copy.deepcopy(col_data)
	col_data_copy = []

	#first remove ?s
	for col in range(len(col_data)):
		new_col = []
		for ex in range(1,len(col_data[col])):
			if col_data[col][ex] != '?':
				col_data[col][ex] = numpy.float(col_data[col][ex])
				new_col.append(numpy.float(col_data[col][ex]))
		col_data_copy.append(new_col)

	#Then calculate the median
	for j in range(len(col_data_copy)):
		medians.append(numpy.float(numpy.median(col_data_copy[j])))

	#Then go through original data and impute
	for i in range(len(col_data)):
		for jj in range(1,len(col_data[i])):
			if col_data[i][jj] == '?' and types[i] == 'numeric':
				col_data[i][jj] = numpy.float(medians[i])
			elif col_data[i][jj] == '?' and types[i] == 'nominal':
				col_data[i][jj] = numpy.float(math.floor(medians[i]))
			else:
				col_data[i][jj] = numpy.float(col_data[i][jj])
	
	col_data = numpy.transpose(col_data)
	col_data_header = col_data[0]
	#convert numpy array to float
	col_data_rest = col_data[1:].astype(numpy.float)

	return header, col_data_rest

def preprocess_dataset(s):
	'''
	Function that takes a csv filename as input and imports it into python using read_csv_data and removes missing values
	using the impute_median function.
	'''
	return impute_median(read_csv_data(s))

def type_of_attr(data): #read in entire dataset
	'''
	Function that takes in a raw dataset and outputs a list of the types of data (nominal or numeric) of each of the predictors
	Determines whether data is nominal or numeric by removing all duplicates and seeing the length of the remaining list
	(A numeric predictor's list without duplicates will be almost just as long as the original length, while a nominal predictor's
	list will drastically reduce in size

	Input: Raw Dataset List of Lists with Header
	Output: List of strings signifying whether a column is a numeric or nominal column
	'''

	#transpose the data so all of the values of a predictor will be in the same row
	col_data = numpy.transpose(data)

	#the number of training examples excluding the output variable
	num_training_ex = len(col_data[-1])

	new_col_data = []
	types = []

	#remove duplicates with set function
	for i in range(len(col_data)):
		new_col_data.append(list(set(col_data[i])))

	#compare length of list without duplicates with original list
	for j in range(len(new_col_data)):
		if len(new_col_data[j]) < math.sqrt(num_training_ex):
			types.append('nominal')
		else:
			types.append('numeric')
	return types

def calc_outcome_distribution(dataset):
	'''
	Function that takes in a list of lists dataset without header, and assuming the outcome variable is in the last column and is binary, 
	returns the distribution of the outcome variable in a tuple as counts.

	Input: List of lists dataset without Header
	Output: Tuple of the distribution of the outcome variable
	'''
	distr = [0,0]
	for example in dataset:
		if(example[-1] == 0):
			distr[0]+=1
		else:
			distr[1]+=1
	return distr

def entropy(curr_distr):
	'''
	Function that calculates the entropy of the current distribution of the binary outcome variable and returns the entropy value.
	Entropy is based on a formula that is the sum over all values the outcome variable can take up of -1 * the prior probability
	of a specific value of the outcome variable * log base 2 of the specific prior probability of a specific value of the outcome
	variable. Log(0) is approximated to be 0 as the best entropy value is 0.

	Input: Tuple signifying the current binary distribution of the outcome variable
	Output: A positive float value signifying the entropy of that distribution
	'''
	denom = curr_distr[0] + curr_distr[1]
	priors = [float(curr_distr[0])/denom, float(curr_distr[1])/denom]
	ent = 0
	for prior in priors:
		if prior != 0:
			ent -= prior * math.log(prior) / math.log(2)
	return abs(ent)

def calc_split_val(sub_train_data):
	'''
	Function that calculates the value at which to split for each predictor given a subset of the dataset.

	Input: Subset of the dataset
	Output: List of values at which to split based on the median assumption for convenience and computational efficiency
	'''
	#transpose the data so all of the values of a predictor will be in the same row
	col_data = numpy.transpose(sub_train_data)

	noise_term = 0.00001

	#array of the median split value that we actually split on
	keys = []
	for attr in col_data:
		attr = numpy.median(list(attr)) + noise_term
		keys.append(attr)
	return keys

def split(sub_train_data, attr_idx):
	'''
	Function that takes in a subset of the dataset and a specific attribute index and outputs the result of a split along that

	Input: Subset of the Dataset, index of a specific attribute (int) based on the header
	Output: Subset of the data that falls on the left of the split value, subset of the data that falls on the right, the split value
	'''
	keys = calc_split_val(sub_train_data)

	#get the split value of the current index
	split_val = keys[attr_idx]
	
	left_data = []
	right_data = []

	#split data into left and right based on whether its value is less than or greater than the chosen split value
	for ex in sub_train_data:
		if(ex[attr_idx] < split_val):
			left_data.append(ex)
		else:
			right_data.append(ex)

	return left_data, right_data, split_val

def info_gain(sub_train_data, attr_idx):
	'''
	Function calculating the information gain at a specific point given a subset of the data, and a split given 
	an attribute index. Information gain is based on the difference in entropy before and after the split. A larger change is better.

	Input: Subset of the training data, and an index of attribute (int) in the header list
	Output: An entropy difference, left subset of the data, right subset of the data, split value, and outcome.

	Outcome is None unless it is a perfect split. This specific situation has a 0 start entropy is created, and we have a leaf node.
	'''

	outcome = None
	start_e = entropy(calc_outcome_distribution(sub_train_data))

	#if the starting entropy is 0, we have a perfect split
	if start_e == 0:
		outcome_distr = calc_outcome_distribution(sub_train_data)
		#print outcome_distr
		if outcome_distr[0] == 0:
			outcome = 1
		else:
			outcome = 0
		return False, False, False, False, outcome

	#retrieve all the data less than the split_val (left_data), greater than or equal to the split_val (right_data), and the split value itself
	left_data,right_data,split_val= split(sub_train_data, attr_idx)

	left_data_size = len(left_data)
	right_data_size = len(right_data)
	total_data_size = left_data_size + right_data_size
	
	#if total data size is 0, the size of both left and right data is 0
	#return just the starting entropy
	if total_data_size == 0:
		return start_e, False, False, False, None

	#calculate the weights of the distribution of data
	left_wt = float(left_data_size)/total_data_size
	right_wt = float(right_data_size)/total_data_size
	
	if left_data_size == 0:
		left_e = 0
		right_e = entropy(calc_outcome_distribution(right_data))
	elif right_data_size == 0:
		left_e = entropy(calc_outcome_distribution(left_data))
		right_e = 0
	else:
		left_e = entropy(calc_outcome_distribution(left_data))
		right_e = entropy(calc_outcome_distribution(right_data))

	#calculate the ending entropy of the split
	end_e = left_wt * left_e + right_wt * right_e

	#return False for left_data or right_data if left_data or right_data is empty, respectively
	if left_data_size == 0:
		return start_e-end_e, False, right_data, split_val, None
	elif right_data_size == 0:
		return start_e-end_e, left_data, False, split_val, None
	else:
		return start_e-end_e, left_data, right_data, split_val, None

def make_tree(data,node_index, dict_tree, header):
	'''
	Function that makes a binary tree stored as a dictionary given some data, a start node_index of 1, an empty dictionary tree, and a header.

	Input: Training Data as a list of lists without a header, node_index initialized to 1 as the first node should have a key value of 1,
		   an empty dictionary signifying an empty tree, and a header list signifying the identities of each of the variables.

	Output: The populated dictionary of the tree with key k as the node of the tree with its left child signified with key 2*k and right
			child signified with key 2*k + 1. The values of the tree are lists. For leaf nodes, those without a split value and perfect
			outcome distributions, the values are [outcome, None, len(data)], where outcome is the decision reached by a node when it
			reaches this node. The second element is a placeholder. The third element signifies the number of training examples that
			reached this node. For non-leaf nodes, the values are [best attribute string, best split value for that attribute, 
			size of training set that reached this node, a tuple signifying the outcome distribution at this point]
	'''

	#variables that will keep track of the best attribute to split on so far and its characteristics
	max_info_gain = 0
	best_attr = None
	best_split_val = None
	chosen_left_data = None
	chosen_right_data = None
	
	#for each attribute
	for attr_idx in range(len(header)-1):

		#return the info gain, left data, right data, split value, and outcome of splitting on this specific attribute
		curr_info_gain, left_data, right_data, split_val, outcome = info_gain(data, attr_idx)

		#if there is no outcome, we are not at a leaf node yet
		if outcome != None:
			dict_tree[node_index] = [outcome, None,len(data)]

		#if there is no update in info_gain, skip checking the rest of this attribute
		if curr_info_gain == False:
			continue

		#if the info gain exists and its greater than the current maximum info gain, update the max variables
		if curr_info_gain > 0 and curr_info_gain > max_info_gain:
			max_info_gain = curr_info_gain
			best_attr = attr_idx

			chosen_left_data = left_data
			chosen_right_data = right_data
			best_split_val = split_val
	
	#add the best attribute to the tree
	if best_attr != None:
		#Name of the attribute, value that we are spliiting on, length of subdata, and outcome distribution
		base_node = [header[best_attr], best_split_val, len(data), calc_outcome_distribution(data)]

		dict_tree[node_index]=base_node
		#repeat on node's children
		if(chosen_left_data != False):
			make_tree(chosen_left_data,2*node_index, dict_tree, header)
		if(chosen_right_data != False):
			make_tree(chosen_right_data,2*node_index+1, dict_tree, header)

		return dict_tree

def dnf(dict_tree):
	'''
	Function that prints a tree in disjunctive normal form.

	Input: Tree represented as a dictionary
	Output: A big string of DNF of the tree with < or > operators and ^ and v operators.
	'''
	paths = []
	curr_path = []
	stack = []
	leaf_ctr = 0
	index = 1
	dnf_output = ""
	stack.append(index)
	#going through stack
	while len(stack) > 0:
		index = stack.pop(len(stack)-1)
		if curr_path.count(index) == 0:
			curr_path.append(index)
		#add the next node to the path
		#check if leaf node
		if index not in dict_tree:
			continue
		if dict_tree[index] == None:
			continue
		if len(dict_tree[index]) == 3:
			if dict_tree[index][0] == 1:
				paths.append(curr_path) 
				leaf_ctr+=1
				if leaf_ctr >= 15: #16 condition
					break
				else:
					curr_path.remove(index)
		if (index*2+1) in dict_tree and dict_tree[index*2+1] != None: #depth first search right child
			stack.append(index * 2  + 1)
		if (index*2) in dict_tree and dict_tree[index*2] != None: #depth first search left child
	 		stack.append(index * 2)

	#prints the output using these paths and the proper operators
	for rows in paths:
		dnf_output += "("
		for j in rows:

			if len(dict_tree[j]) == 4:
				
				dnf_output += str(dict_tree[j][0])
				if j%2:
					dnf_output += " < "
					dnf_output += str(float(dict_tree[j][1]))

				else:
					dnf_output += " > "
					dnf_output += str(float(dict_tree[j][1]))
				dnf_output += " ^ "
		dnf_output += ") v "


	dnf_output=dnf_output.replace('^ )',')')
	if dnf_output[-2:]=='v ':
		dnf_output=dnf_output[:len(dnf_output)-2]

	return dnf_output

def calc_num_splits(tree):
	'''
	Function that calculates the number of splits that a tree has
	'''
	splits = 0
	for key in tree:
		if key not in tree:
			continue
		elif tree[key] == None:
			continue
		elif len(tree[key]) == 4:
			splits += 1
	return splits

def decision_tree(trainfile_name, validfile_name, testfile_name):
	'''
	Main function that takes in strings of the trainfile, validation file, and test file if they are in the same folder/directory in csv format and builds two trees,
	one with pruning and one without. It then prints the sizes of the two trees, and the DNF's of the two trees. Next it prints the validation percentages of the trees
	based on the validation set provided. It then imputes predictions for both of the trees into csv files named unpruned_test_set.csv and pruned_test_set.csv. It then
	calculates data points for a learning curve for each of the types, an unpruned learning curve and a pruned learning curve based on these strategies. The data points
	are put into a csv file named noprunelc.csv and prunelc.csv. These will be plotted and shown using another function, the learningcurve.R function.
	'''

	#first import the data from the files
	header, training_data = preprocess_dataset('btrain.csv')
	header_two, validation_data = preprocess_dataset('bvalidate.csv')
	header_three, test_data = preprocess_dataset('btest.csv')
	dict_tree={}
	
	#make the full tree based on the training data
	full_tree = make_tree(training_data, 1, dict_tree, header)
	full_tree_size = calc_tree_size(full_tree)
	full_tree_splits = calc_num_splits(full_tree)
	#prune that tree such that the best tree is chosen based on our pruning strategy
	pruned_tree = size_pruning(full_tree, len(training_data), validation_data, header_two)
	pruned_tree_size = calc_tree_size(pruned_tree)
	pruned_tree_splits = calc_num_splits(pruned_tree)

	print "Size of Full Tree: " + str(full_tree_size)
	print "Size of Pruned Tree: " + str(pruned_tree_size)

	print "Full Tree Splits: " + str(full_tree_splits)
	print "Pruned Tree Splits: " + str(pruned_tree_splits)

	print "DNF of Full Tree: " + str(dnf(full_tree))
	print "DNF of Pruned Tree: " + str(dnf(pruned_tree))

	#calculate the performance of each tree on the validation set
	best_pruned_perc = validation_percent(validation_data, pruned_tree, header_two)
	unpruned_perc = validation_percent(validation_data, full_tree, header_two)
	
	print "Validation Percentage of Full Tree: " + str(unpruned_perc)
	print "Validation Percentage of Pruned Tree: " + str(best_pruned_perc)
	
	#impute the predicted values for the test set for each of the trees
	pruned_test_data = impute_predictions(test_data, pruned_tree, header_three)
	unpruned_test_data = impute_predictions(test_data,full_tree,header_three)

	#write that data to csv files
	print "Inputing Predictions to Test Set using Full Tree"
	write_test_to_csv('unpruned_test_set.csv', unpruned_test_data, header_three)
	print "Imputation Finished"

	print "Imputing Predictions to Test Set using Pruned Tree"
	write_test_to_csv('pruned_test_set.csv', pruned_test_data, header_three)
	print "Imputation Finished"

	#calculate the learning curves for each type of tree, full and pruned
	#put the data associated with the learning curves into a csv file to be plotted and analyzed in R
	print "Learning Curve Data for Full Trees is being Generated, and put into a csv"
	plot_learning_curve(training_data, validation_data, header, False, 'noprunelc.csv')
	print "Learning Curve Data for Full Tree Finished and put into csv named 'noprunelc.csv'"
	
	print "Learning Curve Data for Pruned Trees is being Generated, and put into a csv"
	print plot_learning_curve(training_data, validation_data, header, True, 'prunelc.csv')
	print "Learning Curve Data for Pruned Trees Finished and put into csv named 'prunelc.csv'"

def write_test_to_csv(filename, data, header):
	'''
	Function that writes the imputed test set into a csv file
	'''

	full_data = copy.deepcopy(data)
	c = csv.writer(open(filename, 'wb'))
	c.writerow(header)
	for row in full_data:
		c.writerow(row)

def write_lc_to_csv(filename, data):
	'''
	Function that writes the learning curve data into a csv file
	'''
	lc_data = copy.deepcopy(data)
	c = csv.writer(open(filename, 'wb'))
	for row in lc_data:
		c.writerow(row)

def calc_tree_size(tree):
	'''
	Simple function that takes in a tree as a dictionary and counts the number of nodes that have not been pruned
	'''
	size = 0
	if len(tree) == 0:
		return 0

	for key in tree:
		if tree[key] != None:
			size += 1

	return size

def validation_percent(valid_set, dict_tree, header):
	'''
	Function that takes in a validation set, tree, and a header and calculates the percent accuracy on that set.
	'''
	curr_attr = None
	curr_split_val = 0
	header_row = header
	curr_attr_idx = -1
	correct = 0
	wrong = 0
	cond_one = 0
	cond_two = 0
	for row in valid_set: #loops through each validation example

		decision_reached = False #decision marker
		key = 1 #starting node of the tree
		tree_path = [] #empty tree path

		while(not decision_reached): #loops until a decision is reached
			if key not in dict_tree or dict_tree[key] == None: #checks if this type of training example has not been seen or if it has been pruned
				cond_one+=1
				new_distr = dict_tree[tree_path[-1]][3] #if so don't append this new node and stay at the old one, and calculate its outcome distribution
				decision = new_distr.index(max(new_distr)) #chose the more abundant outcome value and use that as the decision
				decision_reached = True #mark decision reached

			elif len(dict_tree[key]) == 3: #checks if node is a leaf node
				cond_two+=1
				tree_path.append(key) #add it to path
				decision = dict_tree[key][0] #acquire the decision
				decision_reached = True #mark decision reached
				
			else: #nonterminal node reached
				tree_path.append(key) #add to path
				curr_attr = dict_tree[key][0] #find the attribute to split on here
				curr_split_val = dict_tree[key][1] #find its split value
				attr_idx = header_row.index(curr_attr)
				if(row[attr_idx] < curr_split_val): #go down the appropriate path
					key = 2*key
				else:
					key = 2*key + 1
		
		if int(row[-1]) == int(decision): #keeps track of number correct
			correct += 1
		if int(row[-1]) != int(decision): #keeps track of number wrong
			wrong += 1
	
	return float(correct) / len(valid_set) #returns the percent accuracy

def size_pruning(dict_tree, training_data_size, valid_set, header):
	'''
	Function that prunes a tree given a tree, size of training data, a validation set, and a header. The pruning strategy is based on the number of training examples. First we make a list of thresholds
	that are incremented by one from the 4th root of the training set size to the square root of the training set size. We then prune all nodes that contain less examples than the threshold in the tree.
	We then check the accuracy of each of these trees on the validation set. The tree with the best percent accuracy is chosen as the final tree and is returned.
	'''

	root_size_end = int(math.floor(math.sqrt(training_data_size)))
	root_size_beginning = int(math.floor(math.sqrt(root_size_end)))

	thresholds = []
	#thresholds starts at the 4th root of the training data size and increments by 1 until the square root of the training set size
	for i in range(root_size_beginning, root_size_end, 1):
		thresholds.append(i)

	#baseline accuracy of whole tree
	baseline_acc = validation_percent(valid_set, dict_tree, header)

	best_threshold = 0
	best_tree = dict_tree
	
	for threshold in thresholds: #loop through each threshold
		pruned_tree = copy.deepcopy(dict_tree)
		for key in pruned_tree: #check each node in the tree
			if pruned_tree[key][2] < threshold: #prune if smaller than threshold
				pruned_tree[key] = None

		val_perc = validation_percent(valid_set, pruned_tree, header) #calculate the validation percentage of this tree

		if(val_perc > baseline_acc): #checks for best valiadation percentage
			baseline_acc = val_perc
			best_threshold = threshold
			best_tree = copy.deepcopy(pruned_tree)
	return best_tree #returns the best tree

def data_subset(data_set, percent):
	'''
	Function that selects a subset of a data set based on a percentage needed
	'''
	elements = int(round(percent * len(data_set)))
	return random.sample(data_set, elements)

def impute_predictions(test_data, dict_tree, header):
	'''
	Function that takes in an unlabeled set of test examples and uses an already made decision tree to make decisions for them. It imputes it into the test set.
	'''
	curr_attr = None
	curr_split_val = 0
	header_row = header
	curr_attr_idx = -1
	new_test_data = []
	for row in test_data: #loop through each instance and calculates a decision the same way validation percentage does
		decision_reached = False
		key = 1
		tree_path = []
		while(not decision_reached):
			if key not in dict_tree or dict_tree[key] == None:
				new_distr = dict_tree[tree_path[-1]][3]
				decision = new_distr.index(max(new_distr))
				decision_reached = True
			
				#decision = random.randint(0,1)
				#decision_reached = True;

			elif len(dict_tree[key]) == 3:
				tree_path.append(key)
				decision = dict_tree[key][0]
				decision_reached = True
				
			else:
				tree_path.append(key)
				curr_attr = dict_tree[key][0]
				curr_split_val = dict_tree[key][1]
				attr_idx = header_row.index(curr_attr)
				if(row[attr_idx] < curr_split_val):
					key = 2*key
				else:
					key = 2*key + 1
		
		row[-1] = int(decision) #adds decision in
		new_test_data.append(row) #adds this new example
	return new_test_data #returns a new testset with predictions
			
def plot_learning_curve(training_data, valid_data, header, do_pruning, filename):
	'''
	Function that calculates the learning curve and outputs a list of tuples to be plotted by a function in R.
	This function takes in training and validation data, as well as a header, a boolean indicating if pruning should be done or not, and a filename to store the list of tuples.
	'''
	data_pts = []
	for i in range(2, 21): #Loops through 2-20
		accuracy = []
		perc = float(i)/20 #percentages 0.1 -> 1.0
		num_iter = int(8/i) + 1 #the smaller percentages will be small enough to do more analysis on by averaging values
		for j in range(num_iter): #loops through this num_iter which is larger for smaller percentages
			tree = {}
			pruned_tree = {}
			subset_of_data = data_subset(training_data, perc) #calculates the subset of the data based on this percentage
			tree = make_tree(subset_of_data, 1, {}, header) #makes a full tree with this subset
			
			if do_pruning: #if pruning
				pruned_tree = size_pruning(tree, len(subset_of_data), valid_data, header) #prune the tree and get the best pruned tree
				accuracy.append(validation_percent(valid_data, pruned_tree, header)) #append the validation percent accuracy of the models
			else: #if not pruning			
				accuracy.append(validation_percent(valid_data, tree, header)) #calculate and append the full tree validation percentage
		pt = [perc, sum(accuracy)/len(accuracy)] #average the percent accuracy depending on the number of times its been done
		data_pts.append(pt) #append this datapoint
	write_lc_to_csv(filename, data_pts) #write all this data to a csv
	return data_pts #return the data_pts list as a formality

decision_tree('btrain.csv', 'bvalidate.csv', 'btest.csv')