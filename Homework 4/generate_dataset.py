import csv, random
def generate_num_4():
	'''
	This function writes a dataset to a csv. It is built to perform very well on a classification task using J48, and not as well using k nearest neighbor classification.
	There is one outcome variable termed 'Y' and 50 predictor values termed X1, ... , X50 which are binary. The dataset is generated following two rules. If X1 = 1 AND X3 = 1 AND X7 = 0, Y must be 1. 
	If X4 = 1 AND X9 = 0 AND X22 = 1, Y must also be 1. Otherwise, Y is 0. These rules can be easily learned by J48, while these rules are almost impossible for a nearest neighbor algorithm
	to learn and thus will have some errors.
	'''
	c = csv.writer(open('prob4.csv', 'wb'))
	header = []
	for i in range(50):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	for i in range(1000):
		row = []
		for j in range(50):
			pred_val = random.randint(0,1)
			row.append(pred_val)
		if (row[0] == 1 and row[2] == 1 and row[6] == 0) or (row[3] == 1 and row[8] == 0 and row[21] == 1):
			row.append(1)
		else:
			row.append(0)
		c.writerow(row)

def generate_ec():
	'''
	This function writes a dataset to a csv. It is built to perform very well on a classification task using multilayer perceptron, and not as well using naive bayes classification.
	There is one outcome variable termed 'Y' and 2 predictor values termed X1 and X2 which are all binary. The dataset is generated following one rule, the XOR rule. If X1 = X2, Y = 0, otherwise Y = 1. 
	This rule can be easily learned by the multilayer perceptron because of its multiple layers, while these rules are almost impossible for a naive bayes algorithm to learn 
	and will have some errors. This is due to the naive bayes algorithm having the assumption that given the class, the predictors are independent of each other. This, in the XOR case,
	is false. Thus the naive bayes will perform poorly.
	'''
	c = csv.writer(open('prob_ec.csv', 'wb'))
	header = []
	for i in range(2):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	for i in range(1000):
		row = []
		for j in range(2):
			row.append(random.randint(0,1))
		if row[0] == row[1]:
			row.append(0)
		else:
			row.append(1)
		c.writerow(row)

generate_ec()
