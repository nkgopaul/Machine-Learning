import random, csv

def generate_a():
	c = csv.writer(open('prob1a.csv', 'wb'))
	header = []
	for i in range(100):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	for i in range(10000):
		row = []
		for j in range(100):
			pred_val = random.randint(0,1)
			row.append(pred_val)
		if (row[0] == 1 and row[2] == 1 and row[6] == 0) or (row[44] == 1 and row[66] == 0 and row[98] == 1):
			row.append(1)
		else:
			row.append(0)
		c.writerow(row)

def generate_b():
	c = csv.writer(open('prob1b.csv', 'wb'))
	header = []
	for i in range(100):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	prob1 = [1]*55 + [0]*45
	prob2 = [0]*55 + [1]*45
	for i in range(100):
		row = []
		y = random.randint(0,1)
		if y == 1:
			for j in range(100):
				row.append(random.choice(prob1))
		else:
			for j in range(100):
				row.append(random.choice(prob2))
		row.append(y)
		c.writerow(row)

def generate_c():
	c = csv.writer(open('prob1c.csv','wb'))
	header = []
	for i in range(100):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	for i in range(10000):
		row = []
		count = 0
		for j in range(100):
			pred_val = random.randint(0,1)
			if pred_val == 1:
				count+=1
			row.append(pred_val)
		if count%2 == 0:
			row.append(0)
		else:
			row.append(1)
		c.writerow(row)

def generate_d():
	c = csv.writer(open('prob1d.csv','wb'))
	header = []
	for i in range(100):
		header.append('X' + str(i+1))
	header.append('Y')
	c.writerow(header)
	for i in range(10000):
		row = []
		for j in range(100):
			row.append(random.randint(0,1))
		if (row[2] + row[4] + row[19] + row[20] - row[21] + row[87]*row[98] - row[6]*row[7]/(row[2] + 0.1)) > 0:
			row.append(1)
		else:
			row.append(0)
		c.writerow(row)

#generate_a()
#generate_b()
#generate_c()
#generate_d()