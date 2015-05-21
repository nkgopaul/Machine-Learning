learningcurve <- function(filename){
	data = read.csv(filename)
	plot(data[[1]], data[[2]], xlab = 'percentage of trainingset used', ylab = '% accuracy', main = 'Learning Curve')
}