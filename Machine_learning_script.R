######### Author : Adeyemi Akinade
######### Place : Clemson University, SC
######### lab : Gopalan's lab
######### Date : August 27, 2025
######### Project : Machine Learning with R


######################## LEARNING MACHINE LEARNING WITH R ##########################################

data(iris)
head(iris)

### to allow reproducibility
set.seed(123)
library(caret)
library(ggplot2)
library(lattice)


### divide data into 70% training and 30% testing 
iris_train_index <- createDataPartition(iris$Species, p = 0.7, list =  F)
iris_train <- iris[iris_train_index, ]
iris_test <- iris[-iris_train_index, ]


### Another way to divide samples
set.seed(1)
n <- nrow(iris)           # total number of rows
train_size <- floor (0.7 * n)  # 70% for training
train_indices <- sample(n, size = train_size)

train_data <- iris[train_indices,]
test_data <- iris[-train_indices,]


### using Random-Forest for classification 










