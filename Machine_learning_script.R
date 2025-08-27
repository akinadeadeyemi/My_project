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
### loading data
library(randomForest)
rf_model <- randomForest(Species ~ ., data = train_data)

#### make predictions
predictions <- predict(rf_model, newdata = test_data)

table(Predicted = predictions, Actual = test_data$Species)


### Working on mtcars (using LASSO regression)
data(mtcars)
library(tibble)

data(mtcars) 

X <- as.matrix(mtcars[, -1]) 

y <- mtcars[, 1]

### Fitting the Lasso regression model
##glmnet(): Fits a regularized linear model; family = "gaussian": Specifies linear regression.
## alpha = 1: Sets the model type to Lasso regression.
## summary(): Displays a summary of the fitted model.

mtcar_model <- glmnet(X, y, family = "gaussian", alpha = 1)
summary(mtcar_model)
plot(mtcar_model, label = TRUE)

## Getting model coefficients
## We extract the coefficients at a specific value of lambda.
## coef(): Retrieves coefficients from the fitted model.
## s: Specifies the lambda value at which to extract them.

coef(mtcar_model, , s = 0.1)

## Making predictions with the model
y_pred <- predict(mtcar_model, X)


### Fitting a Lasso model with cross-validation
fit <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
summary(fit)

plot(fit)
