# PLS and PLSDA simulation

# PLS data simulation ----------

train <- 1:40
test <- 41:50

jeu <- create.data(n=50)

X <- jeu$X[train,]
Y <- jeu$Y[train,]
X.test <- jeu$X[test,]
Y.test <- jeu$Y[test,]

result.pls <- PLS(X = X, Y = Y, ncomp = 10, mode = "regression")
Y.chap <- predict.PLS(result.pls, newdata = X.test)

head(Y.test)
head(Y.chap$predict)
sum((Y.test - Y.chap$predict[,,10])^2)
