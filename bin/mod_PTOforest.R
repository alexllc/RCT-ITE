## MODIFIED FROM POWERS

cv.PTOforest <- function(x, tx, y, 
                        validation_fold = NULL, 
                        pscore = rep(.5, nrow(x)), 
                        num.trees = NULL, 
                        mtry = NULL, 
                        min.node.size = NULL,
                        postprocess = TRUE, verbose = FALSE) {

    # Input sanitization

    x = as.matrix(x)

    if (nrow(x) != length(tx)) {
        stop('nrow(x) does not match length(tx)')

    } else if (nrow(x) != length(y)) {
        stop('nrow(x) does not match length(y)')

    } else if (!is.numeric(x)) {
        stop('x must be numeric matrix')

    } else if (!is.numeric(y)) {
        stop('y must be numeric (use 0/1 for binary response)')

    } else if (!is.numeric(tx) | length(setdiff(tx, 0:1)) > 0) {
        stop('tx must be vector of 0s and 1s')
    }

    valID <- sample(1:dim(x)[1], floor(dim(x)[1] / validation_fold))
    train_x <- x[!(1:dim(x)[1] %in% valID), ]
    train_y <- y[!(1:dim(x)[1] %in% valID) ]
    train_tx <- tx[!(1:dim(x)[1] %in% valID) ]
    train_pscore <- pscore[!(1:dim(x)[1] %in% valID) ]

    ho_x <- x[valID, ]
    ho_y <- y[valID]
    ho_tx <- tx[valID]
    ho_pscore <- pscore[valID]

    colnames(x) = paste('x', 1:ncol(x), sep = '')
    fit = list(x = train_x, pscore = train_pscore, postprocess = postprocess)

    train_z = train_tx * train_y / train_pscore - (1 - train_tx) * train_y / (1 - train_pscore)


    if (verbose) cat('fitting IPW treatment forest\n')

    train_data = data.frame(y = train_z, x = train_x)
    colnames(train_data) = c('y', colnames(train_x))
    TOfit = ranger::ranger(data = train_data, 
                    dependent.variable.name = 'y',
                    num.trees = num.trees, 
                    min.node.size = min.node.size, 
                    mtry = mtry,
                    write.forest = TRUE)

    cv_z_pred <- stats::predict(TOfit, ho_x)[["predictions"]]
    ho_z <- ho_tx * ho_y / ho_pscore - (1 - ho_tx) * ho_y / (1 - ho_pscore)

    cv_error <- mean((cv_z_pred - ho_z)^2)
  
    return(cv_error)
}

#'  Pollinate a fitted ranger random forest model
#'
#' @param object a fitted \code{ranger} object
#' @param x matrix of covariates
#' @param y vector of response values
#'
#' @return an object of class \code{pollinated.ranger} which is a \code{ranger}
#'  object that has been pollinated with the data in (x, y)

pollinated.ranger = function(object, x, y) {

  forest = object$forest
  num.trees = forest$num.trees
  which.list = as.list(seq(num.trees))
  split.values = forest$split.values
  split.varIDs = forest$split.varIDs

  for (i in 1:num.trees) {
    which = match(split.varIDs[[i]], 0, FALSE)
    split.values[[i]][which > 0] = seq(sum(which))
    which.list[[i]] = which
  }

  forest$split.values = split.values
  object$forest = forest
  preds = stats::predict(object, x, predict.all = TRUE)$predictions

  ### Get list of means indexed by the unique terminal node values
  pmeans = apply(preds, 2, function(f, y) tapply(y, f, mean), y)

  ### Now populate these terminal nodes with these values
  for (i in 1:num.trees) {
    which = which.list[[i]]
    repvec = rep(NA, sum(which))
    pmean = pmeans[[i]]
    ids = as.integer(names(pmean))
    repvec[ids] = pmean
    split.values[[i]][which > 0] = repvec
  }

  forest$split.values = split.values
  object$forest = forest
  object$mean = mean(y)
  class(object) = c('pollinated.ranger', 'ranger')
  object
}