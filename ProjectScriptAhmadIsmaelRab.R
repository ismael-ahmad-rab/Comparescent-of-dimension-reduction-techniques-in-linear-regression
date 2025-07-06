
library(dplyr) # For dataframe manipulation
library(pls)
library(spcr)
library(ggplot2)
library(MASS)
library(glmnet)
library(caret)

######################################################################################################

pcr_fit = function(formula, X.train, y.train, X.test, y.test, s, response){
    
    # Scale X.train and y.train
    train_means = colMeans(X.train)
    train_sd = apply(X.train, 2, sd)
    X.train.scaled = scale(X.train, center = TRUE, scale = TRUE)
    y.train.scaled = scale(y.train, center = TRUE, scale = TRUE)
    
    # Join again
    df.train.scaled = data.frame(cbind(X.train.scaled, y.train.scaled))
    colnames(df.train.scaled) = c(colnames(X.train), response)
    
    # Scale X.test and y.test accordingly
    X.test.scaled = scale(X.test, center = train_means, scale = train_sd)
    y.test.scaled = scale(y.test, center = mean(y.train), scale = sd(y.train))
    
    # Join again
    df.test.scaled = data.frame(cbind(X.test.scaled, y.test.scaled))
    colnames(df.test.scaled) = c(colnames(X.test), response)
    
    # Perform cross-validation
    formula = as.formula(formula)
    set.seed(s)
    mod.pcr.cv <- pls::pcr(formula, data = df.train.scaled, center = TRUE, scale = TRUE, validation = "CV")
    
    # Get results of cross-validation and get the best number of PCs to use
    cv.results = mod.pcr.cv$validation$PRESS
    dim.pcs = which.min(cv.results)
    
    # Train the final model
    mod.pcr <- pls::pcr(formula, data = df.train.scaled, center = FALSE, scale = FALSE, validation = "none")
    
    # Predict
    y.pred.scaled = predict(mod.pcr, newdata = df.test.scaled, ncomp = dim.pcs)
    y.pred = (sd(y.train)*y.pred.scaled) + mean(y.train)
    pcr.error = mean((y.test - y.pred)^2)
    
    return(list(dim.pcs, mod.pcr, pcr.error))
    
}

######################################################################################################

plsr_fit = function(formula, X.train, y.train, X.test, y.test, s, response){
    
    # Scale X.train and y.train
    train_means = colMeans(X.train)
    train_sd = apply(X.train, 2, sd)
    X.train.scaled = scale(X.train, center = TRUE, scale = TRUE)
    y.train.scaled = scale(y.train, center = TRUE, scale = TRUE)
    
    # Join again
    df.train.scaled = data.frame(cbind(X.train.scaled, y.train.scaled))
    colnames(df.train.scaled) = c(colnames(X.train), response)
    
    # Scale X.test and y.test accordingly
    X.test.scaled = scale(X.test, center = train_means, scale = train_sd)
    y.test.scaled = scale(y.test, center = mean(y.train), scale = sd(y.train))
    
    # Join again
    df.test.scaled = data.frame(cbind(X.test.scaled, y.test.scaled))
    colnames(df.test.scaled) = c(colnames(X.test), response)
    
    # Perform cross-validation
    formula = as.formula(formula)
    set.seed(s)
    mod.plsr.cv <- pls::plsr(formula, data = df.train.scaled, center = TRUE, scale = TRUE, validation = "CV")
    
    # Get results of cross-validation and get the best number of PCs to use
    cv.results = mod.plsr.cv$validation$PRESS
    dim.pcs = which.min(cv.results)
    
    # Train the final model
    mod.plsr <- pls::plsr(formula, data = df.train.scaled, center = FALSE, scale = FALSE, validation = "none")
    
    # Predict
    y.pred.scaled = predict(mod.plsr, newdata = df.test.scaled, ncomp = dim.pcs)
    y.pred = (sd(y.train)*y.pred.scaled) + mean(y.train)
    plsr.error = mean((y.test - y.pred)^2)
    
    return(list(dim.pcs, mod.plsr, plsr.error))
    
}

######################################################################################################

spcr_cross_val = function(X, y, n.folds, ks, lambda.betas, lambda.gammas, s){
    
    # Create grid
    grid = list(
        kss = ks,
        lambda.Bs = lambda.betas,
        lambda.gammass = lambda.gammas
    )
    
    # Initialise vectors of average errors and values
    average.errors = c()
    hp.com = c()
    
    # Set seed and create folds
    set.seed(s)
    flds = createFolds(y = y, k = n.folds, list = TRUE, returnTrain = FALSE)
    
    # Go through each combination of hyper-parameters
    for (k in grid$kss) {
        for (lambda.B in grid$lambda.Bs){
            for (lambda.gamma in grid$lambda.gammass){
                
                error = 0
                com = list(k = k, lambda.B = lambda.B, lambda.gamma = lambda.gamma)
                hp.com = c(hp.com, list(com))
                
                # Go through the folds
                for (i in n.folds){
                    
                    # Get training and validation data
                    train_index = flds[[i]]
                    X.train = X[train_index, ]
                    X.val = X[-train_index, ]
                    y.train = y[train_index]
                    y.val = y[-train_index]
                    
                    # Scale data
                    train_means = colMeans(X.train)
                    train_sd = apply(X.train, 2, sd)
                    
                    X.train.scaled = scale(X.train, center = TRUE, scale = TRUE)
                    X.val.scaled = scale(X.val, center = train_means, scale = train_sd)
                    
                    y.train.scaled = scale(y.train, center = TRUE, scale = TRUE)
                    y.train.scaled = as.vector(y.train.scaled)
                    y.val.scaled = scale(y.val, center = mean(y.train), scale = sd(y.train))
                    y.val.scaled = as.vector(y.val.scaled)
                    
                    # Train the model
                    mod.spcr = spcr(x = X.train.scaled, y = y.train.scaled, k = k,  lambda.B = lambda.B,
                                    lambda.gamma = lambda.gamma, scale = FALSE, center = FALSE)
                    
                    # Get the prediction and error
                    b = mod.spcr$loadings.B
                    coeffs = mod.spcr$gamma
                    coef0 = mod.spcr$gamma0
                    betas = as.matrix(c(coef0, coeffs))
                    
                    X.val.pc = X.val.scaled %*% b
                    X.val.pc = cbind(1, X.val.pc)
                    
                    y.pred.scaled = array(X.val.pc %*% betas)
                    y.pred = (sd(y.train)*y.pred.scaled) + mean(y.train)
                    spcr.error = mean((y.val - y.pred)^2)
                    
                    error = error + spcr.error
                    
                } 
                
                average.errors = c(average.errors, (error/n.folds))
                
            }
        }
    }
    
    # Get best hyperparameters
    index.best.hps = which.min(average.errors)
    best.hyperparameters = hp.com[index.best.hps][[1]]
    
    return(best.hyperparameters)
    
}

######################################################################################################

spcr_fit = function(X.train, X.test, y.train, y.test, k, lambda.B, lambda.gamma){
    
    # Scale data
    train_means = colMeans(X.train)
    train_sd = apply(X.train, 2, sd)
    
    X.train.scaled = scale(X.train, center = TRUE, scale = TRUE)
    X.test.scaled = scale(X.test, center = train_means, scale = train_sd)
    
    y.train.scaled = scale(y.train, center = TRUE, scale = TRUE)
    y.train.scaled = as.vector(y.train.scaled)
    y.test.scaled = scale(y.test, center = mean(y.train), scale = sd(y.train))
    y.test.scaled = as.vector(y.test.scaled)
    
    mod.spcr = spcr(x = X.train.scaled, y = y.train.scaled, k = k,  lambda.B = lambda.B,
                    lambda.gamma = lambda.gamma, scale = FALSE, center = FALSE)
    
    b = mod.spcr$loadings.B
    coeffs = mod.spcr$gamma
    coef0 = mod.spcr$gamma0
    betas = as.matrix(c(coef0, coeffs))
    
    X.test.pc = X.test.scaled %*% b
    X.test.pc = cbind(1, X.test.pc)
    
    y.pred.scaled = array(X.test.pc %*% betas)
    y.pred = (sd(y.train)*y.pred.scaled) + mean(y.train)
    spcr.error = mean((y.test - y.pred)^2)
    
    return(list(mod.spcr, spcr.error))
    
}

######################################################################################################

spcr_coeffs = function(mod, X.train, y.train){
    
    # Retrieve loadings
    b = mod$loadings.B
    
    # Retrieve the coefficients in terms of the PCs
    pc.coeffs = as.matrix(mod$gamma)
    pc.coef0 = mod$gamma0
    
    # Retrieve the number of original variables
    p = nrow(b)
    
    # Retrieve the number of PCs
    k = ncol(b)
    
    # Retrieve the PCs coefficients in terms of the original scaled variables
    or.scaled.coeffs = b%*%pc.coeffs
    
    # Transform to get in the original scale
    pre.or.coeffs = or.scaled.coeffs/apply(X.train, 2, sd)
    or.coeffs = pre.or.coeffs*sd(y.train)
    
    # Also get the original gamma
    pre.or.intercept = pc.coef0 - sum(pre.or.coeffs*apply(X.train, 2, mean))
    or.intercept = (pre.or.intercept*sd(y.train)) + mean(y.train)
    
    # Format into a dataframe
    coeffs.or = data.frame(Coefficient.Value = c(or.intercept, or.coeffs))
    rownames(coeffs.or) = c("Intercept", colnames(X.train))
    
    # Do the same for the loadings
    pc.names = paste("PC.", c(1:k), sep = "")
    df.loadings = data.frame(b)
    rownames(df.loadings) = colnames(X.train)
    colnames(df.loadings) = pc.names
    
    # Return dataframes
    return(list(df.loadings, coeffs.or))
    
}

######################################################################################################

pls_library_coeffs = function(mod, X.train, y.train, ncomp){
    
    # Retrieve the number of original variables
    p = ncol(X.train)
    
    # Retrieve the loadings
    lods = loadings(mod)
    
    # Get what is of interest
    lods = as.matrix(lods)
    lods = lods[1:p, 1:ncomp]
    
    # Format into a dataframe
    pc.names = paste("PC.", c(1:ncomp), sep = "")
    df.loadings = data.frame(lods)
    rownames(df.loadings) = colnames(X.train)
    colnames(df.loadings) = pc.names
    
    # Retrieve the coefficients in terms of the original scaled variables and then the original variables
    or.scaled.coeffs = as.vector(mod$coefficients[, , ncomp])
    pre.or.coeffs = or.scaled.coeffs/apply(X.train, 2, sd)
    or.coeffs = pre.or.coeffs*sd(y.train)
    
    # Compute the intercept
    pre.or.intercept = 0 - sum(pre.or.coeffs*apply(X.train, 2, mean))
    or.intercept = (pre.or.intercept*sd(y.train)) + mean(y.train)
    
    # Return as a dataframe
    coeffs.or = data.frame(Coefficient.Value = c(or.intercept, or.coeffs))
    rownames(coeffs.or) = c("Intercept", colnames(X.train))
    
    # Return dataframes
    return(list(df.loadings, coeffs.or))
    
}

######################################################################################################

sim_data = function(s, mu, matr, n, true.coeffs, err.sd){
    
    # Retrieve the number of features to use
    p = ncol(matr)
    
    # Simulate the features and error
    set.seed(s)
    X = mvrnorm(n, mu = mu, Sigma = matr)
    eps = rnorm(n, mean = 0, sd = err.sd)
    
    # Get the complete feature matrix and compute the response
    X.comp = cbind(1, X)
    y = X.comp%*%true.coeffs + eps
    
    # Return the dataset as a dataframe
    dataset = cbind(X, y)
    df = data.frame(dataset)
    colnames(df) = c(paste("X", 1:p, sep = ""), "y")
    
    return(df)
    
}

######################################################################################################

matr = c(0.5, 0.1, -0.1, 0.1, 0.8, 0.2, 0.1, 0.1, 0.1, 0.1,
         0, 0.5, -0.1, 0.1, 0.2, 0.6, 0.1, 0.1, 0.1, 0.1,
         0, 0, 0.5, -0.2, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
         0, 0, 0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
         0, 0, 0, 0, 0.5, 0.1, 0.1, 0.1, 0.1, 0.1,
         0, 0, 0, 0, 0, 0.5, 0.1, 0.1, 0.1, 0.1,
         0, 0, 0, 0, 0, 0, 0.5, 0.1, 0.1, 0.1,
         0, 0, 0, 0, 0, 0, 0, 0.5, 0.1, 0.1,
         0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0.9,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5)

######################################################################################################

matr = matrix(matr, nrow = 10, ncol = 10)
matr = matr + t(matr)

path = "C:/Users/ismae/Documents/Universities/Carlos III/Master in Statistics for Data Science/Resources/Taken/Term 4/Advanced Regression and Prediction/CW/final/3/covariance_matrix.csv"
matr.save = matr
colnames(matr.save) = paste("X", 1:p, sep = "")
rownames(matr.save) = paste("X", 1:p, sep = "")
#write.csv(matr.save, path)

use.vars = 4
p = ncol(matr)
mu = rep(0, p)
n = 100

set.seed(0)
true.coeffs = runif(use.vars+1, -10, 10)
true.coeffs = c(true.coeffs, rep(0, (p-use.vars)))

path = "C:/Users/ismae/Documents/Universities/Carlos III/Master in Statistics for Data Science/Resources/Taken/Term 4/Advanced Regression and Prediction/CW/final/3/coefficients.csv"
true.coeffs.save = true.coeffs
true.coeffs.save = matrix(true.coeffs.save, ncol = p+1)
colnames(true.coeffs.save) = c("Intercept", paste("X", 1:p, sep = ""))
#write.csv(true.coeffs.save, path)

err.sd = 1.5

######################################################################################################

response = "y"
formula = "y~."

ks = c(1:p)
lambda.betas = seq(from = 0.92, to = 94.887685, length.out = 10)
lambda.gammas = seq(from = 0.92, to = 94.887685, length.out = 10)
s =  0
n.folds = 5

t_prop = 0.8

top.seed = 100
errors = data.frame(PCR = rep(NA, top.seed), 
                    PLRS = rep(NA, top.seed), 
                    SPCR = rep(NA, top.seed))

coeffs.pcr = matrix(NA, nrow = top.seed, ncol = p+1)
coeffs.plsr = matrix(NA, nrow = top.seed, ncol = p+1)
coeffs.spcr = matrix(NA, nrow = top.seed, ncol = p+1)

comps.pcr = rep(NA, top.seed)
comps.plsr = rep(NA, top.seed)
comps.spcr = rep(NA, top.seed)

for (i in 1:top.seed){
    
    df = sim_data(i, mu, matr, n, true.coeffs, err.sd)
    set.seed(i)
    train_rows = round(t_prop*nrow(df))
    train_index = sample(1:train_rows, replace = F)
    
    df.train = df[train_index, ]
    df.test = df[-train_index, ]
    
    X.train = model.matrix(y~.-1, data = df.train)
    X.train = as.matrix(X.train)
    y.train = df.train$y
    
    X.test = model.matrix(y~.-1, data = df.test)
    X.test = as.matrix(X.test)
    y.test = df.test$y
    
    output = pcr_fit(formula, X.train, y.train, X.test, y.test, s, response)
    pcr.ncomp = output[[1]]
    comps.pcr[i] = pcr.ncomp
    mod.pcr = output[[2]]
    errors$PCR[i] = output[[3]]
    output = pls_library_coeffs(mod.pcr, X.train, y.train, pcr.ncomp)
    coeffs.pcr[i, ] = output[[2]]$Coefficient.Value
    
    output = plsr_fit(formula, X.train, y.train, X.test, y.test, s, response)
    plsr.ncomp = output[[1]]
    comps.plsr[i] = plsr.ncomp
    mod.plsr = output[[2]]
    errors$PLRS[i] = output[[3]]
    output = pls_library_coeffs(mod.plsr, X.train, y.train, plsr.ncomp)
    coeffs.plsr[i, ] = output[[2]]$Coefficient.Value
    
    best.hps = spcr_cross_val(X.train, y.train, n.folds, ks, lambda.betas, lambda.gammas, s)
    k = best.hps$k
    lambda.B = best.hps$lambda.B
    lambda.gamma = best.hps$lambda.gamma
    comps.spcr[i] = k
    output = spcr_fit(X.train, X.test, y.train, y.test, k, lambda.B, lambda.gamma)
    mod.spcr = output[[1]]
    errors$SPCR[i] <- output[[2]]
    output = spcr_coeffs(mod.spcr, X.train, y.train)
    coeffs.spcr[i, ] = output[[2]]$Coefficient.Value
    
    print(i)
    
}

coeffs.pcr = data.frame(coeffs.pcr)
colnames(coeffs.pcr) = c("Intercept", paste("X", 1:p, sep = ""))

coeffs.plsr = data.frame(coeffs.plsr)
colnames(coeffs.plsr) = c("Intercept", paste("X", 1:p, sep = ""))

coeffs.spcr = data.frame(coeffs.spcr)
colnames(coeffs.spcr) = c("Intercept", paste("X", 1:p, sep = ""))

boxplot(errors, main = "Distribution of the MSE")

ncomps.tot = data.frame(PCR = comps.pcr,
                        PLSR = comps.plsr,
                        SPCR = comps.spcr)

boxplot(ncomps.tot, main = "Distribution of the number of components used")

######################################################################################################

# Combine all dataframes into one
combined_df <- bind_rows(coeffs.pcr[, 6:ncol(coeffs.pcr)], .id = "id")

# Reshape data to long format
combined_df_long <- combined_df %>% 
    tidyr::pivot_longer(cols = -id, names_to = "Column")

# Create density plot
ggplot(combined_df_long, aes(x = value, fill = id)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Column, scales = "free_y") +
    labs(title = "Density Plots for Coefficient Estimates - PCR", x = "Value", y = "Density") +
    theme_minimal()

######################################################################################################

# Combine all dataframes into one
combined_df <- bind_rows(coeffs.plsr[, 6:ncol(coeffs.plsr)], .id = "id")

# Reshape data to long format
combined_df_long <- combined_df %>% 
    tidyr::pivot_longer(cols = -id, names_to = "Column")

# Create density plot
ggplot(combined_df_long, aes(x = value, fill = id)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Column, scales = "free_y") +
    labs(title = "Density Plots for Coefficient Estimates - PLSR", x = "Value", y = "Density") +
    theme_minimal()

######################################################################################################

# Combine all dataframes into one
combined_df <- bind_rows(coeffs.spcr[, 6:ncol(coeffs.spcr)], .id = "id")

# Reshape data to long format
combined_df_long <- combined_df %>% 
    tidyr::pivot_longer(cols = -id, names_to = "Column")

# Create density plot
ggplot(combined_df_long, aes(x = value, fill = id)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~Column, scales = "free_y") +
    labs(title = "Density Plots for Coefficient Estimates - SPCR", x = "Value", y = "Density") +
    theme_minimal()

######################################################################################################

data.pcr= rbind(apply(coeffs.pcr, 2, mean), apply(coeffs.pcr, 2, sd))
data.pcr = round(data.pcr, 2)
#df.mean.var.pcr = data.frame(data.pcr)
colnames(data.pcr) = c("Intercept", paste("X", 1:p, sep = ""))
rownames(data.pcr) = c("Mean", "Satandard Deviation")
path = "C:/Users/ismae/Documents/Universities/Carlos III/Master in Statistics for Data Science/Resources/Taken/Term 4/Advanced Regression and Prediction/CW/final/3/coeff_info_pcr.csv"
write.csv(data.pcr, path)

data.plsr= rbind(apply(coeffs.plsr, 2, mean), apply(coeffs.plsr, 2, sd))
data.plsr = round(data.plsr, 2)
#df.mean.var.plsr = data.frame(data.plsr)
colnames(data.plsr) = c("Intercept", paste("X", 1:p, sep = ""))
rownames(data.plsr) = c("Mean", "Satandard Deviation")
path = "C:/Users/ismae/Documents/Universities/Carlos III/Master in Statistics for Data Science/Resources/Taken/Term 4/Advanced Regression and Prediction/CW/final/3/coeff_info_plsr.csv"
write.csv(data.plsr, path)

data.spcr= rbind(apply(coeffs.spcr, 2, mean), apply(coeffs.spcr, 2, sd))
data.spcr = round(data.spcr, 2)
#df.mean.var.spcr = data.frame(data.spcr)
colnames(data.spcr) = c("Intercept", paste("X", 1:p, sep = ""))
rownames(data.spcr) = c("Mean", "Satandard Deviation")
path = "C:/Users/ismae/Documents/Universities/Carlos III/Master in Statistics for Data Science/Resources/Taken/Term 4/Advanced Regression and Prediction/CW/final/3/coeff_info_spcr.csv"
write.csv(data.spcr, path)


