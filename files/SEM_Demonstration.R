#############################################################
# SEM                                                       #
#                                                           #
# This is an "0educational deep dive" designed to showcase  #
#   the Matrix Algebra behind Structural Equation Modeling. #
# The model we recreate is "Political Democracy",           #
#   often used as a tutorial model.                         #
#   See https://lavaan.ugent.be/tutorial/sem.html           #
# It presumes at least some understanding of factor analysis#
#  and some knowledge of R                                  #
#                                                           #
# Created By Jeremy Rappel & Shahryar Ebrahimi, 2022        #
# Any mistakes are the fault of my own (Jeremy)             #
#    and please let me know if you spot any                 #
#     (jeremy.rappel@mail.mcgill.ca)                        #
# Last updated December 2025                                #
#############################################################

#First package we're gonna load is lavaan;
#  we use it both because it contains the "Political Democracy"
#  dataset (Bollen, 1989), and to check our results at the end
#install.packages("lavaan")
library(lavaan)

#Second is the "Matrix" package; this contains the
#  bdiag() function, which will help us create sparse matrices
#install.packages("Matrix")
library(Matrix) 

#Here we're renaming the PoliticalDemoncracy dataset to something shorter.
dat <- PoliticalDemocracy

# The model we'll recreate consists of three latent variables:
#  "Democracy 1960" (dem60), consisting of 4 observed variables (y1-y4);
#  "Democracy 1965" (dem65), also consisting of 4 observed variables (y5-y8);
#  "Industry 1960" (ind60), consisting of 3 observed variables (x1-x3)
# "dem60" is regressed onto "ind60", while "dem65" is regressed onto both "dem60" and "ind60"
# "dem60" and "dem65" are considered endogenous variables (they are regressed onto another variable;
#   they have paths leading to them)
# "ind60" is an exogenous variable (not regressed onto another variable;
#    There are no paths leading to it)
# Several of the indicator residuals are also correlated (y1-y5; y2-y4; y2-6; y3-y7; y4-8, and y6-y8), for  total of 6 correlations
# Reference the path diagram at https://lavaan.ugent.be/tutorial/sem.html for a visualization

#On the backend, an SEM model works with 11 matrices:
#   The endogenous lambda (Λ_y) matrix, which contains the loadings for each endogenous latent variable
#   The exogenous lambda (Λ_x) matrix, which contains the loadings for each exogenous latent variable
#   The beta (β) matrix, which contains the regression coefficients between latent endogenous variables
#   The gamma (Γ) matrix, which contains the regression coefficients between latent exogenous variables
#   The theta (Θ) matrix, which contains our estimate error varainces and covariances
#     Theta is broken down into theta_epsilon (Θ_ε) for endogenous variables and theta_delta (Θ_δ) for exogenous variables
#   The Phi (Φ) matrix, which is the covariance matrix between exogenous variables
#   The Psi (Ψ) matrix, which is the covariance matrix between endogenous latent variables
#   The observed covariance/correlation (S) matrix of our indicator variables
#   The implied covariance/correlation (Σ) matrix  between our indicator
#     This last one is the end goal; it's calculated from other matrices (excluding S)
#     Basically it asks "if our model perfectly explained our data, what would the covariance matrix (S) look like?"
#     We determine how well our model fits the data by how closely the Σ and S matrices line up; this is what (absolute) fit indices tell us


# The values in these matrices are iteratively adjusted until our model fit stops improving
# Note we don't model mean structures (which are necessary for more advanced applications like latent growth modeling)

# Because we have 11 indicator variables that load onto latent variables, 6 correlations, 3 regression paths, 11 error variances and 3 latent variable variances
#   this model will estimate 34 separate values in the above matrices
#   (But really we'll only estimate 31, for reasons described below)

# SEM models need starting values for each estimate; here we give 1.5 as the starting value for most of these
start <- rep(1.5,34)
# The order of theses variables is: 
#   1-11 are the latent variable loadings (for the Λ matrix)
#   12 is the regression paths originating from endogenous variables (dem65 ~ dem60; for the β matrix) 
#   13-14 are the regression paths originating from exogenous variables (dem65 ~ ind60, dem60~ind60; for the Γ matrix)
#   15-25 are the indicator variable variances (Diagonals of the Θ matrix)
#   26-31 are the indicator variable covariances (In the Θ matrix)
#   32-34 are the covariances between latent variables (Φ and Ψ matrices)

#There's a matrix algebra quirk that if covariances are equal to variances, a matrix will not be invertible; this will interfere with our estimation
#   Here we fix this by changing the starting values for the variances and covariances
start[26:31] <- 0 
start[32:34] <- .5

#Here we create our custom function. It takes as arguments: 
#   a list of starting values (start)
#   a dataset to analyse (dat)
#   Whether to use population covariances (wish=F) or sample covariances (wish=T)
#      Using sample covariance matrices is called "Wishart estimation", thus the naming convention
#   An estimator to use (estimator)
#      "ML" stands for maximum likelihood
#      We could also use some other options; "ULS" for Unweighted Least Squares, or "GLS" for Generalized Least Squares
#      A discussion of the differences between them is outside the scope of this project
#   Whether we use "variance" or "loading" scaling (scaling)
#      SEM needs to fix some values to 1 in order to estimate things correctly
#      This can be either fixing the first loading of each latent variable to 1 ("loading scaling")
#      Or by fixing the latent variable variances to 1 ("variance" scaling)
#   A list of how many variables are estimated in which matrix (parameterNumb)
#      Here, 11 variables in Λ, 1 in β, 2 in Γ, 11 for the diagonals (variances) of Θ, 6 for the covariances in Θ, 2 in Ψ, and 1 in Φ

semCustom <- function(start, dat, wish = FALSE, estimator = "ML", scaling="variance", parameterNumb = c(11, 1, 2, 11, 6,2,1)){

# functions -----------------------------------
  if (scaling=="variance") {
    # If we're doing variance scaling, we fix factor variances to 1
    start[32:34] = 1 
  } else if (scaling=="loading"){ 
    #Alternatively, if we're doing loading scaling, we fix the first loading of each latent variable to 1
    start[1] = 1 
    start[5] = 1
    start[9] = 1
  }
  
  # calculating cost function
  setFunction <- function(start, dat, wishart = wish, scale=scaling,fit = estimator){
    
    #Λ is our loadings matrix. 
    #We break lambda out for exogenous (lambda_x) and endogenous (lambda_y) variables
    if (scaling=="variance") {
      lambda_y <- bdiag(start[1:4],start[5:8])
      lambda_x <- bdiag(start[9:11])
    } else if (scaling=="loading"){
      # for "loading" scaling, the first indicator of each latent variable is fixed to 1
      lambda_y <- bdiag(c(1,start[2:4]),c(1,start[6:8]))
      lambda_x <- bdiag(c(1,start[10:11]))
    }
    # The Beta matrix contains the regression coefficients between endogenous variables.
    beta <- matrix(c(0,start[12],0,0),nrow=2,ncol=2)
    # The layout here is a bit unintuitive; it creases a 2x2 matrix where only the bottom left cell has a value
    # We can think of it as:
    #        dem60 | dem65
    # dem60 |      |
    # dem65 | VALUE|
    #
    # Since out model has dem65 ~ dem60 as a path (and not dem60 ~ dem60, dem65~dme65, or dem60 ~ dem65), only the bottom left cell has a value
    
    # The Gamma matrix contains the regression coefficients from exogenous to endogenous variables. 
    # The t() function transposes a given matrix
    gamma <- t(t(start[13:14])) #exogenous regression coefficients
    # We can think of this one as
    #        ind60 
    # dem60 | VALUE     
    # dem65 | VALUE     
    # Since the dem60 ~ ind60 and dem60 ~ ind65 paths are both in our model, both these cells will have values in them
    

    # The theta matrix represents our (estimated) error variances and covariances. 
    # Theta_epsilon is for endogenous variables, theta_delta is for exogenous
    theta <- diag(start[15:25],ncol=ncol(PoliticalDemocracy),nrow=ncol(PoliticalDemocracy))
    theta[5,1] <- start[26]
    theta[4,2] <- start[27]
    theta[6,2] <- start[28]
    theta[7,3] <- start[29]
    theta[8,4] <- start[30]
    theta[8,6] <- start[31]
    #The upper and lower triangles of the thera matrix are identical
    theta[upper.tri(theta)] <- t(theta)[upper.tri(theta)]
    #Note how only 6 some non-diagonal cells have values (but that they're mirrored in the upper & lower triangles); these are the cells that correspond to the correlations between y-values in the model
    

    
    
    # The Phi matrix is a covariance matrix between exogenous latent variables
    # The Psi matrix is a covariance matrix between endogenous latent variables
    if (scaling=="variance") {
      # With variance scaling, the latent variable variances are fixed to 1
    psi <- diag(1,2)
    phi <- diag(1)
    } else if (scaling=="loading"){
      # Otherwise, the latent variable variances are freely estimated
      psi <- diag(start[32:33],2)
      phi <- diag(start[34], nrow=length(start[34]))
    }

    # Get the total number of factors
    nfactors <- ncol(lambda_x)+ncol(lambda_y) #1 exogenous plus 2 endogenous, 3 total
    # P is the total number of observed variables; here, 11
    p <- ncol(dat)
    # N is the sample size; here, 75
    N <- nrow(dat)
    
    # Here we break theta into endogenous and exogenous bits
    theta_epsilon <- theta[1:nrow(lambda_y),1:nrow(lambda_y)]
    theta_delta <- theta[(p-nrow(lambda_x)+1):p,(p-nrow(lambda_x)+1):p]
    #Σ is the model-implied covariance matrix; it's a supermatrix of 4 smaller matrices
    #      Σ_YY   Σ_XY
    # Σ =  Σ_YX   Σ_XX
    
    # The %*% indicates matrix multiplication
    # the solve() function inverts a given matrix; indicated mathematically as ()^-1
    # The t() function transposes a given matrix; indicated mathematically as '
    # The diag() function is used here to create Identity matrices
    #   An identity matrix is a matrix with 1 on the diagonals and zero in the other cells
    
    #Σ_yy = Λ_y(I-β)^-1((ΓΦΓ'+Ψ)(I-β)^-1)Λ_y)+ Θ_ε
    sigma_YY <- lambda_y %*% solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta) %*%((gamma%*%phi%*%t(gamma)+psi)%*% t(solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta))%*%t(lambda_y))+theta_epsilon
    
    #Σ_yy = Λ_xΦΓ'(I-β)^-1'Λ_y
    sigma_XY <- lambda_x %*% phi %*% t(gamma) %*% t(solve(diag(1,nrow=nrow(beta),ncol=ncol(beta))-beta)) %*% t(lambda_y)
    #We transpose sigma_XY to get sigma_YX
    
    #Σ_yy = Λ_xΦΛ_x' + Θ_δ
    sigma_XX <- lambda_x %*% phi %*% t(lambda_x) + theta_delta
    
    #Combining the 4 into one supermatrix
    sigma <- cbind(rbind(sigma_YY,sigma_XY),rbind(t(sigma_XY), sigma_XX))
    
    if (wishart==TRUE){
      # If we're using wishart estimation, use the sample covariance matrix
      S <- cov(dat)
    } else{
      #Otherwise, use the population covariance matrix
      S <- cov(dat)*(nrow(dat)-1)/(nrow(dat))
    }

    
    #(Negative) Log likelihood
    tmp <- matrix(array(S%*%solve(sigma)), p ,p)
    LL0 <- (-N*p/2)*log(2*pi) - (N/2)*log(det(sigma)) - (N/2)*(sum(diag(tmp)))
    LL1 <- (-N*p/2)*log(2*pi) - (N/2)*log(det(S)) - (N/2)*p
    #nlminb seems to give an error when we use the negative log likelihood (what stats programs use), but works with the positive log likelihood. Weird.
    
    LL0 <- -LL0
    LL1 <- -LL1
    
    
    if (fit == "ML"){
    #Maximum Likelihood Fit Function
    fit <- log(det(sigma)) - log(det(S)) + sum(diag(solve(sigma)%*%S)) - p
    } else if (fit == "ULS") {
    #Unweighted Least Squares Fit Function
    tmp <- S-sigma
    fit <- 0.5*sum(tmp^2)
    } else if (fit == "GLS") {
    #Generalized Least Squares Fit Function
    tmp <- (S-sigma)%*%solve(sigma)
    fit <- 0.5*sum(tmp^2)
    }else if (fit == "LL0") {
      fit <- LL0
    }
    
    #Standard errors can be calculated by taking the expected (or sometimes observed) information matrix
    #   inverting that matrix, taking the diagonal values, multiplying them by the sample size,
    #   then square rooting those values
    # However, calculating information matrices is outside the scope of this project
    # If you're curious, you can calculate them from lavaan output via:
    #   se <-sqrt(diag(solve(lavInspect(fit, "information.expected")*75)))
    
    return(fit)
  
  }
  
  # calculating model degrees of freedom
  DF <- function(start, dat){    
    
    # P is the total number of observed variables (here, 11)
    p <- ncol(dat)
    
    Basedf <- p*(p+1)/2 # Baseline degress of freedom is equal to the total number of variances and covariances between the observed/indicator variables
    Nulldf <- Basedf-p # The null model is where where all the observed/indicator variables are uncorrelated, and there's no factor structure
    Modeldf <- Basedf-length(start)+sum(start==1) #We have 34 starting values, but the scaling we do fixes three of them (1 for each factor)
    # sum(start==1) is a (lazy) way to calculate the number of latent variables (see lines 22-28)
    
    return(c(Basedf, Nulldf, Modeldf))
    
  }
  
  # Calculating chi-square tests; this will vary a bit depending on if we do Wishart estimation or not
  ChiTest <- function(start, dat, fit, wishart = wish){
    
    if (wishart==TRUE){
      S <- cov(dat)
      chisqmodel <- fit*(nrow(dat)-1)
      chisqnull <- (log(det(diag(diag(S))))-log(det(S))+sum(ncol(dat))-ncol(dat))*(nrow(dat)-1)
    } else{
      S <- cov(dat)*(nrow(dat)-1)/(nrow(dat))
      chisqmodel <- fit*nrow(dat)
      chisqnull <- (log(det(diag(diag(S))))-log(det(S))+sum(ncol(dat))-ncol(dat))*nrow(dat)
    }
    
    return(c(chisqmodel, chisqnull))
    
  }
  
  # calculating fit indices
  Findices <- function(start, dat, fit, wishart = wish){

    chiout      <- ChiTest(start, dat, fit, wishart = wishart)
    dfout       <- DF(start, dat)
    
    chisqmodel <- chiout[1]
    chisqnull  <- chiout[2]
    
    Basedf     <- dfout[1]
    Nulldf     <- dfout[2]
    Modeldf    <- dfout[3]
    
    N          <- nrow(dat)
    
    LL0        <- setFunction(start, dat, wishart = wishart, fit = "LL0")
      
    #Fit indices
    # There are many fit indices, and all the variants can seem overwhelming
    # While norms differ according to journal and field, CFI, TLI, RMSEA, and SRMR are the most common
    # Generally, Absolute fix indices indicate better fit the close they are to 1
    
    #Technically, the Confirmatory factor Index is a correction to Relative Noncentrality Index; CFI avoids bias when sample sizes are small
    CFI <- 1- (chisqmodel-Modeldf)/(chisqnull-Nulldf)
    # Tucker-Lewis Index
    #    Also called the Non-Normed Fit Index
    TLI <- ((chisqnull/Nulldf)-(chisqmodel/Modeldf))/((chisqnull/Nulldf)-1)
    # Bentler-Bonnet Normed Fit Index
    BBNFI <- 1 - (chisqmodel/chisqnull)
    # Parsimony Normed Fit Index
    PNFI <- (Modeldf/Nulldf)*BBNFI
    # Bollen's Relative Fit Index
    BRFI <- ((chisqnull/Nulldf)-(chisqmodel/Modeldf))/(chisqnull/Nulldf)
    # Bollen's Incremental Fit Index
    BIFI <- (chisqnull-chisqmodel)/(chisqnull-Modeldf)
    #Relative Noncentrality Index
    RNI <- ((chisqnull-Nulldf) - (chisqmodel-Modeldf))/(chisqnull-Nulldf)
    
    #Root Mean Squared Error of Approximation
    RMSEA <- sqrt((chisqmodel-Modeldf)/(Modeldf*N))
    
    #Akaike Information Criteria
    AIC <- (-2*(-LL0))+2*Nulldf
    #Bayesian Information Criteria
    BIC <- (-2*(-LL0))+Nulldf*log(N)
    
    #AIC and BIC and relative fit indices; they are used to compare models, with lower values indicating better modle fit (correcting for extra variables added to the model)
    
    #Sample-Size Adjusted BIC
    SSABIC <- (-2*(-LL0))+Nulldf*log((N+2)/24)
    
    #Expected Cross Validation Index
    ECVI <- fit +(2*(Basedf-Modeldf))/N
    
    #Hoelter's Critical N, alpha=.05
    HCN05 <- qchisq(.05, Modeldf, lower.tail=FALSE)/fit+1
    
    #Hoelter's Critical N, alpha=.01
    HCN01 <- qchisq(.01, Modeldf, lower.tail=FALSE)/fit+1
    
    #McDonald's Fit Index
    MFI <- exp(-.5*((chisqmodel-Modeldf)/N))
    
    #Another common fit index is Standardized Root Mean Squared Residual (SRMR)
    #  There are a few variants for calculating it, so numbers will often differ across stats programs
    
    return(c(SSABIC, BIC, AIC, RMSEA, RNI, BIFI, BRFI, PNFI, BBNFI, TLI, CFI, ECVI, HCN05, HCN01, MFI))
    
  }
  
  # printing result
  printResult <- function(output.estimate, output.df, output.chi, output.indices){
    
    parameterName = c("loading", "Beta", "Gamma", "Variance", "Covariance","Psi","Phi")
    IndicesName   = c("SSABIC", "BIC", "AIC", "RMSEA", "RNI", "BIFI", "BRFI", "PNFI", "BBNFI", "TLI", "CFI", "ECVI", "HCN05", "HCN01", "MFI")
    
    writeLines(sprintf(""))
    writeLines(sprintf("--------------------------------"))
    
    for(rep in 1:length(parameterNumb)){
      
      if (rep==1){
        rs <- 1
        re <- parameterNumb[1]
      }else{
        rs <- re+1
        re <- re+parameterNumb[rep]
      }
    
        writeLines(sprintf("%s %s = %f", parameterName[rep], 1:parameterNumb[rep], output.estimate$par[rs:re]))
        writeLines(sprintf("--------------------------------"))
        
    }  
    
    writeLines(sprintf("Base  DF = %f", output.df[1]))
    writeLines(sprintf("Null  DF = %f", output.df[2]))
    writeLines(sprintf("Model DF = %f", output.df[3]))
    
    writeLines(sprintf("--------------------------------"))
    
    writeLines(sprintf("Model Chi-Sq = %f", output.chi[1]))
    writeLines(sprintf("Null  Chi-Sq = %f", output.chi[2]))
    
    writeLines(sprintf("--------------------------------"))
    
    for (rep in 1:length(IndicesName)){
      writeLines(sprintf("%s  = %f", IndicesName[rep], output.indices[rep]))
    }
    
    writeLines(sprintf("--------------------------------"))
    writeLines(sprintf(""))
    
  } 

# Outputs -------------------------------------
# 
  
  # saving outputs
  output.estimate <- nlminb(start, setFunction, dat=PoliticalDemocracy)
  # We use the nlminb() function for estimation
  output.df       <- DF(start, dat)
  output.chi      <- ChiTest(output.estimate$par, dat, fit = output.estimate$objective, wishart = wish)
  output.indices  <- Findices(output.estimate$par, dat, fit = output.estimate$objective, wishart = wish)
  
  printResult(output.estimate, output.df, output.chi, output.indices)
    
  return(list(output.estimate, output.df, output.chi, output.indices))
  
}


res <- semCustom(start, PoliticalDemocracy, wish = FALSE, estimator = "ML", scaling="loading",parameterNumb = c(11, 1, 2, 11, 6,2,1))


#Comparing the results to lavaan
model <- '
  # measurement model
    ind60 =~ x1 + x2 + x3
    dem60 =~ y1 + y2 + y3 + y4
    dem65 =~ y5 + y6 + y7 + y8
  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60
  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'

fit <- sem(model, data = PoliticalDemocracy)
summary(fit, standardized = TRUE)


