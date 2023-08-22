# shrink ----------------------------------------------------------------
## Shrink function from the Rover and Zimmermann Camera Trapping book
# the function returns a shrinked matrix collapsed using nday; 
# if necessary, an X number of columns filled with NA are added
# to adjust the size of the shrinked matrix; 
# be careful that nday makes sense for the size of the matrix, so that
# not many columns of NA are added 
shrink <- function(matrice, nday) {
  dy <- nday
  while (dy < ncol(matrice)) {
    dy <- dy + nday
  }
  addcol <- dy - ncol(matrice)
  if (addcol != 0) {
    matNA <- matrix(NA, nrow = nrow(matrice), ncol = addcol)
    matrice <- data.frame(matrice, matNA)
  }
  
  period <- ncol(matrice)/nday
  newday <- rep(1:period, each = nday)
  
  shr <- function(vec) {
    nav <- is.na(vec)
    dom <- all(nav == T)
    if (dom == T) {
      y <- NA
    } else {
      abb <- sum(vec, na.rm = T)
      y <- ifelse(abb == 0, 0, 1)
    }
    return(y)
  }
  
  matday <- data.frame(newday, t(matrice))
  shrmat <- t(aggregate(matday[, -1], list(matday$newday), shr))
  
  return(shrmat[-1, ])
}

# MASON: You need random starts for every parameter. Here is the function I
# normally modify, you will need to modify it to fit your code,
#  but it'll serve you well! To use it you will need to know;
# 1. All the parameters you need to initialize.
# 2. The dimensions of those parameters (e.g,. are they in a matrix
#    or a vector).
# 3. Plug and chug from there, just name everything and give an
#    appropriate random number (e.g., positive for variance terms)

my_inits <- function(chain){
  gen_list <- function(chain = chain){
    list(
      pred_z = array(1, dim= c(data_list$nsite, data_list$npred, data_list$nseason)),
      pred_beta = matrix(
        rnorm(data_list$npred * data_list$nparm_pred_psi),
        ncol = data_list$nparm_pred_psi,
        nrow = data_list$npred
      ),
      pred_theta = rnorm(data_list$npred),
      pred_alpha = matrix(
        rnorm(data_list$npred * data_list$nparm_pred_rho),
        ncol = data_list$nparm_pred_rho,
        nrow = data_list$npred
      ),
      inxs_b0 = matrix(
        rnorm(data_list$npred * data_list$nprey),
        nrow = data_list$nprey,
        ncol = data_list$npred
      ),
      inxs_beta = array(
        rnorm(data_list$nprey * data_list$npred * (data_list$nparm_prey_psi - 1)),
        dim = c(data_list$nprey, data_list$npred , (data_list$nparm_prey_psi - 1))
      ),
      prey_theta = rnorm(data_list$nprey),
      prey_beta = matrix(
        rnorm(data_list$nprey * data_list$nparm_prey_psi),
        nrow = data_list$nprey,
        ncol = data_list$nparm_prey_psi
      ),
      prey_alpha = matrix(
        rnorm(data_list$nprey * data_list$nparm_pred_rho),
        nrow = data_list$nprey,
        ncol = data_list$nparm_prey_rho
      ),
      prey_z = array(1, dim= c(data_list$nsite, data_list$nprey, data_list$nseason)),
      .RNG.name = switch(chain,
                         "1" = "base::Wichmann-Hill",
                         "2" = "base::Wichmann-Hill",
                         "3" = "base::Super-Duper",
                         "4" = "base::Mersenne-Twister",
                         "5" = "base::Wichmann-Hill",
                         "6" = "base::Marsaglia-Multicarry",
                         "7" = "base::Super-Duper",
                         "8" = "base::Mersenne-Twister"),
      .RNG.seed = sample(1:1e+06, 1)
    )
  }
  return(switch(chain,
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
}


# function to split up mcsamp into a list with correctly shaped arrays
split_mcmc <- function(x){
  # get parameter names
  pars <- colnames(x)
  # unique parameters
  unq_pars <- unique(gsub("\\[.*\\]", "" , pars))
  # make list object to store arrays in
  result_list <- vector(
    "list",
    length = length(unq_pars)
  )
  names(result_list) <- unq_pars
  # fill in the arrays
  for(i in 1:length(result_list)){
    # get just the parameters
    tmp <- pars[grep(
      paste0(
        "^",unq_pars[i], "\\["
      ),
      pars
    )]
    
    if(length(tmp) == 0){
      tmp <- pars[grep(
        paste0("^",unq_pars[i],"$"),
        pars
      )]
    }
    
    # and then the array dimensions
    arr_dim <- gsub(
      paste0(
        unq_pars[i],"\\[|\\]"
      ),
      "",
      tmp
    )
    
    arr_dim <- strsplit(
      arr_dim,
      ","
    )
    
    ndim <- length(arr_dim[[1]])
    npar <- length(arr_dim)
    # make a matrix
    arr_ind <- suppressWarnings(
      matrix(
        as.numeric(
          unlist(arr_dim)
        ),
        ncol = ndim,
        nrow = npar,
        byrow = TRUE
      )
    )
    
    if(nrow(arr_ind) == 1 & ncol(arr_ind) == 1){
      arr_ind[1,1] <- 1
    }
    # get max index for each
    max_ind <- apply(arr_ind, 2, max)
    # and then fill in the array
    result_list[[i]] <- array(
      x[,tmp],
      dim = c(nrow(x), max_ind)
    )
  }
  return(result_list)
}

library(ggtext)
library(rlang)

element_textbox_highlight <- function(..., hi.labels = NULL, hi.fill = NULL,
                                      hi.col = NULL, hi.box.col = NULL, hi.family = NULL) {
  structure(
    c(element_textbox(...),
      list(hi.labels = hi.labels, hi.fill = hi.fill, hi.col = hi.col, hi.box.col = hi.box.col, hi.family = hi.family)
    ),
    class = c("element_textbox_highlight", "element_textbox", "element_text", "element")
  )
}

element_grob.element_textbox_highlight <- function(element, label = "", ...) {
  if (label %in% element$hi.labels) {
    element$fill <- element$hi.fill %||% element$fill
    element$colour <- element$hi.col %||% element$colour
    element$box.colour <- element$hi.box.col %||% element$box.colour
    element$family <- element$hi.family %||% element$family
  }
  NextMethod()
}