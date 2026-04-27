#' Functions to run other methods: for comparisons with TESSERA.
#' @author Florica Constantine
#' Dependencies in file: Matrix, MASS, mgcv, brms, CARBayes, INLA, BRISC.

require(Matrix)
require(MASS)
require(mgcv)
require(brms)
require(CARBayes)
require(INLA)
require(BRISC)

#' Wrapper to fit basic Poisson GLM.
#'
#' This function theoretically works for multiple areas by stacking everything
#' together and ignoring spatial correlation.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from glm, plus a time field for how long the
#'    function ran for.
#'
#' @importFrom stats glm
#' @importFrom stats poisson
#' @importFrom stats offset
#' @export
glm_wrapper <- function(z_list, X_list, library_size_list = NULL) {
  # Start clock
  t0_glm <- Sys.time()

  # Stack into a single vector/matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
    # z_vec <- z_vec / lib_vec

    # Fit GLM
    glm_out <- stats::glm(z_vec ~ 0 + X_mat + stats::offset(log(lib_vec)), family = stats::poisson())

  } else {
    # Fit GLM
    glm_out <- stats::glm(z_vec ~ 0 + X_mat, family = stats::poisson())
  }

  # Stop clock
  t1_glm <- Sys.time()
  glm_out[["time"]] <- difftime(t1_glm, t0_glm)

  return(glm_out)
}

#'  Wrapper to fit basic Negative Binomial GLM.
#'
#'  This function theoretically works for multiple areas by stacking everything
#'  together and ignoring spatial correlation.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from glm.nb, plus a time field for how long the
#'    function ran for.
#'
#' @note Depends on MASS.
#'
#' @importFrom MASS glm.nb
#' @importFrom stats offset
#' @export
glm_nb_wrapper <- function(z_list, X_list, library_size_list = NULL) {
  # Start clock
  t0_glm <- Sys.time()

  # Stack into a single vector/matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
    # z_vec <- z_vec / lib_vec

    # Fit GLM
    glm_out <- MASS::glm.nb(z_vec ~ 0 + X_mat + stats::offset(log(lib_vec)))
  } else {
    # Fit GLM
    glm_out <- MASS::glm.nb(z_vec ~ 0 + X_mat)
  }

  # Stop clock
  t1_glm <- Sys.time()
  glm_out[["time"]] <- difftime(t1_glm, t0_glm)

  return(glm_out)
}

#' Wrapper to fit a spatial spline model.
#'
#' This function theoretically works for multiple areas.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param coords_list List of coordinate matrices (x, y)---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param model_family "poisson" or "gaussian".
#'    Which model family to fit. Note that Poisson models don't scale well in the
#'    multi-area setting---the spline is fit on a per-area basis.
#' @param z_offset In the "gaussian" case, the counts z are transformed as
#'    log(z + z_offset).
#' @param spline_k Dimension of spline basis.
#' @param spline_basis Type of basis---must be valid for mgcv::s.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from mgcv::gam, plus a time field for how long the
#'    function ran for.
#'
#' @note Requires the mgcv library.
#'
#' @importFrom mgcv gam
#' @importFrom mgcv bam
#' @importFrom mgcv s
#' @importFrom stats as.formula
#' @importFrom stats gaussian
#' @importFrom stats poisson
#' @importFrom stats offset
#' @export
mgcv_gam_wrapper <- function(z_list,
                             X_list,
                             coords_list,
                             model_family = "poisson",
                             z_offset = 0.5,
                             spline_k = -1,
                             spline_basis = "tp",
                             library_size_list = NULL) {
  # Start clock
  t0_mgcv <- Sys.time()

  # Create a single set of vectors/covariate matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))
  coords <- as.matrix(Reduce(rbind, coords_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
    # z_vec <- z_vec / lib_vec
  } else {
    lib_vec <- rep(1, length(z_vec))
  }

  xc <- coords[, 1]
  yc <- coords[, 2]

  # Create an indicator for each area
  area_id <- rep(NA, length(z_vec))
  start_idx <- 1
  for (area_idx in 1:length(z_list)) {
    area_id[start_idx:(start_idx + length(z_list[[area_idx]]) - 1)] <- area_idx
    start_idx <- start_idx + length(z_list[[area_idx]])
  }

  # Create data frame
  beta_dim <- ncol(X_mat)
  X_df <- cbind(data.frame(X_mat), z_vec, xc, yc, area_id, lib_vec)
  # Memory
  rm(z_vec, X_mat, xc, yc, area_id, coords)

  # Name columns
  for (idx in 1:beta_dim) {
    colnames(X_df)[idx] <- paste0("X", idx)
  }
  colnames(X_df)[1 + beta_dim] <- "z"
  colnames(X_df)[2 + beta_dim] <- "x"
  colnames(X_df)[3 + beta_dim] <- "y"
  colnames(X_df)[4 + beta_dim] <- "area"
  colnames(X_df)[5 + beta_dim] <- "lib_size"
  # Ensure area is a factor so we can group by it
  X_df$area <- as.factor(X_df$area)

  # Create formula
  spline_str <- paste0(" + s(x, y, by=area, k=",
                       spline_k,
                       ", bs='",
                       spline_basis,
                       "') + offset(log(lib_size))")
  mgcv_formula <- stats::as.formula(paste0("z ~ 0 + ", paste(colnames(X_df)[1:beta_dim], collapse = " + "), spline_str))

  # Fit model
  if ("poisson" == model_family) {
    mgcv_out <- mgcv::bam(mgcv_formula,
                          family = stats::poisson(),
                          data = X_df)
  } else if ("gaussian" == model_family) {
    X_df$z <- log(X_df$z + z_offset)
    mgcv_out <- mgcv::bam(mgcv_formula,
                          family = stats::gaussian(),
                          data = X_df)
  } else {
    stop("Invalid model family.")
  }

  # Stop clock and return
  t1_mgcv <- Sys.time()
  mgcv_out[["time"]] <- difftime(t1_mgcv, t0_mgcv)
  return(mgcv_out)
}

#' Wrapper to fit basic linear model.
#'
#' This function theoretically works for multiple areas by stacking everything
#' together and ignoring spatial correlation.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param transform_z Boolean: log-transform z or not.
#' @param z_offset If transform_z, the counts z are transformed as
#'    log(z + z_offset).
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from lm, plus a time field for how long the
#'    function ran for.
#'
#' @importFrom stats lm
#' @importFrom stats offset
#' @export
lm_wrapper <- function(z_list,
                       X_list,
                       transform_z = TRUE,
                       z_offset = 0.5,
                       library_size_list = NULL) {
  # Start clock
  t0_glm <- Sys.time()

  # Stack into a single vector/matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # Transform z
  if (transform_z) {
    z_vec <- log(z_vec + z_offset)
  }

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))

    # Fit LM
    lm_out <- stats::lm(z_vec ~ 0 + X_mat + stats::offset(log(lib_vec)))
  } else {
    # Fit LM
    lm_out <- stats::lm(z_vec ~ 0 + X_mat)
  }

  # Fit LM
  lm_out <- stats::lm(z_vec ~ 0 + X_mat)

  # Stop clock
  t1_glm <- Sys.time()
  lm_out[["time"]] <- difftime(t1_glm, t0_glm)

  return(lm_out)
}

#' Wrapper to fit a CAR or SAR model with MCMC (stan).
#'
#' This function theoretically works for multiple areas by creating a
#'  block diagonal adjacency matrix.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param W_list List of spatial adjacency matrices---one matrix per area.
#' @param model_type "CAR" or "SAR"
#' @param model_family "poisson" or "gaussian".
#'    Which model family to fit. Note that Poisson models don't scale well in the
#'    multi-area setting---the spline is fit on a per-area basis.
#'    Note: at the time of writing, "SAR" + "poisson" does not work--
#'    this is not implemented.
#' @param z_offset In the "gaussian" case, the counts z are transformed as
#'    log(z + z_offset).
#' @param chains Number of MCMC chains to run.
#' @param iter Number of total iterations per MCMC chain.
#' @param warmup Burn-in for MCMC.
#' @param thin Thinning rate for MCMC.
#' @param cores Number of cores to use for MCMC.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from brms::brm, plus a time field for how long the
#'    function ran for.
#'
#' @note Requires the Matrix library.
#' @note Requires the brms library.
#'
#' @import Matrix
#' @importFrom Matrix bdiag
#' @importFrom brms brm
#' @importFrom stats as.formula
#' @importFrom stats gaussian
#' @importFrom stats poisson
#'
#' @export
BRMS_CAR_SAR_wrapper <- function(z_list,
                                 X_list,
                                 W_list,
                                 model_type = "CAR",
                                 model_family = "poisson",
                                 z_offset = 0.5,
                                 chains = 4,
                                 iter = 2000,
                                 warmup = 200,
                                 thin = 2,
                                 cores = 1,
                                 library_size_list = NULL)  {
  # Start clock
  t0_brm <- Sys.time()

  # Create a single set of vectors/covariate matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
  } else {
    lib_vec <- rep(1, length(z_vec))
  }

  beta_dim <- ncol(X_mat)
  # Create a single adjacency matrix
  W <- Matrix::bdiag(W_list)

  if ("CAR" == model_type) {
    form_str <- " + car(W, type='escar')"
  } else if ("SAR" == model_type) {
    form_str <- " + sar(W, type='lag')"
  } else {
    stop("Invalid model_type.")
  }

  # Create dataframe
  X_df <- cbind(data.frame(X_mat), z_vec, log(lib_vec))
  # Memory
  rm(z_vec, X_mat, lib_vec)
  for (idx in 1:beta_dim) {
    colnames(X_df)[idx] <- paste0("X", idx)
  }
  colnames(X_df)[1 + beta_dim] <- "z"
  colnames(X_df)[2 + beta_dim] <- "log_lib"
  # Add in offset to formula
  form_str <- paste0(form_str, " + offset(log_lib)")

  # Create formula
  BRM_formula <- stats::as.formula(paste0("z ~ 0 + ", paste(colnames(X_df)[1:beta_dim], collapse = " + "), form_str))

  # Fit model
  if ("poisson" == model_family) {
    brms_out <- brms::brm(
      BRM_formula,
      family = stats::poisson(),
      data = X_df,
      data2 = list(W = W),
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      cores = cores
    )
  } else if ("gaussian" == model_family) {
    X_df$z <- log(X_df$z + z_offset)
    brms_out <- brms::brm(
      BRM_formula,
      family = stats::gaussian(),
      data = X_df,
      data2 = list(W = W),
      chains = chains,
      iter = iter,
      warmup = warmup,
      thin = thin,
      cores = cores
    )
  } else {
    stop("Invalid model family.")
  }

  # Stop clock
  t1_brm <- Sys.time()
  brms_out[["time"]] <- difftime(t1_brm, t0_brm)

  return(brms_out)
}

#' Wrapper to fit a Leroux model with MCMC (CARBayes).
#'
#' This function theoretically works for multiple areas by creating a
#'  block diagonal adjacency matrix. However, the CARBayes package
#'  needs a dense, matrix-type object, so this WILL NOT SCALE!
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param W_list List of spatial adjacency matrices---one matrix per area.
#' @param model_family "poisson" or "gaussian".
#'    Which model family to fit. Note that Poisson models don't scale well in the
#'    multi-area setting---the spline is fit on a per-area basis.
#' @param z_offset In the "gaussian" case, the counts z are transformed as
#'    log(z + z_offset).
#' @param chains Number of MCMC chains to run.
#' @param iter Number of total iterations per MCMC chain.
#' @param warmup Burn-in for MCMC.
#' @param thin Thinning rate for MCMC.
#' @param cores Number of cores to use for MCMC.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from CARBayes::S.CARleroux, plus a time field
#'  for how long the function ran for.
#'
#' @note Requires the Matrix library.
#' @note Requires the CARBayes library.
#'
#' @import Matrix
#' @importFrom Matrix bdiag
#' @importFrom CARBayes S.CARleroux
#'
#' @export
CARBayes_Leroux_wrapper <- function(z_list,
                                    X_list,
                                    W_list,
                                    model_family = "poisson",
                                    z_offset = 0.5,
                                    chains = 4,
                                    iter = 2000,
                                    warmup = 200,
                                    thin = 2,
                                    cores = 1,
                                    library_size_list = NULL)  {
  # Start clock
  t0_cb <- Sys.time()

  # Create a single set of vectors/covariate matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
  } else {
    lib_vec <- rep(1, length(z_vec))
  }

  # Create a single adjacency matrix
  # WARNING: THIS HAS TO BE A MATRIX IN R
  # SO, THIS IS TOTALLY DENSE!!
  W <- as.matrix(Matrix::bdiag(W_list))

  # Fit model
  if ("poisson" == model_family) {
    # Content here as a placeholder
  } else if ("gaussian" == model_family) {
    z_vec <- log(z_vec + z_offset)
  } else {
    stop("Invalid model family.")
  }
  cb_out <- CARBayes::S.CARleroux(
    z_vec ~ 0 + X_mat + offset(log(lib_vec)),
    family = model_family,
    burnin = warmup,
    n.sample = iter,
    W = W,
    thin = thin,
    n.chains = chains,
    n.cores = cores,
    verbose = TRUE
  )

  # Stop clock
  t1_cb <- Sys.time()
  cb_out[["time"]] <- difftime(t1_cb, t0_cb)

  return(cb_out)
}

#'  Wrapper to fit sparse NN GP model via BRISC.
#'
#'  This function theoretically works for multiple areas by stacking everything
#'  together and ignoring differences in coordinates.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param coords_list List of coordinate matrices (x, y)---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param k Number of neighbors.
#' @param cov_type String for covariance model type.
#'    "Exp", "Sph", "Gau", and "Mat" are supported.
#' @param transform_z Boolean: log-transform z or not.
#' @param z_offset If transform_z, the counts z are transformed as
#'    log(z + z_offset).
#' @param verbose Whether to print output as BRISC runs.
#'
#' @returns The output list from BRISC, plus a time field for how long the
#'    function ran for.
#'
#' @note Requires the BRISC library.
#'
#' @importFrom BRISC BRISC_estimation
#' @export
BRISC_wrapper <- function(z_list,
                          X_list,
                          coords_list,
                          k = 15,
                          cov_type = "Exp",
                          transform_z = TRUE,
                          z_offset = 0.5,
                          verbose = FALSE) {
  # Start clock
  t0_glm <- Sys.time()

  # Stack into a single vector/matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))
  coords <- as.matrix(Reduce(rbind, coords_list))

  # Transform z
  if (transform_z) {
    z_vec <- log(z_vec + z_offset)
  }

  if ("Exp" == cov_type) {
    brisc_cov = "exponential"
  } else if ("Sph" == cov_type) {
    brisc_cov = "spherical"
  } else if ("Mat" == cov_type) {
    brisc_cov = "matern"
  } else if ("Gau" == cov_type) {
    brisc_cov = "gaussian"
  }
  # Fit LM
  b_out <- BRISC::BRISC_estimation(
    coords,
    z_vec,
    X_mat,
    n.neighbors = k,
    cov.model = brisc_cov,
    verbose = verbose
  )

  # Stop clock
  t1_glm <- Sys.time()
  b_out[["time"]] <- difftime(t1_glm, t0_glm)

  return(b_out)
}

#' Wrapper to fit a multi-sample BYM2 model with INLA (INLA).
#'
#' This function works for multiple areas, but is very slow.
#'
#' @author Florica J Constantine, florica AT berkeley.edu
#'
#' @param z_list List of vectors of observed counts---one vector per area.
#' @param X_list List of design or covariate matrices---one matrix per area.
#'    Same length and ordering as z_list.
#'    Matrices with number of rows equal to length of corresponding vector in
#'    z_list.
#' @param W_list List of spatial adjacency matrices---one matrix per area.
#' @param model_family "poisson" or "gaussian".
#'    Which model family to fit. Note that Poisson models don't scale well in the
#'    multi-area setting---the spline is fit on a per-area basis.
#' @param z_offset In the "gaussian" case, the counts z are transformed as
#'    log(z + z_offset).
#' @param num_threads Number threads to run.
#' @param compute_dic INLA parameter.
#' @param compute_waic INLA parameter.
#' @param compute_cpo INLA parameter.
#' @param verbose Whether INLA is verbose.
#' @param library_size_list A list with a vector of library sizes for each sample.
#'    Total counts for each location if there are more than one gene/measurement;
#'    if NULL, set to 1.
#'
#' @returns The output list from INLA::inla, plus a time field
#'  for how long the function ran for.
#'
#' @note Requires the INLA library.
#'
#' @import Matrix
#' @import INLA
#'
#' @importFrom INLA inla.read.graph
#' @importFrom INLA inla.getOption
#' @importFrom INLA inla.setOption
#' @importFrom INLA inla
#'
#' @export
INLA_BYM2_wrapper <- function(z_list,
                              X_list,
                              W_list,
                              model_family = "poisson",
                              z_offset = 0.5,
                              num_threads = 1,
                              compute_dic = TRUE,
                              compute_waic = TRUE,
                              compute_cpo = FALSE,
                              verbose = TRUE,
                              library_size_list = NULL) {
  t0 <- Sys.time()

  n_blocks <- length(z_list)
  if (!(length(X_list) == n_blocks && length(W_list) == n_blocks)) {
    stop("z_list, X_list, and W_list must have the same length/order.")
  }

  # Create a single set of vectors/covariate matrix
  z_vec <- as.vector(Reduce(c, z_list))
  X_mat <- as.matrix(Reduce(rbind, X_list))

  # If present, apply the library size
  if (!is.null(library_size_list)) {
    lib_vec <- as.vector(Reduce(c, library_size_list))
  } else {
    lib_vec <- rep(1, length(z_vec))
  }
  n_total <- length(z_vec)
  beta_dim <- ncol(X_mat)

  # Apply Gaussian transform if requested
  if ("gaussian" == model_family) {
    z_vec <- log(z_vec + z_offset)
    family_inla <- "gaussian"
  } else if ("poisson" == model_family) {
    family_inla <- "poisson"
  } else {
    stop("Unsupported model_family. Use 'poisson' or 'gaussian'.")
  }

  # Build data frame
  df <- as.data.frame(X_mat)
  names(df) <- paste0("X", seq_len(beta_dim))
  df$z <- z_vec
  df$lib_size <- lib_vec

  # Build block-specific region indices and graphs
  region_idx_list <- vector("list", n_blocks)
  graph_list <- vector("list", n_blocks)
  start <- 1L

  for (b in seq_len(n_blocks)) {
    nb <- length(z_list[[b]])
    idx_vec <- rep(NA_integer_, n_total)
    idx_vec[start:(start + nb - 1L)] <- seq.int(1L, nb)
    region_idx_list[[b]] <- idx_vec
    start <- start + nb

    # Convert adjacency to INLA graph
    Wb <- W_list[[b]]
    if (!is.matrix(Wb) && !inherits(Wb, "Matrix")) {
      stop(sprintf("W_list[[%d]] must be a matrix-like adjacency.", b))
    }
    graph_list[[b]] <- INLA::inla.read.graph(as.matrix(Wb))
  }

  # Add region indices to data frame
  for (b in seq_len(n_blocks)) {
    df[[paste0("region_b", b)]] <- region_idx_list[[b]]
  }

  # Fixed effects
  fixed_part <- paste0(paste0("X", seq_len(beta_dim)), collapse = " + ")

  # BYM2 terms per block
  spatial_terms <- paste0(
    "f(region_b",
    seq_len(n_blocks),
    ", model = 'bym2', graph = graph_list[[",
    seq_len(n_blocks),
    "]], scale.model = TRUE)"
  )

  # Combine into formula string
  formula_str <- paste0(
    "z ~ 0 + ",
    fixed_part,
    " + offset(log(lib_size)) + ",
    paste(spatial_terms, collapse = " + ")
  )
  formula_inla <- stats::as.formula(formula_str)

  # Evaluation environment for INLA to find graph_list
  eval_env <- new.env(parent = environment())
  assign("graph_list", graph_list, envir = eval_env)

  # Set INLA threads
  old_threads <- INLA::inla.getOption("num.threads")
  on.exit({
    INLA::inla.setOption(num.threads = old_threads)
  }, add = TRUE)
  INLA::inla.setOption(num.threads = num_threads)

  # Fit model
  if (verbose) {
    message("Fitting INLA BYM2 model with ", n_blocks, " blocks...")
  }
  fit <- eval(call("with", eval_env, {
    INLA::inla(
      formula_inla,
      family = family_inla,
      data = df,
      control.compute = list(
        dic = compute_dic,
        waic = compute_waic,
        cpo = compute_cpo
      ),
      control.fixed = list(correlation.matrix = TRUE),
      verbose = verbose
    )
  }))

  t1 <- Sys.time()
  fit$time <- difftime(t1, t0, units = "secs")

  # Store useful outputs
  fit$predicted <- as.vector(fit$summary.fitted.values$mean)
  fit$beta_hat <- as.vector(fit$summary.fixed$mean)
  names(fit$beta_hat) <- colnames(X_list[[1]])
  # Parametrize in terms of other models for consistency
  fit$gamma_hat <- fit$summary.hyperpar$mean[grepl("Phi", rownames(fit$summary.hyperpar))]
  fit$tau2_hat <- 1.0 / fit$summary.hyperpar$mean[grepl("Precision", rownames(fit$summary.hyperpar))]
  fit$beta_cov <- fit$misc$lincomb.derived.covariance.matrix
  rownames(fit$beta_cov) <- colnames(X_list[[1]])
  colnames(fit$beta_cov) <- colnames(X_list[[1]])

  return(fit)
}
