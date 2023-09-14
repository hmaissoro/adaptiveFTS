#' Estimate the the standard deviation of the observation error
#'
#' @inheritParams .format_data
#' @param t \code{vector (numeric)}. Observation points at which we want to estimate the standard deviation of the error.
#'
#' @return A data.table with two columns: \code{t} and \code{sig} corresponding to the estimated standard deviation.
#' @export
#'
#' @import data.table
#'
estimate_sigma <- function(data, idcol = NULL, tcol = "tobs", ycol = "X", t = 1/2) {
  # Format data
  data <- .format_data(data = data, idcol = idcol, tcol = tcol, ycol = ycol)

  dt_sig <- data.table::rbindlist(lapply(data[, unique(id_curve)], function(idx, data, t) {
    # Compute get the index the two T_{n,k} closest to time for each X_n
    ## Get the rank i_t and j_t
    data_idx <- data[id_curve == idx]
    dif_mat <- outer(X = t, Y = data_idx[, tobs], function(ti, tobsi) abs(ti - tobsi))
    rank_mat <- apply(X = dif_mat, MARGIN = 1, function(r) rank(r))

    ## Compute [Y_{n,i_t} - Y_{n,j_t}]^2 / 2 for each curve X_n and each t
    Zn <- apply(X = rank_mat, MARGIN = 2, function(c, dt){
      first <- which(c == 1)
      second <- which(c == 2)
      Y1 <- dt[first][, X]
      Y2 <- dt[second][, X]
      Zn <- 1 / 2 * (Y1 - Y2) ** 2
    }, dt = data_idx)
    dt_res <- data.table::data.table("t" = t, "Zn" = Zn)
    rm(data_idx, Zn)
    return(dt_res)
  }, data = data, t = t))
  # Estimate sigma
  dt_sig[, sig := sqrt(mean(Zn)), by = t]
  dt_sig <- unique(dt_sig[, list(t, sig)])

  return(dt_sig)
}

