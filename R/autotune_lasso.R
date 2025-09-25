#' Autotune LASSO: Fitting a linear model with fast and automatic lasso regularization
#' @name autotune_lasso
#' @param xin Matrix of predictors, with one column for each predictor.
#' Dimension will be \code{nobs} \eqn{\times} \code{nvars}; so each row is a new observation
#' and \code{nvars} > 1.
#' @param yin Vector of responses.
#' @param alpha Default 0.01, significance level of sequential F-tests 
#' used for estimation of support set.
#' @param lambda0 Numeric input, default (and recommended) value None.
#' If user provides a \eqn{\lambda_0}, autotune lasso starts finding pathwise solutions
#' from that \eqn{\lambda_0}, else sets the starting \eqn{\lambda_0 = \|x^T y\|_{\infty}/}(\code{nobs} \eqn{\times} \eqn{Var(y)})
#' @param tolerance Numeric input for an additional stopping criteria on the
#' coordinate descent when sigma is updating. Stops that coordinate
#' descent when the relative change between successive iterates of
#' coefficients' estimates is less than the \code{tolerance}.
#' @param beta_tolerance Numeric input for the stopping criteria on the
#' coordinate descent when sigma is not updating. Stops that coordinate
#' descent when the relative change between successive iterates of
#' coefficients' estimates is less than the \code{beta_tolerance}.
#' @param iter_max Maximum number of iterations of coordinate descent
#' allowed when sigma is updating
#' @param beta_iter_max Maximum number of iterations of coordinate 
#' descent allowed when sigma is not updating
#' @param standardize Logical flag for standardization of all variables in x, prior to
#' fitting the model sequence. The coefficients are always returned on the
#' original scale.
#' @param standardize_response Logical flag for demeaning the reponse variable y.
#' @param intercept Should intercept(s) be fitted (default=TRUE) or set to zero
#' (FALSE).
#' @param active Should active set selection be used (TRUE) or avoided
#' (default = FALSE). Autotune fast convergence often alleviates the speed gains
#' coming from active set selection. 
#' @param trace_it Logical input, default \code{FALSE}; if \code{TRUE}
#' prints out the iteration number while running, useful for big datasets that 
#' take a long time to fit.
#' @param l2_partial_residuals Logical flag to whether use the l2 norm of partial 
#' residuals for ordering them instead of the default l1 norm. 
#'
#' 
#' @return A list with various intermediate and final outputs produced by Autotune LASSO in its regularization path.
#' \item{sigma2_seq}{ A \code{no_of_iter} length sequence of estimates
#' noise variance \eqn{\sigma^2}. }
#' \item{beta}{ Final estimates of regression coefficients. }
#' \item{a0}{ Intercept of the fit. }
#' \item{lambda}{Final thresholding value \eqn{\lambda =
#' \lambda_0\hat\sigma^2} used in the coordinate descent after noise
#' variance estimate \eqn{\hat\sigma^2} has converged.}
#' \item{CD.path.details}{ A list of additionals details about the coordinate descent path taken by
#' Autotune lasso:}
#' \itemize{
#'    \item{\code{sorted_predictors:}}{       Decreasing ordering of predictors in terms of
#' their contribution to predicting the response values.}
#'     \item{\code{lambda0:}}{        Value of \eqn{\lambda_0} used in the Autotune LASSO.
#' Refer to the original paper for details.} 
#'    \item{\code{support_set:}}{       Final set of predictors included in the support set for
#'     sigma estimation by autotune lasso.}
#'     \item{\code{no_of_iter_before_lambda_conv:}}{      Number of iterations of coordinate descent performe
#' before noise variance estimate \eqn{\hat{\sigma}^2} converged.}
#' \item{\code{no_of_iter_after_lambda_conv:}}{     After noise variance estimate \eqn{\hat{\sigma}^2} has
#' converged, it is the number of iterations of coordinate descent required
#' for coefficients \eqn{\hat\beta} to converge.}
#' \item{\code{no_of_iterations:}}{     Total of iterations of coordinate descent implemented by Autotune Lasso}
#' \item{\code{count_sig_beta:}}{(\code{no_of_iter})-length vector containing the
#' support set sizes across the coordinate descent iterations while the
#' noise variance estimate is being updated.}
#'   }
#' 
#' 
#' 
#' 
#' @description
#' Fits a linear model via alternative optimization of penalized gaussian maximum likelihood 
#' which is a biconvex function of regression coefficients \eqn{\beta} and noise variance \eqn{\sigma^2}. 
#' Autotune LASSO's regularization path quickly picks out a good lambda for LASSO and then
#' returns the corresponding linear fit along with various attributes related to the fit.
#' 
#' 
#' @export
#'
#' @examples
#' set.seed(10)
#' n <- 80
#' p <- 400
#' s <- 5
#' snr <- 4
#' betatrue <- c(rep(1,s), rep(0, p - s))
#' x <- matrix(rnorm(n * p), ncol = p)
#' error.sd <- sqrt((betatrue %*% betatrue)/snr)
#' err <- rnorm(n, sd = error.sd)
#' y <- x %*% betatrue + err
#' y <- y - mean(y)
#' ans <- autotune_lasso(x, y, trace_it = TRUE)
#' b <- betatrue
#' # The Predictors which are actually significant:
#' which(b != 0)
#' # The Predictors which had nonzero estmated coefficients:
#' which(ans$beta != 0)
#' # Top 10 predictors X_i's in the ranking of X_i's given by autotune:
#' ans$sorted_predictors[1:10]
#' # No of significant predictors in each CD iteration when sigma_hat is allowed to vary:
#' ans$count_sig_beta
#' # Sigma estimates in each CD iteration:
#' ans$sigma2_seq
#' # Empirical noise variance:
#' var(err)
#' 
#' 
autotune_lasso

