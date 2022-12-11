#timer functions tic() and toc()
#' Two clock functions.
#' Place tic() at the line in the code where you want to start timing.
#' Place toc() at the position in the code where you want to stop timing and report.
#'
#' @return
#' @export
#'
#' @examples
#' tic()
#' toc()
tic = function() {assign("timer", Sys.time(), envir=.GlobalEnv)}
toc = function() print(round(Sys.time()-timer, 2))
