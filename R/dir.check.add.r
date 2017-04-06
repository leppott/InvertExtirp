#' A function to check for a directory and add if it doesn't exist.
#'
#' Used in several functions before saving output.  Ensures that the code doesn't fail by making sure the directory exists.
#'
#' @param Dir.main Main/root directory; default is getwd()
#' @param Dir.sub sub directory to check for and add if necessary.
#' @return Nothing to return
#' @keywords directory, check, add, create
#' @examples
#' Dir.main <- getwd()
#' Dir.sub <- "output"
#' dir.check.add(Dir.main, Dir.sub)
#' @export
dir.check.add <- function(Dir.main=getwd(), Dir.sub) {##FUNCTION.dir.check.add.START
  #
  Dir2Check <- file.path(Dir.main, Dir.sub)
  #
  if(!file.exists(Dir2Check)){
    dir.create(Dir2Check)
  }
  #
}##FUNCTION.dir.check.add.END
