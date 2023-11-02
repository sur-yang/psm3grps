#' Extending the Classical Propensity Score Matching Algorithm for Three-Group Datasets
#' @description The psm3grps R package enhances the classic propensity score matching (PSM) for datasets with three groups, aiding in balancing and comparing these groups in observational studies.
#' @param dt A data frame that includes the grouping variable along with other variables requiring balancing.
#' @param group_var A character string specifying the grouping variable. If not provided, it defaults to "treat".
#' @param times A numeric value indicating the number of iterations for the algorithm. If the outcome is not adequately balanced, executing the algorithm multiple times (more than once) is recommended.
#' @param id a character string for variable name of patient id in the data frame.
#' @importFrom dplyr %>%
#'
#' @return a data frame of balanced variables,or a list of data frames if the times > 1
#' @export
#' @examples
#'
#' rm(list=ls())
#' data("exampleData",package="psm3grps")
#'
#' # check unbalanced dataset
#' compareGroups::descrTable(treat~.,method = 4,data = exampleData)
#'
#' # run the code
#' df.psm<-psm3grps(dataset)
#'
#' # check the result
#' compareGroups::descrTable(treat~.,method = 4,data = df.psm)


psm3grps <- function(dt,group_var = "treat",times = 1,id="id") {
  dataset = dt
  group_var = group_var
  cat(paste0("We will perform PSM for ",times," time\n"))
  if (times == 1) {
    # cat(paste0("We will perform PSM for ",times," time\n"))
    dataset.final = psm3grps_v2(dataset, group_var = group_var,id=id)
    cat(paste0("PSM finished! ",nrow(dataset.final)," samples were included"))
    return(dataset.final)
  } else if (times > 1) {
    psm.list <- list()
    for (i in 1:times) {

      if (i == 1) {
        cat(paste0("PSM round ",i," starting \n"))
        psm.list[[i]] <- psm3grps_v2(dataset, group_var = group_var,id=id)
        names(psm.list)[i] = paste0("psm_", i)
        cat(paste0("PSM round ",i," finished! ",nrow(psm.list[[i]])," samples were included\n"))
      }
      if (i > 1) {
        cat(paste0("PSM round ",i," starting \n"))
        psm.list[[i]] <- psm3grps_v2(psm.list[[i - 1]], group_var = group_var,id=id)
        names(psm.list)[i] = paste0("psm_", i)
        cat(paste0("PSM round ",i," finished! ",nrow(psm.list[[i]])," samples were included\n"))
      }
    }
    return(psm.list)
  }
}
