#' Extending the Classical Propensity Score Matching Algorithm for Three-Group Datasets
#' @description This is the foundational function for the primary function "psm3grps" within the package.
#' @param dt A data frame that includes the grouping variable along with other variables requiring balancing.
#' @param group_var A character string specifying the grouping variable. If not provided, it defaults to "treat".
#' @param id a character string for variable name of patient id in the data frame.
#'
#' @return a data frame of balanced variables.
#' @export
#'
#' @examples df.final<-psm3grps_v2(df)
psm3grps_v2<-function(dt,group_var="treat",id="id"){

  dataset = dt
  group_var = group_var

  ############################### prepare data #################################

  # change name
  names(dataset) <- sub(group_var,"treat",names(dataset),fixed = T )
  names(dataset) <- sub(id,"ID",names(dataset),fixed = T )
  # change variable names
  dataset$treat<-sub("group","gggg",dataset$treat,fixed = T)
  # set minimum to group1
  tab<-table(dataset$treat)%>%
    as.data.frame()%>%
    dplyr::arrange(Freq)%>%
    dplyr::mutate(Var1= as.character(Var1))%>%
    dplyr::pull(Var1)
  dataset<-dataset%>%
    dplyr::mutate(treat=sub(tab[1],"group1",treat,fixed = T))%>%
    dplyr::mutate(treat=sub(tab[2],"group2",treat,fixed = T))%>%
    dplyr::mutate(treat=sub(tab[3],"group3",treat,fixed = T))

  table(dataset$treat)

  ############################### psm process  #################################
  # 01 calculate gps
  fors <-as.formula(paste0("treat", " ~ ",
                           paste(names(dataset)[!names(dataset) %in% c("treat", "ID")],
                                 collapse = " + ")))
  fit <- multinom(fors, data = dataset, trace = F)
  Rx <- fitted(fit)
  colnames(Rx) <- paste0("GPS_", colnames(Rx))
  data <- dataset %>%
    cbind(., Rx)

  # 02 remove outliers
  min.max.Ps <- data %>%
    dplyr::group_by(treat) %>%
    dplyr::summarise_at(vars(starts_with("GPS")), list(min = min, max = max))
  Eligible <-
    (
      data$GPS_group1 >= max(min.max.Ps$GPS_group1_min) &
        data$GPS_group1 <= min(min.max.Ps$GPS_group1_max)
    ) &
    (
      data$GPS_group2 >= max(min.max.Ps$GPS_group2_min) &
        data$GPS_group2 <= min(min.max.Ps$GPS_group2_max)
    ) &
    (
      data$GPS_group3 >= max(min.max.Ps$GPS_group3_min) &
        data$GPS_group3 <= min(min.max.Ps$GPS_group3_max)
    )
  data <- dplyr::filter(data, Eligible)


  # 03 Re-GPS
  fit.E <- multinom(fors, data = data, trace = F)
  Rx.E <- fitted(fit.E)
  data[, stringr::str_detect(names(data), "GPS")] <- Rx.E

  # 04 matching
  set.seed(123)
  clustnum <- 3

  temp_group3 <- kmeans(logit(data$GPS_group3), clustnum) # group3
  data$Quint_group3 <- temp_group3$cluster


  temp_group2 <- kmeans(logit(data$GPS_group2), clustnum)
  data$Quint_group2 <- temp_group2$cluster

  group1_group2 <- filter(data, treat != "group3")
  group1_group3 <- filter(data, treat != "group2")


  match12 <-
    Matchby(
      Y = NULL,
      Tr = group1_group2$treat == "group1",
      X = logit(group1_group2$GPS_group1),
      by = group1_group2$Quint_group3,
      caliper = 0.5,
      replace = T,
      print.level = 0
    )

  match13 <- Matchby(
    Y = NULL,
    Tr = group1_group3$treat == "group1",
    X = logit(group1_group3$GPS_group1),
    by = group1_group3$Quint_group2,
    caliper = 0.5,
    replace = T,
    print.level = 0
  )

  # 05 final dataset
  m12t <- group1_group2$ID[match12$index.treated]
  m12c <- group1_group2$ID[match12$index.control]

  m13t <- group1_group3$ID[match13$index.treated]
  m13c <- group1_group3$ID[match13$index.control]

  both <- intersect(m12t, m13t)
  all <- c(both, m12c[match(both, m12t)], m13c[match(both, m13t)])

  data.matched <- data[match(all, data$ID),]
  data.unique <- data[data$ID %in% all,]

  data.index<-data.unique%>%
    dplyr::select(ID)
  names(data.index)[1]<-id

  suppressMessages(dataset.final<-raw%>%
    dplyr::inner_join(data.index))
  return(dataset.final)
}
