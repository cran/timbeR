#' Simulate log extraction using a Kozak (2004) variable-form taper equation that describes the taper of the tree.
#'
#' Simulate the extraction of logs from a tree from its measurements, taper function (Kozak (2004) variable-form taper equation ), trunk quality characteristics and harvest parameters such as stump height and assortments.
#'
#' @param dbh tree diameter at breast height, in centimeters.
#' @param h total tree height, in meters.
#' @param coef numerical vector containing nine coefficients of the Kozak taper equation.
#' @param p numerical value representing the first inflection point calculated in the segmented model of Max and Burkhart (1976).
#' @param assortments a data.frame with five columns and n rows, where n is the number of different wood assortments to be obtained from the tree stem. The first column must contain the names of the assortments, the second, numerical, contains the minimum diameters at the small end of the logs, in centimeters. The third column, numerical, contains the minimum lengths of the logs, in meters. The fourth column, numerical, contains the maximum lengths of the logs, in meters. The fifth column, numerical, contains the values in centimeters referring to the loss of wood due to cutting logs. The algorithm prioritizes the extraction of assortments along the stem in the order presented in the data.frame, starting from the first line, to the last.
#' @param stump_height tree cutting height, in meters. Default is 0.
#' @param downgrade if TRUE, the algorithm,from the defect_height onwards, simulates log extraction only for the last assortment in the assortments data.frame. Default is FALSE.
#' @param broken if TRUE, the algorithm will simulate the extraction of logs only up to the defect_height. Default is FALSE.
#' @param defect_height the height, in meters, from which the logs will be downgraded (if downgrade is TRUE) or log extraction simulation will be stopped (if broken is TRUE). Default is 0 for downgrade = TRUE (the whole tree is downgraded) and h * 0.5 for broken = TRUE (the tree is broken from half its original/estimated total height).
#' @param eliminate if TRUE, the algorithm does not get logs for any assortment present in the assortments table. All will be zero. Default is FALSE.
#' @param total_volume if TRUE, it adds an additional column to the results data.frame with the estimate of the total volume of the tree, from the ground height to h if broken argument is FALSE, or to defect_height if broken is TRUE. Default is FALSE.
#' @param only_vol if TRUE returns only volumes (does not return the number of logs). Default is FALSE.
#'
#' @return a list of two data.frames, the first (volumes) with the calculated volumes per assortment, and the second (logs) with the number of logs per assortment.
#'
#' @details when the `broken` and `downgrade` arguments are set to TRUE, the `defect_height` value is considered as the break height of the tree, and the entire tree is downgraded.
#'
#' @examples
#'
#' library(dplyr)
#' library(minpack.lm)
#' library(timbeR)
#'
#' tree_scaling <- tree_scaling %>%
#' mutate(did = di/dbh,
#'        hih = hi/h)
#'
#' kozak <- nlsLM(di ~ taper_kozak(dbh, h, hih, b0, b1, b2, b3, b4, b5, b6, b7, b8, p),
#'                start=list(b0=1.00,b1=.97,b2=.03,b3=.49,b4=-
#'                             0.87,b5=0.50,b6=3.88,b7=0.03,b8=-0.19, p = .1),
#'                data = tree_scaling,
#'                control = nls.lm.control(maxiter = 1000, maxfev = 2000)
#' )
#'
#' coef_kozak <- coef(kozak)[-10]
#' p_kozak <- coef(kozak)[10]
#'
#' h <- 20
#' dbh <- 25
#'
#' assortments <- data.frame(
#'   NAME = c('15-25','4-15'),
#'   SED = c(15,4),
#'   MINLENGTH = c(2.65,2),
#'   MAXLENGTH = c(2.65,4.2),
#'   LOSS = c(5,5)
#' )
#'
#' kozak_logs(dbh, h, coef_kozak, p_kozak, assortments)
#'
#' @export
kozak_logs <-
  function(dbh,
           h,
           coef,
           p,
           assortments,
           stump_height,
           downgrade,
           broken,
           defect_height,
           eliminate,
           total_volume,
           only_vol) {

    Assortment <- SED <- NULL

    if(missing(only_vol)){
      only_vol <- FALSE
    }

    if (missing(stump_height)) {
      stump_height <- 0
    }
    if (missing(downgrade)) {
      downgrade <- FALSE
    }
    if (missing(eliminate)) {
      eliminate <- FALSE
    }
    if (missing(broken)) {
      broken <- FALSE
    }
    if(!missing(defect_height) & !downgrade & !broken){
      message('The `broken` and `downgrade` arguments are FALSE. The value entered for `defect_height` will be discarded')
    }

    if(broken & downgrade){
      message('The `broken` and `downgrade` arguments are TRUE. The whole tree will be downgraded from `stump_height` to `defect_height`.')
      if(!missing(defect_height)){
        break_height <- defect_height
        defect_height <- 0}
    }

    if (missing(defect_height)) {
      if(downgrade & !broken){
        defect_height <- 0
        message('No value defined for `defect_height`. The whole stem of the tree will be downgraded.')}

      if(broken & !downgrade){
        defect_height <- h * 0.5
        break_height <- defect_height
        message('No value defined for `defect_height`. h * 0.5 will be considered as the break height of the tree.')}

      if(broken & downgrade){
        defect_height <- 0
        break_height <- h * 0.5
        message('No value defined for `defect_height`. h * 0.5 will be considered as the break height of the tree and the whole stem will be downgraded.')
      }

      if(!broken & !downgrade){
        defect_height <- h
      }

    }

    if(!exists('break_height') & !missing(defect_height) & broken){
      break_height <- defect_height
      defect_height <- 0
    }

    if(!exists('break_height')){break_height <- h}

    if (missing(total_volume)) {
      total_volume <- F
    }

    h0 <- stump_height

    colnames(assortments) <- c("Assortment", "SED", "Min Length",
                               "Max Length", "Loss")

    tab_sort <- (assortments %>% dplyr::select(Assortment, SED) %>%
                   tidyr::pivot_wider(names_from = Assortment, values_from = SED))[0,
                   ]

    tab_sort_n <- (assortments %>% dplyr::select(Assortment,
                                                 SED) %>% tidyr::pivot_wider(names_from = Assortment,
                                                                             values_from = SED))[0, ]

    if (any(c(is.na(dbh), is.na(h), eliminate))) {
      tab_sort[1, ] <- 0
      if(total_volume) tab_sort <- tab_sort %>% tibble::add_column(Total = 0)
      tab_sort_n[1, ] <- 0
    }
    else {
      for (i in seq_along(assortments$Assortment)) {
        dsort <- assortments[[i, 2]]
        cminsort <- assortments[[i, 3]]
        cmaxsort <- assortments[[i, 4]]
        psort <- assortments[[i, 4]]/100
        harv_dsort <- timbeR::kozak_hi(dbh, h, dsort, coef, p)
        if ((downgrade & !broken & i < nrow(assortments)) &
            harv_dsort > defect_height) {
          harv_dsort <- defect_height
        }

        if(broken & !downgrade & harv_dsort > break_height){
          harv_dsort <- break_height
        }

        if(broken & downgrade){
          if(i < nrow(assortments)){
            harv_dsort <- defect_height
          }
          else{
            harv_dsort <- break_height
          }
        }

        nlogs <- 0
        vsort <- 0
        while (h0 <= harv_dsort - cminsort) {
          nlogs <- nlogs + 1

          csort <- ifelse(h0 <= harv_dsort - cmaxsort, cmaxsort, harv_dsort-h0)

          h0 <- h0 + csort + psort
          vsort <- vsort + timbeR::kozak_vol(dbh, h, coef, p, h0 - psort, h0 - (psort + csort))
        }
        tab_sort[1, i] <- vsort
        tab_sort_n[1, i] <- nlogs
      }

      if (total_volume) {
        tab_sort <- tab_sort %>% tibble::add_column(Total = ifelse(broken,
                                                                   timbeR::kozak_vol(dbh, h, coef, p, break_height),
                                                                   timbeR::kozak_vol(dbh, h, coef, p, h)))
      }
    }

    if(only_vol){
      return(volumes = tab_sort)
    }
    else{
      return(list(volumes = tab_sort,
                  logs = tab_sort_n))
    }
  }

