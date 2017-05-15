#
# female_fd_path <- "D:/Dropbox/1 Projects/Academic/2016 genetic covariance invasions/data/fecundity_dispersal_raw.csv"
# female_f_path  <- "D:/Dropbox/1 Projects/Academic/2016 genetic covariance invasions/data/fecundity_no_dispersal_raw.csv"
# dates_path <- "D:/Dropbox/1 Projects/Academic/2016 genetic covariance invasions/data/experiment_dates.csv"
# male_d_path <- "D:/Dropbox/1 Projects/Academic/2016 genetic covariance invasions/data/male_dispersal.csv"

#' Wrangle raw data into a form amenable to analysis.
#'
#' @param female_fd_path A string that provides the path to the file that
#'     contains data for female fecundity and dispersal.
#' @param female_f_path A string that provides the path to the file that
#'     contains data for female fecundity data ONLY. These beetles did not
#'     disperse.
#' @param male_d_path A string that provides the path to the file that contains
#'     data for male disperal. This file also contains information on dispersal
#'     arrays.
#' @param dates_path A string that provides the path to the file that contains
#'     data for mating, dispersal, and freezing dates.
#' @return A data frame of all the data, cleaned and compiled.
#' @export
wrangle_beetle_data = function(female_fd_path, female_f_path,
                             male_d_path, dates_path) {
  # load dates of mating, dispersal, and freezing
  dates <- read.csv(dates_path)

  # load dataset that has density-dependent growth AND dispersal
  # (***only contains data for females***)
  female_fd <- read.csv(female_fd_path)
  colnames(female_fd) <- c("sire", "dam", "m_date",
                           "id", "dist", "beans", "d_date", "f", "m",
                           "f_date")

  # load dataset that has density-dependent growth only
  # (***only contains data for females***)
  female_f <- read.csv(female_f_path)
  female_f$id     <- NA
  female_f$dist   <- NA
  female_f$m_date <- NA
  female_f$d_date <- NA
  female_f$f_date <- NA

  # merge female datasets
  female_data <- rbind(female_fd, female_f)

  # load male dispersal data (also contains dispersal array data)
  male_d <- read.csv(male_d_path)

  # reshape male data from wide to long form
  male_data <-  tidyr::gather(male_d, id, dist, -c(sire, dam, array, start, stop))

  # format male ID to be numeric
  male_data$id <- as.numeric(gsub("m", "", male_data$id))

  # unique sire/dam combinations for all families in the experiment
  all_fams  <- unique(female_data[, c("sire", "dam")])

  # unique sire/dam combinations for all families that dispersed
  disp_fams <- unique(male_d[,c("sire", "dam", "array", "start", "stop")])

  # array data for all families in the experiment
  fam_data <- dplyr::left_join(dplyr::left_join(all_fams,
                                                disp_fams,
                                                by = c("sire", "dam")),
                               dates, by = c("sire", "dam"))

  # order fam_data
  fam_data <- fam_data[ order(fam_data$sire, fam_data$dam), ]

  # initiate clean_data
  clean_data <- data.frame(sire   = NA,
                           dam    = NA,
                           sex    = NA,
                           id     = NA,
                           dist   = NA,
                           beans  = NA,
                           f      = NA,
                           m      = NA,
                           t      = NA,
                           m_date = NA,
                           d_date = NA,
                           f_date = NA)[-1,]

  # update missing data that results from combining datasets
  for (i in 1:nrow(fam_data)) {
    # subset data
    tmp_f <- female_data[female_data$sire == fam_data$sire[i] &
                           female_data$dam  == fam_data$dam[i], ]

    tmp_m <- male_data[male_data$sire == fam_data$sire[i] &
                         male_data$dam  == fam_data$dam[i], ]

    # set female ids if missing
    if (length(na.omit(tmp_f$id))) {
      # continue numbering ids from the highest numbered id...
      new_id_range <- c(max(tmp_f$id, na.rm = T) + 1,
                        max(tmp_f$id, na.rm = T) + length(which(is.na(tmp_f$id))))
      tmp_f$id[is.na(tmp_f$id)] <- new_id_range[1]:new_id_range[2]
    } else {
      # ... or start ids from 1
      tmp_f$id[is.na(tmp_f$id)] <- which(is.na(tmp_f$id))
    }

    # clean female dataframe
    clean_f <- data.frame(sire  = fam_data$sire[i],
                          dam   = fam_data$dam[i],
                          sex   = "f",
                          id    = tmp_f$id,
                          dist  = tmp_f$dist,
                          beans = tmp_f$beans,
                          f     = tmp_f$f,
                          m     = tmp_f$m,
                          t     = tmp_f$f + tmp_f$m,
                          m_date = fam_data$m_date[i],
                          d_date = fam_data$d_date[i],
                          f_date = fam_data$f_date[i])
    if (nrow(tmp_m) > 0) {

      #clean male dataframe
      clean_m <- data.frame(sire  = fam_data$sire[i],
                            dam   = fam_data$dam[i],
                            sex   = "m",
                            id    = tmp_m$id,
                            dist  = tmp_m$dist,
                            beans = NA,
                            f     = NA,
                            m     = NA,
                            t     = NA,
                            m_date = fam_data$m_date[i],
                            d_date = fam_data$d_date[i],
                            f_date = fam_data$f_date[i])
    } else {
      clean_m <- clean_f[1, ][-1, ]
    }

    # append tmp to clean_data
    clean_data <- rbind(clean_data, rbind(clean_f, clean_m))
  }

  # each array was split into 3 sub-arrays, so dispersal distances need to be
  # centered on 0 relative to the sub-array they were in
  sub_low <- which(-6 <= clean_data$dist & clean_data$dist <= 41)
  sub_mid <- which(42 <= clean_data$dist & clean_data$dist <= 89)
  sub_hi  <- which(90 <= clean_data$dist & clean_data$dist <= 137)

  clean_data[sub_low, ]$dist <- clean_data[sub_low, ]$dist - 18
  clean_data[sub_mid, ]$dist <- clean_data[sub_mid, ]$dist - 66
  clean_data[sub_hi,  ]$dist <- clean_data[sub_hi,  ]$dist - 114

  return(clean_data)
}
