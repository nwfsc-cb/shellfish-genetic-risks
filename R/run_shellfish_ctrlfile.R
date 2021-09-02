#' Run shellfishrisk from control file
#'
#' @param batch the name of the batch
#' @param reps the number of reps to run
#' @param coreid no id
#' @param ctrl_file_path the file path to the control file to be used
#'
#' @return nothing
#' @export
#'
run_shellfishrisk_ctrlfile <-
  function(batch = "demo",
           reps = 1,
           coreid = 1,
           ctrl_file_path) {
    # parse control file

    # batch <- "demo"

    # ctrl_file_path <- file.path("inst", "ctrl.csv")

    ctrl_file <- readr::read_csv(ctrl_file_path,show_col_types =FALSE) # using read_csv since read.csv was reading in weird things on windows

    parsefoo <- function(x) {
      tmp <- ctrl_file$value[ctrl_file$variable == x]

      qan <- purrr::quietly(as.numeric)

      check <- qan(tmp)

      if (!is.na(check$result)) {
        out <- as.numeric(tmp)
      } else {
        out <- eval(parse(text = tmp))

      }

      return(out)
    }

    ctrl <- lapply(ctrl_file$variable, parsefoo)

    names(ctrl) <- stringr::str_trim(ctrl_file$variable)
    # run shellfish risks with controlfile parameters


    shellfishrisk(
      batch = batch,
      reps = reps,
      coreid = coreid,
      freq = ctrl$freq,
      pre_farm_years = ctrl$pre_farm_years,
      farm_years  = ctrl$farm_years ,
      post_farm_years = ctrl$post_farm_years,
      wild_N_init = ctrl$wild_N_init,
      rec_a = ctrl$rec_a,
      sd_recruit = ctrl$sd_recruit,
      numWildOffspring_par = ctrl$numWildOffspring_par,
      wild_mig_rate_L = ctrl$wild_mig_rate_L,
      wild_mig_rate_J = ctrl$wild_mig_rate_J,
      wild_mig_rate_A = ctrl$wild_mig_rate_A,
      seed_batch_size = ctrl$seed_batch_size,
      numFarmOffspring_par = ctrl$numFarmOffspring_par,
      sd_seed = ctrl$sd_seed,
      gamEsc_rec = ctrl$gamEsc_rec,
      source_idx = ctrl$source_idx,
      local_wild_idx = ctrl$local_wild_idx,
      L_escape_rate = ctrl$L_escape_rate,
      J_escape_rate = ctrl$J_escape_rate,
      A_escape_rate = ctrl$A_escape_rate,
      farm_nonB_recruit_factor = ctrl$farm_nonB_recruit_factor,
      numGameteEscapeOffspring_par = ctrl$numGameteEscapeOffspring_par,
      prob_repro_by_month = ctrl$prob_repro_by_month,
      prob_prodseed_by_month = ctrl$prob_prodseed_by_month,
      prob_L_G_escape_by_month = ctrl$prob_L_G_escape_by_month,
      prob_J_A_escape_by_month = ctrl$prob_J_A_escape_by_month,
      num_aloci_stage = ctrl$num_aloci_stage,
      num_nloci_ibd = ctrl$num_nloci_ibd,
      num_nloci_fst = ctrl$num_nloci_fst,
      good_start_AF = ctrl$good_start_AF,
      s = ctrl$s,
      z = ctrl$z,
      settle_age = ctrl$settle_age,
      repro_age = ctrl$repro_age,
      ageX = ctrl$ageX,
      max_age = ctrl$max_age,
      harvest_age = ctrl$harvest_age,
      max_harvest_age = ctrl$max_harvest_age,
      harvest_rate = ctrl$harvest_rate,
      kill_used_bstock = ctrl$kill_used_bstock,
      adult_annual_M = ctrl$adult_annual_M,
      lar_annual_M = ctrl$lar_annual_M,
      juv_annual_M = ctrl$juv_annual_M,
      var_in_mort = ctrl$var_in_mort,
      num_broodstock = ctrl$num_broodstock,
      eq_sex = ctrl$eq_sex,
      deg_add = ctrl$deg_add,
      farm_reduced_mort = ctrl$farm_reduced_mort,
      xyrs_new_bstock = ctrl$xyrs_new_bstock
    )


    # save controlfile to results folder for posterity


    file.copy(from = ctrl_file_path,
              to = file.path("results", batch, "ctrl.csv"))

    out <- paste0("Results are in ",file.path("results",batch))

  }
