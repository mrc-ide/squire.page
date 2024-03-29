generate_diff_func <- function(){
  warning("For now please ensure all deriv files in the odin script are single lines!")
  odin_file_loc <- file.path("inst", "odin")
  #find all R files
  files <- list.files(odin_file_loc)
  files <- files[stringr::str_detect(files, ".json", negate = TRUE)]
  #check which ones are difference models already
  difference <- files[stringr::str_detect(files, "_diff.R")]
  files <- setdiff(files, difference)
  if(length(difference) > 0){
    #check for poorly named files
    poorly_named <- difference[
      !stringr::str_remove(difference, "_diff") %in% files
    ]
    if(length(poorly_named) > 0 & !all(is.na(poorly_named))){
      warning(paste0(
        "The following odin files are named as difference models (with _diff.R) but do not have a corresponding ode model:\n",
        paste0(poorly_named, collapse = ", "),
        ".\nThese will be ignored."
      ))
    }
    difference <- setdiff(difference, poorly_named)
    message(
      "\nRemoving existing difference models"
    )
    unlink(file.path(odin_file_loc, difference))
  }
  #generate difference model code
  for(file in files){
    ode_script <- readLines(file.path(odin_file_loc, file))
    #replace derived() with updates
    diff_script <- purrr::map(ode_script, function(line){
      if(stringr::str_detect(line, "deriv") & stringr::str_sub(stringr::str_trim(line), 1, 1) != "#"){
        compartment <- stringr::str_split(line, "<-")[[1]][1]
        compartment <- stringr::str_split(compartment, "\\(")[[1]][2]
        compartment <- stringr::str_remove(compartment, "\\) ")
        #add indexes to compartment
        index_compartment <- stringr::str_split(compartment, "[\\[\\],]")[[1]]
        index_compartment[3] <- " j"
        index_compartment[2] <- "i"
        index_compartment <- paste0(index_compartment[1], "[", index_compartment[2], ",", index_compartment[3], "]")
        change_terms <- stringr::str_split(line, "<-")[[1]][2]
        line <- paste0("update(", compartment, ") <- ", index_compartment, " + dt*(", change_terms, ")")
      } else if(stringr::str_detect(line, "interpolate") & stringr::str_sub(stringr::str_trim(line), 1, 1) != "#"){
        #increase time frame of interpolations
        tt_parameter <- stringr::str_split(line, "interpolate\\(")[[1]][2]
        tt_parameter <- stringr::str_split(tt_parameter, ",")[[1]][1]
        tt_parameter <- stringr::str_trim(tt_parameter)
        extra_line <- paste0(
          tt_parameter, "_dt[] <- ", tt_parameter, "[i] / dt"
        )
        dim_line <- paste0(
          "dim(", tt_parameter, "_dt) <- length(", tt_parameter, ")"
        )
        new_line <- stringr::str_replace(line, tt_parameter, paste0(tt_parameter, "_dt"))
        line <- c(extra_line, dim_line, new_line)
      }
      line
    }) %>% unlist()
    #sort out dt and time
    time_lines <- which(stringr::str_detect(diff_script, "time") & stringr::str_detect(diff_script, "#", negate = TRUE))
    dt_line <- "dt <- user()"
    #replace first line and remove the rest
    diff_script[time_lines[1]] <- dt_line
    diff_script <- diff_script[-time_lines[-1]]
    #add header
    diff_script <- c("#### DO NOT EDIT BY HAND!!! ####",
                     "#### File generated by inst/r_gen/create_difference_model.R ####",
                     "#### After editing the other odin files, recreate this with source('inst/r_gen/create_difference_model.R' ####",
                     diff_script)
    #get new name
    new_file <- stringr::str_remove(file, ".R") %>%
      paste0("_diff.R")
    writeLines(diff_script, file.path(odin_file_loc, new_file))
  }
}
generate_diff_func()
#generate odin code
odin::odin_package(here::here())
