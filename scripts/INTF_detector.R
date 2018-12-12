#!/usr/bin/env Rscript

# INTF_detector.R 
# COMMON MODULE - Detection of constant and continuous interferences in RAW weather radar data
# P. Altube Vazquez - March 2018

###################################################################################################

# Functions
# library(here) # retrieve current script path
# source(paste(here(), fns_detector.R, sep="/"))
source("/home/pav/progR/SunINTFCal/Scripts/functions/fns_detector.R")

hdr_line_num <- 15
acc_period <- "1 year"

# command line arguments
# opt <- parse_cl()

# TESTING
opt <- list(date=as.Date("2018-02-26"),
            rad="CDV", root="/home/pav/progR/SunINTFCal", set="SunINTFCal.ini")

# User settings from settings file
set <- get_settings(paste(opt$root, opt$set, sep="/"))

# Paths
path_sun <- create_paths(root=opt$root, radar=opt$rad, log="LOG", module="SOLAR", out="INTF")
path_emi <- create_paths(root=opt$root, radar=opt$rad, log="LOG", module="EMITTER", out="INTF")

# OUTPUT FILES
log_file <- paste("SunINTFCal_", format(opt$date, format="%y%m%d"), "_", 
                  opt$rad, "_LOG.txt", sep="")
ofile_emi <- paste("INTF_EMITTER_", opt$rad, ".txt", sep="")
ofile_sun <- paste("INTF_SOLAR_", opt$rad, ".txt", sep="")
ofile_day_emi <- paste("INTF_EMITTER_", format(opt$date, format="%y%m%d"), "_", 
                  opt$rad, ".txt", sep="")
ofile_day_sun <- paste("INTF_SOLAR_", format(opt$date, format="%y%m%d"), "_", 
                      opt$rad, ".txt", sep="")

# Input folder and raw files
in_path <- paste(opt$root, "tmp", sep="/")

# Accumulation files
acc_file_sun <- paste(path_sun$module, ofile_sun, sep="/") 
acc_file_emi <- paste(path_emi$module, ofile_emi, sep="/") 

# LOG file messages
log_msg <- list("ini"=c(paste("DATE:", opt$date), paste("RADAR:", opt$rad)),
                "exe"=c(paste("EXECUTION: INTF_detector.R"), as.character(Sys.time())),
                "set"=set$msg,
                "out"="NO OUTPUT FILE(S) GENERATED", 
                "inp"="NO INPUT FILE(S) FOUND")

###################################################################################################

# Read accumulation files and remove oldest entries (defined in acc_period)
old_date <- seq(opt$date, length=2, by=paste0("-", acc_period))[2]
intfs_acc_emi <- rm_intfs_by_date(date=old_date, file=acc_file_emi, 
                                  hdr_lines=hdr_line_num, type="EMITTER")
intfs_acc_sun <- rm_intfs_by_date(date=old_date, file=acc_file_sun, hdr_lines=hdr_line_num)

# Find all raw files in folder
raw_files <- find_files(in_path, pattern="*.RAW*")

if (length(raw_files)!=0){
  
  intfs_all <- init_intfs()
  # Update log messages
  log_msg$out <- paste("OUTPUT file(s): ", paste(path_sun$out, ofile_day_sun, sep="/"), 
                       paste(path_emi$out, ofile_day_emi, sep="/"), sep="\n")
  log_msg$inp <- "INTERFERENCES DETECTED:"

  for (rw in raw_files){

    # Read file and find interferences
    intfs <- intf_detector(raw_file=rw, radar=opt$rad, 
                           min_range_clt=set$min_range_clt, min_range_pow=set$min_range_pow, 
                           min_gate_frac=set$min_gate_frac, max_sd_pow=set$max_sd_pow, 
                           max_delta_az=set$max_delta_az, max_delta_el=set$max_delta_el, 
                           d_sun=set$d_sun, gas_att_2way=set$gas_att_2way, atm_eq_height=set$atm_eq_height,
                           k_cnt=set$k_cnt, refractivity_0=set$refractivity_0, r_earth=set$r_earth)
    
    intfs_all <- rbind(intfs_all, intfs$intfs)
    
    # Files for which any interferences have been detected are indicated in LOG message/file
    if (nrow(intfs$intfs)!=0){
      # Update log message
      log_msg$inp <- c(log_msg$inp, paste(gsub(pattern=in_path, replacement="", x=rw), 
                                          ": ", nrow(intfs$intfs), " interference(s) found", sep=""))
    }

  }
  
  # Power scaling factor due to scanning
  scan_factor <- round(pow_scaling(bw=(intfs$meta$beam_width_h+intfs$meta$beam_width_v)/2, 
                                   d_ray=intfs$meta$radial_res, d_sun=set$d_sun), 3)
  # Output file headers
  hdr_sun_day <- gen_intf_hdr(opt$rad, metadata=intfs$meta, title="SOLAR INTERFERENCES", 
                              xtra_hdr=c(paste("pow_sc_fact:", scan_factor, "[dB]"), set$msg))
  hdr_emi_day <- gen_intf_hdr(opt$rad, metadata=intfs$meta, title="EMITTER INTERFERENCES", 
                              xtra_hdr=c(set$msg))
  
  # Format output
  intfs_all <- format_df_num(intfs_all, ndec=5)

  intfs_day_sun <- intfs_all[intfs_all$label=="SOLAR", !names(intfs_all) %in% c("label")]
  intfs_day_emi <- intfs_all[intfs_all$label=="EMITTER", 
                             !names(intfs_all) %in% c("x", "y", "pow_corr", "l_gas", "l_scan", "label")]

  # Write daily output files
  write_df(intfs_day_sun, ofile_day_sun, path_sun$out, hdr=hdr_sun_day, create=TRUE)
  write_df(intfs_day_emi, ofile_day_emi, path_emi$out, hdr=hdr_emi_day, create=TRUE)
  
  # Add detected interferences to accumulation data
  if (nrow(intfs_day_sun)!=0){
    rbind(intfs_acc_emi, intfs_day_emi)
  }
  if (nrow(intfs_day_sun)!=0){
    rbind(intfs_acc_sun, intfs_day_sun)
  }
  
}

hdr_sun_acc <- gen_intf_hdr(opt$rad, title="SOLAR INTERFERENCES (dynamic accumulation)")
hdr_emi_acc <- gen_intf_hdr(opt$rad, title="EMITTER INTERFERENCES (dynamic accumulation)")

# Write LOG file and EMITTER interference accumulation file
write_log(log_msg, log_file, path_sun$log)
write_df(intfs_acc_emi, ofile_emi, path_emi$module, hdr=hdr_emi_acc, create=TRUE, append=FALSE)
write_df(intfs_acc_sun, ofile_sun, path_sun$module, hdr=hdr_sun_acc, create=TRUE, append=FALSE)