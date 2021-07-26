#' Runs MAC
#' @details Creates a batch file on the given host in the destination directory which, once executed, 
#' runs MAC/GMC on all of the input files in the destination directory. Assumes the destination directory 
#' contains the MAC/GMC executable. 
#' @export
runMAC <- function(src_dir, dest_dir, base_fn, user, host = "cryoterm1", mac="mac4z-3_8_TR", n_jobs=25){
  #system("plink jstuckne@cryoterm1 find")
  #system("find")
  src_dir <- normalizePath(path = src_dir)
  src_files <- list.files(src_dir)
  src_mac_files <- src_files[grep(src_files, pattern = "*.mac")]
  
  ## Generate the bash file
  n_files <- length(src_mac_files)
  # print(c("n_files", n_files))
  # cpu_avail <- seq(0, n_jobs, 5)
  # cpu_avail[1] = 1
  # print(c("cpu_avail", cpu_avail))
  # res <- n_files/cpu_avail
  # n_cpus <- cpu_avail[max(which(res == floor(res)))]
  # print(c("n_cpus", n_cpus))
  
  offset <- 0L
  if (n_jobs == n_files){
    mac_lines <- sprintf("./%s %s_%s", mac, base_fn, sprintf("$((SLURM_ARRAY_TASK_ID))"))
  } else { # this might mess up if n_files < n_cpus
    indices <- seq(0, n_files-1, n_jobs)
    #print(c("indices", indices))
    mac_lines <- sprintf('j=$(printf "%s_%%09d" %s); ./%s $j', base_fn, sprintf("$((SLURM_ARRAY_TASK_ID + %d))", indices), mac)
    #mac_lines <- sprintf("./%s %s_%s", mac, base_fn, sprintf("%%09d $((SLURM_ARRAY_TASK_ID + %d))", indices))
  }
  if (length(mac_lines) == 0){ stop("Stopping; no mac files found.") }
  sbatch_lines <- c(
    "#!/bin/bash",
    "#SBATCH -p acp",
    "# Submit an array job with tasks 1-n_jobs, allow n_jobs to run simultaneously:",
    sprintf("#SBATCH --array=1-%i%%%i", n_jobs, n_jobs),
            #"#SBATCH --nodelist=acp3",
    "#SBATCH -o slurm-job-%j.output",
    "#SBATCH -e slurm-job-%j.error",
    "#SBATCH --mem_bind=local",
    "#SBATCH --export=ALL",
    "echo $SLURM_ARRAY_TASK_ID > slurm_task_id_save",
    "ulimit -s unlimited",
    "module load intel", 
    mac_lines, 
    "module rm intel",
    "false")
  readr::write_lines(x = sbatch_lines, path = file.path(src_dir, "sbatch_script.sh"))

  
  
  ## Remove any current MAC or OUT files 
  message("Deleting old files...")
  rm_files <- "find %s -maxdepth 1 -type f -name '%s' -delete"
  dquote <- function(x){ paste0('"', x, '"') }
  squote <- function(x){ paste0("'", x, "'") }
  #print(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.mac")))
  system(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.mac")), ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  #print(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.out")))
  system(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.out")), ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  #print(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.data")))
  system(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "*.data")), ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  #print(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "slurm*")), ignore.stdout = FALSE, ignore.stderr = FALSE)
  system(sprintf("plink %s@%s %s", user, host, sprintf(rm_files, dest_dir, "slurm*")), ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  
  
  ## Copy all the files over 
  message("Copying files over...")
  scp_command <- paste0("pscp ", 
                        dquote(file.path(normalizePath(src_dir), "*")), 
                        sprintf(" %s@%s:%s", user, host, dest_dir))
  #print(scp_command)
  system(scp_command, ignore.stdout = TRUE, ignore.stderr = TRUE) ## copies the batch file + the MAC files 
  
  ## Run the sbatch script
  message("Running MAC/GMC...")
  #print(sprintf("plink %s@%s %s", user, host, sprintf("cd %s && sbatch sbatch_script.sh", dest_dir)))
  system(sprintf("plink %s@%s %s", user, host, sprintf("cd %s && sbatch sbatch_script.sh", dest_dir)), ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  ## Periodically check the status
  status_check_cmd <- sprintf("plink %s@%s 'squeue'", user, host)
  status_check_cmd2 <- sprintf('plink %s@%s "ls %s"', user, host, file.path(dest_dir, "*.out"))
  pb <- txtProgressBar(min = 0, max = length(src_mac_files), style = 3)
  status_check <- function(){
    #print(status_check_cmd2)
    status_lines <- system(status_check_cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)
    #print(status_check_cmd2)
    out_files <- system(status_check_cmd2, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = FALSE)
    Sys.sleep(0.5)
    #print(length(out_files))
    status <- do.call(rbind, strsplit(status_lines, split = "\\s+"))
    setTxtProgressBar(pb, value = length(out_files))
    (!user %in% status[, 5]) && length(out_files) == length(src_mac_files)
  }
  
  ## Wait for MAC to finish 
  while(!status_check()){ Sys.sleep(time = 2) }
  close(pb)
  
  ## Retrieve the out files 
  src_dir = normalizePath(src_dir, winslash="/")
  retrieve_data_base <- sprintf("pscp %s@%s:%s %s", user, host, file.path(dest_dir, "*.data"), dquote(src_dir))
  retrieve_out_base <- sprintf("pscp %s@%s:%s %s", user, host, file.path(dest_dir, "*.out"), dquote(src_dir))
  
  message(paste0("Retrieving data files with: ", retrieve_data_base))
  #print(sprintf(retrieve_data_base))
  system(sprintf(retrieve_data_base), ignore.stdout = FALSE, ignore.stderr = FALSE)
  
  message(paste0("Retrieving data files with: ", retrieve_data_base))
  #print(sprintf(retrieve_out_base))
  system(sprintf(retrieve_out_base), ignore.stdout = TRUE, ignore.stderr = TRUE)
}
# runMAC <- function(src_dir, dest_dir, base_fn, user, pass, host = "cryoterm1", mac="mac4z-3_8_TR"){
#   
#   src_dir <- normalizePath(path = src_dir)
#   src_files <- list.files(src_dir)
#   src_mac_files <- src_files[grep(src_files, pattern = "*.mac")]
#   
#   ## Generate the bash file
#   n_files <- length(src_mac_files)
#   cpu_avail <- c(1, 5, 10, 15, 20, 25)
#   res <- n_files/c(1, 5, 10, 15, 20, 25)
#   n_cpus <- cpu_avail[max(which(res == floor(res)))]
#   offset <- 0L
#   if (n_cpus == n_files){
#     mac_lines <- sprintf("./%s %s_%s", mac, base_fn, sprintf("$((SLURM_ARRAY_TASK_ID))"))
#   } else {
#     indices <- cumsum(rep(n_cpus, (n_files/n_cpus))) - n_cpus
#     mac_lines <- sprintf("./%s %s_%s", mac, base_fn, sprintf("$((SLURM_ARRAY_TASK_ID + %d))", indices))
#   }
#   if (length(mac_lines) == 0){ stop("Stopping; no mac files found.") }
#   sbatch_lines <- c(
#     "#!/bin/bash",
#     "#SBATCH -p acp",
#     "# Submit an array job with tasks 1-25, allow 25 to run simultaneously:",
#     "#SBATCH --array=1-25%25",
#     "#SBATCH --nodelist=acp3",
#     "#SBATCH -o slurm-job-%j.output",
#     "#SBATCH -e slurm-job-%j.error",
#     "#SBATCH --mem_bind=local",
#     "#SBATCH --export=ALL",
#     "echo $SLURM_ARRAY_TASK_ID > slurm_task_id_save",
#     "ulimit -s unlimited",
#     "module load intel", 
#     mac_lines, 
#     "module rm intel",
#     "false")
#   readr::write_lines(x = sbatch_lines, path = file.path(src_dir, "sbatch_script.sh"))
#   sec_check <- paste0("sshpass -p '", eval(pass), "'")
#   
#   ## Remove any current MAC or OUT files 
#   rm_files <- "find %s -maxdepth 1 -type f -name '%s' -delete"
#   dquote <- function(x){ paste0("\"", x, "\"") }
#   squote <- function(x){ paste0("'", x, "'") }
#   system(sprintf("%s ssh %s@%s %s", sec_check, user, host, dquote(sprintf(rm_files, dest_dir, "*.mac"))), ignore.stdout = FALSE, ignore.stderr = FALSE)
#   system(sprintf("%s ssh %s@%s %s", sec_check, user, host, dquote(sprintf(rm_files, dest_dir, "*.out"))), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   system(sprintf("%s ssh %s@%s %s", sec_check, user, host, dquote(sprintf(rm_files, dest_dir, "*.data"))), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   system(sprintf("%s ssh %s@%s %s", sec_check, user, host, dquote(sprintf(rm_files, dest_dir, "slurm*"))), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   #system(sprintf("%s ssh %s@%s %s", sec_check, user, host, dquote(sprintf(rm_files, dest_dir, "*.mac"))), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   
#   
#   ## Copy all the files over 
#   message("Copying files over...")
#   scp_command <- paste0(sec_check, 
#                         " scp ", 
#                         file.path(normalizePath(src_dir), "*"), 
#                         sprintf(" %s@%s:%s", user, host, dest_dir))
#   system(scp_command, ignore.stdout = TRUE, ignore.stderr = TRUE) ## copies the batch file + the MAC files 
#   
#   ## Run the sbatch script
#   message("Running MAC/GMC...")
#   system(sprintf("%s ssh %s@%s %s", sec_check, user, host, squote(sprintf("cd %s && sbatch sbatch_script.sh", dest_dir))), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   
#   ## Periodically check the status
#   status_check_cmd <- sprintf("%s ssh %s@%s 'squeue'", sec_check, user, host)
#   status_check_cmd2 <- sprintf("%s ssh %s@%s 'ls %s'", sec_check, user, host, file.path(dest_dir, "*.out"))
#   pb <- txtProgressBar(min = 0, max = length(src_mac_files), style = 3)
#   status_check <- function(){
#     status_lines <- system(status_check_cmd, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
#     out_files <- system(status_check_cmd2, intern = TRUE, ignore.stdout = FALSE, ignore.stderr = TRUE)
#     Sys.sleep(0.5)
#     status <- do.call(rbind, strsplit(status_lines, split = "\\s+"))
#     setTxtProgressBar(pb, value = length(out_files))
#     (!user %in% status[, 5]) && length(out_files) == length(src_mac_files)
#   }
#   
#   ## Wait for MAC to finish 
#   while(!status_check()){ Sys.sleep(time = 2) }
#   close(pb)
#   
#   ## Retrieve the out files 
#   retrieve_data_base <- sprintf("scp %s@%s:%s %s", user, host, file.path(dest_dir, "*.data"), src_dir)
#   retrieve_out_base <- sprintf("scp %s@%s:%s %s", user, host, file.path(dest_dir, "*.out"), src_dir)
# 
#   message(paste0("Retrieving data files with: ", retrieve_data_base))
#   system(sprintf("%s %s", sec_check, retrieve_data_base), ignore.stdout = TRUE, ignore.stderr = TRUE)
#   
#   message(paste0("Retrieving data files with: ", retrieve_data_base))
#   system(sprintf("%s %s", sec_check, retrieve_out_base), ignore.stdout = TRUE, ignore.stderr = TRUE)
# }


