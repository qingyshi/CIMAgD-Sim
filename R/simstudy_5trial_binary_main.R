# Appendix B simulation driver: 5trial_binary_main.
# Default: smoke test (2 iterations, scenarios 1 and 16).
# Full run: Rscript simstudy_5trial_binary_main.R --full --workers=4
# See SIMULATION_RUN_GUIDE.txt; keep simulation_common.R in this directory.

.sim_source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.sim_script <- if (length(.sim_source_files)) tail(.sim_source_files, 1L)[[1L]] else {
  .sim_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(.sim_arg)) stop("Use source() or Rscript to run this driver.")
  sub("^--file=", "", .sim_arg[1L])
}
.sim_script_dir <- dirname(normalizePath(.sim_script))
source(file.path(.sim_script_dir, "simulation_common.R"), local = TRUE)

# Sourcing defines the setting and functions without starting a simulation.
setting <- "5trial_binary_main"
scenarios <- sim_grid(setting)
if (sys.nframe() == 0L) sim_run(setting, .sim_script_dir)
