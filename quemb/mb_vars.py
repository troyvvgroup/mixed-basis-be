# Set relevant variables for the mixed-basis embedding code

# Set the QuEmb scratch path
# This will be used as the WorkDir.from_environment(user_defined_root=scratch_path)
# This sets the root of the scratch, while WorkDir will generate the directory
# name from Slurm or the PID
# Else: jobs will be run in the default, /tmp, and ERIs deleted at end automatically
scratch_path=''

# Verbose levels: Quiet (0), Error (1), Warn (2), Note (3), Info (4), Debug (5)
log_verbose_level=0
