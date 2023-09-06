#!/bin/bash -l

# Set SCC project
#$ -P fpkmc3d

# Specify hard time limit for the job.
#   The job will be aborted if it runs longer than this time.
#   The default time is 12 hours
#$ -l h_rt=3:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m a

# Give job a name
#$ -N testsurrogate

# Combine output and error files into a single file
#$ -j y

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="


/project/fpkmc3d/julia/bin/julia $1 $2 $3 $4 $5 $6
