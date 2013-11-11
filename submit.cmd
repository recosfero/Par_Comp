#!/bin/bash
##
## optional: energy policy tags
##
#
#@ job_type = parallel
#@ class = fat
#@ node = 4
###@ island_count= not needed for
#### class general
##@ total_tasks= 156
## other example
#@ tasks_per_node = 1
#@ wall_clock_limit = 1:20:30
##                    1 h 20 min 30 secs
#@ job_name = mytest
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/
#@ output = job$(jobid).out
#@ error = job$(jobid).err
#@ notification=always
#@ notify_user=youremailaddress@yomama.com
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh

./mmult > testrun.txt
