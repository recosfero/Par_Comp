#!/bin/bash
##
## optional: energy policy tags
##
#
#@ job_type = parallel
#@ class = test
#@ node = 1
#@ tasks_per_node = 1
#@ wall_clock_limit = 1:20:30
##                    1 h 20 min 30 secs
#@ job_name = mmult
#@ network.MPI = sn_all,not_shared,us
#@ initialdir = $(home)/
#@ output = mmult$(jobid).out
#@ error = mmult$(jobid).err
#@ notification=always
#@ notify_user=youremailaddress@yomama.com
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh

./mmult > testrun.txt
