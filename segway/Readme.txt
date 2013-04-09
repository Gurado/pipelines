This pipeline is for running segway v1.1.0 on the cluster and utilizes Wiggler to get all data tracks onto the same page.

run_segway.sh comes in 6 run modes:
* -1 step 1 - collect the bam data from gagri
* -2 step 2 - groom the bam data into bedGraph format
* -3 step 3 - put the data as tracks into a genomedata archive
* -4 step 4 - train the segway model
* -5 step 5 - predict using the segway model
* -6 step 6 - evaluate the output using segtools

to perform any of these steps you need a config script that holds all relevant
parametes and locations. An example is provided in config_example.sh.

Important note:
run_segway.sh generates scripts that can be submited to the cluster.
These scripts are not actually submitted/run in default mode. For that to
happen you need to start segway in armed mode  (i.e. called with the parameter -a).

Segway needs to run on the head node as it submits jobs itself and the compute
nodes are not permitted to submit new jobs. Therefore it is advised to run
every of the above scripts in a screen. 
To name and get back a screen of name segway, type:
screen -R segway

