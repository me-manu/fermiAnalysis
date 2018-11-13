#!/bin/bash

script="/u/gl/mmeyer/projects/VarFSRQ/varFSRQ/manage.py"
yaml=(  "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/srcconf_weekly_daily_orbit.yaml"
        "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/energy_bins/300-1000MeV.yaml"
        "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/energy_bins/100-300MeV.yaml")
yaml=(  "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/energy_bins/300-1000MeV.yaml" )
#yaml=(  "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/energy_bins/above1000MeV.yaml" )
#yaml=(  "/u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/srcconf_weekly_daily_orbit.yaml" )
src="none"
#src="3C454.3"
concur=10
option=" --conf /u/gl/mmeyer/projects/VarFSRQ/config/brightest_fsrqs/template.yaml --interval orbit --state avgspec --dry 0 --lightcurve 1 --concurrent $concur --select_src $src --srcconf "

i=0
maxrepeat=10
sleeptime=15
initsleep=600
maxjobsremain=100

kinit -r 1d -l 1d


# initial submission
for y in "${yaml[@]}"
    do
        echo "python $script $option $y"
        python $script $option $y
        sleep 0.5
    done
echo "Sleeping for $initsleep"
sleep $initsleep

br=`bjobs -r | grep RUN | wc -l`
bp=`bjobs -p | grep PEN | wc -l`

# repitions
while [ $i -lt $maxrepeat ]
    do
        echo $i
        while [ $(($bp + $br)) -gt $maxjobsremain ]
            do

                date
                echo "RUN $br PEND $bp"
                echo "Sleeping for $sleeptime seconds"
                sleep $sleeptime
                br=`bjobs -r | grep RUN | wc -l`
                bp=`bjobs -p | grep PEN | wc -l`
            done

        for y in "${yaml[@]}"
            do
                echo "python $script $option $y"
                python $script $option $y
                sleep 0.5
            done
        i=$[$i+1]
        kinit -R
        aklog
    done
