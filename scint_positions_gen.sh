#!/bin/bash
#

umask 0002
verbose=T

# Grid run job iterator
NEWPROCESS=`echo "($JOBSUBJOBSECTION - 1) " | bc`

# Arrays with the (x, y, z) for the scintillation points 
#------------------------------------------------------------
x_=( 10. 30. 50. 70. 90. 110. 130. 150. 170. 190. 195.)
y_=( 40. 50. 60. 90. 100. 110. 170. 180. 190. )
z_=( 285. 345. 405. 465. )
# Some variable definitions
counter=0
N=1
N_photons_per_job=100000
# loop over our chosen points
for (( i = 0; i < 11; i++ )) ### loop in x coordinate ###
do
   
    for (( j = 0 ; j < 9; j++ )) ### loop in y coordinate ###
    do
	
        for (( k = 0 ; k < 4; k++ )) ### loop in z coordinate ###
	do
	    for (( l = 0; l < $N; l++ )) ### we repeat xN the number of jobs with the same (x,y,z) ###
	    do            
		x=${x_[$i]}
		y=${y_[$j]}
		z=${z_[$k]}
	    
		if [ "$counter" = "${NEWPROCESS}" ]
		then
		    echo "  x      y     z      t      dx      dy      dz      dt      p         dp        n" >> myLightSourceSteering_$counter.txt
		    echo "  $x    $y     $z     0.0     0.0     0.0     0.0     0.0     9.69     0.25     ${N_photons_per_job}" >> myLightSourceSteering_$counter.txt
		fi
		counter=`expr $counter + 1`
	    done
	done
	
    done
    
    echo "" #### print the new line ###
done

#------------------------------------------------------------

chmod 755 myLightSourceSteering_${NEWPROCESS}.txt
mv myLightSourceSteering_${NEWPROCESS}.txt myLightSourceSteering.txt

echo $FHICL_FILE_PATH "fhicl file path"
echo "services.TFileService.fileName: \"${NEWPROCESS}_scint.root\"" >> $FCL
