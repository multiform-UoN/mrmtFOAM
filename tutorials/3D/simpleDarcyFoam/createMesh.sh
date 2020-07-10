#!/bin/bash


 ofprog="createmesh"

    filelog=$(echo $ofprog".log")
    filerr=$(echo $ofprog".err")

 echo "Starting calculation with $ofprog"

 blockMesh  > $filelog  2> $filerr
 mirrorMesh -dict mirrorMeshDict.x -overwrite  >> $filelog  2>> $filerr
 mirrorMesh -dict mirrorMeshDict.y -overwrite  >> $filelog  2>> $filerr
