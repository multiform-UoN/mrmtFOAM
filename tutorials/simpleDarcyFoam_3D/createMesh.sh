#!/bin/bash

# Source tutorial run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

 runApplication blockMesh  
 runApplication mirrorMesh -dict mirrorMeshDict.x -overwrite  
 mv log.mirrorMesh log.mirrorMesh.x
 runApplication mirrorMesh -dict mirrorMeshDict.y -overwrite  
 mv log.mirrorMesh log.mirrorMesh.y
