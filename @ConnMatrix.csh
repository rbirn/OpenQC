#!/bin/tcsh

#--------------------------------------------------------
# Script to compute connectivity matrix from fMRI dataset
#--------------------------------------------------------
# Author: R.M. Birn
# Date:   2016-11-15
#--------------------------------------------------------

if ($#argv < 3) then
   echo "Usage: "
   echo "   @ConnMatrix.csh <file> <SeedMask> <OutputPath> <OutputFile>"
   echo "      <file>       = input 3d+t dataset (AFNI or NIfTI format)"
   echo "      <SeedMask>   = mask of ROIs (unique value for each ROI)"
   echo "      <OutputFile> = prefix for output file"
   echo "" 
   exit()
endif

set file     = $argv[1]
set SeedMask = $argv[2]
set file_out = $argv[3]

#--- determine paths and filenames ---
set filename = `basename $file`
set InPath   = `dirname $file`
set output   = `basename $file_out`
set OutPath  = `dirname $file_out`

set TMP = $OutPath/tmpX

MATRIX:
   echo "Extracting seeds from $file to $TMP.$filename.1D ..."
   3dROIstats -mask $SeedMask -mask_f2short -quiet $file > $TMP.$filename.1D

   echo "Computing connectivity matrix ..."
   1ddot -terse -okzero -dem $TMP.$filename.1D > $TMP.cc.$filename.1D
   
   set NS = `wc -l < $TMP.cc.$filename.1D`
   @ NS1 = $NS - 1

   #--- convert to single vector ---
   echo "Converting to vector ..."
   rm -f $TMP.vec.cc.$filename.1D
   foreach ns (`count -digits 3 0 $NS1`)
      1dcat $TMP.cc.$filename.1D'['$ns']' >> $TMP.vec.cc.$filename.1D
   end

   #--- compute indexes ---
   if (! -e $TMP.indexes.$NS.1D) then
   echo "Computing indexes ..."
   3dUndump -dimen $NS $NS 1 -prefix $TMP.indexes.$NS 
   3dmaskdump $TMP.indexes.$NS+orig | awk '{print $1 " " $2 " " $3}' > $TMP.indexes.$NS.1D
   rm -f $TMP.indexes.$NS+orig.*
   endif

   #--- concat indexes and connectivity vals ---
   1dcat $TMP.indexes.$NS.1D $TMP.vec.cc.$filename.1D > $TMP.indx.cc.$filename.1D

   #--- Create AFNI dataset for matrix ---
   echo "Creating connectivity matrix ..." 
   3dUndump -dimen $NS $NS 1 -datum float -prefix $OutPath/$output $TMP.indx.cc.$filename.1D

CLEANUP:
   echo "Cleaning temporary files ..." 
   rm -f $TMP.$filename.1D 
   rm -f $TMP.cc.$filename.1D
   rm -f $TMP.vec.cc.$filename.1D
   rm -f $TMP.indexes.$NS.1D
   rm -f $TMP.indx.cc.$filename.1D

