#!/bin/tcsh

#---
# QC stats at group level
#---
# Author: R.M. Birn
#---

#---------------------------
# Paths:
#---------------------------
set RootPath   = ~/src/OpenQC
set ScriptPath = $RootPath/scripts
set MergePath  = $RootPath/MERGE
set RoiPath    = $RootPath/ROIs 
set ParcPath   = $RootPath/ROIs/Parcels

#---------------------------
# Parameters:
#---------------------------
set ign             = 4   #number of time points to ignore
set dxyz            = 2   #resampled resolution (in template space)
set fwhm            = 8   #full width at half maximum - target for smoothing
set mot_thresh      = 0.4 #motion censoring threshold
set anaticor_fwhm   = 40  #FWHM for anaticor (not yet implemented)

#---------------------------
# Options:
#---------------------------
set CENSOR_PREV     = 1		#censor previous time point as well as current
set GRAYPLOT        = 0		#produce Grayplot
set DO_PEARSON      = 1		#compute Pearson's correlation
set CREATE_DIRFILE  = 0#1	#create directory file (include only dirs where output exists)
set START_CLEAN     = 0#1	#clean temporary files before running

set DIR_FILE = $ScriptPath/dirs.ok.txt

if ($CENSOR_PREV) then
   set cens_prev_label  = "_wPrev"
   set cens_prev_string = "-censor_prev_TR"
else
   set cens_prev_label  = ""
   set cens_prev_string = ""
endif

#--- create list of "good" directories ---
if ($CREATE_DIRFILE) then
   set n = 0
   set DIR_FILE_good = $ScriptPath/dirs.ok.good_${mot_thresh}.txt
   rm -f $DIR_FILE_good
   foreach dir (`cat $DIR_FILE`)
      set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
      set file     = ${dir}_ses-01_task-rest_run-01_bold
      set fileF    = $file.xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni
      set CCmatrix = CCmatrix.GordonCS.$fileF.nii.gz 
      if (! -e $DataPath/$CCmatrix) then
	 echo "Missing $DataPath/$CCmatrix. Skipping."
	 continue
      endif
      @ n ++
      echo "$dir" >> $DIR_FILE_good
   end
   exit()
endif

#--- Uncomment to select the directory file to use ---
#set DIR_FILE_use = $ScriptPath/dirs.ok.good_0.2.txt
#set group_label  = N121

set DIR_FILE_use = $ScriptPath/dirs.ok.good_0.4.txt
set group_label  = N136

#set DIR_FILE_use = $ScriptPath/dirs.ok.good_1.0.txt
#set group_label  = N137

#--- Uncomment to jump ahead ---
#----------------------
#goto SIMILARITY_DISTANCE_TO_MEAN #SIMILARITY_TO_MEAN #DISTANCE_DEPENDENCE #
#----------------------

set cov_enorm = $MergePath/cov.enorm.$group_label.${mot_thresh}.1D
rm -f $cov_enorm

set ListCC = $MergePath/list.CCmatrix.$group_label.${mot_thresh}.txt
rm -f $ListCC


set n = 0
foreach dir (`cat $DIR_FILE_use`) #($argv) #

   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
   
   cd $DataPath/
   
   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold
   
   set anat0    = $anat
   set anatF    = anatQQ.$dir.2mm
   set anatA    = anatQQ.$dir.aff
   set anatW    = anatQQ.${dir}_WARP
   
   set file0    = $file
   set file1    = $file.xtoru.aw
   set fileF    = $file.xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni  #<-- with GSR, filter
   #set fileF    = $file.xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni  #<-- with GSR, no filter
   #set fileF    = $file.xtoru.aw.Nwcdm12fXz${mot_thresh}S8.mni  #<-- no GSR, filter
   
   set CCmatrix = CCmatrix.GordonCS.$fileF.Xk.nii.gz 	#<-- Xk=with censored time points removed
   #set CCmatrix = CCmatrix.GordonCS.$fileF.nii.gz 	
   
   set motfile          = mot.$file0.1D
   set CensorFilePrefix = mot_censor.$file0.X${mot_thresh}${cens_prev_label}
   set CensorFile       = ${CensorFilePrefix}_censor.1D
   set EnormFile        = ${CensorFilePrefix}_enorm.1D
   
   if (! -e $DataPath/$CCmatrix) then
      echo "Missing $DataPath/$CCmatrix. Skipping."
      continue
   endif
   
   @ n ++
   
   set mean_enorm = `3dTstat -mean -prefix - $EnormFile\'`
   echo $mean_enorm >> $cov_enorm
   
   #--- Fisher-Z transform ---
   if (! -e Z.$CCmatrix) then
   3dcalc \
      -datum float \
      -a $DataPath/$CCmatrix \
      -expr "log((1+a)/(1-a))" \
      -prefix $DataPath/Z.$CCmatrix
   endif
   
   echo "$DataPath/Z.$CCmatrix" >> $ListCC
   
end #dir

echo "N: $n"

cd $MergePath

#========================================
QC_FC:
#----------------------------------------
# Correlation between quality control (mean_enorm) and functional connectivity
#----------------------------------------
set corF = xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni.Xk	    #<-- with GSR, filter
#set corF = xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni.Xk	    #<-- with GSR, no filter
#set corF = xtoru.aw.Nwcdm12fXz${mot_thresh}S8.mni.Xk	    #<-- no GSR, filter

set CCout = $group_label.Z.CCmatrix.GordonCS.$corF.nii.gz

if ($START_CLEAN) then
   rm -f CAT.$CCout
   rm -f CC_pearson.mean_enorm.$CCout
   rm -f CC_spearman.mean_enorm.$CCout
endif

if (! -e CAT.$CCout) then
   3dTcat -prefix CAT.$CCout `cat $ListCC`
endif

if (! -e CC_pearson.mean_enorm.$CCout) then
3dTcorr1D \
   -prefix CC_pearson.mean_enorm.$CCout \
   CAT.$CCout \
   $cov_enorm
endif

if (! -e CC_spearman.mean_enorm.$CCout) then
3dTcorr1D \
   -spearman \
   -prefix CC_spearman.mean_enorm.$CCout \
   CAT.$CCout \
   $cov_enorm
endif


#========================================
HISTOGRAMS:
#----------------------------------------
# Compute histograms of QC_FC 
#----------------------------------------

set covar_label = mean_enorm
set cc_matrix   = $CCout
set CatFile     = CAT.$CCout

cd $MergePath/

if ($START_CLEAN) then
   rm -f hist_enorm.$group_label.1D
   rm -f histCCpearson.$covar_label.$cc_matrix.1D
   rm -f histCCspearman.$covar_label.$cc_matrix.1D
endif

#--- histogram of motion ---
if (! -e hist_enorm.$group_label.1D) then
   echo "Computing histogram for enorm"
   3dhistog \
      -nbin 20 \
      -pdf \
      -dind 0 \
      -min 0 \
      -max 1.3 \
      -prefix hist_enorm.$group_label \
      $cov_enorm
endif

#--- histogram of cc-values w/ motion ---
if (! -e histCCpearson.$covar_label.$cc_matrix.1D) then
echo "Computing histogram for CC.$covar_label.$cc_matrix"
3dhistog \
   -nbin 100 \
   -pdf \
   -dind 0 \
   -min -1 \
   -max 1 \
   -prefix histCCpearson.$covar_label.$cc_matrix \
   CC_pearson.$covar_label.$cc_matrix
endif

 #--- histogram of cc-values w/ motion ---
if (! -e histCCspearman.$covar_label.$cc_matrix.1D) then
echo "Computing histogram for CC.$covar_label.$cc_matrix"
3dhistog \
   -nbin 100 \
   -pdf \
   -dind 0 \
   -min -1 \
   -max 1 \
   -prefix histCCspearman.$covar_label.$cc_matrix \
   CC_spearman.$covar_label.$cc_matrix
endif

exit()

#========================================
PERMUTATION:
#----------------------------------------
#Permutations to determine null distribution 
#----------------------------------------

set NS = $n
set NI = 10 #1000

set ListCCs = $MergePath/list.CCs.txt
rm -f $ListCCs

foreach ni (`count -digits 3 1 $NI`)
   
   if ($START_CLEAN) then
      rm -f $MergePath/CC_pearson.$group_label.rnd_${ni}.$cc_matrix
      rm -f $MergePath/CC_spearman.$group_label.rnd_${ni}.$cc_matrix
      rm -f histCCpearson.$group_label.rnd_${ni}.$cc_matrix.1D
      rm -f histCCspearman.$group_label.rnd_${ni}.$cc_matrix.1D
   endif
    
   set cov_random = $MergePath/cov.enorm.$group_label.random.$NS.${ni}.1D
   1deval -num $NS -expr "gran(0,1)" > tmp.$group_label.random.$NS.${ni}.1D
   1dcat tmp.$group_label.random.$NS.${ni}.1D $cov_enorm  | sort -n | awk '{print $2}' > $cov_random
   rm -f tmp.$group_label.random.$NS.${ni}.1D	

   if ($DO_PEARSON) then
   if (! -e $MergePath/CC_pearson.$group_label.rnd_${ni}.$cc_matrix) then
   echo "Computing across subject Pearson correlation"
   3dTcorr1D \
      -prefix $MergePath/CC_pearson.$group_label.rnd_${ni}.$cc_matrix \
      $MergePath/$CatFile \
      $cov_random
   endif
   endif

   if (! -e $MergePath/CC_spearman.$group_label.rnd_${ni}.$cc_matrix) then
   echo "Computing across subject Spearman correlation"
   3dTcorr1D \
      -spearman \
      -prefix $MergePath/CC_spearman.$group_label.rnd_${ni}.$cc_matrix \
      $MergePath/$CatFile \
      $cov_random
   endif

   #--- save for later combine ---
   echo "$MergePath/CC_spearman.$group_label.rnd_${ni}.$cc_matrix" >> $ListCCs

   #--- Histograms ---
   cd $MergePath/

   #--- histogram of cc-values w/ motion ---
   if ($DO_PEARSON) then
   if (! -e histCCpearson.$group_label.rnd_${ni}.$cc_matrix.1D) then
   echo "Computing histogram for CC.$group_label.rnd_${ni}.$cc_matrix"
   3dhistog \
      -nbin 100 \
      -pdf \
      -dind 0 \
      -min -1 \
      -max 1 \
      -prefix histCCpearson.$group_label.rnd_${ni}.$cc_matrix \
      CC_pearson.$group_label.rnd_${ni}.$cc_matrix
   endif
   endif

    #--- histogram of cc-values w/ motion ---
   if (! -e histCCspearman.$group_label.rnd_${ni}.$cc_matrix.1D) then
   echo "Computing histogram for CC.$group_label.rnd_${ni}.$cc_matrix"
   3dhistog \
      -nbin 100 \
      -pdf \
      -dind 0 \
      -min -1 \
      -max 1 \
      -prefix histCCspearman.$group_label.rnd_${ni}.$cc_matrix \
      CC_spearman.$group_label.rnd_${ni}.$cc_matrix
   endif
  
  #--- cleanup ---
  mkdir -p HISTOGRAMS
  mkdir -p TMP
  mv histCC* HISTOGRAMS/
  mv CC*rnd* TMP/
  mv cov*random* TMP/
   
end #ni
exit()

#========================================
SIMILARITY_TO_MEAN:
#----------------------------------------
#similarity of connectivity matrix to group mean 
#----------------------------------------

#--- final corrections ---
set corF = xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni  #<--with GSR, filter
#set corF = xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni  #<--with GSR, no filter
#set corF = xtoru.aw.Nwcdm12fXz${mot_thresh}S8.mni  #<--no GSR, filter

#- corF for group mean (to compare to)
#(N.B. only works if AVG was previously computed)
#set corFgroup = xtoru.aw.Ngwcdm12fXz0.2S8.mni	#<-- with GSR, filter, mot<0.2mm
#set corFgroup = xtoru.aw.Ngwcdm12Xz0.2S8.mni 	#<-- with GSR, no filter, mot<0.2mm
#set corFgroup = xtoru.aw.Ngwcdm12fXz0.4S8.mni	#<-- with GSR, filter, mot<0.4mm
#set corFgroup = xtoru.aw.Ngwcdm12Xz0.4S8.mni 	#<-- with GSR, no filter, mot<0.4mm
set corFgroup = $corF				#<-- group corrections same as individual

set out_label = .ToItself #.To04f #.To02f #""

set group_label_base  = N136 #N121 #
set DIR_FILE_use_base = $ScriptPath/dirs.ok.good_0.4.txt #$ScriptPath/dirs.ok.good_0.2.txt #

if ($START_CLEAN) then
rm -f AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz
endif

#--- compute group mean CC matrix ---
if (! -e $MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz) then
   
   echo "Computing average: AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz"
   
   set ListCC = $MergePath/list.Z.CCmatrix.GordonCS.$corFgroup.txt
   rm -f $ListCC
   foreach dir (`cat $DIR_FILE_use_base`) #($argv) #

      set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
      set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
      set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed

      cd $DataPath/

      set anat     = ${dir}_ses-01_run-01_T1w
      set file     = ${dir}_ses-01_task-rest_run-01_bold

      set anat0    = $anat
      set anatF    = anatQQ.$dir.2mm
      set anatA    = anatQQ.$dir.aff
      set anatW    = anatQQ.${dir}_WARP

      set file0    = $file
      set file1    = $file.xtoru.aw
      set fileF    = $file.$corF
      set fileFgroup = $file.$corFgroup


      set CCmatrix = Z.CCmatrix.GordonCS.$fileFgroup.nii.gz 
      echo "$DataPath/$CCmatrix" >> $ListCC

   end #dir
   3dmerge -prefix $MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz `cat $ListCC`
else
   echo "$MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz exists"
endif

set out_cc = $MergePath/similarity.cc.${group_label}${out_label}.$corF.txt
rm -f $out_cc
foreach dir (`cat $DIR_FILE_use`) #($argv) #
   
   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed

   cd $DataPath/

   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold

   set anat0    = $anat
   set anatF    = anatQQ.$dir.2mm
   set anatA    = anatQQ.$dir.aff
   set anatW    = anatQQ.${dir}_WARP

   set file0    = $file
   set file1    = $file.xtoru.aw
   #set fileF    = $file.xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni  #<-- with GSR, filter
   #set fileF    = $file.xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni   #<-- with GSR, no-filter
   set fileF    = $file.$corF					  #<-- corrections as above
   

   set CCmatrix = Z.CCmatrix.GordonCS.$fileF.nii.gz 
   set AvgCCmatrix = AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz
   
   set cc = `3ddot -demean $MergePath/$AvgCCmatrix $DataPath/$CCmatrix`
   echo "$dir $cc" >> $out_cc

end #dir
exit()

#========================================
SIMILARITY_DISTANCE_TO_MEAN:
#----------------------------------------
# similarity of connectivity matrix to group mean 
# here similarity is computed as the Euclidean distance
#----------------------------------------

#--- final corrections ---
set corF = xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni	#<-- with GSR, filter
#set corF = xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni 	#<-- with GSR, no-filter
#set corF = xtoru.aw.Nwcdm12fXz${mot_thresh}S8.mni   	#<-- without GSR, filter

#- corF for group mean (to compare to)
#(N.B. only works if AVG was previously computed)
#set corFgroup = xtoru.aw.Ngwcdm12fXz0.2S8.mni	#<-- with GSR, filter, mot<0.2mm
#set corFgroup = xtoru.aw.Ngwcdm12Xz0.2S8.mni 	#<-- with GSR, no filter, mot<0.2mm
#set corFgroup = xtoru.aw.Ngwcdm12fXz0.4S8.mni	#<-- with GSR, filter, mot<0.4mm
#set corFgroup = xtoru.aw.Ngwcdm12Xz0.4S8.mni 	#<-- with GSR, no filter, mot<0.4mm
set corFgroup = $corF				#<-- group corrections same as individual

set out_label = .ToItself #.To04 #.To04f #.To02f #""

#set group_label_base  = N121 
#set DIR_FILE_use_base = $ScriptPath/dirs.ok.good_0.2.txt #

set group_label_base  = N136 
set DIR_FILE_use_base = $ScriptPath/dirs.ok.good_0.4.txt 

#set group_label_base  = N137 
#set DIR_FILE_use_base = $ScriptPath/dirs.ok.good_1.0.txt #

if ($START_CLEAN) then
rm -f AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz
endif

#--- compute group mean CC matrix ---
if (! -e $MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz) then
   
   echo "Computing average: AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz"
   
   set ListCC = $MergePath/list.Z.CCmatrix.GordonCS.$corFgroup.txt
   rm -f $ListCC
   foreach dir (`cat $DIR_FILE_use_base`) #($argv) #

      set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
      set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
      set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed

      cd $DataPath/

      set anat     = ${dir}_ses-01_run-01_T1w
      set file     = ${dir}_ses-01_task-rest_run-01_bold

      set anat0    = $anat
      set anatF    = anatQQ.$dir.2mm
      set anatA    = anatQQ.$dir.aff
      set anatW    = anatQQ.${dir}_WARP

      set file0    = $file
      set file1    = $file.xtoru.aw
      set fileF    = $file.$corF
      set fileFgroup = $file.$corFgroup


      set CCmatrix = Z.CCmatrix.GordonCS.$fileFgroup.nii.gz 
      echo "$DataPath/$CCmatrix" >> $ListCC

   end #dir
   3dmerge -prefix $MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz `cat $ListCC`
else
   echo "$MergePath/AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz exists"
endif

set out_cc = $MergePath/similarity.cc_dist.${group_label}${out_label}.$corF.txt
rm -f $out_cc
foreach dir (`cat $DIR_FILE_use`) #($argv) #
   
   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed

   cd $DataPath/

   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold

   set anat0    = $anat
   set anatF    = anatQQ.$dir.2mm
   set anatA    = anatQQ.$dir.aff
   set anatW    = anatQQ.${dir}_WARP

   set file0    = $file
   set file1    = $file.xtoru.aw
   #set fileF    = $file.xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni  #<-- with GSR, filter
   #set fileF    = $file.xtoru.aw.Ngwcdm12Xz${mot_thresh}S8.mni   #<-- with GSR, no-filter
   set fileF    = $file.$corF					  #<-- corrections as above
   
   echo "Computing Euclidean distance of CC matrix to mean for $dir ($corF)..."

   set CCmatrix = Z.CCmatrix.GordonCS.$fileF.nii.gz 
   set AvgCCmatrix = AVG.$group_label_base.CCmatrix.GordonCS.$corFgroup.nii.gz
   
   #--- create mask of ccmatrix (upper triangle) ---
   set mask = Mask.UpperTriangle.$AvgCCmatrix
   if (! -e $MergePath/$mask) then
      3dcalc -a $MergePath/$AvgCCmatrix -expr "step(i-j)" -prefix $MergePath/$mask
   endif
   
   rm -f $DataPath/tmp.SqDiff.$CCmatrix
   if (! -e $DataPath/tmp.SqDiff.$CCmatrix) then
   3dcalc -a $MergePath/$AvgCCmatrix -b $DataPath/$CCmatrix -expr "(a-b)*(a-b)" -prefix $DataPath/tmp.SqDiff.$CCmatrix
   endif
   
   #--- Average of Squares (scaled version of sum-of-squares) ---
   set cc2 = `3dROIstats -mask $MergePath/$mask -mask_f2short -quiet $DataPath/tmp.SqDiff.$CCmatrix`
   echo "$dir $cc2" >> $out_cc

end #dir
exit()


#========================================
DISTANCE_DEPENDENCE:
#----------------------------------------
# Compute distance dependence of QC_FC
#----------------------------------------
cd $MergePath/

#--- final corrections ---
#set corF = xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni.Xk
set corF = xtoru.aw.Nwcdm12fXz${mot_thresh}S8.mni.Xk

set dist_matrix = Distance.GordonCS.nii.gz
set qcfc_matrix = CC_spearman.mean_enorm.$group_label.Z.CCmatrix.GordonCS.$corF.nii.gz
set out_label = 

#--- create triangular mask ---
set mask_triangle = Mask.Triangle.GordonCS.nii.gz
3dcalc -a $MergePath/$qcfc_matrix -expr "step(i-j)" -prefix $mask_triangle 

#--- dump distance and qcfc ---
3dmaskdump -mask $mask_triangle -noijk $MergePath/$qcfc_matrix > $MergePath/qcfc.$qcfc_matrix.1D
3dmaskdump -mask $mask_triangle -noijk $ParcPath/$dist_matrix  > $MergePath/dist.$qcfc_matrix.1D

1dCorrelate $MergePath/qcfc.$qcfc_matrix.1D $MergePath/dist.$qcfc_matrix.1D 


