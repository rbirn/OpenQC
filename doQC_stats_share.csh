#!/bin/tcsh

#---
# individual subject stats for OpenQC project
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
set CENSOR_PREV     = 1	#censor previous time point as well as current
set GRAYPLOT        = 0	#produce Grayplot

#--- file containing list of seed region ---
set SeedFile = $RoiPath/list.ROIs.QC.txt  

set DIR_FILE = $ScriptPath/dirs.txt
set OUTPUT   = $MergePath/stats.th0.2.txt
set MISSING  = $ScriptPath/missing.txt

set OUTPUT_dvars = $MergePath/dvars_summary.txt

echo "missing final output for the following dirs:" > $MISSING
echo "dir TSNRo TSNRf TR NT num_censored num_good pcnt_censored good_time dof mean_enorm max_enorm mean_dvars fcs reho fwhmx fwhmy fwhmz fwhmc nvox_brain nvox_gm nvox_wm nvox_csf sigma_w_x sigma_w_y sigma_w_z" > $OUTPUT

echo "dir mean_dvars" > $OUTPUT_dvars

if ($CENSOR_PREV) then
   set cens_prev_label  = "_wPrev"
   set cens_prev_string = "-censor_prev_TR"
else
   set cens_prev_label  = ""
   set cens_prev_string = ""
endif

#--- jump ahead if desired ---
#goto CAT
#goto CAT_FIM
#goto CHECK
#goto ORIG

#=============================
# Compute various QC measures 
#-----------------------------
foreach dir (`cat $DIR_FILE`) #($argv) #

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
   set fileF    = $file.xtoru.aw.Ngwcdm12fXz${mot_thresh}S8.mni
   
   set MaskGM   = Mask.GM.anatQQ.$dir.nii.gz
   set MaskWM   = Mask.WM.anatQQ.$dir.nii.gz
   set MaskCSF  = Mask.CSF.anatQQ.$dir.nii.gz
   set MaskBrain = Mask.Brain.anatQQ.$dir.2mm.nii.gz
   
   set WarpFile = ${anatW}_WARP.nii
   
   set motfile          = mot.$file0.1D
   set CensorFilePrefix = mot_censor.$file0.X${mot_thresh}${cens_prev_label}
   set CensorFile       = ${CensorFilePrefix}_censor.1D
   set EnormFile        = ${CensorFilePrefix}_enorm.1D
   
   set dvars_file       = dvars_s.$file0.1D
   
   set mask0 = Mask.Brain.$file.nii.gz
   set maskF = Mask.Brain.$file1.mni.nii.gz
   set maskA = Mask.Brain.$anatA.nii.gz
   
   if (! -e $DataPath/$fileF.nii.gz) then
      echo "$dir" >> $MISSING
   endif
   
   #--- motion ---
   set num_good = `1dsum $CensorFile`
   set NT = `wc -l < $CensorFile`
   set TR = `3dinfo -tr $OrigPath/$file.nii.gz`
   
   set num_censored = `ccalc -form int "$NT - $num_good"`
   set pcnt_censored = `ccalc "$num_censored * 100 / $NT"`
   set good_time = `ccalc "($NT - $num_censored) * $TR"`
   set mean_enorm = `3dTstat -mean -prefix - $EnormFile\'`
   set max_enorm = `sort -n $EnormFile | tail -1`
   
   if ("$num_good" == "") set num_good = 0
   if ("$NT" == "") set NT = 0
   if ("$TR" == "") set TR = 0
   if ("$num_censored" == "") set num_censored = -1
   if ("$pcnt_censored" == "") set pcnt_censored = -1
   if ("$good_time" == "") set good_time = 0
   if ("$mean_enorm" == "") set mean_enorm = 0
   if ("$max_enorm" == "") set max_enorm = 0
   
   #--- dvars ---
   set mean_dvars = `3dTstat -mean -prefix - $dvars_file\'`
   
   echo "$dir $mean_dvars" >> $OUTPUT_dvars
   
   #--- degrees of freedom ---
   1dBport -band 0.01 0.1 -nodata $NT $TR > $DataPath/bport.1D
   set bp_ort =  `head -1 $DataPath/bport.1D | awk '{print NF}'`
   set dof = `ccalc -form int -expr "$num_good - $bp_ort - 3 - 6 - 12"` #3=polort, 6=gwcd, 12=mot
   if ("$dof" == "") set dof = NA
   
   #--- determine censor string ---
   set uncensor_string = `1d_tool.py -infile $motfile -censor_motion $mot_thresh tmp.$CensorFilePrefix $cens_prev_string -show_trs_uncensored encoded -overwrite`
   
   #--- volumetrics (note: GM/WM/CSF=1mm, Brain=2mm) ---
   set nvox_brain = `3dROIstats -quiet -nzvoxels -mask $MaskBrain $MaskBrain | awk '{print $2}'`
   set nvox_gm    = `3dROIstats -quiet -nzvoxels -mask $MaskGM $MaskGM | awk '{print $2}'`
   set nvox_wm    = `3dROIstats -quiet -nzvoxels -mask $MaskWM $MaskWM | awk '{print $2}'`
   set nvox_csf   = `3dROIstats -quiet -nzvoxels -mask $MaskCSF $MaskCSF | awk '{print $2}'`
   
   if ("$nvox_brain" == "") set nvox_brain = 0
   if ("$nvox_gm" == "") set nvox_gm = 0
   if ("$nvox_wm" == "") set nvox_wm = 0
   if ("$nvox_csf" == "") set nvox_csf = 0
   
   #--- warp (amount of deviation from template) ---
   # recompute affine aligned anat
   if (! -e $anatA.nii.gz) then
   3dAllineate -1Dmatrix_apply anatQQ.$dir.aff12.1D -prefix $anatA.nii.gz -master $anatW.nii -input anatSS.${dir}.nii.gz
   endif
   if (! -e $maskA) then
   3dAutomask -prefix $maskA $anatA.nii.gz
   endif
   set sigma_warps = (`3dROIstats -quiet -sigma -mask $maskA $anatW.nii | awk '{print $2}'`)
   if ("$sigma_warps" == "") then
      set sigma_warps = (0 0 0)
   endif
   
   #--- whole brain masks ---
   if (! -e $mask0) then
   3dAutomask -prefix $mask0 $OrigPath/$file.nii.gz'['$ign']'
   endif
   
   if (! -e $maskF) then
      if (! -e Mask.Brain.$file1.nii.gz) then
      3dAutomask -prefix Mask.Brain.$file1.nii.gz $file1.nii.gz
      endif
      3dresample -master $RoiPath/MNI_avg152T1+tlrc -input Mask.Brain.$file1.nii.gz -prefix $maskF
   endif
   
   #--- TSNR ---
   if (! -e Stats.$file.nii.gz) then
   3dTstat -mean -stdev -cvarinv -prefix Stats.$file.nii.gz $OrigPath/$file.nii.gz'['$ign'-$]'
   endif
   
   rm -f Stats.$fileF.nii.gz #TMP
   if (! -e Stats.$fileF.nii.gz) then
   3dTstat -mean -stdev -cvarinv -prefix Stats.$fileF.nii.gz $fileF.nii.gz'['$uncensor_string']'
   endif
   
   set TSNRo = `3dROIstats -mask $mask0 -quiet Stats.$file.nii.gz'[2]'`
   set TSNRf = `3dROIstats -mask $maskF -quiet Stats.$fileF.nii.gz'[2]'`
   
   if ("$TSNRo" == "") set TSNRo = 0
   if ("$TSNRf" == "") set TSNRf = 0
    
   #-- FWHM ---
   set fwhm_file = $DataPath/fwhm.$file.txt
   if (! -e $fwhm_file || -z $fwhm_file) then
      set fwhm = (`3dFWHMx -mask $mask0 -ShowMeClassicFWHM -detrend -input $OrigPath/$file.nii.gz'['$ign'-$]' | head -1`)
      echo "$fwhm" |& tee $fwhm_file
   endif
   if (-e $fwhm_file) then
      set fwhm = (`cat $fwhm_file`)
   else
      set fwhm = (0 0 0)
   endif
   
   #--- FCS ---
   #- Scale by stdev -
   if (! -e $fileF.sc.nii.gz) then
   3dcalc -datum float -a $fileF.nii.gz -b Stats.$fileF.nii.gz'[1]' -expr "a/b" -prefix $fileF.sc.nii.gz
   endif
   #- Compute (global) mean for each time point -
   if (! -e ts.avg.$fileF.1D || -z ts.avg.$fileF.1D) then
   3dROIstats -mask $maskF -quiet $fileF.sc.nii.gz > ts.avg.$fileF.1D
   endif
   #- Compute Correlation -
   if (! -e FCS.$fileF.nii.gz) then
   3dTcorr1D -prefix FCS.$fileF.nii.gz -mask $maskF $fileF.sc.nii.gz ts.avg.$fileF.1D
   endif
   set fcs = `3dROIstats -mask $maskF -quiet FCS.$fileF.nii.gz`       
   if ("$fcs" == "") set fcs = 0
   
   #--- ReHo ---
   if (! -e ReHo.$fileF.nii.gz) then
   3dReHo -prefix ReHo.$fileF.nii.gz -inset $fileF.nii.gz -mask $maskF
   endif
   set reho = `3dROIstats -mask $maskF -quiet ReHo.$fileF.nii.gz`
   if ("$reho" == "") set reho = 0
   
   #--- OUTPUT ---
   echo "$dir $TSNRo $TSNRf $TR $NT $num_censored $num_good $pcnt_censored $good_time $dof $mean_enorm $max_enorm $mean_dvars $fcs $reho $fwhm $nvox_brain $nvox_gm $nvox_wm $nvox_csf $sigma_warps[1] $sigma_warps[2] $sigma_warps[3]" >> $OUTPUT
   
   #--- Grayplots ---
   if ($GRAYPLOT) then
      
      #--- resample ---
      if (! -e Mask.WM.$anatF.mni.nii.gz) then
      3dresample -master $RoiPath/MNI_avg152T1+tlrc -input Mask.WM.$anatF.nii.gz -prefix Mask.WM.$anatF.mni.nii.gz
      3dresample -master $RoiPath/MNI_avg152T1+tlrc -input Mask.GM.$anatF.nii.gz -prefix Mask.GM.$anatF.mni.nii.gz
      3dresample -master $RoiPath/MNI_avg152T1+tlrc -input Mask.CSF.$anatF.nii.gz -prefix Mask.CSF.$anatF.mni.nii.gz
      endif
      
      set MaskWM  = Mask.WM.$anatF.mni.nii.gz
      set MaskGM  = Mask.GM.$anatF.mni.nii.gz
      set MaskCSF = Mask.CSF.$anatF.mni.nii.gz
      
      
      if (! -e MaskSeg.$anatF.mni.nii.gz) then
      3dcalc -a $MaskGM -b $MaskCSF -c $MaskWM -prefix MaskSeg.$anatF.mni.nii.gz -expr 'a+(2*b)+(3*c)' -short
      endif
      
      if (! -e 3dGp.$fileF.png) then
      3dGrayplot -mask MaskSeg.$anatF.mni.nii.gz -peelorder -prefix 3dGp.$fileF.png $fileF.nii.gz
      endif
      
   endif
   
   #exit() #TMP
   
end

exit()

#===========
CAT:
#===========
# Concatenate measures across subjects 
#-----------
set ListA = $MergePath/list.Anat.txt
set ListE = $MergePath/list.EPI.txt
set ListR = $MergePath/list.ReHo.txt
set ListF = $MergePath/list.FCS.txt
set ListS = $MergePath/list.StDev.txt
set ListT = $MergePath/list.TSNR.txt
rm -f $ListA
rm -f $ListE
rm -f $ListR
rm -f $ListF
rm -f $ListS
rm -f $ListT

foreach dir (`cat $DIR_FILE`)
   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
   
   cd $DataPath/
   
   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold
   
   set anat0    = $anat
   set anatF    = anatQQ.$dir.2mm
   
   set cor      = aw_aZ
   
   set file0    = $file
   set file1    = $file.xtoru.$cor.1
   set fileF    = $file.xtoru.aw.Ngwcdm12fXz0.4S8.mni
   
   set f_reho   = ReHo.$fileF
   set f_fcs    = FCS.$fileF
   set f_stats  = Stats.$fileF
   
   if (-e $DataPath/$anatF.nii.gz && -e $DataPath/$file1.nii.gz) then
   echo "$DataPath/$anatF.nii.gz" >> $ListA
   echo "$DataPath/$file1.nii.gz" >> $ListE
   endif
   
   if (-e $DataPath/$f_reho.nii.gz) then
   echo $DataPath/$f_reho.nii.gz >> $ListR
   endif
   
   if (-e $DataPath/$f_fcs.nii.gz) then
   echo $DataPath/$f_fcs.nii.gz >> $ListF
   endif
   
   if (-e $DataPath/$f_stats.nii.gz) then
   echo $DataPath/$f_stats.nii.gz'[1]' >> $ListS #1=StDev
   endif
   
   if (-e $DataPath/$f_stats.nii.gz) then
   echo $DataPath/$f_stats.nii.gz'[2]' >> $ListS #2=TSNR
   endif
   
end

if (! -e $MergePath/CAT.Anat.nii.gz) then
3dTcat -prefix $MergePath/CAT.Anat.nii.gz `cat $ListA`
endif
if (! -e $MergePath/CAT.EPI.$cor.nii.gz) then
3dTcat -prefix $MergePath/CAT.EPI.$cor.nii.gz `cat $ListE`
endif
if (! -e $MergePath/CAT.ReHo.nii.gz) then
3dTcat -prefix $MergePath/CAT.ReHo.nii.gz `cat $ListR`
endif
if (! -e $MergePath/CAT.FCS.nii.gz) then
3dTcat -prefix $MergePath/CAT.FCS.nii.gz `cat $ListF`
endif

if (! -e $MergePath/CAT.StDev.nii.gz) then
3dTcat -prefix $MergePath/CAT.StDev.nii.gz `cat $ListS`
endif

exit()

#===========
CAT_FIM:
#===========
# Concatenate seed-based Functional connectivity Images
#-----------
echo "Running CAT_FIM ..."
set MISSING  = $ScriptPath/missing.fc.txt
rm -f $MISSING

foreach seed (`cat $SeedFile`)
   set ListF = $MergePath/list.FC.$seed.txt
   set ListA = $MergePath/list.AnatFC.txt
   rm -f $ListF
   rm -f $ListA

   foreach dir (`cat $DIR_FILE`)
      set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
      set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
      set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed

      cd $DataPath/

      set file     = ${dir}_ses-01_task-rest_run-01_bold
      set fileF    = $file.xtoru.aw.Ngwcdm12fXz0.4S8.mni
      
      set anat     = ${dir}_ses-01_run-01_T1w
      set anatF    = anatQQ.$dir.2mm

      set FimFile = tmp.Z.Fim.$seed.$fileF.nii.gz
      
      if (-e $DataPath/$FimFile) then
         echo "Adding $dir ($seed) ..."
	 echo $DataPath/$FimFile >> $ListF
	 echo "$DataPath/$anatF.nii.gz" >> $ListA
      else
         echo "($dir) Missing $FimFile" |& tee -a $MISSING
      endif
      
      
   end
   
   if (! -e $MergePath/CAT.FC.$seed.nii.gz) then
   echo "Concatenating $seed ..."
   3dTcat -prefix $MergePath/CAT.FC.$seed.nii.gz `cat $ListF`
   endif
   
   if (! -e $MergePath/CAT.Anat.fc.nii.gz) then
   3dTcat -prefix $MergePath/CAT.Anat.fc.nii.gz `cat $ListA`
   endif

end
exit()


#===========
CHECK:
#===========
# Create links for more easily checking the results in AFNI
#-----------

foreach dir (`cat $DIR_FILE`)
   set OrigPath  = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath  = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath  = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
   set CheckPath = $RootPath/fmri-open-qc-rest/CHECK/$dir
   
   cd $DataPath/
   
   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold
   
   set anat0    = $anat
   set anat1    = ${dir}_ses-01_run-01_T1w.o    
   set anatF    = anatQQ.$dir.2mm
   
   set file0    = $file
   set file1    = $file.xtoru.1
   set file2    = $file.xtoru.aw.1
   set file3    = $file.xtoru.aw_aZ.1
   set fileF    = $file.xtoru.aw.Ngwcdm12fXz0.4S8.mni
   set f_reho   = ReHo.$fileF
   set f_fcs    = FCS.$fileF
   set f_reho_o = ReHo.$file.x
   set f_fcs_o  = FCS.$file.x
   set f_stats_o = Stats.$file
   
   if (! -e $DataPath/$anat1.nii.gz) then
      echo "missing anat for $dir"
   endif
   if (! -e $DataPath/$file1.nii.gz) then
      echo "missing epi for $dir"
   endif
   
   if (-e $DataPath/$anat1.nii.gz && -e $DataPath/$file1.nii.gz) then
      if (! -e $CheckPath) then
         mkdir $CheckPath
      endif
      #echo "Linking $dir ..."
      if (! -e $CheckPath/$file0.nii.gz) ln -s $OrigPath/$file0.nii.gz $CheckPath/$file0.nii.gz
      if (! -e $CheckPath/$anat1.nii.gz) ln -s $DataPath/$anat1.nii.gz $CheckPath/$anat1.nii.gz
      if (! -e $CheckPath/$file1.nii.gz) ln -s $DataPath/$file1.nii.gz $CheckPath/$file1.nii.gz 
      if (! -e $CheckPath/$file2.nii.gz) ln -s $DataPath/$file1.nii.gz $CheckPath/$file2.nii.gz 
      if (! -e $CheckPath/$file3.nii.gz) ln -s $DataPath/$file1.nii.gz $CheckPath/$file3.nii.gz 
   endif
   
   if (-e $DataPath/$f_reho.nii.gz) then
      if (! -e $CheckPath) then
         mkdir $CheckPath
      endif
      echo "Linking ReHo $dir ..."
      rm -f $CheckPath/$f_reho.nii.gz
      if (! -e $CheckPath/$f_reho.nii.gz) ln -s $DataPath/$f_reho.nii.gz $CheckPath/$f_reho.nii.gz
   endif
   
   if (-e $DataPath/$f_fcs.nii.gz) then
      if (! -e $CheckPath) then
         mkdir $CheckPath
      endif
      echo "Linking fcs $dir ..."
      rm -f $CheckPath/$f_fcs.nii.gz
      if (! -e $CheckPath/$f_fcs.nii.gz) ln -s $DataPath/$f_fcs.nii.gz $CheckPath/$f_fcs.nii.gz
   endif
   
   if (-e $DataPath/$f_reho.nii.gz) then
      echo "Linking ReHo_o $dir ..."
      if (! -e $CheckPath/$f_reho_o.nii.gz) ln -s $DataPath/$f_reho_o.nii.gz $CheckPath/$f_reho_o.nii.gz
   endif
   if (-e $DataPath/$f_fcs_o.nii.gz) then
      echo "Linking fcs_o $dir ..."
      if (! -e $CheckPath/$f_fcs_o.nii.gz) ln -s $DataPath/$f_fcs_o.nii.gz $CheckPath/$f_fcs_o.nii.gz
   endif
   
   if (-e $DataPath/$f_stats_o.nii.gz) then
      echo "Linking stats $dir ..."
      if (! -e $CheckPath/$f_stats_o.nii.gz) ln -s $DataPath/$f_stats_o.nii.gz $CheckPath/$f_stats_o.nii.gz
   endif
   
end
exit()

#===========
ORIG:
#===========
# Compute stats from original (unprocessed) files 
#-----------
foreach dir ($argv) #(`cat $DIR_FILE`) #

   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
   
   cd $DataPath/
   
   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold
   
   set file0 = $file
   
   #--- remove first ---
   if (! -e $file.x.nii.gz) then
   3dTcat -prefix $file.x.nii.gz $OrigPath/$file.nii.gz'['$ign'-$]'
   endif
   set file = $file.x
   
   #--- mask ---
   set mask = Mask.Brain.$file.nii.gz
   if (! -e $mask) then
   3dAutomask -prefix $mask $file.nii.gz'[0]'
   endif
   
   
   #--- TSNR ---
   if (! -e Stats.$file.nii.gz) then
   3dTstat -mean -stdev -cvarinv -prefix Stats.$file.nii.gz $file.nii.gz
   endif
   
   #--- FCS ---
   #- Scale by stdev -
   if (! -e $file.sc.nii.gz) then
   3dcalc -datum float -a $file.nii.gz -b Stats.$file.nii.gz'[1]' -expr "a/b" -prefix $file.sc.nii.gz
   endif
   #- Compute (global) mean for each time point -
   if (! -e ts.avg.$file.1D || -z ts.avg.$file.1D) then
   3dROIstats -mask $mask -quiet $file.sc.nii.gz > ts.avg.$file.1D
   endif
   #- Compute Correlation -
   if (! -e FCS.$file.nii.gz) then
   3dTcorr1D -prefix FCS.$file.nii.gz -mask $mask $file.sc.nii.gz ts.avg.$file.1D
   endif
   set fcs = `3dROIstats -mask $mask -quiet FCS.$file.nii.gz`       
   if ("$fcs" == "") set fcs = 0
   
   #--- ReHo ---
   if (! -e ReHo.$file.nii.gz) then
   3dReHo -prefix ReHo.$file.nii.gz -inset $file.nii.gz -mask $mask
   endif
   set reho = `3dROIstats -mask $mask -quiet ReHo.$file.nii.gz`
   if ("$reho" == "") set reho = 0
   
   #--- Cleanup ---
   rm -f $DataPath/$file0.x.nii.gz
   
end
