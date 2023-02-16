#!/bin/tcsh

#---
# Preprocessing for OpenQC project
#---
# Author: R.M. Birn
#---

#---------------------------
# Paths:
#---------------------------
set RootPath     = ~/src/OpenQC
set ScriptPath   = $RootPath/scripts
set MergePath    = $RootPath/MERGE
set RoiPath      = $RootPath/ROIs 
set TemplatePath = $RootPath/mni_icbm152_nlin_asym_09c

#--- file containing list of seed region ---
set SeedFile = $RoiPath/list.ROIs.QC.txt #$ScriptPath/list.ROIs.Use.txt 

#---------------------------
# Parameters:
#---------------------------
set ign             = 4		#number of time points to ignore
set dxyz            = 2		#resampled resolution (in template space)
set fwhm            = 8		#full width at half maximum - target for smoothing
set mot_thresh      = 0.4 #1.0 #5.0 #0.4 #0.2 #motion censoring threshold
set anaticor_fwhm   = 40	#FWHM for anaticor (not yet implemented)
      
set template   = mni_icbm152_t1_tal_nlin_asym_09c.nii 

set BANDPASS       = 1  #Perform bandpass filtering (0.01-0.1Hz)
set CLEANUP        = 1  #Cleanup files at end. Not Yet Implemented.
set DEBUG          = 1  #output debugging info
set CENSOR_PREV    = 1  #censor previous time point as well as current
set USE_BFC        = 1  #use bias field corrected files
set USE_NLWARP     = 1  #use non-linear warp (auto_warp.py)
set USE_SSWARP     = 1  #use @SSwarper (instead of 3dSkullStrip and auto_warp.py) overrides USE_NLWARP
set MOT            = 12 #number of motion regressors (0, 6, 12, 24, 36, motsim)
set WM_CSF         = 1  #regress out white matter and CSF (eroded)
set GLOBAL         = 1  #regress out global signal
set CENSOR         = 1  #censor high-motion time points

set SKIP_AFTER_ALIGN = 0 #skip processing steps after alignment (useful if debugging alignment between EPI and T1)

if ($CENSOR_PREV) then
   set cens_prev_label  = "_wPrev"
   set cens_prev_string = "-censor_prev_TR"
else
   set cens_prev_label  = ""
   set cens_prev_string = ""
endif

foreach dir ($argv) 
  
   set OrigPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/func
   set AnatPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/anat
   set DataPath = $RootPath/fmri-open-qc-rest/$dir/ses-01/processed
   set TmpPath  = /scratch/rbirn/$dir
   
   set anat     = ${dir}_ses-01_run-01_T1w
   set file     = ${dir}_ses-01_task-rest_run-01_bold
   
   set anat0        = $anat
   set file0        = $file
   set file_base    = ${dir}_ses-01_task-rest_run-01_bold.xto 
   set file1        = $file0.xtoru.aw
   
   if ($GLOBAL) then
      set corF1     = Ngwcdm12
   else
      set corF1     = Nwcdm12
   endif
   
   if ($BANDPASS) then
      set fileF     = $file0.xtoru.aw.${corF1}fXz${mot_thresh}S${fwhm}.mni
   else
      set fileF     = $file0.xtoru.aw.${corF1}Xz${mot_thresh}S${fwhm}.mni
   endif
   
   set motfile          = mot.$file0.1D
   set CensorFilePrefix = mot_censor.$file0.X${mot_thresh}${cens_prev_label}
   set CensorFile       = ${CensorFilePrefix}_censor.1D
   
   if (! -e $DataPath) then
      mkdir $DataPath
   endif
   
   if (! -e $TmpPath) then
      mkdir $TmpPath
   endif
   
   set LOG_FILE = $DataPath/log.proc.txt
   echo "--------------------------------" >> $LOG_FILE
   date                                    >> $LOG_FILE
   echo "--------------------------------" >> $LOG_FILE
   
   cd $DataPath/
   
   if (! -e $OrigPath/$file0.nii.gz) then
      echo "$OrigPath/$file0.nii.gz does not exist. Skipping"
      continue
   endif
   
   #--- Options to jump ahead when re-running script ---
   #TMP
   #if (-e $fileF.nii.gz) then
   #   echo "$fileF.nii.gz already exists. Aborting."
   #   continue
   #endif
   
   if (0) then
   if (-e $fileF.nii.gz) then
      echo "$fileF.nii.gz already exists. Jumping ahead."
      set file = $fileF
      set MaskTotal = Mask.Brain.AE.${file0}.xtoru.aw.mni.nii.gz
      goto CONNECTIVITY_MATRIX
   endif
   endif #0
   
   #if (0) then
   if (-e $TmpPath/$file1.nii.gz) then
      echo "$TmpPath/$file1.nii.gz already exists. Jumping ahead."
      set file = $file1
      set MaskTotal = Mask.Brain.AE.${file1}.nii.gz
      set MaskCSF = Mask.CSF.anatQQ.${dir}.2mm
      set MaskGM  = Mask.GM.anatQQ.${dir}.2mmm
      set MaskWM  = Mask.WM.anatQQ.${dir}.2mm
   
      goto NUISANCE
   endif
   #endif
   
   
   #====================================================================================
   DEOBLIQUE_ANAT:
   #------------------------------------------------------------------------------------
   #--- DeOblique Anat ---
   set anat_out = $anat.o
   if (! -e $anat_out.nii.gz) then
      if ($DEBUG) echo "Running Deoblique for anat..."
      3dWarp -deoblique -prefix $DataPath/$anat_out.nii.gz $AnatPath/$anat.nii.gz |& tee -a $LOG_FILE
   endif
   set anat = $anat_out
   
   #====================================================================================
   SKULLSTRIP:
   #------------------------------------------------------------------------------------
   #--- Skull Strip ---
   if (! $USE_SSWARP) then
      set anat_out = $anat.ns
      if (! -e $anat_out.nii.gz) then
	 if ($DEBUG) echo "Running Skullstrip..."
	 3dSkullStrip -input $anat.nii.gz -prefix $anat_out.nii.gz |& tee -a $LOG_FILE
      endif
      set anat    = $anat_out
   endif
   
   #====================================================================================
   ALIGN_ANAT_TO_TEMPLATE:
   #------------------------------------------------------------------------------------
   #--- Warp Anatomical to Template (NONLINEAR) ---
   if ($USE_SSWARP) then
      if (! -e anatQQ.$dir.nii.gz) then
	 echo "---@SSwarper---"
	 @SSwarper \
            -base MNI152_2009_template_SSW.nii.gz \
	    -input $anat.nii.gz -subid $dir 

	 echo "Gzipping anatSS.$dir.nii" |& tee -a $LOG_FILE
	 gzip anatSS.$dir.nii
	 gzip anatQQ.$dir.nii
      endif
      
      set anat_ns = anatSS.$dir
      set anat_warp   = anatQQ.${dir}_WARP.nii
      set anat_affine = anatQQ.${dir}.aff12.1D 
      set anat_al     = anatQQ.${dir}
   else
      if (! -e $anat.aw.nii.gz) then
	 echo "Running auto_warp.py (Nonlinear warp)"
	 auto_warp.py -base $template -input $anat.nii.gz -suffix .aw |& tee -a $LOG_FILE

	 if (-e awpy/$anat.aw.nii && ! -e $anat.aw.nii) then
	 echo "Moving $anat.aw.nii out of awpy" |& tee -a $LOG_FILE
	 mv ./awpy/$anat.aw.nii ./
	 endif

	 echo "Gzipping $anat.aw.nii" |& tee -a $LOG_FILE
	 gzip $anat.aw.nii
	 
      endif
      
      #--- link to affine aligned anat ---
      if (! -e $anat.un.aff.nii) then
         ln -s ./awpy/anat.un.aff.nii $anat.un.aff.nii
      endif
      if (! -e anat.un.aff.nii.Xaff12.1D) then
         ln -s ./awpy/anat.un.aff.nii.Xaff12.1D $anat.un.aff.nii.Xaff12.1D
      endif
      if (! -e anat.un.aff.qw_WARP.nii) then
         ln -s ./awpy/anat.un.aff.qw_WARP.nii $anat.un.aff.qw_WARP.nii
      endif
      
      
      set anat_ns = $anat
      set anat_warp   = $anat.un.aff.qw_WARP.nii
      set anat_affine = $anat.un.aff.nii.Xaff12.1D
      
      if ($USE_NLWARP) then
         set anat_al     = $anat.aw
      else
         set anat_al = $anat.un.aff
      endif
   endif
   # do not set anat to anat.aw yet. We need the dataset prior to warp for Anat-EPI Alignment.
   
   #====================================================================================
   REMOVE_FIRST_TIMEPOINTS:
   #------------------------------------------------------------------------------------
   set file_out = ${file}.x
   if (! -e $file_out.nii.gz) then
   echo "Removing first $ign time points"
   3dcalc -a $OrigPath/$file.nii.gz'['$ign'-$]' -expr 'a' -prefix $TmpPath/$file_out.nii.gz |& tee -a $LOG_FILE
   endif
   set file = $file_out
   
   #====================================================================================
   TSHIFT:
   #------------------------------------------------------------------------------------
   #--- Slice timing correction ---
   # (NOTE: This assumes dicom->nifty preserved the slice timing information correctly)
   set file_out = ${file}t
   if (! -e $TmpPath/$file_out.nii.gz) then
   3dTshift -tzero 0 -prefix $TmpPath/$file_out.nii.gz $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
   endif
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) rm -f $TmpPath/$file.nii.gz
   set file = $file_out
   
   #====================================================================================
   DEOBLIQUE_EPI:
   #------------------------------------------------------------------------------------
   set file_out = ${file}o
   if (! -e $TmpPath/$file_out.nii.gz) then
      3dWarp -deoblique -prefix $TmpPath/$file_out.nii.gz $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
   endif
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) rm -f $TmpPath/$file.nii.gz
   set file = $file_out
   
   #====================================================================================
   VOLREG:
   #------------------------------------------------------------------------------------
   #--- motion correction (& chop off first 3 volumes) ---
   set file_out = ${file}r
   if (! -e $TmpPath/$file_out.nii.gz) then
   3dvolreg -base $TmpPath/$file_base.nii.gz'[0]' -1Dfile $DataPath/$motfile -prefix $TmpPath/$file_out.nii.gz $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
   endif
   #--- if failed, try 3dAllineate (e.g. if grids don't match)---
   if (! -e $TmpPath/$file_out.nii.gz) then
   3dAllineate -base $TmpPath/$file_base.nii.gz -master $TmpPath/$file_base.nii.gz -warp shr -1Dparam_save $DataPath/$motfile -prefix $TmpPath/$file_out.nii.gz -input $TmpPath/$file.nii.gz 
   endif
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) rm -f $TmpPath/$file.nii.gz
   set file = $file_out
   
   #====================================================================================
   #--- motion censor file ---
   if (! -e ${CensorFilePrefix}_censor.1D) then
      echo "---Motion Censoring---"  |& tee -a $LOG_FILE
      1d_tool.py -infile $motfile -censor_motion $mot_thresh $CensorFilePrefix $cens_prev_string
   endif
      
   #====================================================================================
   BIAS_FIELD_CORRECTION:
   #------------------------------------------------------------------------------------
   if ($USE_BFC) then
      set file_out = ${file}u
      set curr_path = `pwd`
      cd $TmpPath/
      if (! -e $file_out.nii.gz) then
         3dcalc -a $file.nii.gz'[0]' -expr 'a' -prefix tmp.1.$file.nii.gz  |& tee -a $LOG_FILE
	 N4BiasFieldCorrection  -d 3 -i tmp.1.$file.nii.gz --output '[corr_ign.nii,bias_field_'$file'.nii]' -s 2 -b '[150]' -c '[200x200,0.0]'  |& tee -a $LOG_FILE
         3dcalc -a $file.nii.gz -b bias_field_$file.nii -expr "a/b" -prefix $file_out.nii.gz  |& tee -a $LOG_FILE
	 rm -f tmp.1.$file.nii.gz
      endif 
      if (-e $file_out.nii.gz && $CLEANUP) rm -f $file.nii.gz   
      cd $curr_path/
      set file = $file_out
   endif
   
   #====================================================================================
   ALIGN_EPI_TO_ANAT:
   #------------------------------------------------------------------------------------
   #--- Align EPI to Anatomical ---
   set cor = aZgc #aZ #label for aligned dataset (a=aligned, Z=lpc+ZZ, g=giant_move, c=cmass)
   
   #--- Extract one volume ---
   if (! -e $DataPath/$file.1.nii.gz) then
   3dcalc -a $TmpPath/$file.nii.gz'[0]' -expr 'a' -prefix $DataPath/$file.1.nii.gz
   endif
   
   #--- Align one volume ---
   if (! -e xform.${file}-${anat}.$cor.aff12.1D || -z xform.${file}-${anat}.$cor.aff12.1D) then
      echo "Aligning ${file} to ${anat} ..." |& tee -a $LOG_FILE
      if (! -e $file.${cor}+orig.HEAD) then
      align_epi_anat.py \
	 -anat $anat.nii.gz \
	 -epi $file.1.nii.gz \
	 -epi_base 0 \
	 -epi2anat \
	 -volreg off \
	 -tshift off \
	 -anat_has_skull no \
	 -epi_strip 3dAutomask \
	 -cost lpc+ZZ \
	 -cmass cmass \
	 -giant_move \
	 -suffix ".${cor}"  |& tee -a $LOG_FILE 
	 
      endif
      #--- Save transforms ---
      echo "Copying ${file}.1.${cor}_mat.aff12.1D to xform.${file}-${anat}.$cor.aff12.1D"
      /bin/cp -f ${file}.1.${cor}_mat.aff12.1D xform.${file}-${anat}.$cor.aff12.1D

      #--- Concatenate transforms ---
      #cat_matvec -ONELINE awpy/anat.un.aff.Xat.1D xform.$file-Anat.aff12.1D > xform.$file-Template.aff12.1D

   else
      echo "xform.${file}-${anat}.aff12.1D exists."
   endif
   
   
   #====================================================================================
   APPLY_NL_WARP:
   #------------------------------------------------------------------------------------
   
   #--- Concatenate affine transforms (EPI->Anat, Anat->Template) ---
   set xform_affine = xform.$file-Template.$cor.aff12.1D
   if (! -e $xform_affine || -z $xform_affine) then
   cat_matvec -ONELINE $anat_affine xform.${file}-${anat}.$cor.aff12.1D > $xform_affine
   endif
   
   #------------------------------------------------------------------------------------
   # now that alignment is done, we can set anat to the template-aligned version
   #------------------------------------------------------------------------------------
   set anat = $anat_al
   #------------------------------------------------------------------------------------
   
   #--- nii vs nii.gz ---
   if (-e $anat.nii.gz) then
      set anat_full = $anat.nii.gz
   else if (-e $anat.nii) then
      set anat_full = $anat.nii
   else
      echo "Could not find $anat.nii or $anat.nii.gz"
      exit()
   endif
   
   #--- resample anatomical to serve as output template ---
   if (! -e $anat.${dxyz}mm.nii.gz) then
   3dresample -input $anat_full -dxyz $dxyz $dxyz $dxyz -prefix $anat.${dxyz}mm.nii.gz
   endif
   
   if ($USE_NLWARP) then
      #--- Apply nonlinear warp to EPI ---
      set file_out = ${file}.aw
      #set file_out = ${file}.aw_${cor}
      
      if (! -e $TmpPath/$file_out.nii.gz) then
      echo "Applying NL warp to $file"
      3dNwarpApply \
	 -nwarp  ${anat_warp} ${xform_affine} \
	 -prefix $TmpPath/$file_out.nii.gz \
	 -master $anat.${dxyz}mm.nii.gz \
	 -source $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
      endif

   else
      
      #--- Apply affine warp to EPI ---
      set file_out = ${file}.a

      if (! -e $file_out.nii.gz) then
      echo "Applying NL warp to $file"
      3dAllineate \
	 -1Dmatrix_apply ${xform_affine} \
	 -prefix $TmpPath/$file_out.nii.gz \
	 -master $anat.${dxyz}mm.nii.gz \
	 $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
      endif

   endif
   if (-e $file_out.nii.gz && $CLEANUP) rm -f $file.nii.gz   
   set file = $file_out
   
   #--- save one for later ---
   if (! -e $file.1.nii.gz) then
   3dcalc -a $TmpPath/$file.nii.gz'[0]' -expr 'a' -prefix $DataPath/$file.1.nii.gz |& tee -a $LOG_FILE
   endif
   
   #TMP
   if ($SKIP_AFTER_ALIGN) then
      echo "Skipping other processing steps after alignment ..." |& tee -a $LOG_FILE
      continue
   endif
   
   #====================================================================================
   SEGMENT:
   #------------------------------------------------------------------------------------
   #--- Segment Anatomical ---
   set MaskCSF = Mask.CSF.$anat
   set MaskGM  = Mask.GM.$anat
   set MaskWM  = Mask.WM.$anat
   if (-e segment_${anat}_seg_0.nii.gz || -e $MaskCSF.nii.gz) then
      echo "FSL segmentation already exists." |& tee -a $LOG_FILE
   else
      echo "starting FSL segmentation..." |& tee -a $LOG_FILE
      fast -t 1 -g -p -o segment_$anat $anat.nii.gz
      echo "...done" |& tee -a $LOG_FILE
   endif
   if (! -e $MaskCSF.nii.gz) then
   3dcopy segment_${anat}_seg_0.nii.gz $MaskCSF.nii.gz |& tee -a $LOG_FILE
   3dcopy segment_${anat}_seg_1.nii.gz $MaskGM.nii.gz |& tee -a $LOG_FILE
   3dcopy segment_${anat}_seg_2.nii.gz $MaskWM.nii.gz |& tee -a $LOG_FILE
   endif
   
   #if (! -e $MaskCSF) exit()
   
   #--- resample masks ---
   if (! -e $MaskCSF.${dxyz}mm.nii.gz) then
      3dresample -input $MaskCSF.nii.gz -dxyz $dxyz $dxyz $dxyz -prefix $MaskCSF.${dxyz}mm.nii.gz |& tee -a $LOG_FILE
      3dresample -input $MaskGM.nii.gz -dxyz $dxyz $dxyz $dxyz -prefix  $MaskGM.${dxyz}mm.nii.gz |& tee -a $LOG_FILE
      3dresample -input $MaskWM.nii.gz -dxyz $dxyz $dxyz $dxyz -prefix  $MaskWM.${dxyz}mm.nii.gz |& tee -a $LOG_FILE
   endif
   set MaskCSF = $MaskCSF.${dxyz}mm
   set MaskGM  = $MaskGM.${dxyz}mm
   set MaskWM  = $MaskWM.${dxyz}mm
   
   #====================================================================================
   MASK:
   #------------------------------------------------------------------------------------
   #--- Create Mask ---
   echo "Creating total brain mask for $file ..."
   set MaskTotal = Mask.Brain.AE.$file.nii.gz
   if (! -e Mask.Brain.$anat.${dxyz}mm.nii.gz) then
   3dAutomask -prefix Mask.Brain.$anat.${dxyz}mm.nii.gz $anat.${dxyz}mm.nii.gz |& tee -a $LOG_FILE
   endif
   if (! -e Mask.Brain.$file.nii.gz) then
   3dAutomask -prefix Mask.Brain.$file.nii.gz $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
   endif
   if (! -e $MaskTotal) then
   3dcalc -a Mask.Brain.$file.nii.gz -b Mask.Brain.$anat.${dxyz}mm.nii.gz -expr "step(a+b)" -prefix $MaskTotal |& tee -a $LOG_FILE
   endif
   
   
   #====================================================================================
   NUISANCE:
   #------------------------------------------------------------------------------------
   #--- Nuisance regression, Temporal Filtering, Spatial Smoothing ---
   
   #----------------------------------------
   # Motion-12
   #----------------------------------------
   if ($DEBUG) echo "...motion..."
   if (! -e ts.1x.$motfile) then
   1dnorm -demean $motfile'[0]' - > ts.1x.$motfile
   1dnorm -demean $motfile'[1]' - > ts.2x.$motfile
   1dnorm -demean $motfile'[2]' - > ts.3x.$motfile
   1dnorm -demean $motfile'[3]' - > ts.4x.$motfile
   1dnorm -demean $motfile'[4]' - > ts.5x.$motfile
   1dnorm -demean $motfile'[5]' - > ts.6x.$motfile
   endif

   if (! -e ts.1x.dx.$motfile) then
   1d_tool.py -derivative -infile ts.1x.$motfile -write - | 1dnorm -demean - - > ts.1x.dx.$motfile
   1d_tool.py -derivative -infile ts.2x.$motfile -write - | 1dnorm -demean - - > ts.2x.dx.$motfile
   1d_tool.py -derivative -infile ts.3x.$motfile -write - | 1dnorm -demean - - > ts.3x.dx.$motfile
   1d_tool.py -derivative -infile ts.4x.$motfile -write - | 1dnorm -demean - - > ts.4x.dx.$motfile
   1d_tool.py -derivative -infile ts.5x.$motfile -write - | 1dnorm -demean - - > ts.5x.dx.$motfile
   1d_tool.py -derivative -infile ts.6x.$motfile -write - | 1dnorm -demean - - > ts.6x.dx.$motfile
   endif
   
   #----------------------------------------
   # Motion-24
   #----------------------------------------
   if ($MOT == "24" || $MOT == "36") then
   if ($DEBUG) echo "...motion24..."
   set NTm = `wc -l < $motfile`
   @ NTm1 = $NTm - 1

   #--- motion at preceding time point ---
   if (! -e ts.s1.$motfile) then
   cat $motfile | head -$NTm1 > tmp.head.$motfile
   cat $motfile | head -1 > tmp.first.$motfile
   cat tmp.first.$motfile tmp.head.$motfile > ts.s1.$motfile
   rm -f tmp.head.$motfile tmp.first.$motfile
   endif

   if (! -e ts.1x.s1.$motfile) then
   1dnorm -demean ts.s1.$motfile'[0]' - > ts.1x.s1.$motfile
   1dnorm -demean ts.s1.$motfile'[1]' - > ts.2x.s1.$motfile
   1dnorm -demean ts.s1.$motfile'[2]' - > ts.3x.s1.$motfile
   1dnorm -demean ts.s1.$motfile'[3]' - > ts.4x.s1.$motfile
   1dnorm -demean ts.s1.$motfile'[4]' - > ts.5x.s1.$motfile
   1dnorm -demean ts.s1.$motfile'[5]' - > ts.6x.s1.$motfile
   endif

   #--- motion squared ---
   if (! -e ts.1x.sq.$motfile) then
   1deval -a $motfile'[0]' -expr "a*a" | 1dnorm -demean - - > ts.1x.sq.$motfile
   1deval -a $motfile'[1]' -expr "a*a" | 1dnorm -demean - - > ts.2x.sq.$motfile
   1deval -a $motfile'[2]' -expr "a*a" | 1dnorm -demean - - > ts.3x.sq.$motfile
   1deval -a $motfile'[3]' -expr "a*a" | 1dnorm -demean - - > ts.4x.sq.$motfile
   1deval -a $motfile'[4]' -expr "a*a" | 1dnorm -demean - - > ts.5x.sq.$motfile
   1deval -a $motfile'[5]' -expr "a*a" | 1dnorm -demean - - > ts.6x.sq.$motfile
   endif

   #--- motion at preceding time point squared ---
   if (! -e ts.1x.s1.sq.$motfile) then
   1deval -a ts.s1.$motfile'[0]' -expr "a*a" | 1dnorm -demean - - > ts.1x.s1.sq.$motfile
   1deval -a ts.s1.$motfile'[1]' -expr "a*a" | 1dnorm -demean - - > ts.2x.s1.sq.$motfile
   1deval -a ts.s1.$motfile'[2]' -expr "a*a" | 1dnorm -demean - - > ts.3x.s1.sq.$motfile
   1deval -a ts.s1.$motfile'[3]' -expr "a*a" | 1dnorm -demean - - > ts.4x.s1.sq.$motfile
   1deval -a ts.s1.$motfile'[4]' -expr "a*a" | 1dnorm -demean - - > ts.5x.s1.sq.$motfile
   1deval -a ts.s1.$motfile'[5]' -expr "a*a" | 1dnorm -demean - - > ts.6x.s1.sq.$motfile
   endif
   endif

   #----------------------------------------
   # Motion-36
   #----------------------------------------
   if ($MOT == "36") then
   if ($DEBUG) echo "...motion36..."

   #--- motion at 2*preceding time point ---
   if (! -e ts.s2.$motfile) then
   cat ts.s1.$motfile | head -$NTm1 > tmp.head2.$motfile
   cat ts.s1.$motfile | head -1 > tmp.first2.$motfile
   cat tmp.first2.$motfile tmp.head2.$motfile > ts.s2.$motfile
   rm -f tmp.head2.$motfile tmp.first2.$motfile
   endif

   if (! -e ts.1x.s2.$motfile) then
   1dnorm -demean ts.s2.$motfile'[0]' - > ts.1x.s2.$motfile
   1dnorm -demean ts.s2.$motfile'[1]' - > ts.2x.s2.$motfile
   1dnorm -demean ts.s2.$motfile'[2]' - > ts.3x.s2.$motfile
   1dnorm -demean ts.s2.$motfile'[3]' - > ts.4x.s2.$motfile
   1dnorm -demean ts.s2.$motfile'[4]' - > ts.5x.s2.$motfile
   1dnorm -demean ts.s2.$motfile'[5]' - > ts.6x.s2.$motfile
   endif

   #--- motion at 2*preceding time point squared ---
   if (! -e ts.1x.s2.sq.$motfile) then
   1deval -a ts.s2.$motfile'[0]' -expr "a*a" | 1dnorm -demean - - > ts.1x.s2.sq.$motfile
   1deval -a ts.s2.$motfile'[1]' -expr "a*a" | 1dnorm -demean - - > ts.2x.s2.sq.$motfile
   1deval -a ts.s2.$motfile'[2]' -expr "a*a" | 1dnorm -demean - - > ts.3x.s2.sq.$motfile
   1deval -a ts.s2.$motfile'[3]' -expr "a*a" | 1dnorm -demean - - > ts.4x.s2.sq.$motfile
   1deval -a ts.s2.$motfile'[4]' -expr "a*a" | 1dnorm -demean - - > ts.5x.s2.sq.$motfile
   1deval -a ts.s2.$motfile'[5]' -expr "a*a" | 1dnorm -demean - - > ts.6x.s2.sq.$motfile
   endif
   endif

   #----------------------------------------
   # MotSim
   #----------------------------------------
   #--- Compute motsim ---
   if ($MOT == "motsim") then
      #--- redo volreg if needed ---
      if (! -e $file0.xr.nii.gz) then
         echo "---Redoing Motion Correction---" 
	 3dvolreg -base $file_base.nii.gz'['$ign']' -prefix $file0.xr.nii.gz $file0.nii'['$ign'-$]'
      endif
      ~birn/bin/@MotSim.csh $file0.xr.nii.gz 0 motsim 12
   endif


   #----------------------------------------
   # Global
   #----------------------------------------
   #--- Brain mask ---
   if ($DEBUG) echo "...global..."
   if (! -e Mask.Brain.$file.nii.gz) then
   3dAutomask -prefix Mask.Brain.$file.nii.gz $TmpPath/$file.nii.gz
   endif
   3dROIstats -mask Mask.Brain.$file.nii.gz -mask_f2short -quiet $TmpPath/$file.nii.gz > tmp.ts.global.$file.1D
   rm -f ts.global.$file.1D
   1dnorm -demean tmp.ts.global.$file.1D ts.global.$file.1D
   rm -f tmp.ts.global.$file.1D
   
   #---------------------------------------------------
   # Derivative of Global
   #---------------------------------------------------
   if ($DEBUG) echo "...global derivative..."
   1d_tool.py -derivative -infile ts.global.$file.1D -write - > tmp.ts.global.d.$file.1D 
   rm -f ts.global.d.$file.1D
   1dnorm -demean tmp.ts.global.d.$file.1D ts.global.d.$file.1D
   rm -f tmp.ts.global.d.$file.1D
   
   #----------------------------------------
   # WM
   #----------------------------------------
   #--- Erode mask ---
   if ($DEBUG) echo "...WM..."
   if (! -e $MaskWM.erode.nii.gz) then
   3dcalc \
      -datum short \
      -a $MaskWM.nii.gz \
      -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
      -expr "a*(1-amongst(0,b,c,d,e,f,g))" \
      -prefix $MaskWM.erode.nii.gz
   endif
   3dROIstats -mask $MaskWM.erode.nii.gz -mask_f2short -quiet $TmpPath/$file.nii.gz > tmp.ts.WMe.$file.1D
   rm -f ts.WMe.$file.1D
   1dnorm -demean tmp.ts.WMe.$file.1D ts.WMe.$file.1D
   rm -f tmp.ts.WMe.$file.1D
   
   #---------------------------------------------------
   # Derivative of WM 
   #---------------------------------------------------
   if ($DEBUG) echo "...WM derivative..."
   1d_tool.py -derivative -infile ts.WMe.$file.1D -write - > tmp.ts.WMe.d.$file.1D 
   rm -f ts.WMe.d.$file.1D
   1dnorm -demean tmp.ts.WMe.d.$file.1D ts.WMe.d.$file.1D
   rm -f tmp.ts.WMe.d.$file.1D
   
   #----------------------------------------
   # CSF 
   #----------------------------------------
   #--- Erode mask ---
   if ($DEBUG) echo "...CSF..."
   if (! -e $MaskCSF.erode.nii.gz) then
   3dcalc \
      -datum short \
      -a $MaskCSF.nii.gz \
      -b a+i -c a-i -d a+j -e a-j -f a+k -g a-k \
      -expr "a*(1-amongst(0,b,c,d,e,f,g))" \
      -prefix $MaskCSF.erode.nii.gz
   endif
   3dROIstats -mask $MaskCSF.erode.nii.gz -mask_f2short -quiet $TmpPath/$file.nii.gz > tmp.ts.CSFe.$file.1D
   rm -f ts.CSFe.$file.1D
   1dnorm -demean tmp.ts.CSFe.$file.1D ts.CSFe.$file.1D
   rm -f tmp.ts.CSFe.$file.1D
   
   #----------------------------------------------------
   # Derivative of CSF 
   #----------------------------------------------------
   if ($DEBUG) echo "...CSF derivative..."
   1d_tool.py -derivative -infile ts.CSFe.$file.1D -write - > tmp.ts.CSFe.d.$file.1D 
   rm -f ts.CSFe.d.$file.1D
   1dnorm -demean tmp.ts.CSFe.d.$file.1D ts.CSFe.d.$file.1D
   rm -f tmp.ts.CSFe.d.$file.1D
   
   #~~~
   set nuis_string = "N"
   set ort_file  = $DataPath/ort_string.txt
   rm -f $ort_file

   if ($GLOBAL) then
      set nuis_string = ${nuis_string}g
      echo "-ort $DataPath/ts.global.$file.1D \
	    -ort $DataPath/ts.global.d.$file.1D " >> $ort_file

   endif

   if ($WM_CSF) then
      set nuis_string = ${nuis_string}wcd
      echo "-ort $DataPath/ts.WMe.$file.1D \
	    -ort $DataPath/ts.WMe.d.$file.1D \
	    -ort $DataPath/ts.CSFe.$file.1D \
	    -ort $DataPath/ts.CSFe.d.$file.1D " >> $ort_file
   endif

   if ($MOT == 6) then
      set nuis_string = ${nuis_string}m6
      echo "-ort $DataPath/ts.1x.$motfile \
	    -ort $DataPath/ts.2x.$motfile \
	    -ort $DataPath/ts.3x.$motfile \
	    -ort $DataPath/ts.4x.$motfile \
	    -ort $DataPath/ts.5x.$motfile \
	    -ort $DataPath/ts.6x.$motfile " >> $ort_file
   else if ($MOT == 12) then
      set nuis_string = ${nuis_string}m12
      echo "-ort $DataPath/ts.1x.$motfile \
	    -ort $DataPath/ts.2x.$motfile \
	    -ort $DataPath/ts.3x.$motfile \
	    -ort $DataPath/ts.4x.$motfile \
	    -ort $DataPath/ts.5x.$motfile \
	    -ort $DataPath/ts.6x.$motfile \
	    -ort $DataPath/ts.1x.dx.$motfile \
	    -ort $DataPath/ts.2x.dx.$motfile \
	    -ort $DataPath/ts.3x.dx.$motfile \
	    -ort $DataPath/ts.4x.dx.$motfile \
	    -ort $DataPath/ts.5x.dx.$motfile \
	    -ort $DataPath/ts.6x.dx.$motfile " >> $ort_file
   else if ($MOT == 24) then
      set nuis_string = ${nuis_string}m24
      echo "-ort $DataPath/ts.1x.$motfile \
	    -ort $DataPath/ts.2x.$motfile \
	    -ort $DataPath/ts.3x.$motfile \
	    -ort $DataPath/ts.4x.$motfile \
	    -ort $DataPath/ts.5x.$motfile \
	    -ort $DataPath/ts.6x.$motfile \
	    -ort $DataPath/ts.1x.s1.$motfile \
	    -ort $DataPath/ts.2x.s1.$motfile \
	    -ort $DataPath/ts.3x.s1.$motfile \
	    -ort $DataPath/ts.4x.s1.$motfile \
	    -ort $DataPath/ts.5x.s1.$motfile \
	    -ort $DataPath/ts.6x.s1.$motfile \
	    -ort $DataPath/ts.1x.sq.$motfile \
	    -ort $DataPath/ts.2x.sq.$motfile \
	    -ort $DataPath/ts.3x.sq.$motfile \
	    -ort $DataPath/ts.4x.sq.$motfile \
	    -ort $DataPath/ts.5x.sq.$motfile \
	    -ort $DataPath/ts.6x.sq.$motfile \
	    -ort $DataPath/ts.1x.s1.sq.$motfile \
	    -ort $DataPath/ts.2x.s1.sq.$motfile \
	    -ort $DataPath/ts.3x.s1.sq.$motfile \
	    -ort $DataPath/ts.4x.s1.sq.$motfile \
	    -ort $DataPath/ts.5x.s1.sq.$motfile \
	    -ort $DataPath/ts.6x.s1.sq.$motfile " >> $ort_file
   else if ($MOT == 36) then
      set nuis_string = ${nuis_string}m36
      echo "-ort $DataPath/ts.1x.$motfile \
	    -ort $DataPath/ts.2x.$motfile \
	    -ort $DataPath/ts.3x.$motfile \
	    -ort $DataPath/ts.4x.$motfile \
	    -ort $DataPath/ts.5x.$motfile \
	    -ort $DataPath/ts.6x.$motfile \
	    -ort $DataPath/ts.1x.sq.$motfile \
	    -ort $DataPath/ts.2x.sq.$motfile \
	    -ort $DataPath/ts.3x.sq.$motfile \
	    -ort $DataPath/ts.4x.sq.$motfile \
	    -ort $DataPath/ts.5x.sq.$motfile \
	    -ort $DataPath/ts.6x.sq.$motfile \
	    -ort $DataPath/ts.1x.s1.$motfile \
	    -ort $DataPath/ts.2x.s1.$motfile \
	    -ort $DataPath/ts.3x.s1.$motfile \
	    -ort $DataPath/ts.4x.s1.$motfile \
	    -ort $DataPath/ts.5x.s1.$motfile \
	    -ort $DataPath/ts.6x.s1.$motfile \
	    -ort $DataPath/ts.1x.s1.sq.$motfile \
	    -ort $DataPath/ts.2x.s1.sq.$motfile \
	    -ort $DataPath/ts.3x.s1.sq.$motfile \
	    -ort $DataPath/ts.4x.s1.sq.$motfile \
	    -ort $DataPath/ts.5x.s1.sq.$motfile \
	    -ort $DataPath/ts.6x.s1.sq.$motfile \
	    -ort $DataPath/ts.1x.s2.$motfile \
	    -ort $DataPath/ts.2x.s2.$motfile \
	    -ort $DataPath/ts.3x.s2.$motfile \
	    -ort $DataPath/ts.4x.s2.$motfile \
	    -ort $DataPath/ts.5x.s2.$motfile \
	    -ort $DataPath/ts.6x.s2.$motfile \
	    -ort $DataPath/ts.1x.s2.sq.$motfile \
	    -ort $DataPath/ts.2x.s2.sq.$motfile \
	    -ort $DataPath/ts.3x.s2.sq.$motfile \
	    -ort $DataPath/ts.4x.s2.sq.$motfile \
	    -ort $DataPath/ts.5x.s2.sq.$motfile \
	    -ort $DataPath/ts.6x.s2.sq.$motfile" >> $ort_file
   else if ($MOT == "motsim") then
      set nuis_string = ${nuis_string}msim12
      echo "-ort $DataPath/motsim.both00.1D \
	    -ort $DataPath/motsim.both01.1D \
	    -ort $DataPath/motsim.both02.1D \
	    -ort $DataPath/motsim.both03.1D \
	    -ort $DataPath/motsim.both04.1D \
	    -ort $DataPath/motsim.both05.1D \
	    -ort $DataPath/motsim.both06.1D \
	    -ort $DataPath/motsim.both07.1D \
	    -ort $DataPath/motsim.both08.1D \
	    -ort $DataPath/motsim.both09.1D \
	    -ort $DataPath/motsim.both10.1D \
	    -ort $DataPath/motsim.both11.1D " >> $ort_file
   else if ($MOT == 0) then
      if ($DEBUG) echo "No motion regressors selected"
   else
      echo "Could not determine type of motion correction"
      continue
   endif

   if($BANDPASS) then
      set nuis_string = ${nuis_string}f
      set bp_cmd   = "-bandpass 0.01 0.1"
   else
      set bp_cmd   = ""
   endif

   if ($CENSOR) then
      set nuis_string = ${nuis_string}Xz${mot_thresh}
      set censor_cmd = "-censor $CensorFile -cenmode ZERO"
   else
      set censor_cmd = ""
   endif  
   
   #--- redo motion censor file in case it doesn't exist ---
   if (! -e ${CensorFilePrefix}_censor.1D) then
      echo "---Motion Censoring---"  |& tee -a $LOG_FILE
      1d_tool.py -infile $motfile -censor_motion $mot_thresh $CensorFilePrefix $cens_prev_string
   endif
   

   #--- save ort command for record ---
   /bin/cp -fp $ort_file $DataPath/ort_cmd.$file.$nuis_string.txt

   set file_out = ${file}.${nuis_string} 

   if (! -e $DataPath/$file_out.nii.gz) then

      #--- Save Mean before nuisance regression ---
      if (! -e Mean.$file.nii.gz) then
      3dTstat -mean -prefix Mean.$file.nii.gz $TmpPath/$file.nii.gz |& tee -a $LOG_FILE
      endif
   
   
      #----------------------------------------------------
      # Nuisance Regression
      #----------------------------------------------------
      3dTproject \
	 -input $TmpPath/$file.nii.gz \
	 -prefix $TmpPath/tmp.$file_out.nii.gz \
	 -mask $MaskTotal \
	 $censor_cmd \
	 #-blur $fwhm \
	 $bp_cmd \
	 `cat $ort_file`  |& tee -a $LOG_FILE

      #--- Add meanback ---
      3dcalc -a Mean.$file.nii.gz -b $TmpPath/tmp.$file_out.nii.gz -expr "a+b" -prefix $TmpPath/$file_out.nii.gz |& tee -a $LOG_FILE
      
   endif #file_out exists 
   
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) then
      rm -f $TmpPath/tmp.$file_out.nii.gz
      #rm -f $TmpPath/$file.nii.gz
   endif
   set file = $file_out
   
   #TMP
   #--- jump ahead ---
   if (-e ${file}S${fwhm}.mni.nii.gz) then
      set file = ${file}S${fwhm}.mni
      set MaskTotalPrefix = $MaskTotal:r:r
      set MaskTotal = $MaskTotalPrefix.mni.nii.gz
      goto CONNECTIVITY_MATRIX
   endif
   
   #====================================================================================
   BLUR_TO_FWHM:
   #------------------------------------------------------------------------------------
   #--- Spatial smoothing ---
   set file_out = ${file}S${fwhm}
   if (! -e $TmpPath/$file_out.nii.gz && ! -e $TmpPath/$file_out.nii) then
   if ($DEBUG) echo "Running Blur To FWHM ($fwhm)..."
   3dBlurToFWHM -input $TmpPath/$file.nii.gz -mask $MaskTotal -FWHM $fwhm -prefix $TmpPath/$file_out.nii.gz |& tee -a $LOG_FILE
   endif
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) rm -f $TmpPath/$file.nii.gz
   set file = $file_out
   
   #TMP
   #continue
   #goto CONNECTIVITY_MATRIX
   
   #====================================================================================
   RESAMPLE_TO_MNI:
   #------------------------------------------------------------------------------------
   #--- Resample to MNI ---
   set file_out = $file.mni
   if (! -e $file_out.nii.gz) then
   3dresample -master $RoiPath/MNI_avg152T1+tlrc \
              -input $TmpPath/$file.nii.gz -prefix $TmpPath/$file_out.nii.gz |& tee -a $LOG_FILE
   endif
   set MaskTotalPrefix = $MaskTotal:r:r
   if (! -e $MaskTotalPrefix.mni.nii.gz) then
   3dresample -master $RoiPath/MNI_avg152T1+tlrc \
              -input $MaskTotalPrefix.nii.gz -prefix $MaskTotalPrefix.mni.nii.gz |& tee -a $LOG_FILE
   endif
   if (-e $TmpPath/$file_out.nii.gz && $CLEANUP) rm -f $TmpPath/$file.nii.gz
   set file = $file_out
   set MaskTotal = $MaskTotalPrefix.mni.nii.gz

   #TMP
   goto CONNECTIVITY_MATRIX
   
   #====================================================================================
   CONNECTIVITY:
   #------------------------------------------------------------------------------------
   #--- Connectivity ---
   foreach seed (`cat $SeedFile`)
      if (! -e Fim.$seed.$file.nii.gz) then
      
         echo "Running connectivity for $seed ($file)..."
      
	 3dmaskave -mask $RoiPath/Mask.$seed+tlrc -quiet $TmpPath/$file.nii.gz > tmp.ts.$seed.$file.1D
	 rm -f ts.$seed.$file.1D
	 1dnorm -demean tmp.ts.$seed.$file.1D ts.$seed.$file.1D
	 rm -f tmp.ts.$seed.$file.1D

	 3dDeconvolve \
            -censor $CensorFile \
            -input $TmpPath/$file.nii.gz \
            -mask $MaskTotal \
            -num_stimts 1 \
            -stim_file 1 ts.$seed.$file.1D -stim_label 1 seed.$seed \
            -fout -tout -rout \
            -bucket Fim.$seed.$file.nii.gz |& tee -a $LOG_FILE
      else
         echo "Fim.$seed.$file.nii.gz already exists"
      endif
   end 
   
   #====================================================================================
   FISHER_Z:
   #------------------------------------------------------------------------------------
   #--- Fisher-Z transform ---
   foreach seed (`cat $SeedFile`)
   
      set FimFile = Fim.$seed.$file.nii.gz
      if (! -e tmp.Z.$FimFile) then
      
      if (! -e tmp.R.$FimFile) then
      3dcalc \
         -a $FimFile'[4]' \
         -b $FimFile'[3]' \
         -expr "sqrt(a)*(b/abs(b))" \
         -prefix tmp.R.$FimFile 
      endif
      if( ! -e tmp.Z.$FimFile) then
      3dcalc \
         -datum float \
         -a tmp.R.$FimFile \
         -expr "log((1+a)/(1-a))" \
         -prefix tmp.Z.$FimFile
      endif
      
      endif
   end
      
   #====================================================================================
   CONNECTIVITY_MATRIX:
   #------------------------------------------------------------------------------------
   set mask = $RoiPath/Parcels/Parcels_MNI_222_CommunitySort.nii 
   set out_label = GordonCS #Gordon Parcellation, sorted by community
   
   #TMP
   #rm -f $DataPath/CCmatrix.$out_label.$file.nii.gz
   if (! -e $DataPath/CCmatrix.$out_label.$file.nii.gz) then
      echo "Computing connectivity matrix ..."
      
      #--- remove censored time points ---
      $ScriptPath/@CensorChop $file.nii.gz $CensorFile $file.Xk.nii.gz
      
      #--- compute connectivity matrix ---
      $ScriptPath/@ConnMatrix.csh $file.Xk.nii.gz $mask $DataPath/CCmatrix.$out_label.$file.Xk.nii.gz
      
      #--- cleanup ---
      rm -f $file.Xk.nii.gz
   else
      echo "$DataPath/CCmatrix.$out_label.$file.nii.gz already exists"
   endif
   
   #--- Fisher-Z transform ---
   set CCmatrix = CCmatrix.$out_label.$file.nii.gz
   rm -f Z.$CCmatrix #redo - start clean
   if (! -e Z.$CCmatrix) then
   3dcalc \
      -datum float \
      -a $DataPath/$CCmatrix \
      -expr "log((1+a)/(1-a))" \
      -prefix $DataPath/Z.$CCmatrix
   endif
   
   
   #====================================================================================
   CLEANUP:
   #------------------------------------------------------------------------------------
   
end #dir
