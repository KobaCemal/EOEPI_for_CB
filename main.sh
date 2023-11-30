#Run fmriprep on docker 
fmriprep-docker rest_blnd_can_raw/ rest_blnd_can_raw/derivatives/ participant --participant-label sub-sccb01 --ignore slicetiming --fs-license-file /home/koba/freesurfer/license.txt --fs-no-reconall

# Resample anat files to 1x1x1, and create the mean image (for possible later use)
cd derivatives/fmriprep
for i in $(find . -type f | grep space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz | sort ); do name=$(echo $i | cut -d . -f 2); 3dresample -prefix ."$name"_1mm.nii.gz -dxyz 1 1 1 -input $i;done
cd ../..
3dbucket -prefix code/allT1s.nii.gz $(find . -type f | grep pace-MNI152NLin2009cAsym_desc-preproc_T1w_1mm.nii.gz | sort ) -overwrite
3dTstat -prefix code/allT1s_mean.nii.gz code/allT1s.nii.gz

# Band-pass filtering
for sub in $(cat code/sublist.txt); do echo "================================================================================== $sub"; 3dBandpass -prefix  derivatives/fmriprep/"$sub"/func/"$sub"_filtered.nii.gz -band 0.01 0.1 -input derivatives/fmriprep/"$sub"/func/"$sub"_task-rest_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz -overwrite; done


# Create the ortvec file
for sub in $(cat code/sublist.txt); do echo "================================================================================== $sub"; awk ' NR==1 { for (i=1; i<=NF; i++) {f[$i] = i }}{ print $(f["csf"]), $(f["white_matter"]), $(f["rot_x"]), $(f["rot_y"]), $(f["rot_z"]), $(f["trans_x"]), $(f["trans_y"]), $(f["trans_z"])} ' derivatives/fmriprep/"$sub"/func/"$sub"_task-rest_desc-confounds_timeseries.tsv > derivatives/fmriprep/"$sub"/func/"$sub"_ortvec_temp.txt; sed 's~n/a~0~g' derivatives/fmriprep/"$sub"/func/"$sub"_ortvec_temp.txt > derivatives/fmriprep/"$sub"/func/"$sub"_ortvec_temp1.txt;tail -n +2 derivatives/fmriprep/"$sub"/func/"$sub"_ortvec_temp1.txt  > derivatives/fmriprep/"$sub"/func/"$sub"_ortvec.txt; rm derivatives/fmriprep/"$sub"/func/"$sub"_ortvec_temp*; done

# Clean the resting state data
for sub in $(cat code/sublist.txt); do 3dDeconvolve -jobs 6 -float -input derivatives/fmriprep/"$sub"/func/"$sub"_filtered.nii.gz -polort -1 -num_stimts 0  \
-ortvec derivatives/fmriprep/"$sub"/func/"$sub"_ortvec.txt \
-noFDR -tout -bout -nofullf_atall -nodmbase -nfirst 0 \
-errts derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts.nii.gz \
-bucket derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_bucket.nii.gz \
-overwrite; done

# Smooth
for sub in $(cat code/sublist.txt); do 3dBlurToFWHM -prefix derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz  -FWHM 6 -input derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts.nii.gz; done

# Resample the eye mask file to func
3dresample -prefix code/eyemask_res.nii.gz -master derivatives/fmriprep/sub-cb01/func/sub-sub-cb01_cleanRS_errts_sm.nii.gz -input code/MNI_1mm_bothEyes.nii.gz -overwrite

# Extract the mean eye signal and convolve
for sub in $(cat code/sublist.txt); do echo "================================================================================== $sub"; 3dmaskave -quiet -mask code/eyemask_res.nii.gz derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz  > derivatives/fmriprep/"$sub"/func/"$sub"_eyes_raw_sm.txt; waver -numout 132 -dt 2.2 -input derivatives/fmriprep/"$sub"/func/"$sub"_eyes_raw_sm.txt -WAV >  derivatives/fmriprep/"$sub"/func/"$sub"_eyes_wav_sm.txt; done

# Find the correlates of eyes 
for sub in $(cat code/sublist.txt); do 3dDeconvolve -jobs 6 -float -input derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz -polort -1 -num_stimts 1  \
-stim_file 1 derivatives/fmriprep/"$sub"/func/"$sub"_eyes_wav_sm.txt -stim_base 1 -stim_label 1 eyes \
-noFDR -tout -bout -nofullf_atall -nodmbase -nfirst 0 \
-errts derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_errts_sm_wav.nii.gz \
-bucket derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_bucket_sm_wav.nii.gz; done

for sub in $(cat code/sublist.txt); do 3dDeconvolve -jobs 6 -float -input derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz -polort -1 -num_stimts 1  \
-stim_file 1 derivatives/fmriprep/"$sub"/func/"$sub"_eyes_raw_sm.txt -stim_base 1 -stim_label 1 eyes \
-noFDR -tout -bout -nofullf_atall -nodmbase -nfirst 0 \
-errts derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_errts_sm_raw.nii.gz \
-bucket derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_bucket_sm_raw.nii.gz; done

#Gather the necessary coefficients into one file
mkdir groupanalyses
mkdir groupanalyses/replication
mkdir groupanalyses/replication/sightedonly
for sub in $(cat code/sublist.txt); do 3dbucket -prefix derivatives/fmriprep/"$sub"/func/"$sub"-eyeCoefWav.nii.gz derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_bucket_sm_wav.nii.gz[0] -overwrite;done
for sub in $(cat code/sublist.txt); do 3dbucket -prefix derivatives/fmriprep/"$sub"/func/"$sub"-eyeTWav.nii.gz derivatives/fmriprep/"$sub"/func/sub-"$sub"_eyes_bucket_sm_wav.nii.gz[1] -overwrite;done
fslmerge -t groupanalyses/replication/eyeCoefsSightedOnly.nii.gz $(find derivatives/fmriprep/ -type f | grep eyeCoefWav | sort  | tail -25)
fslmerge -t groupanalyses/replication/eyeCoefsBlindSighted.nii.gz $(find derivatives/fmriprep/ -type f | grep eyeCoefWav | sort)
fslmerge -t groupanalyses/replication/eyeCoefsBlindSighted.nii.gz $(find derivatives/fmriprep/ -type f | grep eyeCoefWav | sort | grep -v sccb01 | grep -v sclb02 | grep -v sub-lb01)

fslmerge -t groupanalyses/replication/eyeCoefsBlindOnlyRaw.nii.gz  $(find derivatives/fmriprep/ -type f | grep eyeCoefraw.ni | sort | grep sub-cb)
# Do the statistical analysis
randomise -i groupanalyses/replication/eyeCoefsSightedOnly.nii.gz -o groupanalyses/replication/sightedonly -1 -n 10000 -D -T -v 4

randomise -i groupanalyses/replication/eyeCoefsBlindSighted_raw.nii.gz -o groupanalyses/anova_nonparam_raw -m derivatives/fmriprep/sub-cb02/func/sub-cb02_task-rest_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz -d code/design.mat -t code/design.con -T -n 10000 -v 4

# Do the analysis with seperate eyes

for sub in $(cat code/sublist.txt);
do
	for side in right left
	do
		echo "================================================================================== $sub     $side " ;
		3dmaskave -quiet -mask code/eyes_"$side"_res.nii.gz derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz  > derivatives/fmriprep/"$sub"/func/"$sub"_eyes_raw_sm_"$side".txt
		waver -numout 132 -dt 2.2 -input derivatives/fmriprep/"$sub"/func/"$sub"_eyes_raw_sm_"$side".txt -WAV >  derivatives/fmriprep/"$sub"/func/"$sub"_eyes_wav_sm_"$side".txt

		for type in raw wav
		do
			 do
			 	3dDeconvolve -jobs 6 -float -input derivatives/fmriprep/"$sub"/func/sub-"$sub"_cleanRS_errts_sm.nii.gz -polort -1 -num_stimts 1  \
				-stim_file 1 derivatives/fmriprep/"$sub"/func/"$sub"_eyes_"$type"_sm_"$side".txt -stim_base 1 -stim_label 1 eyes \
				-noFDR -tout -bout -nofullf_atall -nodmbase -nfirst 0 \
				-errts derivatives/fmriprep/"$sub"/func/"$sub"_eyes_"$type"_sm_"$side"_errts.nii.gz \
				-bucket derivatives/fmriprep/"$sub"/func/"$sub"_eyes_"$type"_sm_"$side"_bucket.nii.gz
		done
	done
done
