# obtain up to date bmrb entry list
rsync -avh rsync://bmrb.io:/bmrb_entries/ bmrb_entries

# soft link all star3 files to separate directory
ln -s bmrb_entries/bmr*/bmr*_3.str star3_files/