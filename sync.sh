# obtain up to date bmrb entry list
rsync -avh rsync://bmrb.io:/bmrb_entries/ bmrb_entries

# soft link all star3 files to separate directory
cp -s bmrb_entries/bmr*/bmr*_3.str star3_files/

# watch, make sure protein is in unbuffered -u
watch -d 'grep -v bmr foo | sort | uniq -c | sort -nrk 1'