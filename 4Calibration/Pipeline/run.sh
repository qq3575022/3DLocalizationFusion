# detect stationary states and save to res.mat
python3 static.py

# input ground truth gravity value
#!/bin/bash
matlab -nodesktop -nosplash -nodisplay -r "run inputG.m"

# get acceleration calibrated
matlab -nodesktop -nosplash -nodisplay -r "run cali.m"


# verify calibrated results
matlab -nodesktop -nosplash -nodisplay -nodesktop -r "run veri.m"