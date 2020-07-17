#gaussian_MP
## Create gaussian input .com from SMILES
put all the molecular SMILES into list\_smi in `smile2com_neu.py`,  `smile2com_pos.py` and `smile2com_neg.py`.  
And execute  
`python smile2com_neg.py`  
`python smile2com_neu.py`  
`python smile2com_pos.py`  
Next, upload all the .com file we generate to our machine.
## Submit the calculation job on our machine
`sh create_file_run.sh` please define the range of for-loop.
## Calculate IE and EA
After finishing calculation,  
`sh calculate_IE_EA`  
you will get a result file.
