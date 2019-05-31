source activate py3
python  a_split.py  > split_cmds.txt
source split_cmds.txt
python b_alignment_rev.py
source ../qsub/b_alignment_rev/_commands.sh
python c6_polish_rev.py
source ../qsub/c6_polish_rev/_commands.sh
python e_newgenotype_rev.py
source ../qsub/e_newgenotype_rev/_commands.sh