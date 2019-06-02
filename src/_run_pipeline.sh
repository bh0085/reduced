source activate py3
python  a_split.py  > split_cmds.txt
source split_cmds.txt
python b_alignment_rev.py
source ../qsubs/b_alignment_rev/_commands.txt
python c6_polish_rev.py
source ../qsubs/c6_polish_rev/_commands.txt
python e_newgenotype_rev.py
source ../qsubs/e_newgenotype_rev/_commands.txt