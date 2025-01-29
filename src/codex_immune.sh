#!/bin/bash

# cd /home/makowlg/Documents/Immune-CCI/src
#conda init
conda activate mkpy

python /home/makowlg/Documents/Immune-CCI/src/scripts/scripts_immune/dge_immune.py

python /home/makowlg/Documents/Immune-CCI/src/scripts/scripts_immune/new_immune.py

python /home/makowlg/Documents/Immune-CCI/src/scripts/scripts_immune/immune_canonical.py

python /home/makowlg/Documents/Immune-CCI/src/scripts/general_mako/excel_merge.py

Rscript /home/makowlg/Documents/Immune-CCI/src/scripts/general_mako/singlec_generate_tstat.r

python /home/makowlg/Documents/Immune-CCI/src/scripts/gsea/new_gsea.py

python /home/makowlg/Documents/Immune-CCI/src/scripts/gsea/gsea_summary.py
