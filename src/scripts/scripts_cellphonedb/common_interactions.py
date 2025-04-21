# This script will take the significant interations from cellphonedb output and its gonna be removing the interations from the uninjered/sham file 
# which will only be left with the injured new ones 
#
# Goncalo Alves Msc Thesis 2025

import pandas as pd
import scanpy as sc
import numpy as np

def load_tables():
 """
    This function will load the txts that contain the significant interactions from the 3 condictions (injured15, injured60 and sham/uninjured)
    Return will be the tables in a df

"""



def remove_rows():
  """
    This function will remove the unwanted rows from each df 
    (we only want the interations that are only significant and unique in the injured condictions)
    
    Ent a minha ideia é tendo em conta o facto que o txt dos pvalues tem todas as interacoes que tem e n tem pvalues decentes
    a ideia seria primeiro remover todas as rows que n tem pelo menos um pvalues menor que 0.05.
    dps podiamos isolar as interacoes itself, ou seja temos uma row com um id, e dps remover todas as colunas que tem o pvalue superior a 0.05.
    dps o proximo step seria criar um df com cada id de interacao por row e as colunas seriam os conjuntos de clusters que sao detetados as interacoes.
   e dps ias ver os ids e clusters que se verificam no uninjured e ver quais os que se repetem os outros e se sim retirar
"""
