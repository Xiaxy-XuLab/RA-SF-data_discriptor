import scrublet as scr
import scanpy as sc
import pandas as pd
from pathlib import Path
import sys

#Load raw counts matrix and gene list
# inx = 'test.h5Seurat'
inx = sys.argv[1]
anaData = sc.read_h5ad(inx)
sample_index = Path(inx).stem 

#Initialize Scrublet object
scrub = scr.Scrublet(anaData.raw.X,
                     expected_doublet_rate=0.1,
                     sim_doublet_ratio=2,
                     n_neighbors = 8)

#Run the default pipeline
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=1, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=25)

# scrub.call_doublets(threshold=0.25)
#pd.DataFrame({'doublets_scores': doublet_scores, 'predicted_doublets': predicted_doublets, index = anaData.obs.index })
pd.DataFrame(anaData.obs.index[predicted_doublets]).to_csv(sample_index + '.csv', header = False, index = False)
