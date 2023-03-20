import numpy as np
import pandas as pd

ens_num = 61
m = 10
typ = 'corr_v10-2_sumrsamp30'

n = 0
i = 0
hop_syn = pd.read_feather('~/syn_forecast_test/output/%s/feather/samp%d/syn-hop_%d.feather' %(typ,n+1,i+1))

hop_new = np.empty((np.shape(hop_syn)[0],np.shape(hop_syn)[1],ens_num))
lam_new = np.empty((np.shape(hop_syn)[0],np.shape(hop_syn)[1],ens_num))
uka_new = np.empty((np.shape(hop_syn)[0],np.shape(hop_syn)[1],ens_num))

for n in range(m):
    for i in range(ens_num):
        hop_syn = pd.read_feather('~/syn_forecast_test/output/%s/feather/samp%d/syn-hop_%d.feather' %(typ,n+1,i+1))
        hop_new[:,:,i] = hop_syn
        lam_syn = pd.read_feather('~/syn_forecast_test/output/%s/feather/samp%d/syn-lam_%d.feather' %(typ,n+1,i+1))
        lam_new[:,:,i] = lam_syn
        uka_syn = pd.read_feather('~/syn_forecast_test/output/%s/feather/samp%d/syn-uka_%d.feather' %(typ,n+1,i+1))
        uka_new[:,:,i] = uka_syn
        
    np.savez_compressed('%s/npz/LkMendoSynHcst_%s_samp%d.npz' %(typ,typ,n+1), syn_hcstLm = lam_new, syn_hcstWf = uka_new, syn_hcstHop = hop_new)
    

#END