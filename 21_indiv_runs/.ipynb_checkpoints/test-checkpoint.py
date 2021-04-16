import fn

folder_name = 'TEST'
vessels_off = False

data = fn.import_and_set_up(folder_name,vessels_off)

# fn.plt_single(data,'Q_norm', '$Q_[norm]$')
fn.plt_single(data,'phi','$\phi$')
# fn.plt_single(data,'phi_min','$\phi_[min]$')
# fn.plt_single(data,'norm R_tot','Normalsied R')
print('singles')
fn.plt_triple(data,'Ap','Vp','Dp','Proportion of pericytes\nin state')
# fn.plt_triple(data,'An','Vn','Dn','Proportion of neurons\nin state')
print('doubles')
# fn.plt_base(data,'c','c',0.909090909)
# fn.plt_base(data,'pt_volume_averaged','$p_[t]$',44.113829801035656)

# fn.plt_norm(data,'c','c norm')
# fn.plt_norm(data,'pt_volume_averaged','$p_[t]$ norm')

# fn.plt_vessels(data,'partial pressure blood(mmHg)','partial pressure blood')
# fn.plt_vessels(data,'partial pressure blood(mmHg)','partial pressure blood',normal=1)


print('done')
