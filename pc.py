#!/home/zrollins/miniconda3/bin/python
# coding: utf-8

# In[ ]:


#L1_90ns
import MDAnalysis as mda
#from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import pca, align
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nglview as nv
import warnings

warnings.filterwarnings('ignore')
#get_ipython().run_line_magic('matplotlib', 'inline')
u = mda.Universe('pc.gro','pc.xtc')
aligner = align.AlignTraj(u, u, select='backbone',
                          in_memory=True).run()
pc = pca.PCA(u, select='backbone',
             align=False, mean=None,
             n_components=None).run()
backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print('There are {} backbone atoms in the analysis'.format(n_bb))
print(pc.p_components.shape)
pc.variance[0]
plt.plot(pc.cumulated_variance[:10], color = 'purple')
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance');
plt.title('DMF5-MART1 Glycosylated 50 ns PCA Cumulative Variance')
plt.savefig('PCA_cv.png')

transformed = pc.transform(backbone, n_components=5)
transformed.shape
df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(5)])
df['Time (ps)'] = df.index * u.trajectory.dt
df.head()

import seaborn as sns

g = sns.PairGrid(df, hue='Time (ps)',
                 palette=sns.color_palette('Purples',
                                           n_colors=len(df)))
g.map(plt.scatter, marker='.')
g.savefig('PCA.png')

pc1 = pc.p_components[:, 0]
trans1 = transformed[:, 0]
projected = np.outer(trans1, pc1) + pc.mean
coordinates = projected.reshape(len(trans1), -1, 3)

proj1 = mda.Merge(backbone)
proj1.load_new(coordinates, order="fac")


pca1 = proj1.atoms
pca1.write('pc1.gro')
pca1.write('pc1.xtc', frames='all')

pc2 = pc.p_components[:, 1]
trans2 = transformed[:, 1]
projected2 = np.outer(trans2, pc2) + pc.mean
coordinates2 = projected2.reshape(len(trans2), -1, 3)

proj2 = mda.Merge(backbone)
proj2.load_new(coordinates2, order="fac")


pca2 = proj2.atoms
pca2.write('pc2.gro')
pca2.write('pc2.xtc', frames='all')

pc3 = pc.p_components[:, 2]
trans3 = transformed[:, 2]
projected3 = np.outer(trans3, pc3) + pc.mean
coordinates3 = projected3.reshape(len(trans3), -1, 3)

proj3 = mda.Merge(backbone)
proj3.load_new(coordinates3, order="fac")


pca3 = proj3.atoms
pca3.write('pc3.gro')
pca3.write('pc3.xtc', frames='all')

pc4 = pc.p_components[:, 3]
trans4 = transformed[:, 3]
projected4 = np.outer(trans4, pc4) + pc.mean
coordinates4 = projected4.reshape(len(trans4), -1, 3)

proj4 = mda.Merge(backbone)
proj4.load_new(coordinates4, order="fac")


pca4 = proj4.atoms
pca4.write('pc4.gro')
pca4.write('pc4.xtc', frames='all')

pc5 = pc.p_components[:, 4]
trans5 = transformed[:, 4]
projected5 = np.outer(trans5, pc5) + pc.mean
coordinates5 = projected5.reshape(len(trans5), -1, 3)

proj5 = mda.Merge(backbone)
proj5.load_new(coordinates5, order="fac")


pca5 = proj5.atoms
pca5.write('pc5.gro')
pca5.write('pc5.xtc', frames='all')

