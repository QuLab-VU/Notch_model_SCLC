from sklearn import datasets
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import sklearn.metrics as sm

#import some data
iris=datasets.load_iris()
x=iris.data
y=iris.target
#X_std = StandardScaler().fit_transform(x)

print(x, y)

# pca = PCA(n_components=2)
# data=pca.fit_transform(x,y)

kmeans = KMeans(n_clusters=3)
data = kmeans.fit_transform(x,y)

# plt.scatter(*zip(*data), c=y)
# plt.figure()
plt.scatter(*zip(*data), c=y,)
plt.show()