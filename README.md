# image-registration
### Image registration, clustering, point matching, and quantification for correlative STORM imaging
#### Combined live- and fixed-animal imaging
Correlative imaging can enable new insights by harnessing the strengths of multiple imaging techniques.  Meaningful  interpretation of correlated images is dependent on their accurate registration (alignment).  This project provides a framework for the registration and subsequent analysis of point cloud images from different sources.  Specifically,  calcium fluorescence localization images of live-organism *Drosophila* larva (corresponding to neural synaptic vesicle  release events) are correlated with STORM images of active zone protein scaffolds.
### Dependencies:
- [STORM image I/O](https://github.com/sjkenny/common) for STORM molecule lists generated with Insight3
- [DBSCAN and flood fill clustering](https://github.com/sjkenny/clustering) for feature extraction

### Required inputs:
Conventional fluorescence images of the same target are required from each microscope for initial registration.  Point clouds may also be binned into pixels and used for initial registration if fluorescence reference images are not available. Each point cloud must have three fields: x,y coordinates, as well as local density (code for local density calculation available [here]).

### Sample workflow:
#### 1. Manual affine image registration by control point selection (img_reg_similarity.m)
#### 2. Identification of active zones in each image by DBSCAN clustering (reg_clustering_affine_v3.m)
#### 3. [Point matching for non-rigid registration (Zheng 2006)](http://ieeexplore.ieee.org/abstract/document/1597120/) (reg_clustering_affine_v3.m)
#### 4. Interpolation of matched clusters by vector field spline fitting (interpolate_matches.m)
#### Additional scripts included for post-registration quantification of matched sites:
- Jarvis march-based neuron outlining (compute_mask.m)
- Cluster alignment and averaging (reg_clustering_affine_v3.m)
- Protein count normalization and quantification (quantify_matches.m)


Fluorescent images of a membrane label are taken on each microscope as references for initial alignment. Hand selection of control points provides an initial rigid transformation.  This is mostly necessary to crop the two images to the same size and shape, in order to improve the correspondence between features extracted in each image.  Exact overlap is unnecessary.

![Affine transformation](https://i.imgur.com/q9T03eA.png)

The resulting transformation is applied to both point cloud images:

![Registered point clouds](https://i.imgur.com/egHURyY.png)

[Clusters are extracted](https://github.com/sjkenny/clustering) from each point cloud using either binning and flood-fill or dbscan algorithms. Cluster centers-of-mass are matched between the images with a method which preserves local neighborhood geometry, but allows for shape variations across longer distances, as described by [Zheng, et. al](http://ieeexplore.ieee.org/abstract/document/1597120/)  Crucially, this method is tolerant to outliers in one direction, which is often the case in our experiments - each STORM cluster has a matching calcium-localization cluster, but not vice versa.

![Clustered point clouds](https://i.imgur.com/9BW2D7N.png)

Since not all STORM clusters have a corresponding match, interpolation is used to map them to regions on the live-animal image.  This enables site quantification even when there are insufficient calcium localizations to be detected as a cluster.

![Vector field](https://i.imgur.com/PzmQrlm.png)
