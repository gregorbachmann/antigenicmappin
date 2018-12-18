# Antigenic Mapping

Repository storing the code used in my project on Antigenic Mappings.

## Requirements and Packages
The code is implemented in R using **RStudio Version 1.1.456**
The following R-packages are needed in order to run the code:
* **ggplot**
* **gridExtra**
* **phytools**
* **vegan**
* **plot3D** and **plot3plot3Drgl**
* **smacof**
* **MASS**
* **clusterCrit**
* **reshape**

## Structure of Code

The main structure of my code consists of files that contain functions and algorithms which get called in the script files.
A quick overview over the files:
* **EmbeddingAlgos**: Contains the multidimensional scaling algorithms
* **HelpFunctions**: Contains smaller functions that are needed to process data
* **Read_In_File**: Function that reads a .csv file and converts the data as needed.
* **ClusterCreator**: Contains all cluster-generating functions needed for the evaluation.
* **EvaluateAlgos**: Contains the functions for the frameworks discussed in my work.
* **EvaluationScript**: Evaluates all algorithms using the frameworks
* **PlottingScript** Produces embeddings for the HIV data.

## Form of data needed
Assume we have IC-50 distances between N viral strains and M antibodies.
I mainly worked with data in form of a list, containing the following:
* **Standard distance matrix**: A matrix of size MxN containing the pairwise distance between antibodies and viral strains
* **Symmetrized distance matrix**: The symmetrized version of size (M+N)x(M+N), with zeros on the diagonal and NAs for strain and antibody pairs.
* **Isolation time**: Vector of size M+N, containing isolation time measured in weeks for each antibody and strain (optional)
* **Labels**: Names of strains and antibodies (optional)
### EmbeddingAlgos
#### isomapvis
Relies heavily on the function **isomap** from **vegan**.
This function calculates the embedding by filling out the incomplete matrix using the neighbourhood graph. Note that **nn** needs to be chosen high enough in order to avoid fragmentation of the data.

**Arguments**:
* **data**: A list in the explained format.
* **dim**: The desired embedding dimension
* **n**: Number of nearest neighbours to consider for graph
* **plot**: Boolean, the plot is displayed if TRUE
* **label**: Boolean, labels are displayed in the plot if TRUE

**Returns**:
* Embedding returned by isomap
* Plot object


#### mdsvis
Relies heavily on the function **smacofSym** from the package **smacof**.
This function calculates the embedding using the classical metric scaling approaches discussed in the project.

**Arguments**:
* **data**: A list in the explained format.
* **dim**: The desired embedding dimension
* **method**: One of "ratio", "interval", "mspline" or "ordinal" as discussed in the project.
* **init**: "torgerson" for a deterministic initial configuration or "random" for a random one.
* **plot**: Boolean, the plot is displayed if TRUE.
* **label**: Boolean, labels are displayed in the plot if TRUE.
* **real**: Boolean, indicating whether data has time and name attribut.

**Returns**:
* Embedding returned by **smacofSym**
* Plot object


#### sammon
Relies heavily on function **sammon** from the package **MASS**.
This function calculates the embedding for the Sammon approach as discussed in the project.
**Arguments**:

* **data**: A list in the explained format.
* **dim**: The desired embedding dimension
* **plot**: Boolean, the plot is displayed if TRUE
* **label**: Boolean, labels are displayed in the plot if TRUE

**Returns**:
* Embedding returned by **sammon**
### HelpFunctions
#### generate_shape_space_coord

Function to get random initial configurations for sammon.

**Arguments**
* **n.v**: Number of viruses
* **n.a**: Number of antibodies
* **k**: Dimension that produced data should have

**Returns**
* Matrix containing the coordinates for antibodies and strains as rows.

#### calc_dist
Function that calculates the standard distance matrix given coordinates of the strains and antibodies.

**Arguments**:
* **coord**: Matrix containing coordinates of antibodies and strains
* **n.v**: Number of viruses
* **n.a**: Number of antibodies

**Returns**:
* Standard distance matrix

#### sym_matrix
Function that calculates the symmetrized distance matrix

**Arguments**:
* **dist_matrix**: Standard distance matrix

**Returns**:
* Symmetrized distance matrix

### Read_In_File
#### read_in
Function to read in the data from a csv file and returning it in the format as discussed above.

**Arguments**:
* **link**: String indicating the link to the csv file on your machine
* **n.v**: Number of viruses
* **n.a**: Number of antibodies
* **ordering** 0 if csv is ordered w.r.t. virus and 1 if ordered w.r.t. antibody

**Returns**:
* List of the form as disussed above

### ClusterCreator
#### create_test_clusters
This function implements the first step in the first framework in my project. Namely clusters are created in a high dimensional hypercube.

**Arguments**
* **numclusters**: Number of clusters desired
* **dim**: Dimension in which clusters should be formed.
* **numpoints**: Vector of length numcluster containing the number of samples belonging to each cluster.
* **plot**: TRUE if plot of the clusters should be produced, only works for dim=2.

**Returns**:
* List consisting of virus and antibody coordinates, cluster centers and memberships
* Plot of clusters if plot=TRUE

#### preprocess
Function to calculate the standard and symmetrized distance matrix for the cluster data produced by create_test_clusters

**Arguments**:
* cluster: List of the form returned by create_test_clusters

**Returns**:
* The two distances matrices as a list
