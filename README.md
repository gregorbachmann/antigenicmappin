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

**Arguments**:
* **data**: A list in the explained format.
* **dim**: The desired embedding dimension
* **n**: Number of nearest neighbours to consider for graph
* **plot**: Boolean, the plot is displayed if TRUE
* **label**: Boolean, labels are displayed in the plot if TRUE

**Returns**:
* Embedding returned by isomap
* Plot object

This function calculates the embedding by filling out the incomplete matrix using the neighbourhood graph. Note that **nn** needs to be chosen high enough in order to avoid fragmentation of the data.

#### mdsvis
Relies heavily on the function **smacofSym** from the package **smacof**.

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
