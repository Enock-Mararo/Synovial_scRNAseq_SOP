# scRNA-seq analysis: Old vs. Young mice synovial tissues

### **Prerequisite: Server & Environment**
* **todata2:** Supports Seurat v5 natively.
* **todata3:** Default environment is locked to Seurat v4.3.
* **Solution:** This SOP uses a **Miniconda** virtual environment to execute Seurat v5 on `todata3`.

---
### 1. Setup Server and Environment
**Objective:** Establish a secure connection to the MVLS computational server (`todata3`) from a local Windows machine and configure a protected environment for Seurat v5.

#### Accessing the todata3 Server via SSH
Because sequencing data must remain secure, ensure you are connected to the University of Glasgow VPN (Cisco Secure Client). Open the Windows PowerShell terminal and initiate a Secure Shell (SSH) connection using your GUID:

```bash
ssh yourGUID@todata3.mvls.gla.ac.uk
```

***


### 2. Activate the Miniconda Environment
Once logged into the server, bypass the default R environment by activating the specific Seurat v5 container:

```bash
# Activate the protected Miniconda environment containing Seurat v5
conda activate seurat5_env
```

> **🔍 Sanity Check:** Your terminal prompt must change from `(base)` to `(seurat5_env)` before proceeding to the next step.

***
### 2.1 Installation of Analytical Dependencies (Internal to Environment)
To ensure research continuity and computational reproducibility without requiring system-wide administrative (`sudo`) permissions, all specialized libraries are installed directly into this protected environment. Run these commands sequentially within your activated terminal:

```bash
# 1. Install the core Seurat v5 suite and essential object handlers
conda install -c conda-forge r-seurat=5.0.0 r-seuratobject -y

# 2. Install performance-enhancing libraries for large matrix operations
conda install -c conda-forge r-glmgampoi r-matrix -y

# 3. Install batch-correction and high-resolution visualization tools
conda install -c conda-forge r-harmony r-patchwork r-dplyr -y
```

> 💡 **Why this is required:** By installing these within the `seurat5_env`, we establish a self-contained ecosystem. This prevents "Namespace" conflicts with the server's default R libraries and allows the student to manage a professional-grade single-cell pipeline independently.

***

### **2.2 Technical Breakdown of Dependencies**
Understanding the specific role of each library is essential for troubleshooting and experimental design. The following dependencies form the backbone of the synovial tissue analysis:

#### **Core Analysis Engine**
* **`Seurat (v5.0.0)`:** The primary analytical framework. Version 5 introduces `Assay5` objects, allowing for more efficient memory handling of the massive synovial count matrices.
* **`SeuratObject`:** Provides the underlying data structures. It ensures that metadata (like age or treatment group) remains synchronized with the gene expression data during complex transformations.

#### **Performance & Mathematical Optimization**
* **`Matrix`:** Handles "sparse matrices." Since single-cell data is mostly zeros (genes not expressed in a specific cell), this library stores only the non-zero values, reducing RAM usage by ~90%.
* **`glmGamPoi`:** A fast function used during normalization. It utilizes a Gamma-Poisson generalized linear model to fit the data, significantly outperforming standard methods in both speed and accuracy for large-scale datasets.

#### **Integration & Batch Correction**
* **`Harmony`:** A high-performance algorithm for batch correction. It projects cells into a shared low-dimensional space and iteratively clusters them to remove technical "batch effects" (like different sequencing runs) while preserving biological variation.
* **`dplyr`:** The industry standard for data manipulation. It allows for fast, intuitive filtering and renaming of cell identities within the Seurat metadata.

#### **Visualization Architecture**
* **`ggplot2`:** The foundational graphics engine. It allows for the layer-by-layer construction of professional, publication-quality figures.
* **`patchwork`:** A utility for layout management. It enables the seamless combining of multiple plots (e.g., side-by-side UMAPs of Aged vs. CIA) using simple mathematical operators like `+` or `/`.

***

### 3. Registering the Jupyter Kernel
To use this environment visually in the browser, we must build a computational "bridge" (the IRkernel). Inside the active environment in PowerShell, start R and register the kernel:

```bash
# Start the R terminal
R
```

```R
# Register the environment to Jupyter
IRkernel::installspec(name = 'seurat5_env', displayname = 'R 4.3 (PRO V5)')

# Quit R
q()
```

### 4. Accessing the Jupyter Web Interface
Now, transition from the PowerShell terminal to your standard web browser (e.g., Chrome, Edge).

1. Navigate to the server address: http://todata3.mvls.gla.ac.uk:8000/
2. You will be prompted to log in. Enter your standard university GUID and password.
3. Once logged into the Jupyter dashboard, create a new notebook. Ensure you select **R 4.3 (PRO V5)** from the kernel dropdown menu located at the top right of the screen.

### 5. High-Memory Capacity Reservation
Inside your new Jupyter Notebook, run this block immediately in the first cell to prevent "memory exhausted" crashes during downstream integration:

```R
options(future.globals.maxSize = 800 * 1024^3)
Sys.setenv("R_MAX_VSIZE" = "900Gb")
gc()
cat(" ✅ High memory reserved successfully. Workspace is clean.")
```

### 6. Load the Data
**Objective:** Securely transfer raw 10X Genomics sequencing matrices from MARS and construct the initial Seurat object.

#### Transferring 10X data from MARS
Return to your PowerShell terminal. We use `rsync` instead of a standard copy command because `rsync` will automatically resume the transfer if the connection drops, preventing corrupted data. 

```bash
rsync -avz /mnt/autofs/data/userdata/project0067/Lab_Data/mouse_ageing/ ~/UofG_SingleCell_Enock/mouse_ageing_data/
```

> **🔍 Sanity Check:** Run `du -sh ~/UofG_SingleCell_Enock/mouse_ageing_data/` to confirm the folder is multiple Gigabytes in size, proving the transfer finished completely.

#### Defining the Root Directory in Jupyter
Moving back to your Jupyter Notebook, set the "Main Address" for your data. We map the specific subfolders to ensure the raw matrices are accurately pulled for both the Aged cohort and the Young CIA controls.

```R
base_path <- "~/UofG_SingleCell_Enock/mouse_ageing_data/"

# Map the exact folder names from the MARS transfer
sample_map <- list(
  "Old1"         = "OLD1/filtered_feature_bc_matrix/",                             # Path to Aged Sample 1
  "Old2"         = "OLD2/filtered_feature_bc_matrix/",                             # Path to Aged Sample 2
  "Old3"         = "OLD3/filtered_feature_bc_matrix/",                             # Path to Aged Sample 3
  "Con1_CD45neg" = "CIA_controls/Con1_CD45neg/filtered_feature_bc_matrix/",         # Path to Young Control 1
  "Con2_CD45neg" = "CIA_controls/Con2_CD45neg/filtered_feature_bc_matrix/",         # Path to Young Control 2 
  "Con3_CD45neg" = "CIA_controls/Con3_CD45_neg/filtered_feature_bc_matrix/",        # Path to Young Control 3 (Note the underscore)
  "Con4_CD45neg" = "CIA_controls/Con4_CD45_neg/filtered_feature_bc_matrix/"         # Path to Young Control 4 (Note the underscore)
)
```

#### Reading and Merging the Matrices
We use an automated loop to read every folder and convert them into Seurat objects simultaneously.

```R
library(Seurat)

seurat_list <- list()
for (sample_name in names(sample_map)) {
  full_path <- paste0(base_path, sample_map[[sample_name]])
  counts <- Read10X(data.dir = full_path)
  seurat_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
}

# Merge the first object with the remaining six to create a single master dataset
combined_mouse <- merge(x = seurat_list[[1]], 
                        y = seurat_list[2:7], 
                        add.cell.ids = names(sample_map))
```

#### 📌 Checkpoint 2: Verify the dimensions of the initial merged cohort
```R
# Output the total number of genes (rows) and total number of cells (columns) in the raw object
dim(combined_mouse)

# Output the initial cell count distribution across the 7 biological replicates
table(combined_mouse\$orig.ident)

# Visually audit the quality control metrics using a Violin Plot
VlnPlot(combined_mouse, features = c("nFeature_RNA", "percent.mt"), ncol = 2)
```

<img width="4200" height="2400" alt="Fig1_Adaptive_QC_Audit" src="https://github.com/user-attachments/assets/2a1746f8-0686-47ee-b3e1-95f9ca6374ce" />

### 7. Quality Control (QC) Filters
**Objective:** Remove technical artifacts (empty droplets, multiplets, and apoptotic cells) and exclude sub-optimal samples to establish a statistically sound 3-vs-3 comparative cohort.

#### Calculating Mitochondrial RNA
During cell preparation, dying or stressed cells often leak cytoplasmic RNA while retaining mitochondrial RNA. We quantify the percentage of mitochondrial transcripts in every cell to identify and remove these stressed populations.

```R
# Calculate the percentage of reads that map to the mitochondrial genome 
# The "^mt-" pattern is specific to murine (mouse) genomic nomenclature
combined_mouse[["percent.mt"]] <- PercentageFeatureSet(combined_mouse, pattern = "^mt-")
```

#### Removing the Con2 Outlier and Applying Cutoffs
Before proceeding to downstream biological analysis, we must apply strict mathematical thresholds to clean the data matrix:

1. **Gene Capture Limits:** Cells must possess between 500 and 6,000 unique genes (`nFeature_RNA`). Counts below 500 indicate ambient RNA debris (empty droplets), while counts above 6,000 strongly suggest technical multiplets (two cells captured in one droplet).
2. **Mitochondrial Limit:** Cells must have `< 10%` mitochondrial RNA to ensure we only analyze healthy, viable cells.
3. **Cohort Balancing (Dropping Con2):** The `Con2_CD45neg` sample displayed critically low transcriptional complexity (a median gene capture of ~415 genes) and elevated stress markers. It is excluded entirely at this stage to prevent it from skewing the high-resolution fibroblast sub-clustering, finalizing a balanced 3-vs-3 experimental design (Young vs. Old).

We apply all of these rules simultaneously using the `subset` function:

```R
# Apply the specific QC parameters and explicitly filter out the Con2 sample
combined_mouse_clean <- subset(combined_mouse, 
                               subset = nFeature_RNA > 500 & 
                                        nFeature_RNA < 6000 & 
                                        percent.mt < 10 & 
                                        orig.ident != "Con2_CD45neg")
```

#### 📌 Checkpoint 3: Quantify Cell "Survival Rate"
Verify the cell counts after applying the QC thresholds to ensure that the remaining samples have sufficient cell depth for batch integration.

```R
# View the cell counts BEFORE filtering
cat("--- Cells BEFORE QC ---\n")
table(combined_mouse\$orig.ident)

# View the cell counts AFTER filtering (Con2 should be completely absent)
cat("\n--- Cells AFTER QC ---\n")
table(combined_mouse_clean\$orig.ident)

# Expected Output Verification:
# Con1  Con3  Con4  Old1  Old2  Old3 
#  411  3638  4813  2941  4534  5670 

```
***
### 8. Save the "Clean" Object
We serialize the filtered, batch-balanced object to a protected directory. This `.rds` (R Data Serialized) file serves as the definitive foundational starting point for the comparative analysis between different integration methods, such as CCA and Harmony.

```R
# Define the path for the protected results storage directory
results_path <- "~/UofG_SingleCell_Enock/results/CCA_vs_Harmony/"

# Ensure the results directory is created if it does not already exist
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# Serialize the refined Seurat object for long-term version control and storage
saveRDS(combined_mouse_clean, file = paste0(results_path, "Filtered_3vs3_Synovial_Object.rds"))

# Output a final confirmation of successful cohort refinement
cat("\n✅ Cohort refined. Con2 omitted and technical noise trimmed.")
```

***
### 9. Normalize and Run PCA
**Objective:** Prepare the filtered count matrices for dimensionality reduction by standardizing variance and identifying the key genes driving biological heterogeneity.

#### Log-Normalization and Finding Variable Genes
Raw sequence counts cannot be directly compared because some cells are sequenced more deeply than others. We apply Log-Normalization to mathematically standardize the gene expression across all cells. Following this, we isolate the top 2,000 highly variable genes (HVGs) that represent the true biological differences between our cell types, ignoring genes that are flat or uninformative across the cohort.

#### Scaling and Running PCA
Before calculating the Principal Components, the data must be scaled. This ensures that genes with naturally high expression (like housekeeping genes) do not mathematically drown out critically important, but low-expressing, transcription factors. Once scaled, we run Principal Component Analysis (PCA) to compress the thousands of variable genes into the top dimensions representing the majority of the variance.

```R
# Apply Log-Normalization to correct for sequencing depth bias
combined_mouse_clean <- NormalizeData(combined_mouse_clean, 
                                      normalization.method = "LogNormalize", 
                                      scale.factor = 10000)

# Identify the top 2,000 highly variable features for downstream clustering
combined_mouse_clean <- FindVariableFeatures(combined_mouse_clean, 
                                             selection.method = "vst", 
                                             nfeatures = 2000)

# Scale the data matrix to give equal weight to all variable features
combined_mouse_clean <- ScaleData(combined_mouse_clean)

# Execute Linear Dimensionality Reduction (PCA)
combined_mouse_clean <- RunPCA(combined_mouse_clean, 
                               features = VariableFeatures(object = combined_mouse_clean),
                               npcs = 50)
```

#### 📌 Checkpoint 5: Generate an Elbow Plot
This plot visualizes the standard deviation of each Principal Component. It is required to statistically justify how many dimensions will be fed into the final integration algorithm.

```R
# Generate the Elbow Plot to determine the dimensionality of the dataset
ElbowPlot(combined_mouse_clean, ndims = 50)
```

<img width="2400" height="1800" alt="Fig2_PCA_Elbow_Audit" src="https://github.com/user-attachments/assets/40a697f8-064b-4f35-a0b4-63d01226522f" />

***

### 10. CCA Integration
**Objective:** Align the distinct experimental batches (Young CIA vs. Aged) to identify shared biological cell states and remove technical batch variance.

#### Why We Use CCA (30 PCs)
Because our samples come from distinct biological groups and batches, we must mathematically integrate them. We utilize Canonical Correlation Analysis (CCA) rather than Harmony for this specific pipeline. CCA is an anchor-based approach that aggressively identifies shared biological structures (Mutual Nearest Neighbors) across the batches. 

> ⚠️ **Note on Dimensionality:** We strictly use 30 Principal Components (PCs) for CCA. Because CCA is mathematically aggressive, using too many PCs (e.g., 50) can result in "over-correction," where distinct biological populations, like subtle fibroblast sub-states, are artificially crushed together.

#### Running the Integration
Using the Seurat v5 framework, we execute the integration directly on the normalized layers.

```R
# Execute Canonical Correlation Analysis (CCA) Integration using 30 PCs
combined_mouse_clean <- IntegrateLayers(
  object = combined_mouse_clean, 
  method = CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  dims = 1:30, 
  verbose = FALSE
)
```

#### Joining Layers (Crucial for Seurat v5)
**CRITICAL STEP:** In Seurat v5, the matrices are physically split by sample (e.g., `counts.Old1`, `counts.Con3`) to allow for integration. However, downstream Differential Expression (DE) algorithms cannot read across split layers. We must physically collapse these individual matrices back into a single unified matrix.

```R
# Re-unify the raw sequence matrices for downstream marker identification
combined_mouse_clean[["RNA"]] <- JoinLayers(combined_mouse_clean[["RNA"]])
```

#### Generating the UMAP and Clusters
With the data integrated, we project the cells into a 2D space (UMAP). We cluster the cells using a Resolution of 0.6, which is the community standard for a cohort of this size (~20,000 cells) to accurately separate major fibroblast and myeloid lineages without over-clustering technical noise.

```R
# Run UMAP using the newly integrated CCA reduction
combined_mouse_clean <- RunUMAP(combined_mouse_clean, reduction = "integrated.cca", dims = 1:30)

# Build the network graph and find clusters
combined_mouse_clean <- FindNeighbors(combined_mouse_clean, reduction = "integrated.cca", dims = 1:30)
combined_mouse_clean <- FindClusters(combined_mouse_clean, resolution = 0.6)
```

#### 📌 Checkpoint 6: Visually Verify the Integration
Generate side-by-side UMAPs to confirm that the Young and Old samples are mixing evenly (left plot) and that distinct biological clusters have formed (right plot).

```R
# Generate the validation plots
p_cca_sample <- DimPlot(combined_mouse_clean, reduction = "umap", group.by = "orig.ident") + 
                labs(title = "CCA: Sample Alignment")
p_cca_cluster <- DimPlot(combined_mouse_clean, reduction = "umap", label = TRUE) + 
                 NoLegend() + labs(title = "CCA: Clusters (Res 0.6)")

# Display plots side-by-side
p_cca_sample | p_cca_cluster
```
<img width="4200" height="1800" alt="Fig3_CCA_Integration_UMAP" src="https://github.com/user-attachments/assets/67e735d7-7d0d-47e0-9e7c-ede8bfaa6ea3" />

***
## **11. Save and Organize the Object**
**Objective:** Archive the fully integrated, batch-corrected cohort to a persistent server directory, establishing a stable starting point for biological interpretation.

### **Organizing the Final .rds File**
Once the Seurat object is serialized into an `.rds` file, it is transferred to a dedicated, clearly labeled project directory. This centralized storage prevents accidental deletion, ensures strict version control, and allows incoming researchers to bypass upstream processing.

```r
# Define the absolute path of the original processed object for internal tracking
# Path: /home/emm44q/UofG_SingleCell_Enock/results/CCA_vs_Harmony/Filtered_3vs3_Synovial_Object.rds
old_file_path <- path.expand("~/UofG_SingleCell_Enock/results/CCA_vs_Harmony/Filtered_3vs3_Synovial_Object.rds")

# Define the standardized student project directory
target_dir <- path.expand("~/Mahony_Lab/scRNA_seq_data/Eidan")

# Create the protected directory structure if it does not exist
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# Define the new target path using the finalized lab naming convention
# Full Absolute Path: /home/emm44q/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds
new_file_path <- paste0(target_dir, "/Synovial_3vs3_CCA_Integrated.rds")

# Execute the file transfer and rename the object for better discoverability
file.copy(from = old_file_path, to = new_file_path, overwrite = TRUE)

# Output the confirmation of successful archiving
cat(" ✅ Data archived successfully to:", new_file_path)
```

> **Quick Start for Students:** If you need to begin biological analysis immediately without re-running the full pipeline, use the command below.

### **⚡ Fast-Track: Load Integrated Data (Read-In-Place)**
This command allows you to bypass the upstream data processing and immediately load the batch-corrected master object into your active session.


# Load the Seurat library
library(Seurat)

# Read the pre-integrated CCA object
```r
# Absolute Path: /home/emm44q/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds
synovial_cohort <- readRDS("~/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds")

# Verify the object was loaded correctly by checking the cluster identity distribution
table(synovial_cohort$seurat_clusters)
```
---

### 12. Data Access Summary Table
Add this table to the bottom of your README to make the metadata immediately discoverable for incoming researchers.


| Resource Type | Server | Absolute System Path |
| :--- | :--- | :--- |
| **Integrated Object** | `todata3` | `/home/emm44q/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds` |
| **Control Raw Data** | MARS | `/mnt/autofs/data/userdata/project0067/Lab_Data/mouse_ageing/CIA_controls/` |
| **Software Env** | `todata3` | `/home/emm44q/.conda/envs/seurat5_env` |

***




