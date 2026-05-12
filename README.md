# Synovial_scRNAseq_SOP
Standard Operating Procedure:







scRNAseq analysis 
(Young (CIA) vs. old mice synovial tissues)
 
Table of contents
Prerequisite Note: The todata2 server natively supports Seurat v5. The todata3 server is locked to Seurat v4.3. This SOP utilizes Miniconda to run Seurat v5 on todata3.
1. Setup Server and Environment
•	Accessing the todata3 server
•	Setting up Miniconda
•	Loading R and Seurat v5
2. Load the Data
•	Transferring 10X data from MARS
•	Reading the matrices into Seurat
•	Merging the 7 samples
3. Quality Control (QC) Filters
•	Calculating mitochondrial RNA
•	Applying cutoffs (500–6000 genes; <10% MT)
•	Removing the Con2 outlier
4. Normalize and Run PCA
•	Log-normalization and finding variable genes
•	Running PCA
5. CCA Integration
•	Why we use CCA (30 PCs)
•	Running the Integration
•	Joining layers (Crucial for Seurat v5)
•	Generating the UMAP
6. Save the Object
•	Saving the final .rds file
•	Server path for the student to access the prepared data






1. Setup Server and Environment
Objective: Establish a secure connection to the MVLS computational server (todata3) from a local Windows machine and configure a protected environment for Seurat v5.
Accessing the todata3 server via SSH Because sequencing data must remain secure, ensure you are connected to the University of Glasgow VPN (Cisco Secure Client). Open the Windows PowerShell terminal and initiate a Secure Shell (SSH) connection using your GUID:
ssh yourGUID@todata3.mvls.gla.ac.uk
Activating the Miniconda Environment The default R environment on todata3 is locked to Seurat v4.3. We bypass this using a Miniconda virtual environment containing Seurat v5. In PowerShell, activate the environment:
conda activate seurat5_env
(Sanity Check: Your command prompt must change from (base) to (seurat5_env).)
Registering the Jupyter Kernel To use this environment visually in the browser, we must build a computational "bridge" (the IRkernel). Inside the active environment in PowerShell, start R and register the kernel:
# Start the R terminal

# Register the environment to Jupyter
IRkernel::installspec(name = 'seurat5_env', displayname = 'R 4.3 (PRO V5)')

# Quit R
q()
Accessing the Jupyter Web Interface Now, transition from the PowerShell terminal to your standard web browser (e.g., Chrome, Edge).
1.	Navigate to the server address: http://todata3.mvls.gla.ac.uk:8000/
2.	You will be prompted to log in. Enter your standard university GUID and password.
3.	Once logged into the Jupyter dashboard, create a new notebook. Ensure you select R 4.3 (PRO V5) from the kernel dropdown menu located at the top right of the screen.
High-Memory Capacity Reservation Inside your new Jupyter Notebook, run this block immediately in the first cell to prevent "memory exhausted" crashes during downstream integration:
options(future.globals.maxSize = 800 * 1024^3)
Sys.setenv("R_MAX_VSIZE" = "900Gb")
gc()
cat(" ✅ High memory reserved successfully. Workspace is clean.")



2. Load the Data
Objective: Securely transfer raw 10X Genomics sequencing matrices from MARS and construct the initial Seurat object.
Transferring 10X data from MARS Return to your PowerShell terminal. We use rsync instead of a standard copy command because rsync will automatically resume the transfer if the connection drops, preventing corrupted data. 
rsync -avz /mnt/autofs/data/userdata/project0067/Lab_Data/mouse_ageing/ ~/UofG_SingleCell_Enock/mouse_ageing_data/
(Sanity Check: Run du -sh ~/UofG_SingleCell_Enock/mouse_ageing_data/ to confirm the folder is multiple Gigabytes in size, proving the transfer finished completely).
Defining the Root Directory in Jupyter Moving back to your Jupyter Notebook, set the "Main Address" for your data. We map the specific subfolders to ensure the raw matrices are accurately pulled for both the Aged cohort and the Young CIA controls.
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
Reading and Merging the Matrices We use an automated loop to read every folder and convert them into Seurat objects simultaneously.
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
•	Checkpoint 2: Verify the dimensions of the initial merged cohort.
# Output the total number of genes (rows) and total number of cells (columns) in the raw object
dim(combined_mouse)
# Output the initial cell count distribution across the 7 biological replicates
table(combined_mouse$orig.ident)

# Visually audit the quality control metrics using a Violin Plot
VlnPlot(combined_mouse, features = c("nFeature_RNA", "percent.mt"), ncol = 2)
<img width="4200" height="2400" alt="Fig1_Adaptive_QC_Audit" src="https://github.com/user-attachments/assets/4f57ca76-b6be-4c2c-91e1-afb678af7ef4" />




3. Quality Control (QC) Filters
Objective: Remove technical artifacts (empty droplets, multiplets, and apoptotic cells) and exclude sub-optimal samples to establish a statistically sound 3-vs-3 comparative cohort.
Calculating mitochondrial RNA During cell preparation, dying or stressed cells often leak cytoplasmic RNA while retaining mitochondrial RNA. We quantify the percentage of mitochondrial transcripts in every cell to identify and remove these stressed populations.
# Calculate the percentage of reads that map to the mitochondrial genome 
# The "^mt-" pattern is specific to murine (mouse) genomic nomenclature
combined_mouse[["percent.mt"]] <- PercentageFeatureSet(combined_mouse, pattern = "^mt-")
Removing the Con2 outlier and Applying cutoffs Before proceeding to downstream biological analysis, we must apply strict mathematical thresholds to clean the data matrix:
1.	Gene Capture Limits: Cells must possess between 500 and 6,000 unique genes (nFeature_RNA). Counts below 500 indicate ambient RNA debris (empty droplets), while counts above 6,000 strongly suggest technical multiplets (two cells captured in one droplet).
2.	Mitochondrial Limit: Cells must have < 10% mitochondrial RNA to ensure we only analyze healthy, viable cells.
3.	Cohort Balancing (Dropping Con2): The Con2_CD45neg sample displayed critically low transcriptional complexity (a median gene capture of ~415 genes) and elevated stress markers. It is excluded entirely at this stage to prevent it from skewing the high-resolution fibroblast sub-clustering, finalizing a balanced 3-vs-3 experimental design (Young vs. Old).
We apply all of these rules simultaneously using the subset function:
# Apply the specific QC parameters and explicitly filter out the Con2 sample
combined_mouse_clean <- subset(combined_mouse, 
                               subset = nFeature_RNA > 500 & 
                                        nFeature_RNA < 6000 & 
                                        percent.mt < 10 & 
                                        orig.ident != "Con2_CD45neg")
•	Checkpoint 3: Quantify the "Survival Rate" of the cells after applying the QC thresholds to verify that the remaining samples have sufficient cell depth for integration.
# View the cell counts BEFORE filtering
cat("--- Cells BEFORE QC ---\n")
table(combined_mouse$orig.ident)

# View the cell counts AFTER filtering (Con2 should be completely absent)
cat("\n--- Cells AFTER QC ---\n")
table(combined_mouse_clean$orig.ident)

4. Normalize and Run PCA
Objective: Prepare the filtered count matrices for dimensionality reduction by standardizing variance and identifying the key genes driving biological heterogeneity.
Log-normalization and finding variable genes Raw sequence counts cannot be directly compared because some cells are sequenced more deeply than others. We apply Log-Normalization to mathematically standardize the gene expression across all cells. Following this, we isolate the top 2,000 highly variable genes (HVGs) that represent the true biological differences between our cell types, ignoring genes that are flat or uninformative across the cohort.
Scaling and Running PCA Before calculating the Principal Components, the data must be scaled. This ensures that genes with naturally high expression (like housekeeping genes) do not mathematically drown out critically important, but low-expressing, transcription factors. Once scaled, we run Principal Component Analysis (PCA) to compress the thousands of variable genes into the top dimensions representing the majority of the variance.

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
•	Checkpoint 4: Generate an Elbow Plot. This plot visualizes the standard deviation of each Principal Component. It is required to statistically justify how many dimensions will be fed into the final integration algorithm.

# Generate the Elbow Plot to determine the dimensionality of the dataset
ElbowPlot(combined_mouse_clean, ndims = 50)
<img width="2400" height="1800" alt="Fig2_PCA_Elbow_Audit" src="https://github.com/user-attachments/assets/5f8903d4-143c-4724-b55b-56e7d0feaac6" />



5. CCA Integration
Objective: Align the distinct experimental batches (Young CIA vs. Aged) to identify shared biological cell states and remove technical batch variance.
Why we use CCA (30 PCs) Because our samples come from distinct biological groups and batches, we must mathematically integrate them. We utilize Canonical Correlation Analysis (CCA) rather than Harmony for this specific pipeline. CCA is an anchor-based approach that aggressively identifies shared biological structures (Mutual Nearest Neighbors) across the batches. Note on Dimensionality: We strictly use 30 Principal Components (PCs) for CCA. Because CCA is mathematically aggressive, using too many PCs (e.g., 50) can result in "over-correction," where distinct biological populations, like subtle fibroblast sub-states, are artificially crushed together.
Running the Integration Using the Seurat v5 framework, we execute the integration directly on the normalized layers.
# Execute Canonical Correlation Analysis (CCA) Integration using 30 PCs
combined_mouse_clean <- IntegrateLayers(
  object = combined_mouse_clean, 
  method = CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  dims = 1:30, 
  verbose = FALSE
)
Joining layers (Crucial for Seurat v5) CRITICAL STEP: In Seurat v5, the matrices are physically split by sample (e.g., counts.Old1, counts.Con3) to allow for integration. However, downstream Differential Expression (DE) algorithms cannot read across split layers. We must physically collapse these individual matrices back into a single unified matrix.
# Re-unify the raw sequence matrices for downstream marker identification
combined_mouse_clean[["RNA"]] <- JoinLayers(combined_mouse_clean[["RNA"]])
Generating the UMAP and Clusters With the data integrated, we project the cells into a 2D space (UMAP). We cluster the cells using a Resolution of 0.6, which is the community standard for a cohort of this size (~20,000 cells) to accurately separate major fibroblast and myeloid lineages without over-clustering technical noise.

# Run UMAP using the newly integrated CCA reduction
combined_mouse_clean <- RunUMAP(combined_mouse_clean, reduction = "integrated.cca", dims = 1:30)

# Build the network graph and find clusters
combined_mouse_clean <- FindNeighbors(combined_mouse_clean, reduction = "integrated.cca", dims = 1:30)
combined_mouse_clean <- FindClusters(combined_mouse_clean, resolution = 0.6)
•	Checkpoint 5: Visually verify the integration. Generate side-by-side UMAPs to confirm that the Young and Old samples are mixing evenly (left plot) and that distinct biological clusters have formed (right plot).

# Generate the validation plots
p_cca_sample <- DimPlot(combined_mouse_clean, reduction = "umap", group.by = "orig.ident") + 
                labs(title = "CCA: Sample Alignment")
p_cca_cluster <- DimPlot(combined_mouse_clean, reduction = "umap", label = TRUE) + 
                 NoLegend() + labs(title = "CCA: Clusters (Res 0.6)")

# Display plots side-by-side
p_cca_sample | p_cca_cluster
<img width="4200" height="1800" alt="Fig3_CCA_Integration_UMAP" src="https://github.com/user-attachments/assets/c86cdf80-082e-4ead-a35e-d8098b9bfe91" />

 
6. Save and Organize the Object
Objective: Archive the fully integrated, batch-corrected cohort to a persistent server directory, establishing a stable starting point for biological interpretation.
Organizing the Final .rds File Once the Seurat object is serialized into an .rds file, we must transfer it to a dedicated, clearly labelled project directory for incoming researchers. This prevents accidental deletion and ensures strict version control. We also utilize standard naming conventions (Synovial_3vs3_CCA_Integrated.rds) so the file contents are immediately clear.
# Define the exact location of the original processed object
old_file_path <- path.expand("~/UofG_SingleCell_Enock/results/CCA_vs_Harmony/Filtered_3vs3_Synovial_Object.rds")

# Define the new student project directory and create it
target_dir <- path.expand("~/Mahony_Lab/scRNA_seq_data/Eidan")
dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

# Define the target path with the finalized naming convention and execute the copy
new_file_path <- paste0(target_dir, "/Synovial_3vs3_CCA_Integrated.rds")
file.copy(from = old_file_path, to = new_file_path, overwrite = TRUE)

cat(" ✅ Data archived successfully to:", new_file_path)

7. Fast-Track Instructions for Downstream Analysis 
Note for incoming students: If you are operating on a strict timeline, you do not need to run the data ingestion and QC infrastructure detailed in Sections 1–5. You can immediately begin identifying biological markers by loading the pre-processed master object.
Lab Recommendation: Option A is strongly preferred. Single-cell .rds files are massive (often exceeding 5–10 Gigabytes). Loading the file directly into memory without copying it conserves critical server storage and guarantees you are working from the exact same version-controlled dataset as the rest of the lab. Option B is provided strictly for cases where you need to perform modifications to the foundational data. 
Option A: Read-In-Place (Preferred)
You do not need to duplicate this file into your personal folders. Duplicating massive .rds files quickly exhausts server storage limits. Instead, simply load the master object directly from the protected lab directory into your active Jupyter Notebook session:
# Load the required Seurat library into the active R session
library(Seurat)

# Read the pre-integrated CCA object directly from the centralized lab directory
synovial_cohort <- readRDS("~/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds")
Option B: Terminal Copy 
To isolate your work and physically duplicate the file into your personal workspace, open your PowerShell terminal, ensure you are connected to todata3, and use the standard Linux copy (cp) command. You do not need rsync for internal server transfers.

# Create a new directory within your personal home space to hold the data
mkdir -p ~/My_scRNA_Project/data/

# Copy the integrated object from the centralized lab directory to your personal folder
cp ~/Mahony_Lab/scRNA_seq_data/Eidan/Synovial_3vs3_CCA_Integrated.rds ~/My_scRNA_Project/data/
Once the copy is complete, open a new Jupyter Notebook and load your personal copy into memory:
# Load the required Seurat library into the active R session
library(Seurat)

# Read the copied CCA object from your personal project directory
synovial_cohort <- readRDS("~/My_scRNA_Project/data/Synovial_3vs3_CCA_Integrated.rds")

