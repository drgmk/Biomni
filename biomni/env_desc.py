# Data lake dictionary with detailed descriptions
data_lake_dict = {
    # --- interaction / drug / screening datasets (commented out for RNA-only use) ---
    # "affinity_capture-ms.parquet": "Affinity capture mass spectrometry protein–protein interactions.",
    # "affinity_capture-rna.parquet": "Affinity capture assays for protein–RNA interactions.",
    # "BindingDB_All_202409.tsv": "Small molecule–protein binding affinities from BindingDB.",
    # "broad_repurposing_hub_molecule_with_smiles.parquet": "Molecules from the Broad Drug Repurposing Hub with SMILES.",
    # "broad_repurposing_hub_phase_moa_target_info.parquet": "Drug phases, mechanisms of action and targets for the repurposing hub.",
    # "co-fractionation.parquet": "Protein–protein interactions from co-fractionation experiments.",

    # --- RNA / single-cell relevant resources (kept active) ---
    "czi_census_datasets_v4.parquet": "Single-cell datasets from the Chan Zuckerberg Initiative Cell Census.",
    # DepMap & DDInter etc. – useful but not needed for pure RNA workflows
    # "DepMap_CRISPRGeneDependency.csv": "CRISPR gene dependency scores for cancer cell lines.",
    # "DepMap_CRISPRGeneEffect.csv": "CRISPR gene effect estimates for cancer cell lines.",
    # "DepMap_Model.csv": "Metadata describing DepMap cancer models.",
    # "DepMap_OmicsExpressionProteinCodingGenesTPMLogp1.csv": "Gene expression (TPM log1p) for DepMap models.",
    # "ddinter_alimentary_tract_metabolism.csv": "Drug–drug interactions for alimentary tract and metabolism drugs (DDInter 2.0).",
    # "ddinter_antineoplastic.csv": "Drug–drug interactions for antineoplastic and immunomodulating agents (DDInter 2.0).",
    # "ddinter_antiparasitic.csv": "Drug–drug interactions for antiparasitic products (DDInter 2.0).",
    # "ddinter_blood_organs.csv": "Drug–drug interactions for blood and blood-forming organ drugs (DDInter 2.0).",
    # "ddinter_dermatological.csv": "Drug–drug interactions for dermatological drugs (DDInter 2.0).",
    # "ddinter_hormonal.csv": "Drug–drug interactions for hormonal preparations (DDInter 2.0).",
    # "ddinter_respiratory.csv": "Drug–drug interactions for respiratory system drugs (DDInter 2.0).",
    # "ddinter_various.csv": "Drug–drug interactions for various other drugs (DDInter 2.0).",
    # "DisGeNET.parquet": "Gene–disease associations from multiple sources (DisGeNET).",
    # "dosage_growth_defect.parquet": "Growth phenotypes linked to gene dosage perturbations.",
    # "enamine_cloud_library_smiles.pkl": "Enamine REAL library compounds with SMILES.",
    # "evebio_assay_table.csv": "Assay metadata from EveBio screening.",
    # "evebio_bundle_table.csv": "EveBio target bundle metadata.",
    # "evebio_compound_table.csv": "EveBio compound metadata and identifiers.",
    # "evebio_control_table.csv": "Control datapoints for EveBio screening plates.",
    # "evebio_detailed_result_table.csv": "Detailed EveBio screening results with curve fit information.",
    # "evebio_observed_points_table.csv": "Raw observed points from EveBio screening and profiling.",
    # "evebio_summary_result_table.csv": "Summary results for each EveBio assay–compound combination.",
    # "evebio_target_table.csv": "Target metadata for EveBio screens.",

    # --- genetics / GWAS / misc (commented for RNA-only) ---
    # "genebass_missense_LC_filtered.pkl": "Filtered missense variants from GeneBass.",
    # "genebass_pLoF_filtered.pkl": "Predicted loss-of-function variants from GeneBass.",
    # "genebass_synonymous_filtered.pkl": "Filtered synonymous variants from GeneBass.",

    # --- core gene / ontology / expression (kept) ---
    "gene_info.parquet": "Basic gene information and identifiers.",
    # genetic interaction network etc.
    "genetic_interaction.parquet": "Genetic interaction scores across pairs of genes.",
    "go-plus.json": "Gene Ontology (GO) in extended JSON format.",
    "gtex_tissue_gene_tpm.parquet": "Bulk RNA-seq expression (TPM) across human tissues from GTEx.",
    # "gwas_catalog.pkl": "GWAS Catalog variant–trait associations.",
    "hp.obo": "Human Phenotype Ontology (HPO) in OBO format.",
    # "kg.csv": "Knowledge graph edges summarising biological relationships.",

    "marker_celltype.parquet": "Marker genes for cell types to aid annotation.",
    "McPAS-TCR.parquet": "TCR sequences and associated pathogens/antigens.",
    # "miRDB_v6.0_results.parquet": "microRNA target predictions from miRDB.",
    # "miRTarBase_microRNA_target_interaction.parquet": "Experimentally validated miRNA–target interactions (MiRTarBase).",
    # "miRTarBase_microRNA_target_interaction_pubmed_abtract.txt": "PubMed abstracts for MiRTarBase interactions.",
    # "miRTarBase_MicroRNA_Target_Sites.parquet": "miRNA target sites from MiRTarBase.",

    # --- MouseMine gene sets (kept) ---
    "mousemine_m1_positional_geneset.parquet": "Mouse positional gene sets from MouseMine (MSigDB M1 style).",
    "mousemine_m2_curated_geneset.parquet": "Curated mouse gene sets from MouseMine (M2).",
    "mousemine_m3_regulatory_target_geneset.parquet": "Regulatory target gene sets (mouse) from MouseMine (M3).",
    "mousemine_m5_ontology_geneset.parquet": "Ontology-based mouse gene sets from MouseMine (M5).",
    "mousemine_m8_celltype_signature_geneset.parquet": "Mouse cell type signature gene sets from MouseMine (M8).",
    "mousemine_mh_hallmark_geneset.parquet": "Mouse hallmark gene sets from MouseMine (MH).",

    # --- MSigDB human gene sets (kept) ---
    "msigdb_human_c1_positional_geneset.parquet": "MSigDB human C1 positional gene sets.",
    "msigdb_human_c2_curated_geneset.parquet": "MSigDB human C2 curated gene sets.",
    "msigdb_human_c3_regulatory_target_geneset.parquet": "MSigDB human C3 regulatory target gene sets.",
    "msigdb_human_c3_subset_transcription_factor_targets_from_GTRD.parquet": "MSigDB/GTRD transcription factor target gene sets.",
    "msigdb_human_c4_computational_geneset.parquet": "MSigDB human C4 computational gene sets.",
    "msigdb_human_c5_ontology_geneset.parquet": "MSigDB human C5 ontology gene sets (GO etc.).",
    "msigdb_human_c6_oncogenic_signature_geneset.parquet": "MSigDB human C6 oncogenic signatures.",
    "msigdb_human_c7_immunologic_signature_geneset.parquet": "MSigDB human C7 immunologic signatures.",
    "msigdb_human_c8_celltype_signature_geneset.parquet": "MSigDB human C8 cell type signature gene sets.",
    "msigdb_human_h_hallmark_geneset.parquet": "MSigDB human H hallmark gene sets.",

    # --- other large resources not needed for basic RNA workflows ---
    # "omim.parquet": "Gene–disease associations from OMIM.",
    # "proteinatlas.tsv": "Protein expression profiles from Human Protein Atlas.",
    # "proximity_label-ms.parquet": "Protein interactions from proximity labeling MS experiments.",
    # "reconstituted_complex.parquet": "Reconstituted protein complex interaction data.",
    # "sgRNA_KO_SP_mouse.txt": "sgRNA knockout screen data (mouse).",
    # "sgRNA_KO_SP_human.txt": "sgRNA knockout screen data (human).",
    # "synthetic_growth_defect.parquet": "Synthetic growth defect genetic interactions.",
    # "synthetic_lethality.parquet": "Synthetic lethal genetic interactions.",
    # "synthetic_rescue.parquet": "Synthetic rescue genetic interactions.",
    # "two-hybrid.parquet": "Yeast two-hybrid protein–protein interactions.",
    # "variant_table.parquet": "Annotated variant table across genes and samples.",
    # "Virus-Host_PPI_P-HIPSTER_2020.parquet": "Virus–host protein–protein interactions (P-HIPSTER).",
    # "txgnn_name_mapping.pkl": "Name mapping metadata for TXGNN.",
    # "txgnn_prediction.pkl": "Prediction outputs from TXGNN.",
}


# Library content dictionary with descriptions of Python, R, and CLI tools
library_content_dict = {
    # === PYTHON PACKAGES (bio / RNA / general DS) ===
    "biopython": "[Python] Tools for biological sequences and file formats.",
    "biom-format": "[Python] Support for BIOM tables used in many omics workflows.",
    "scanpy": "[Python] Single-cell RNA-seq analysis toolkit using AnnData.",
    "scikit-bio": "[Python] General bioinformatics data structures and algorithms.",
    "anndata": "[Python] Annotated matrix object for single-cell and omics data.",
    "mudata": "[Python] Multimodal extension of AnnData for multiple omics layers.",

    # These are more niche / structural / chemoinformatics → commented out
    # "pyliftover": "[Python] Genomic coordinate liftover between assemblies.",
    # "biopandas": "[Python] Molecular structure handling in pandas.",
    # "biotite": "[Python] Structural bioinformatics and macromolecular analysis.",
    # "lazyslide": "[Python] Whole-slide image handling for digital pathology.",

    "gget": "[Python] Convenience access to common genomic databases and annotations.",
    "lifelines": "[Python] Survival analysis (Kaplan–Meier, Cox models, etc.).",
    "gseapy": "[Python] Gene set enrichment analysis and visualisation.",
    "scrublet": "[Python] Doublet detection for single-cell RNA-seq.",
    "cellxgene-census": "[Python] Interface to the CellxGene Census single-cell resource.",
    "hyperopt": "[Python] Hyperparameter optimisation library (useful for ML workflows).",
    "scvelo": "[Python] RNA velocity analysis for single-cell data.",
    "pysam": "[Python] SAM/BAM/VCF IO and alignment utilities.",
    "pyfaidx": "[Python] Efficient random access to FASTA sequences.",
    "pyranges": "[Python] Genomic interval operations with a pandas-like interface.",
    "pybedtools": "[Python] Python wrapper for BEDTools genomic interval operations.",

    # gmk added
    "celltypist": "[Python] Automated cell type annotation for single-cell RNA-seq data.",
    "decoupler": "[Python] Gene regulatory network analysis, activity inference, cell annotation.",
    "pydeseq2": "[Python] Differential expression analysis using Python port of DESeq2",
    "cellphonedb": "[Python] Database of ligand–receptor interactions for cell–cell communication analysis.",

    # Chemoinformatics / docking / protein ML – not needed for RNA-only
    # "rdkit": "[Python] Chemoinformatics and molecule handling.",
    # "deeppurpose": "[Python] Deep learning for drug–target interaction and ADMET.",
    # "pyscreener": "[Python] Virtual screening utilities.",
    # "openbabel": "[Python] Chemical toolbox for SMILES/structure conversions.",
    # "descriptastorus": "[Python] Molecular descriptors and featurisation.",
    # "openmm": "[Python] Molecular dynamics simulation toolkit.",
    # "pytdc": "[Python] Therapeutics Data Commons access library.",

    # Core data science stack
    "pandas": "[Python] DataFrames and general data analysis.",
    "numpy": "[Python] N-dimensional arrays and numerical computing.",
    "scipy": "[Python] Scientific computing utilities and statistics.",
    "scikit-learn": "[Python] Machine learning algorithms and model utilities.",
    "matplotlib": "[Python] Plotting library for static visualisations.",
    "seaborn": "[Python] High-level statistical plotting on top of matplotlib.",
    "statsmodels": "[Python] Statistical models and tests.",
    # "pymc3": "[Python] Probabilistic programming and Bayesian inference.",

    "umap-learn": "[Python] UMAP dimensionality reduction.",
    "faiss-cpu": "[Python] Fast nearest-neighbour search and clustering.",
    "harmony-pytorch": "[Python] PyTorch implementation of Harmony batch correction.",

    "tiledb": "[Python] Array storage engine; useful for large omics matrices.",
    "tiledbsoma": "[Python] SOMA/TileDB tools for large single-cell datasets.",
    "h5py": "[Python] HDF5 file reading and writing.",
    "tqdm": "[Python] Progress bars for loops and pipelines.",
    "joblib": "[Python] Simple parallelism and caching utilities.",

    # Extra utilities mostly not needed for this RNA-focused setup
    # "opencv-python": "[Python] Computer vision and image processing.",
    # "PyPDF2": "[Python] PDF manipulation and parsing.",
    # "googlesearch-python": "[Python] Simple Google search wrapper.",
    # "scikit-image": "[Python] Image processing algorithms.",
    # "pymed": "[Python] PubMed / MEDLINE search interface.",
    # "arxiv": "[Python] ArXiv metadata and paper search.",
    # "scholarly": "[Python] Google Scholar scraping utilities.",
    # "cryosparc-tools": "[Python] CryoSPARC API client.",
    # "mageck": "[Python] MAGeCK CRISPR screen analysis tools.",
    # "igraph": "[Python] Graph analysis; can be installed separately if needed.",
    "pyscenic": "[Python] Gene regulatory network inference for single-cell RNA-seq.",
    # "cooler": "[Python] Hi-C contact matrix handling.",
    # "trackpy": "[Python] Particle tracking in videos.",
    # "nnunet": "[Python] Deep learning segmentation (medical images).",
    # "cellpose": "[Python] Generalist cell segmentation models.",
    # "viennarna": "[Python] RNA secondary structure prediction and thermodynamics.",

    # Metabolomics / modelling / flow / simulation – commented out
    # "PyMassSpec": "[Python] Mass spectrometry utilities.",
    # "python-libsbml": "[Python] SBML model handling.",
    # "cobra": "[Python] Constraint-based metabolic modelling.",
    # "reportlab": "[Python] PDF generation library.",
    # "flowkit": "[Python] Flow cytometry analysis tools.",
    # "hmmlearn": "[Python] Hidden Markov models.",
    # "msprime": "[Python] Coalescent simulation.",
    # "tskit": "[Python] Tree sequence storage and analysis.",
    # "cyvcf2": "[Python] Fast VCF parsing.",
    # "pykalman": "[Python] Kalman filtering tools.",
    # "fanc": "[Python] Hi-C/3D genome analysis.",
    # "loompy": "[Python] Loom file utilities for single-cell data.",
    # "pyBigWig": "[Python] BigWig track reading.",
    # "pymzml": "[Python] mzML mass-spec file parsing.",
    # "optlang": "[Python] Optimisation modelling language.",
    # "FlowIO": "[Python] Flow cytometry file IO.",
    # "FlowUtils": "[Python] Flow cytometry utilities.",

    "arboreto": "[Python] Gene regulatory network inference (GRN) from expression data.",
    # "pdbfixer": "[Python] Fixing and preparing PDB protein structures.",

    # === R PACKAGES (kept) ===
    "ggplot2": "[R] Grammar of Graphics plotting for R.",
    "dplyr": "[R] Grammar of data manipulation.",
    "tidyr": "[R] Tidy data reshaping.",
    "readr": "[R] Fast rectangular data import.",
    "stringr": "[R] String manipulation helpers.",
    "Matrix": "[R] Sparse and dense matrix support.",
    "DESeq2": "[R] Differential expression for bulk RNA-seq.",
    "clusterProfiler": "[R] Functional enrichment and gene set analysis.",
    "edgeR": "[R] Differential expression for count data.",
    "limma": "[R] Linear models and voom for expression data.",
    "harmony": "[R] Batch correction and integration (Harmony).",
    "WGCNA": "[R] Weighted gene co-expression network analysis.",

    # gmk added
    "nichenetr": "[R] Ligand–target regulatory network analysis for cell–cell communication.",
    "multinichenetr": "[R] Extension of NicheNet for multi-sample differential analysis between conditions.",

    # === CLI TOOLS (subset kept) ===
    "samtools": "[CLI] BAM/CRAM/VCF utilities (alignment and variant files).",
    "bowtie2": "[CLI] Short-read aligner to reference genomes.",
    "bwa": "[CLI] Burrows–Wheeler aligner for NGS reads.",
    "bedtools": "[CLI] Genomic interval arithmetic and set operations.",
    # "macs2": "[CLI] Peak calling for ChIP/ATAC-seq.",
    "fastqc": "[CLI] QC for high-throughput sequencing reads.",
    "trimmomatic": "[CLI] Read trimming and adapter removal.",
    "mafft": "[CLI] Multiple sequence alignment.",

    # Other CLI tools not needed for this RNA-only setup
    # "Homer": "[CLI] Motif discovery and next-gen sequencing analysis.",
    # "FastTree": "[CLI] Approximate maximum-likelihood phylogenetic trees.",
    # "muscle": "[CLI] Multiple sequence alignment.",
    # "plink": "[CLI] GWAS and population genetics toolkit.",
    # "plink2": "[CLI] Updated PLINK GWAS toolkit.",
    # "gcta64": "[CLI] Genome-wide Complex Trait Analysis.",
    # "iqtree2": "[CLI] Phylogenetic inference with ultrafast bootstrapping.",
    # "ADFR": "[CLI] AutoDock for Receptors docking suite.",
    # "diamond": "[CLI] Fast protein sequence aligner.",
    # "fcsparser": "[CLI] Flow cytometry (FCS) file parser.",
    # "plannotate": "[CLI] Plasmid annotation.",
    # "vina": "[CLI] AutoDock Vina docking.",
    # "autosite": "[CLI] Binding site detection for protein structures.",
    # "PyLabRobot": "[Python] Control of liquid-handling robots and lab automation.",
}
