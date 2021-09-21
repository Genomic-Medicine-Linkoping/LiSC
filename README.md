# LiSC
LiSC is a [Seurat (v4.0.4)](https://satijalab.org/seurat/ 'Seurat main page') wrapper for scRNA-seq data analysis. LiSC supports integrated or non-integrated data analysis with or without SCTransform (please refer to the Seurat documentation). LiSC implements [scCATCH (v2.1)](https://github.com/ZJUFanLab/scCATCH 'scCATCH') to perform cell type prediction.

<br>

## Input files
LiSC takes in input a [CellRanger (v6.1.1)](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation 'Cell Ranger Installation') folder or a list of paths and you can pass them through a .csv file with some other information such as the Sample name and Condition. You don't need to assign any column name to this file. Please find templates to run both integrated (info_integrate.csv) and non-integrated analysis (info_single.csv) [here](https://github.com/pyrevo/LiSC/tree/main/info 'csv templates'). Basically, supposing that you have to run a Seurat analysis integrating data from different samples and conditions, you have to compile a file like this:

```bash=
id1,cond1,path/to/cellranger/id1
id2,cond1,path/to/cellranger/id2
id3,cond2,path/to/cellranger/id3
...
```

<br>

## Output Files
LiSC gives 6 .pdf files and 1 R object in output. The .pdf files show the results of the different steps of the analysis while the R object is meant to load the entire analysis in the R console (or in [RStudio](https://www.rstudio.com/)) Here is a brief description of the output files content:

| Output         | Description |
| -------------- | ------------------------------------------------------------------ |
| 1_preQC.pdf    | quality control stats before filtering                             |
| 2_postQC.pdf   | quality control after filtering                                    |
| 3_PCA.pdf      | results of the linear dimensional reduction                        |
| 4_UMAP.pdf     | results of the non-linear dimensional reduction                    |
| 5_Markers.pdf  | markers that define different clusters via differential expression |
| 6_annoUMAP.pdf | results of the cell type annotation from scCATCH                   |
| id.rds         | R object that can be loaded in Seurat for further analyses         |

Actually, LiSC performs the search for markers of different cell populations with the function *FindAllMarkers()* in Seurat that performs an All vs All comparison. However, it does not take into account the different conditions that you might have. For this reason, you can load the R object (which inherits the project id) and perform further analyses. After having imported the Seurat library, you can load this object like this:

```r
library(Seurat)
mySeuratObj <- readRDS("path/to/results/id.rds")
```

<br>

## Run LiSC
You require just one line to run the single basic analysis or you can provide several samples but then just the first line will be used for the analysis.

This is how the help screen looks like. Required arguments are the algorithms (you can choose between single or integrated analysis), id of the project, the .csv file with the information described above and the path in which you want to store the results. A folder with the id of the project will be created at this path storing all the outputs generated within the analysis. In addition, you can tune some parameters to drive the filtering and clustering steps and chose to apply SCTransform or not (disabled by default):

```bash
LiSC is a Seurat (v4.0.4) wrapper and uses scCATCH (v2.1) for cell-type annotation.

Usage:
  LiSC.R (single | integrate) <id> <csv> <out> [--minGPC=<n>, --maxGPC=<n>, -npc=<n>, --res=<n>, --species=<c>, --tissue=<c>, --SCT]
  LiSC.R (-h | --help)
  LiSC.R --version

Options:
  -h --help         Show this screen.
  --version         Show version.
  --SCT             Use SCTransform algorithm for normalization.
  --minGPC=<n>      Min number of genes per cell [default: 200].
  --maxGPC=<n>      Max number of genes per cell [default: 2500].
  --percMT=<n>      Max percentage of MT [default: 5].
  -n --npc=<n>      Max number of Principal Components [default: 30].
  -r --res=<n>      Cluster resolution [default: 0.5].
  -s --species=<c>  Species for scCATCH [default: Human].
  -t --tissue=<c>   Tissue for scCATCH [default: Blood].
```

Here are some examples of commands to run LiSC with the different algorithms provided:

```bash=
# Run single analysis without SCT
Rscript LiSC.R single testSingle_vst info_single.csv ./
# Run single analysis with SCT
Rscript LiSC.R single testSingle_sct info_single.csv ./ --SCT
# Run integrated analysis without SCT
Rscript LiSC.R integrate testIntegrate_vst info_integrate.csv ./
# Run integrated analysis with SCT
Rscript LiSC.R integrate testIntegrate_sct info_integrate.csv ./ --SCT
```

<br>

## scCATCH
[scCATCH](https://github.com/ZJUFanLab/scCATCH 'scCATCH') allows automatic annotation on cell types of clusters found in Seurat but you need to provide a species name for the database (Human and Mouse are currently supported) and the tissue. You can pass this information through the parameters `--species` (Human by default) and `--tissue` (Blood by default). You can find the complete list of the supported tissues below:


```bash=
#1.1 For Human tissue, tissue types are listed as follows:

Adipose tissue-related: Abdominal adipose tissue; Adipose tissue; Brown adipose tissue; Fat pad; Subcutaneous adipose tissue; Visceral adipose tissue; White adipose tissue.

Bladder-related: Bladder; Urine.

Blood-related: Blood; Peripheral blood; Plasma; Serum; Umbilical cord blood; Venous blood.

Bone-related: Anterior cruciate ligament; Bone; Bone marrow; Cartilage; Intervertebral disc; Meniscus; Nucleus pulposus; Osteoarthritic cartilage; Periosteum; Skeletal muscle; Spinal cord; Synovial fluid; Synovium.

Brain-related: Brain; Dorsolateral prefrontal cortex; Embryonic brain; Embryonic prefrontal cortex; Fetal brain; Hippocampus; Inferior colliculus; Midbrain; Sympathetic ganglion.

Breast-related: Breast; Mammary epithelium.

Embryo-related: Embryo; Embryoid body; Embryonic brain; Embryonic prefrontal cortex; Embryonic stem cell; Germ; Primitive streak.

Esophagus-related: Esophagus.

Eye-related: Cornea; Corneal endothelium; Corneal epithelium; Eye; Lacrimal gland; Limbal epithelium; Optic nerve; Retina; Retinal pigment epithelium; Sclerocorneal tissue.

Fetus-related: Amniotic fluid; Amniotic membrane; Fetal brain; Fetal gonad; Fetal kidney; Fetal liver; Fetal pancreas; Placenta; Umbilical cord; Umbilical cord blood; Umbilical vein.

Gonad-related: Corpus luteum; Fetal gonad; Foreskin; Gonad; Ovarian cortex; Ovarian follicle; Ovary; Seminal plasma; Testis.

Hair-related: Chorionic villus; Hair follicle; Scalp.

Heart-related: Heart; Myocardium.

Intestine-related: Colon; Colorectum; Gastrointestinal tract; Gut; Intestinal crypt; Intestine; Jejunum; Large intestine; Small intestinal crypt; Small intestine.

Kidney-related: Adrenal gland; Fetal kidney; Kidney; Renal glomerulus.

Liver-related: Fetal liver; Liver.

Lung-related:Airway epithelium; Alveolus; Bronchoalveolar system; Lung.

Lymph-related: Lymph; Lymph node; Lymphoid tissue.

Muscle-related: Muscle; Skeletal muscle.

Nose-related: Nasal concha; Nasal epithelium; Sinonasal mucosa.

Oral cavity-related: Laryngeal squamous epithelium; Oral mucosa; Salivary gland; Sputum; Submandibular gland; Thyroid; Tonsil; Vocal fold.

Ovary-related: Corpus luteum; Ovarian cortex; Ovarian follicle; Ovary; Oviduct.

Pancreas-related: Fetal pancreas; Pancreas; Pancreatic acinar tissue; Pancreatic islet.

Prostate-related: Prostate.

Skin-related: Dermis; Skin.

Spleen-related: Spleen; Splenic red pulp.

Stomach-related: Gastric corpus; Gastric epithelium; Gastric gland; Gastrointestinal tract; Pyloric gland; Stomach.

Testis-related: Testis.

Tooth-related: Deciduous tooth; Dental pulp; Gingiva; Molar; Periodontal ligament; Premolar; Tooth.

Uterus-related: Endometrium; Endometrium stroma; Myometrium; Uterus; Vagina.

Vessel-related: Adventitia; Antecubital vein; Artery; Blood vessel; Umbilical vein.

Others: Ascites; Epithelium; Ligament; Pluripotent stem cell; Thymus; Whartons jelly.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#1.2 For Human tissue about cancer, cancer types and the corresponding tissue types are listed as follows:

Acute Myelogenous Leukemia: Blood.

Acute Myeloid Leukemia: Bone marrow.

Adenoid Cystic Carcinoma: Salivary gland.

Alveolar Cell Carcinoma: Serum.

Astrocytoma: Brain.

B-Cell Lymphoma: Lymphoid tissue.

Bladder Cancer: Bladder.

Brain Cancer: Blood vessel; Brain.

Breast Cancer: Blood: Breast; Mammary gland.

Cholangiocarcinoma: Liver; Serum.

Chronic Lymphocytic Leukemia: Blood.

Chronic Myeloid Leukemia: Blood.

CNS Primitive Neuroectodermal Tumor: Brain.

Colon Cancer: Blood; Colon; Serum.

Colorectal Cancer: Blood; Colon; Colorectum; Gastrointestinal tract; Intestine; Liver; Lung; Venous blood.

Cutaneous Squamous Cell Carcinoma: Skin.

Endometrial Cancer: Endometrium.

Ependymoma: Brain.

Esophageal Adenocarcinoma: Esophagus.

Fibroid: Myometrium.

Follicular Lymphoma: Lymph node.

Gallbladder Cancer: Gall bladder; Gastrointestinal tract.

Gastric Cancer: Blood; Peripheral blood; Serum; Stomach.

Glioblastoma: Blood; Brain.

Glioma: Blood vessel; Brain.

Gonadoblastoma: Embryo.

Head and Neck Cancer: Blood; Brain; Oral cavity.

Hepatoblastoma: Liver.

Hepatocellular Cancer: Blood; Bone marrow; Embryo; Liver.

High-grade glioma: Brain.

Infantile Hemangiomas: Placenta.

Intestinal Cancer: Gastrointestinal tract.

Intracranial Aneurysm: Brain.

Kaposi s Sarcoma: Lymph node.

Larynx Cancer: Larynx.

Leukemia: Bone marrow; Peripheral blood.

Lipoma: Adipose tissue.

Liver Cancer: Blood; Liver.

Lung Adenocarcinoma: Lung.

Lung Cancer: Blood; Lung.

Lung Squamous Cell Carcinoma: Lung.

Lymphoma: Blood; Brain; Kidney; Liver; Lymph; Lymph node.

Malignant Insulinoma: Pancreas.

Malignant Mesothelioma: Lung; Pleura.

Malignant Peripheral Nerve Sheath Tumor: Brain.

Medulloblastoma: Brain.

Melanoma: Blood; Peripheral blood; Skin.

Mucoepidermoid Carcinoma: Salivary gland.

Multiple Myeloma: Bone marrow; Peripheral blood.

Myeloma: Bone marrow.

Natural Killer Cell Lymphoma: Lymph node.

Nephroblastoma: Kidney.

Non-Small Cell Lung Cancer: Blood; Lung; Peripheral blood.

Oesophageal Cancer: Blood.

Oligodendroglioma: Brain.

Oral Cancer: Oral cavity.

Oral Squamous Cell Carcinoma: Oral cavity; Peripheral blood.

Osteosarcoma: Bone.

Ovarian Cancer: Ascites; Ovarian cortex; Ovary; Peripheral blood.

Pancreatic Cancer: Blood vessel; Pancreas.

Pancreatic Ductal Adenocarcinomas: Pancreas.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Blood; Peripheral blood; Prostate.

Renal Cell Carcinoma: Kidney; Serum.

Renal Clear Cell Carcinoma: Lymph node.

Retinoblastoma: Eye.

Salivary Gland Tumor: Parotid gland; Salivary gland.

Sarcoma: Muscle.

Small Cell Lung Cancer: Lung.

Testicular Germ Cell Tumor: Peripheral blood; Testis.

Thyroid Cancer: Thyroid.

Tongue Cancer: Tongue.

Uterine Leiomyoma: Uterus.

Vascular Tumour: Lymph node.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#2.1 For Mouse tissue, tissue types are listed as follows:

Adipose tissue-related: Adipose tissue; White adipose tissue.

Bladder-related: Bladder.

Blood-related: Blood; Peripheral blood; Serum; Umbilical cord blood.

Bone-related: Bone; Bone marrow; Meniscus; Skeletal muscle; Spinal cord.

Brain-related: Brain; Cerebellum; Fetal brain; Hippocampus; Neural tube.

Breast-related: Mammary epithelium; Mammary gland.

Calvaria-related: Neonatal calvaria.

Ear-related: Cochlea; Inner Ear.

Embryo-related: Embryo; Embryoid body; Embryonic heart; Embryonic stem cell.

Esophagus-related: Esophagus.

Eye-related: Corneal epithelium; Eye; Ganglion cell layer of retina; Inner nuclear layer of retina; Lacrimal gland; Retina.

Fetus-related: Fetal brain; Fetal intestine; Fetal liver; Fetal lung; Fetal stomach; Placenta; Umbilical cord; Umbilical cord blood.

Gonad-related: Gonad; Ovary; Testis; Yolk sac.

Hair-related: Hair follicle.

Heart-related: Embryonic heart; Heart; Heart muscle; Neonatal heart.

Intestine-related: Colon; Colon epithelium; Fetal intestine; Gastrointestinal tract; Ileum; Intestinal crypt; Intestine; Mesenteric lymph node; Small intestine.

Kidney-related: Kidney; Mesonephros.

Liver-related: Fetal liver; Liver.

Lung-related: Bronchiole; Fetal lung; Lung; Trachea.

Lymph-related: Lymph node; Lymphoid tissue; Mesenteric lymph node; Peyer patch.

Muscle-related: Heart muscle; Muscle; Neonatal muscle; Skeletal muscle.

Neonate-related: Neonatal calvaria; Neonatal heart; Neonatal muscle; Neonatal pancreas; Neonatal rib; Neonatal skin.

Oral cavity-related: Submandibular gland; Taste bud.

Ovary-related: Ovary; Yolk sac.

Pancreas-related: Neonatal pancreas; Pancreas; Pancreatic islet.

Prostate-related: Prostate.

Skin-related: Dermis; Epidermis; Neonatal skin; Skin.

Spleen-related: Spleen.

Stomach-related: Fetal stomach; Gastrointestinal tract; Stomach.

Testis-related: Testis.

Uterus-related: Uterus.

Vessel-related: Aorta; Artery; Blood vessel; Carotid artery.

Others: Basilar membrane; Epithelium; Peritoneal cavity; Thymus.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#2.2 For Mouse tissue about cancer, cancer types and the corresponding tissue types are listed as follows:

Breast Cancer: Lymph node; Breast; Lung.

Chronic Myeloid Leukemia: Blood.

Colon Cancer: Colon.

Colorectal Cancer: Lymph node; Colon; Colorectum.

Cutaneous Squamous Cell Carcinoma: Skin.

Hepatocellular Cancer: Blood; Liver.

Liver Cancer: Liver.

Lung Cancer: Lung.

Melanoma: Lung.

Pancreatic Cancer: Blood.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Prostate.

Renal Cell Carcinoma: Kidney.

Supratentorial Primitive Neuroectodermal Tumor: Brain. https://github.com/ZJUFanLab/scCATCH1.1 For Human tissue, tissue types are listed as follows:

Adipose tissue-related: Abdominal adipose tissue; Adipose tissue; Brown adipose tissue; Fat pad; Subcutaneous adipose tissue; Visceral adipose tissue; White adipose tissue.

Bladder-related: Bladder; Urine.

Blood-related: Blood; Peripheral blood; Plasma; Serum; Umbilical cord blood; Venous blood.

Bone-related: Anterior cruciate ligament; Bone; Bone marrow; Cartilage; Intervertebral disc; Meniscus; Nucleus pulposus; Osteoarthritic cartilage; Periosteum; Skeletal muscle; Spinal cord; Synovial fluid; Synovium.

Brain-related: Brain; Dorsolateral prefrontal cortex; Embryonic brain; Embryonic prefrontal cortex; Fetal brain; Hippocampus; Inferior colliculus; Midbrain; Sympathetic ganglion.

Breast-related: Breast; Mammary epithelium.

Embryo-related: Embryo; Embryoid body; Embryonic brain; Embryonic prefrontal cortex; Embryonic stem cell; Germ; Primitive streak.

Esophagus-related: Esophagus.

Eye-related: Cornea; Corneal endothelium; Corneal epithelium; Eye; Lacrimal gland; Limbal epithelium; Optic nerve; Retina; Retinal pigment epithelium; Sclerocorneal tissue.

Fetus-related: Amniotic fluid; Amniotic membrane; Fetal brain; Fetal gonad; Fetal kidney; Fetal liver; Fetal pancreas; Placenta; Umbilical cord; Umbilical cord blood; Umbilical vein.

Gonad-related: Corpus luteum; Fetal gonad; Foreskin; Gonad; Ovarian cortex; Ovarian follicle; Ovary; Seminal plasma; Testis.

Hair-related: Chorionic villus; Hair follicle; Scalp.

Heart-related: Heart; Myocardium.

Intestine-related: Colon; Colorectum; Gastrointestinal tract; Gut; Intestinal crypt; Intestine; Jejunum; Large intestine; Small intestinal crypt; Small intestine.

Kidney-related: Adrenal gland; Fetal kidney; Kidney; Renal glomerulus.

Liver-related: Fetal liver; Liver.

Lung-related:Airway epithelium; Alveolus; Bronchoalveolar system; Lung.

Lymph-related: Lymph; Lymph node; Lymphoid tissue.

Muscle-related: Muscle; Skeletal muscle.

Nose-related: Nasal concha; Nasal epithelium; Sinonasal mucosa.

Oral cavity-related: Laryngeal squamous epithelium; Oral mucosa; Salivary gland; Sputum; Submandibular gland; Thyroid; Tonsil; Vocal fold.

Ovary-related: Corpus luteum; Ovarian cortex; Ovarian follicle; Ovary; Oviduct.

Pancreas-related: Fetal pancreas; Pancreas; Pancreatic acinar tissue; Pancreatic islet.

Prostate-related: Prostate.

Skin-related: Dermis; Skin.

Spleen-related: Spleen; Splenic red pulp.

Stomach-related: Gastric corpus; Gastric epithelium; Gastric gland; Gastrointestinal tract; Pyloric gland; Stomach.

Testis-related: Testis.

Tooth-related: Deciduous tooth; Dental pulp; Gingiva; Molar; Periodontal ligament; Premolar; Tooth.

Uterus-related: Endometrium; Endometrium stroma; Myometrium; Uterus; Vagina.

Vessel-related: Adventitia; Antecubital vein; Artery; Blood vessel; Umbilical vein.

Others: Ascites; Epithelium; Ligament; Pluripotent stem cell; Thymus; Whartons jelly.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#1.2 For Human tissue about cancer, cancer types and the corresponding tissue types are listed as follows:

Acute Myelogenous Leukemia: Blood.

Acute Myeloid Leukemia: Bone marrow.

Adenoid Cystic Carcinoma: Salivary gland.

Alveolar Cell Carcinoma: Serum.

Astrocytoma: Brain.

B-Cell Lymphoma: Lymphoid tissue.

Bladder Cancer: Bladder.

Brain Cancer: Blood vessel; Brain.

Breast Cancer: Blood: Breast; Mammary gland.

Cholangiocarcinoma: Liver; Serum.

Chronic Lymphocytic Leukemia: Blood.

Chronic Myeloid Leukemia: Blood.

CNS Primitive Neuroectodermal Tumor: Brain.

Colon Cancer: Blood; Colon; Serum.

Colorectal Cancer: Blood; Colon; Colorectum; Gastrointestinal tract; Intestine; Liver; Lung; Venous blood.

Cutaneous Squamous Cell Carcinoma: Skin.

Endometrial Cancer: Endometrium.

Ependymoma: Brain.

Esophageal Adenocarcinoma: Esophagus.

Fibroid: Myometrium.

Follicular Lymphoma: Lymph node.

Gallbladder Cancer: Gall bladder; Gastrointestinal tract.

Gastric Cancer: Blood; Peripheral blood; Serum; Stomach.

Glioblastoma: Blood; Brain.

Glioma: Blood vessel; Brain.

Gonadoblastoma: Embryo.

Head and Neck Cancer: Blood; Brain; Oral cavity.

Hepatoblastoma: Liver.

Hepatocellular Cancer: Blood; Bone marrow; Embryo; Liver.

High-grade glioma: Brain.

Infantile Hemangiomas: Placenta.

Intestinal Cancer: Gastrointestinal tract.

Intracranial Aneurysm: Brain.

Kaposi s Sarcoma: Lymph node.

Larynx Cancer: Larynx.

Leukemia: Bone marrow; Peripheral blood.

Lipoma: Adipose tissue.

Liver Cancer: Blood; Liver.

Lung Adenocarcinoma: Lung.

Lung Cancer: Blood; Lung.

Lung Squamous Cell Carcinoma: Lung.

Lymphoma: Blood; Brain; Kidney; Liver; Lymph; Lymph node.

Malignant Insulinoma: Pancreas.

Malignant Mesothelioma: Lung; Pleura.

Malignant Peripheral Nerve Sheath Tumor: Brain.

Medulloblastoma: Brain.

Melanoma: Blood; Peripheral blood; Skin.

Mucoepidermoid Carcinoma: Salivary gland.

Multiple Myeloma: Bone marrow; Peripheral blood.

Myeloma: Bone marrow.

Natural Killer Cell Lymphoma: Lymph node.

Nephroblastoma: Kidney.

Non-Small Cell Lung Cancer: Blood; Lung; Peripheral blood.

Oesophageal Cancer: Blood.

Oligodendroglioma: Brain.

Oral Cancer: Oral cavity.

Oral Squamous Cell Carcinoma: Oral cavity; Peripheral blood.

Osteosarcoma: Bone.

Ovarian Cancer: Ascites; Ovarian cortex; Ovary; Peripheral blood.

Pancreatic Cancer: Blood vessel; Pancreas.

Pancreatic Ductal Adenocarcinomas: Pancreas.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Blood; Peripheral blood; Prostate.

Renal Cell Carcinoma: Kidney; Serum.

Renal Clear Cell Carcinoma: Lymph node.

Retinoblastoma: Eye.

Salivary Gland Tumor: Parotid gland; Salivary gland.

Sarcoma: Muscle.

Small Cell Lung Cancer: Lung.

Testicular Germ Cell Tumor: Peripheral blood; Testis.

Thyroid Cancer: Thyroid.

Tongue Cancer: Tongue.

Uterine Leiomyoma: Uterus.

Vascular Tumour: Lymph node.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#2.1 For Mouse tissue, tissue types are listed as follows:

Adipose tissue-related: Adipose tissue; White adipose tissue.

Bladder-related: Bladder.

Blood-related: Blood; Peripheral blood; Serum; Umbilical cord blood.

Bone-related: Bone; Bone marrow; Meniscus; Skeletal muscle; Spinal cord.

Brain-related: Brain; Cerebellum; Fetal brain; Hippocampus; Neural tube.

Breast-related: Mammary epithelium; Mammary gland.

Calvaria-related: Neonatal calvaria.

Ear-related: Cochlea; Inner Ear.

Embryo-related: Embryo; Embryoid body; Embryonic heart; Embryonic stem cell.

Esophagus-related: Esophagus.

Eye-related: Corneal epithelium; Eye; Ganglion cell layer of retina; Inner nuclear layer of retina; Lacrimal gland; Retina.

Fetus-related: Fetal brain; Fetal intestine; Fetal liver; Fetal lung; Fetal stomach; Placenta; Umbilical cord; Umbilical cord blood.

Gonad-related: Gonad; Ovary; Testis; Yolk sac.

Hair-related: Hair follicle.

Heart-related: Embryonic heart; Heart; Heart muscle; Neonatal heart.

Intestine-related: Colon; Colon epithelium; Fetal intestine; Gastrointestinal tract; Ileum; Intestinal crypt; Intestine; Mesenteric lymph node; Small intestine.

Kidney-related: Kidney; Mesonephros.

Liver-related: Fetal liver; Liver.

Lung-related: Bronchiole; Fetal lung; Lung; Trachea.

Lymph-related: Lymph node; Lymphoid tissue; Mesenteric lymph node; Peyer patch.

Muscle-related: Heart muscle; Muscle; Neonatal muscle; Skeletal muscle.

Neonate-related: Neonatal calvaria; Neonatal heart; Neonatal muscle; Neonatal pancreas; Neonatal rib; Neonatal skin.

Oral cavity-related: Submandibular gland; Taste bud.

Ovary-related: Ovary; Yolk sac.

Pancreas-related: Neonatal pancreas; Pancreas; Pancreatic islet.

Prostate-related: Prostate.

Skin-related: Dermis; Epidermis; Neonatal skin; Skin.

Spleen-related: Spleen.

Stomach-related: Fetal stomach; Gastrointestinal tract; Stomach.

Testis-related: Testis.

Uterus-related: Uterus.

Vessel-related: Aorta; Artery; Blood vessel; Carotid artery.

Others: Basilar membrane; Epithelium; Peritoneal cavity; Thymus.

#———————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————————

#2.2 For Mouse tissue about cancer, cancer types and the corresponding tissue types are listed as follows:

Breast Cancer: Lymph node; Breast; Lung.

Chronic Myeloid Leukemia: Blood.

Colon Cancer: Colon.

Colorectal Cancer: Lymph node; Colon; Colorectum.

Cutaneous Squamous Cell Carcinoma: Skin.

Hepatocellular Cancer: Blood; Liver.

Liver Cancer: Liver.

Lung Cancer: Lung.

Melanoma: Lung.

Pancreatic Cancer: Blood.

Papillary Thyroid Carcinoma: Thyroid.

Prostate Cancer: Prostate.

Renal Cell Carcinoma: Kidney.

Supratentorial Primitive Neuroectodermal Tumor: Brain. 
```
