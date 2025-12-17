## *VastDB* **case example: Identification and characterization of conserved** ***Srrm4*** **targets in mammals**

This document is a companion of the publication *Computational analysis of alternative splicing using vast-tools and the VastDB framework*. There, we show how the tools and resources in the vastDB framework (see figure below) can be used to (i) quantify AS and identify differentially spliced AS events using RNA-seq data ([vast-tools](https://github.com/vastgroup/vast-tools)), (ii) perform multiple genomic and sequence analyses for sets of AS events ([Matt](https://gitlab.com/aghr/matt)), (iii) identify AS events with genomic and regulatory conservation among species ([ExOrthist](https://github.com/biocorecrg/ExOrthist)), and (iv) help with the biological interpretation of the results, and, ultimately, with the identification of interesting AS events to design wet-lab experiments ([VastDB](https://vastdb.crg.eu/) and [PastDB](https://pastdb.crg.eu/)).  

![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastdb_framework.png)
*The VastDB framework: schematic representation of the different modules and commands of vast-tools, Matt and ExOrthist, as well as the different layers of information covered in VastDB and PastDB.*  

This tutorial leads the reader through the code example of Section 6: **The identification of microexons as the main targets of the splicing factor *Srrm4***.

### Table of contents
* [1. *VastDB*: *Srrm4* expression across different cell and tissue types](#1-vastdb-srrm4-expression-across-different-cell-and-tissue-types)  
* [2. *vast-tools*: *Srrm4* splicing quantification and identification of *Srrm4* regulated exons](#2-vast-tools-srrm4-splicing-quantification-and-identification-of-srrm4-regulated-exons)  
* [3. GO-term enrichment analysis using gene-id lists provided by *vast-tools*](#3-go-term-enrichment-analysis-using-gene-id-lists-provided-by-vast-tools)  
* [4. *Matt*: Identifying potential genomic and sequence features associated with Srrm4 regulation](#4-matt-identifying-potential-genomic-and-sequence-features-associated-with-srrm4-regulation)  
* [5. Use of *VastDB* resources to identify tissue specific AS events](#5-use-of-vastdb-resources-to-identify-tissue-specific-as-events)  
* [6. Use of *ExOrthist* for conservation analysis of exons](#6-use-of-ExOrthist-for-conservation-analysis-of-exons)  

### 1. *VastDB*: *Srrm4* expression across different cell and tissue types

The *VastDB* [Gene View Page of *Srrm4* in mouse](https://vastdb.crg.eu/gene/ENSMUSG00000063919@mm10) confirms the neural-specific expression of *Srrm4* across mouse cell and tissue types.

![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastdb_srrm4.png)
*Expression of mouse Srrm4 across cell and tissue types form the main VastDB Gene view sample panel plot, using the cRPKM metrics.*  

We see more expression profiles of interest in the *VastDB* Special Datasets section of the *Srrm4* Gene View Page, e.g. *Neural differentiation time course* showing the temporal dynamics of *Srrm4* expression during neuronal differentiation.
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastdb_srrm4_special.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics.* 

### 2. *vast-tools*: *Srrm4* splicing quantification and identification of *Srrm4* regulated exons

To quantify *Srrm4* splicing and identify *Srrm4* regulated exons with *vast-tools*, we download RNA-seq data for a knockdown (KD) of *Srrm4* in neuroblastoma N2A cells (SRP041656), including two replicates for KD and control, using *Matt*:
```bash
matt retr_rnaseq accessions_mouse.txt
```
where [accessions_mouse.txt](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/vast-tools_files/accessions_mouse.txt) contains the GEO SRA IDs and sample names.

We use vast-tools align to process each RNA-seq sample separately with respect to the mouse mm10 transcriptome and specifying `vast_out/mm10/MMB` as the central output directory
```bash
vast-tools align CL_N2A_Srrm4_Cont_a_R1.fq.gz CL_N2A_Srrm4_Cont_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm4_Cont_a 

vast-tools align CL_N2A_Srrm4_Cont_b_R1.fq.gz CL_N2A_Srrm4_Cont_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm4_Cont_b

vast-tools align CL_N2A_Srrm4_KD_a_R1.fq.gz CL_N2A_Srrm4_KD_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm4_KD_a

vast-tools align CL_N2A_Srrm4_KD_b_R1.fq.gz CL_N2A_Srrm4_KD_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm4_KD_b 
```

Afterwards, we combine the results of the four datasets and generate a single results table, called INCLUSION table, with PSIs for each sample
```bash
vast-tools combine -sp mm10 -o vast_out/mm10/MMB
```

Inside the central output directory (`vast_out/mm10/MMB`) we apply ```vast-tools compare``` on the INCLUSION table to extract the differentially regulated AS events, here with thresholds |ΔPSI| ≥ 25 and minimum range ≥ 5. Plus, the flags `--GO` and `--print_sets` enable the generation of lists of gene IDs to perform Gene Ontology analyses and of control AS events for later Matt analyses, respectively:
```bash
cd vast_out/mm10/MMB

vast-tools compare INCLUSION_LEVELS_FULL-mm10-4.tab \ 
   -a CL_N2A_Srrm4_Cont_a,CL_N2A_Srrm4_Cont_b \
   -b CL_N2A_Srrm4_KD_a,CL_N2A_Srrm4_KD_b \ 
   --print_dPSI --GO -sp mm10 --print_sets \ 
   --min_dPSI 25 --min_range 5 \
   -name_A Control -name_B Srrm4_KD
```

The summary output confirms a clear tendency: most regulated AS events are cassette exons, especially microexons with length ≤ 27 nt that show higher inclusion in the control (108 vs. 0 for microexons and 52 vs. 9 for longer exons), consistent with the known role of Srrm4 enhancing inclusion of very short exons.  

<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vasttools_compare_output.png" width=700 height=550 />  
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vasttools_compare_output.png) -->

In addition to the summary statistics, the run of ```vast-tools compare``` with the `--GO` and `--print_sets` flags generates the following useful files:  
1. **DiffAS**: The set of differentially regulated AS events with average ΔPSI in the last column
```bash
DiffAS-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab
```
2. **AS_NC**: AS events (10<PSI<90 in at least one group) that do not change (|ΔPSI| ≤ 5)
```bash
AS_NC-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI-Max_dPSI5.tab
```
3. **CS**: Constitutive exons (PSI > 95 in both groups) 
```bash
CS-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab
```
4. **CR**: Cryptic exons (PSI < 5 in both groups)
```bash
CR-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab
```
5. **GeneIDs** of cassette exons (AltEx) and background (BG) genes for GO analysis
```bash
AltEx-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.txt
BG-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.txt
```

### 3. GO-term enrichment analysis using gene-id lists provided by *vast-tools*
To perform a GO enrichment analysis, we upload the two Gene-ID lists generated by ```vast-tools compare``` to [DAVID](https://davidbioinformatics.nih.gov/summary_new.jsp), download the chart results and plot the p-values (-log10) of the resulting categories as histograms. The GO terms associated with genes containing cassettes exons reveal enrichments in gene functions related to GTPase regulation, synaptic organization and cytoskeleton (as previously described for *Srrm4*).  

<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/GOTEA_srrm4.png" width=600 height=700 />  
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/GOTEA_srrm4.png) -->  

*Gene Ontology enrichment analyses using DAVID for genes containing exons that are differentially regulated upon Srrm4 KD. * 

### 4. *Matt*: Identifying potential genomic and sequence features associated with Srrm4 regulation

We then use *Matt* to identify potential genomic and sequence features associated with Srrm4-regulated exons. The [mm10.gtf.gz](http://vastdb.crg.eu/FRAMEWORK/mm10.gtf.gz) and [mm10.fasta.gz](http://vastdb.crg.eu/FRAMEWORK/mm10.fasta.gz) files required by this analysis can be downloaded from the relative link. The files need to be uncompressed in order to be used as *Matt*'s input.  
We first apply *Matt*'s table manipulation commands (`add_val`, `rand_rows`, `add_rows` and `get_rows`) to prepare the input table containing all the exons to be compared together with the relative group ID (reported in the column GROUP). We exploit ` grep -P "(MmuEX|EVENT)"` to selectively extract exon AS events from the *vast-tools* table, and we randomly down-sample 1000 non-changing, constitutive, and cryptic exons in order to preserve a reasonable dataset size.  

```bash
matt add_val AS_NC-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI-Max_dPSI5.tab GROUP AS_NC \
    | grep -P "(MmuEX|EVENT)" \
    | matt rand_rows - 1000 > tmp.tab

matt add_val CR-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab GROUP CR \
    | grep -P "(MmuEX|EVENT)" \
    | matt rand_rows - 1000\
    | matt add_rows tmp.tab -

matt add_val CS-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab GROUP CS \
    | grep -P "(MmuEX|EVENT)" \
    | matt rand_rows - 1000  \
    | matt add_rows tmp.tab -

matt get_rows DiffAS-mm10-4-dPSI25-range5-min_ALT_use25-upreg_ALT_Control-vs-Srrm4_KD-with_dPSI.tab dPSI[-100,-15] \
    | matt add_val - GROUP Srrm4_DOWN \
    | grep -P "(MmuEX|EVENT)" \
    | matt add_rows tmp.tab -
```

The table `tmp.tab` contains exons (i.e. MmuEX) downregulated by Srrm4 KD, and up to 1000 random constitutive exons, non-regulated alternative exons and cryptic exons. 

Next, we use `matt get_vast` to extract the relevant genomic information (chromosome, start and end coordinates and strand) from the vast-tools formatted table: 
```bash
matt get_vast tmp.tab COORD FullCO COMPLEX LENGTH -gtf mm10.gtf > Matt_input_Srrm4_ex.tab
```
Eventually, we run `matt cmpr_exons` on the resulting input table (Matt_input_Srrm4_ex.tab) to extract exon related features and compare them among the four exon groups: 
```bash
matt cmpr_exons Matt_input_Srrm4_ex.tab START END SCAFFOLD STRAND GENEID \
       mm10.gtf mm10.fasta Mmus 150 GROUP[Srrm4_DOWN,CR,AS_NC,CS] \
       Matt_Srrm4_KD -notrbts -colors:red,white,lightgray,darkgray
``` 
The automatically generated PDF summary shows the main regulatory features known to be associated with Srrm4-regulated exons: weak 3′ splice sites but strong 5′ splice sites, as well as much smaller exon lengths.  
<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/matt_exon_length.png" width=600 height=500 />  
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/matt_exon_length.png) -->
*Example of a boxplot from matt cmpr_exons visualizing the much shorter length distribution of Srrm4-regulated exons compared to other exon sets.*

Next, we want to investigate the enrichment of motifs associated with Srrm4-regulated exons. Since UGC motifs are known to be enriched in the upstream intron, we use `matt rna_maps` to confirm this association. First, we need to create a table with motifs to be included in this analysis containing only the UGC motif.
```bash
echo -e 'TYPE\tNAME\tEXPR_FILE\tTHRESH\tBGMODEL\nREGEXP\tUGC\tTGC\tNA\tNA' > ugc_motif.tab
```
Then, we generate the RNA map with Matt:
```bash
matt rna_maps Matt_input_Srrm4_ex.tab UPSTRM_EX_BORDER START END \
      DOSTRM_EX_BORDER SCAFFOLD STRAND GROUP 15 50 150 mm10.fasta \
      ugc_motif.tab TYPE NAME EXPR_FILE THRESH BGMODEL \
      -d UGC_map_Matt_Srrm4_KD
```
where 15 is the length of the sliding window, 50 and 150 are the number of exon and intron positions to be considered, respectively.
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/ugc_rna_map.png)
*RNA map for the UGC motif. The motif is strongly enriched in a window from approx. -30 to -5 upstream of Srrm4-regulated exons, but not in other exons.*

### 5. Use of *VastDB* resources to identify tissue specific AS events
We download two files from *VastDB*: 
1. [PSI_TABLE-mm10.tab.gz](https://vastdb.crg.eu/downloads/mm10/PSI_TABLE-mm10.tab.gz) containing the PSIs across the main sample panel for each AS event.  
2. [PROT_DISORDER-mm10.tab.gz](https://vastdb.crg.eu/downloads/mm10/PROT_DISORDER-mm10.tab.gz) with the overlap with disordered regions for each AS event.  
We first calculate the ΔPSI between neural and non-neural samples. We use the utility script [Get_Tissue_Specific_AS.pl](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Scripts/Get_Tissue_Specific_AS.pl). This script takes as input a vast-tools INCLUSION table and a config file (a tab-separated table with the tissue group for each sample) and identifies tissue-specific AS events. If the option `--test_tis` is provided as below, it will generate a table with the average ΔPSI between neural and non-neural samples (as specified in the config file) for all AS events with sufficient read coverage for the comparison (at least five replicates in each group; -min_rep 5).
```bash
perl Get_Tissue_Specific_AS.pl PSI_TABLE-mm10.tab.gz \ 
      -g Config_Neural.txt -N 5 \ 
      -test_tis Neural -min_rep 2
```
With this information, we plot the ΔPSI per type of exon, which shows a very strong tendency for neurally upregulated exons among Srrm4-regulated exons.  
<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_A.png" width="500" height="500" />  
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_A.png) -->

Next, we use the disorder information downloaded from *VastDB* and plot the average percentage of disorder residues for the alternative as well as the upstream (C1) and downstream (C2) exons for each exon set. As expected for tissue-specific exons, *Srrm4*-regulated exons more often encode disorder regions.  

<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_B.png" width="600" height="500" />  
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_B.png) -->
However, some Srrm4-regulated microexons are also known to be inserted within structured domains, as exemplified by a 15-nt exon in *Vav2* with VastID [MmuEX0051282](https://vastdb.crg.eu/event/MmuEX0051282@mm10)  

<!-- -->
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_C.png)

### 6. Use of *ExOrthist* for conservation analysis of exons
Finally, we perform a conservation analysis using *ExOrthist*. We run *ExOrthist* `main.nf` for mouse and human, using default parameters for the short evolutionary distance range, adding all vast-tools exons, and using Ensembl 1-to-1 orthologs as gene orthogroups. 
```bash
nextflow main.nf
```
`main.nf` needs two config files in the pipeline execution directory. First a [nextflow.config](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/nextflow.config) file, which is also automatically downloaded during ExOrthist installation and contains the pipeline configuration properties. Second, a [params.config](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/params.config) file specifying the required parameters for the run, the location of the input files as well as of the output folder (`hg38_mm10_output`, for this example). Among other files of interest, the output folder will include the [EX_clusters.tab](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/EX_clusters.tab.gz) file, which will contain all exon orthogroups within the provided human and mouse gene orthogroups. The following files are necessary to replicate the run, and need to be organized into a `hg38_mm10_inputs` input folder and relative subfolders as follows (and specified in the params.config):  
* [evodists.txt](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/evodists.txt): evolutionary distance range between human and mouse.  
* [hg38_mm10_annotated.tab.gz](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/hg38_mm10_annotated.tab.gz): gene orthogroups derived from human and mouse Ensembl 1:1 orthologs (i.e. all orthogroups include one gene per species).  
* **GTF** folder: it has to contain human and mouse annotation files ([hg38.gtf.gz](http://vastdb.crg.eu/FRAMEWORK/hg38.gtf.gz), [mm10.gtf.gz](http://vastdb.crg.eu/FRAMEWORK/mm10.gtf.gz)).  
* **GENOME** folder: it has to contain human and mouse genome fasta files ([hg38.fasta.gz](http://vastdb.crg.eu/FRAMEWORK/hg38.fasta.gz), [mm10.fasta.gz](http://vastdb.crg.eu/FRAMEWORK/mm10.fasta.gz)).  
* **EXTRAEXONS** folder: it has to contains human and mouse vast-tools identified exons ([hg38_extraexons.tab](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/hg38_extraexons.tab.gz), [mm10_extraexons.tab](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/mm10_extraexons.tab.gz)).  


Next, we perform two types of analysis with `compare_exon_sets.pl`: 
1. genome conservation of all identified mouse exon sets (DiffAS, AS_NC, CR, CS) in human.  
2. genome and regulatory conservation of Srrm4-regulated (DiffAS) exons between human and mouse.  

From the first comparison, using a single list of Srrm4-regulated exons, we obtain:
```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 
     -exon_list_sp1 Exons_mm10-Srrm4_KD.txt \ 
     -main_folder hg38_mm10_output/
```
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab1.png)
Among other statistics, ExOrthist highlights that the percentage of mouse Srrm4-regulated exons hosted by genes with 1:1 human orthologs with an exon ortholog in human is 89%.   

We then run the same analysis for each control exon set (AS_NC, CR, CS), and plot the percentage of genome conservation in human, showing higher conservation of Srrm4-regulated exons compared to non-regulated alternative and cryptic exons.

```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-AS_NC.txt -main_folder hg38_mm10_output/

perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-CS.txt -main_folder hg38_mm10_output/

perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-CR.txt -main_folder hg38_mm10_output/
```  

<img align="middle" src="https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_fig1.png" width="600" height="500" />    
<!-- ![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_fig1.png) -->

*Genome conservation between human and mouse for various exon sets. *  

To assess the regulatory conservation of *Srrm4*-regulated exons between mouse and human, and identify orthologs exons regulated in both species, we utilized compare_exon_sets.pl for two exon lists. The second list, *SRRM4*-regulated exons in human, was obtained by analyzing with vast-tools a RNA-seq dataset (**SRP149913**) in which *SRRM4* or GFP (as control) were ectopically expressed in HEK293 cells.  

```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-4-dPSI25.txt -exon_list_sp2 \
      Exons_hg38-2-dPSI50.txt -main_folder hg38_mm10_output/ -print_out
```
This provides us a richer output about the percent of conservation at the gene level:
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab2_1.png)
and at the exon level:
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab2_2.png)

Moreover, with the option `-print_out`, `compare_exon_sets.pl` generates a text file with all the exon orthogroups (i.e. clusters of conserved exons between human and mouse) containing Srrm4-regulated exons for either one or both species. The format of the output file (`Conserved_exons-mm10-hg38.tab`) is the following:  

![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_output.png)
