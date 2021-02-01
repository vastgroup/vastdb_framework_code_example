## *VastDB* **case example: Identification and characterization of conserved** ***Srrm4*** **targets in mammals**

This document is a companion of the publication *Computational analysis of alternative splicing using vast-tools and the VastDB framework* and leads the reader through the code example of Section 6: **The identification of microexons as the main targets of the splicing factor *Srrm4***.



### 1. *VastDB*: *Srrm4* expression across different cell and tissue types

The *VastDB* [Gene View Page of *Srrm4* in mouse](https://vastdb.crg.eu/gene/ENSMUSG00000063919@mm10) confirms the neural-specific expression of *Srrm4* across mouse cell and tissue types.

![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/vastdb_srrm4.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics.*

We see more expression profiles of interest in the *VastDB* Special Datasets section of the *Srrm4* Gene View Page, e.g. *Neural differentiation time course* showing the temporal dynamics of *Srrm4* expression during neuronal differentiation.
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/Figures/vastdb_srrm4_special.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics.* 

### 2. *vast-tools*: *Srrm4* splicing quantification and identification of *Srrm4* regulated exons

To quantify *Srrm4* splicing and identify *Srrm4* regulated exons with *vast-tools*, we download RNA-seq data for a knockdown (KD) of *Srrm4* in neuroblastoma N2A cells (SRP041656), including two replicates for KD and control, using *Matt*:
```bash
matt retr_rnaseq accessions_mouse.txt
```
where [accessions_mouse.txt](https://github.com/vastgroup/molbio2021_code_companion/blob/main/vast-tool_files/accessions_mouse.txt) contains the GEO SRA IDs and sample names.

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

Inside the central output directory we apply ```vast-tools compare``` to the INCLUSION table to extract the differentially regulated AS events, here with thresholds |ΔPSI| ≥ 25 and minimum range ≥ 5, as well as lists of gene IDs to perform Gene Ontology analyses and of control AS events for later Matt analyses:
```bash
cd vast_out/mm10/MMB

vast-tools compare INCLUSION_LEVELS_FULL-mm10-4.tab \ 
   -a CL_N2A_Srrm4_Cont_a,CL_N2A_Srrm4_Cont_b \
   -b CL_N2A_Srrm4_KD_a,CL_N2A_Srrm4_KD_b \ 
   --print_dPSI --GO -sp mm10 --print_sets \ 
   --min_dPSI 25 --min_range 5 \
   -name_A Control -name_B Srrm4_KD \
```

The summary output confirms a clear tendency: most regulated AS events are exons, especially microexons with length ≤ 27 nt that show higher inclusion in the control (108 vs. 0 for microexons and 52 vs. 9 for longer exons), consistent with the known role of Srrm4 enhancing inclusion of very short exons.
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/vasttools_compare_output.png)

In addition to the descriptive tables, the run of ```vast-tools compare``` generates several useful output files, among them:
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
To perform a GO enrichment analysis, we upload the two Gene-ID lists to [DAVID](https://david.ncifcrf.gov), download the chart results and plot the p-values (-log10) of the resulting categories as histograms. The found GO terms reveal enrichment in gene functions associated with GTPase regulation, synaptic organization and cytoskeleton, as previously described for *Srrm4*.
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/GOTEA_srrm4.png)
*Gene Ontology enrichment analysis using DAVID.* 

### 4. *Matt*: Identifying potential genomic and sequence features associated with Srrm4 regulation

To identify potential genomic and sequence features associated with Srrm4-regulated exons we will use *Matt*. 
To prepare the input table where we would need to gather all exon test to be compared together with a column GROUP with exon group IDs, we apply Matt's table manipulation commands (add_val, rand_rows, add_rows and get_rows). To extract from the *vast-tools* table exon AS events only, we exploit ` grep -P "(MmuEX|EVENT)"` and to keep the data set at a reasonable size, we randomly down-sample 1000 non-changing, constitutive, and cryptic exons.

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

Next, we use `matt get_vast` to extract from the vast-tools formatted table the relevant genomic information (chromosome, start and end coordinates and strand):
```bash
matt get_vast tmp.tab COORD FullCO COMPLEX LENGTH -gtf mm10.gtf > Matt_input_Srrm4_ex.tab
```
Eventually, to extract exon related features and compare them among the four exon groups, we run `matt cmpr_exons` on the resulting input table (Matt_input_Srrm4_ex.tab):
```bash
matt cmpr_exons Matt_input_Srrm4_ex.tab START END SCAFFOLD STRAND GENEID \
       mm10.gtf mm10.fasta Mmus 150 GROUP[Srrm4_DOWN,CR,AS_NC,CS] \
       Matt_Srrm4_KD -notrbts -colors:red,white,lightgray,darkgray
``` 
The automatically generated PDF summary shows the main regulatory features known to be associated with Srrm4-regulated exons: weak 3′ splice sites but strong 5′ splice sites, as well as much smaller exon lengths.
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/matt_exon_length.png)
*Example of a boxplot from matt cmpr_exons visualizing the much shorter length distribution of Srrm4-regulated exons compared to other exon sets.*

Next, we want to investigate the enrichment of motifs associated with Srrm4-regulated exons. Since UGC motifs are known to be enriched in the upstream intron we use `matt rna_maps` to confirm this association. First, we need to create a table with motifs to be included in this analysis containing only the UGC motif.
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
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/Figures/ugc_rna_map.png)
*RNA map for the UGC motif. The motif is strongly enriched in a window from approx. -30 to -5 upstream of Srrm4-regulated exons, but not in other exons.*

### 5. Use of *VastDB* resources to identify tissue specific AS events
We download from *VastDB* two files: 
1. [PSI_TABLE-mm10.tab.gz](https://vastdb.crg.eu/downloads/mm10/PSI_TABLE-mm10.tab.gz) containing the PSIs across the main sample panel for each AS event
2. [PROT_DISORDER-mm10.tab.gz](https://vastdb.crg.eu/downloads/mm10/PROT_DISORDER-mm10.tab.gz) with the overlap with disordered regions for each AS event
We first calculate the ΔPSI between neural and non-neural samples. We use the utility script [Get_Tissue_Specific_AS.pl](https://github.com/vastdb-pastdb/pastdb/blob/master/bin/Get_Tissue_Specific_AS.pl) from  [PastDB](https://github.com/vastdb-pastdb/pastdb). This script takes as input a vast-tools INCLUSION table and a config file (a tab-separated table with the tissue groups for each sample) and identifies tissue-specific AS events. If the option `--test_tis` is provided, it will generate a table with the average ΔPSI between neural and non-neural samples (as specified in the config file) for all AS events with sufficient read coverage for the comparison (at least five replicates in each group; -min_rep 5).
```bash
perl Get_Tissue_Specific_AS.pl PSI_TABLE-mm10.tab.gz \ 
      -g Config_Neural.txt -min_N 2 \ 
      -test_tis Neural -min_rep 5
```
With this information, we plot the ΔPSI per type of exon, which shows a very strong tendency for neurally upregulated exons among Srrm4-regulated exons.
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_A.png)

Next, we use the disorder information downloaded from *VastDB* and plot the average percentage of disorder residues for the alternative as well as the upstream (C1) and downstream (C2) exons for each exon set. As expected for tissue-specific exons, *Srrm4*-regulated exons more often encode disorder regions.
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_B.png)
However, some *Srrm4*-regulated microexons are also known to be inserted within structured domains, as exemplified by a 15-nt exon in *Vav2* with VastID [MmuEX0051282](https://vastdb.crg.eu/event/MmuEX0051282@mm10).
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/vastbd_resources_C.png)

### 6. Use of *ExOrthist* for conservation analysis of exons
Finally, we perform a conservation analysis using *ExOrthist*. We run *ExOrthist* `main.nf` for mouse and human, using default parameters for the short evolutionary distance range, adding all vast-tools exons, and using Ensembl 1-to-1 orthologs as gene orthogroups. 
```bash
nextflow main.nf
```
`main.nf` needs a Nextflow config file [params.config](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/params.config) specifying the location of the input files as well as the output folder; `hg38_mm10_output` in this example. Among other output files of interest, this folder contains the file [EX_clusters.tab](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/ExOrthist_files/EX_clusters.tab.gz) which contains all exon orthogroups between human and mouse.

Next, we perform two types of analysis with `compare_exon_sets.pl`: 
1. genome conservation of each mouse exon set in human
2. evolutionary conservation of human and mouse Srrm4-regulated exons

From the first comparison, using a single list, we obtain:
```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 
     -exon_list_sp1 Exons_mm10-Srrm4_KD.txt \ 
     -main_folder hg38_mm10_output/
```
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab1.png)
Among other statistics, the percentage of exons with an exon ortholog in human, for those genes with 1-to-1 orthologs is 89%. 

We then run similar analyses for each exon set (non-changing, constitutive, cryptic), and plot the percentage of genome conservation in human, showing a much higher conservation of Srrm4-regulated exons compared to non-regulated alternative exons.
```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-AS_NC.txt -main_folder hg38_mm10_output/

perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-CS.txt -main_folder hg38_mm10_output/

perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-CR.txt -main_folder hg38_mm10_output/
```
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_fig1.png)
*Genome conservation between human and mouse for various exon sets.*

To assess the regulatory conservation of *Srrm4*-regulated exons between mouse and human, and identify ortholog exons regulated in both species, we utilize `compare_exon_sets.pl` for two lists. The second list, *Srmm4*-regulated exons, was obtained by analzying with vast-tools an RNA-seq dataset with SRRM4 or GFP (as control) ectopically expressed in HEK293 cells.
```bash
perl ~/ExOrthist/bin/compare_exon_sets.pl -sp1 mm10 -sp2 hg38 \
      -exon_list_sp1 Exons_mm10-4-dPSI25.txt -exon_list_sp2 \
      Exons_hg38-2-dPSI50.txt -main_folder hg38_mm10_output/ -print_out
```
This provides us a richer output about the percent of conservation on the gene level:
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab2_1.png)
and on the exon level:
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_tab2_2.png)

Moreover, with option `-print_out` *ExOrthist*  generates a list of orthologous exons, containing *Srrm4*-regulated exons in one or the two species (Conserved_exons-mm10-hg38.tab), with the following format:
![](https://github.com/vastgroup/vastdb_framework_code_example/blob/main/Figures/exorthist_output.png)
