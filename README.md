## *VastDB* **case example: Identification and characterization of conserved** ***Srrm4*** **targets in mammals**

This document is a companion of the publication *Computational analysis of alternative splicing using vast-tools and the VastDB framework* and leads the reader through the code example of Section 6: **The identification of microexons as the main targets of the splicing factor *Srrm4***.



### 1. *VastDB*: *Srrm4* expression across different cell and tissue types

The *VastDB* [Gene View Page of *Srrm4* in mouse](https://vastdb.crg.eu/gene/ENSMUSG00000063919@mm10) confirms the neural-specific expression of *Srrm4* across mouse cell and tissue types.

![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/vastdb_srrm4.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics.*

We see more expression profiles of interest in the *VastDB* Special Datasets section of the *Srrm4* Gene View Page, e.g. *Neural differentiation time course* showing the temporal dynamics of *Srrm4* expression during neuronal differentiation.
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/vastdb_srrm4_special.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics.* 

### 2. *vast-tools*: *Srrm4* splicing quantification and identification of *Srrm4* regulated exons

To quantify *Srrm4* splicing and identify *Srrm4* regulated exons with *vast-tools*, we download RNA-seq data for a knockdown (KD) of *Srrm4* in neuroblastoma N2A cells (SRP041656), including two replicates for KD and control, using *Matt*:
```bash
matt retr_rnaseq accessions_mouse.txt
```
where [accessions_mouse.txt](https://github.com/vastgroup/molbio2021_code_companion/blob/main/accessions_mouse.txt) contains the GEO SRA IDs and sample names.

We use vast-tools align to process each RNA-seq sample separately with respect to the mouse mm10 transcriptome and specifying `vast_out/mm10/MMB` as the central output directory
```bash
vast-tools align CL_N2A_Cont_a_R1.fq.gz CL_N2A_Cont_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Cont_a 

vast-tools align CL_N2A_Cont_b_R1.fq.gz CL_N2A_Cont_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Cont_b 

vast-tools align CL_N2A_Srrm34_KD_a_R1.fq.gz CL_N2A_Srrm34_KD_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm34_KD_a 

vast-tools align CL_N2A_Srrm34_KD_b_R1-153.fq.gz CL_N2A_Srrm34_KD_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm34_KD_b 
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
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/vasttools_compare_output.png)

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
![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/GOTEA_srrm4.png)
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

