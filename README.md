#### *VastDB* **case example: Identification and characterization of conserved** ***Srrm4*** **targets in mammals**

This document is a companion of the publication *Computational analysis of alternative splicing using vast-tools and the VastDB framework* and leads the reader through the code example of Section 6: **The identification of microexons as the main targets of the splicing factor *Srrm4***.



##### 1. *VastDB*: *Srrm4* expression across different cell and tissue types

The *VastDB* [Gene View Page of *Srrm4* in mouse](https://vastdb.crg.eu/gene/ENSMUSG00000063919@mm10) confirms the neural-specific expression of *Srrm4* across mouse cell and tissue types.

![](https://github.com/vastgroup/molbio2021_code_companion/blob/main/vastdb_srrm4.png)
*Expression of mouse Srrm4 during differentiation of embryonic stem cells (ESC) to glutamatergic neurons, using the cRPKM metrics*

#### 2. *vast-tools*: *Srrm4* splicing quantification and identification of *Srrm4* regulated exons

To quantify *Srrm4* splicing and identify *Srrm4* regulated exons with *vast-tools*, we download RNA-seq data for a knockdown (KD) of *Srrm4* in neuroblastoma N2A cells (SRP041656), including two replicates for KD and control, using matt retr_rnaseq:
```bash
matt retr_rnaseq accessions_mouse.txt
```
where [accessions_mouse.txt](https://github.com/vastgroup/molbio2021_code_companion/blob/main/accessions_mouse.txt) contains the GEO SRA IDs and sample names.

We use vast-tools align to process each RNA-seq sample separately.
```bash
vast-tools align CL_N2A_Cont_a_R1.fq.gz CL_N2A_Cont_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Cont_a 

vast-tools align CL_N2A_Cont_b_R1.fq.gz CL_N2A_Cont_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Cont_b 

vast-tools align CL_N2A_Srrm34_KD_a_R1.fq.gz CL_N2A_Srrm34_KD_a_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm34_KD_a 

vast-tools align CL_N2A_Srrm34_KD_b_R1-153.fq.gz CL_N2A_Srrm34_KD_b_R2.fq.gz  -sp mm10 -o vast_out/mm10/MMB --expr  --IR_version 2 -c 8 -n CL_N2A_Srrm34_KD_b 
```
