# Smart-SNPer
An automated pipeline that efficiently designs crRNA and isothermal primers specifically optimized for single-tube, Cas12a-mediated SNP detection. This end-to-end solution rapidly converts raw SNP information into ready-to-use assay components within hours, thereby eliminating the need for manual multistep processing and ensuring optimal performance in the final one-pot RPA-Cas12a detection system.version 2.1, updated on May 9, 2025.

# Usage
Download sequence.fasta, reference genome, GCF_000001405, and Smart-SNPerv4.exe. Follow the instructions provided within Smart-SNPerv4.exe for further steps.

Branch A only needs to download the reference sequence. Branch B requires the download of Chrome browser and reference sequence. C branches require sequence.fasta, reference genome and GCF_000001405. Suggested B branch for personal use.

During operation, it is important not to blindly use the default genome. First, test with a smaller-scale background genome to ensure that the computer configuration can support the execution of Smart-SNPer. If the default genome is required, ensure that the computer has at least 32GB of RAM and a CPU with performance superior to the i5-12400F

The rs121913237 folder contains the results of rs121913237 run at branch B using the default genome, which on a personal computer takes about 7h per crRNA.

# sequence.fasta, reference genome and GCF_000001405

These three files are larger, so please access them via the links below:

https://drive.google.com/file/d/1-kvRtvxbfcFFYzZOrCPKtA3eqKVPHJVG/view?usp=drive_link

https://drive.google.com/file/d/1Iu-GHBVHntzsv-l2RjizKoa0EePZw5Jx/view?usp=drive_link

https://drive.google.com/file/d/1xu4a70-9X01b9WWlxJxw_vFFuFMZY2aG/view?usp=sharing

# Dependencies：
Python v3.8.7

biopython v1.81

selenium v4.23.1

requests v2.31.0

bio v1.6.0

ViennaRNA v2.6.4

pandas v2.0.3

beautifulsoup4 v4.12.3

# Citations:
Lorenz, R., Bernhart, S. H., Höner zu Siederdissen, C., Tafer, H., Flamm, C., Stadler, P. F., & Hofacker, I. L. (2011). ViennaRNA Package 2.0. Algorithms for molecular biology, 6, 1-14.https://doi.org/10.1186/1748-7188-6-26

Asgar H Ansari, Manoj Kumar, Sajal Sarkar, Souvik Maiti, Debojyoti Chakraborty (2023) CriSNPr: a single interface for the curated and de novo design of gRNAs for CRISPR diagnostics using diverse Cas systems eLife 12:e77976 https://doi.org/10.7554/eLife.77976

Lin, C., Chen, F., Huang, D., Li, W., He, C., Tang, Y., ... & Huang, L. (2023). A universal all-in-one RPA-Cas12a strategy with de novo autodesigner and its application in on-site ultrasensitive detection of DNA and RNA viruses. Biosensors and Bioelectronics, 239, 115609.https://doi.org/10.1016/j.bios.2023.115609

