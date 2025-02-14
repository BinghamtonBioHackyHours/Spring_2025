1. What is the full path to your home directory?

```bash
robda@LAPTOP-5VM2NF6A MINGW64 ~
$ pwd
/c/Users/robda
```

2. Go up one directory?
    - How many files does it contain?

```bash
robda@LAPTOP-5VM2NF6A MINGW64 ~
$ ls | wc -l
32
```

    - How many directories?

```bash

robda@LAPTOP-5VM2NF6A MINGW64 ~
$ ls -F | grep '/' | wc -l
15
 ```


3. Make a directory called `hackyhour_problemset` in your home directory.

```bash

robda@LAPTOP-5VM2NF6A MINGW64 ~
$ mkdir hackyhour_problemset

```

4. Navigate into this new directory called `hackyhour_problemset`. Verify that you are in the correct directory by using `pwd`.

```bash
robda@LAPTOP-5VM2NF6A MINGW64 ~
$ cd hackyhour_problemset/

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ pwd
/c/Users/robda/hackyhour_problemset
```

5. Use `curl` or `wget` to copy <https://raw.githubusercontent.com/BinghamtonBioHackyHours/Spring_2025/refs/heads/main/data/week_03/sequences.nt.fa> from the web into your problemsets directory. If `wget` is not available on your system, use `curl -O` as an alternative.

```bash
robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ curl -O https://raw.githubusercontent.com/BinghamtonBioHackyHours/Spring_2025/refs/heads/main/data/week_03/sequences.nt.fa
```

6. Without using a text editor calculate or report these qualities for the file `sequences.nt.fa`.
  This file can be found here <https://raw.githubusercontent.com/BinghamtonBioHackyHours/Spring_2025/refs/heads/main/data/week_03/sequences.nt.fa>
      - How many lines does this file contain?   
      - How many characters?    (Hint: check out the options of wc)
      - What is the first line of this file?    (Hint: read the man page of head)
      - What are the last 3 lines?    (Hint: read the man page of tail)
      - How many sequences are in the file?    (Hint: use grep) (Note: The start of a sequence is indicated by a `>` character.)    

```bash
# - How many lines does this file contain?

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ wc -l sequences.nt.fa
148 sequences.nt.fa


# - How many characters?    (Hint: check out the options of wc)

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ wc -c sequences.nt.fa
10267 sequences.nt.fa

# - What is the first line of this file?    (Hint: read the man page of head)

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ head -n 1 sequences.nt.fa
>NM_001126114.2 Homo sapiens tumor protein p53 (TP53), transcript variant 3, mRNA

# - What are the last 3 lines?    (Hint: read the man page of tail)

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ tail -n 3 sequences.nt.fa
TTTATAGCTGTTGGAAGGACTAGGTCTTCCCTAGCCCCCCCAGTGTGCAAGGGCAGTGAAGACTTGATTG
TACAAAATACGTTTTGTAAATGTTGTGCTGTTAACACTGCAAATAAACTTGGTAGCAAACACTTCCAAAA
AAAAAAAAAAAAAA

# - How many sequences are in the file?    (Hint: use grep) (Note: The start of a sequence is indicated by a `>` character.)

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ grep '>' sequences.nt.fa | wc -l
2

```


7. Rename `sequences.nt.fa` to `cancer_genes.fasta`. (Hint: read the man page for mv)

```bash
robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ mv sequences.nt.fa cancer_genes.fasta
```

8. Copy/download cuffdiff.txt, to your problemset directory using `curl` or `wget`. Here is the url you can use: https://raw.githubusercontent.com/BinghamtonBioHackyHours/Spring_2025/refs/heads/main/data/week_03/cuffdiff.txt. Use -O is the letter O not the number zero 0.

```
curl -O https://raw.githubusercontent.com/BinghamtonBioHackyHours/Spring_2025/refs/heads/main/data/week_03/cuffdiff.txt
```


9. Do the following to `cuffdiff.txt`. The descriptions of each column in the file are in the table below.
    - Look at the first few lines of the file
    - Sort the file by log fold change 'log2(fold_change)', from highest to lowest, and save in a new file in your directory called sorted.cuffdiff.out
    - Sort the file (log fold change highest to lowest) then print out only the first 100 lines. Save in a file called `top100.sorted.cuffdiff.out`.
    - Sort the file by log fold change, print out the top 100, print only first column. This will be a list of the genes with the largest change in expression. Make sure your list is sorted by gene name and is unique. Save this curated list in a file called `differentially.expressed.genes.txt`.

```bash
# - Look at the first few lines of the file
robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ head cuffdiff.txt

# - Sort the file by log fold change 'log2(fold_change)', from highest to lowest, and save in a new file in your directory called sorted.cuffdiff.out

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ sort -k 10 -g -r cuffdiff.txt >sorted.cuffdiff.out

# - Sort the file (log fold change highest to lowest) then print out only the first 100 lines. Save in a file called `top100.sorted.cuffdiff.out`.

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ sort -k 10 -g -r cuffdiff.txt | head -n 100 >top100.sorted.cuffdiff.out

# - Sort the file by log fold change, print out the top 100, print only first column. This will be a list of the genes with the largest change in expression. Make sure your list is sorted by gene name and is unique. Save this curated list in a file called `differentially.expressed.genes.txt`.

robda@LAPTOP-5VM2NF6A MINGW64 ~/hackyhour_problemset
$ sort -k 10 -g -r cuffdiff.txt | head -n 100 | cut -f 1 | sort | uniq  >differentially.expressed.genes.up.txt

```
