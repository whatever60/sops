# Demultiplexing Novogene sequencing results

whatever60 2024.12.05

## Contents
1. [Get a Linux machine](#get-a-linux-machine)
2. [Install `bcl-convert` on the Linux machine](#install-bcl-convert-on-the-linux-machine)
3. [Get data from Novogene](#get-data-from-novogene)
4. [Prepare sample sheet](#prepare-sample-sheet)
4.1. [Name your experiment](#name-your-experiment)
4.2. [Know your instrument](#know-your-instrument)
4.3. [Figure out which lane is used](#figure-out-which-lane-is-used)
4.4. [Decide cycle configuration](#decide-cycle-configuration)
4.5. [Another example](#another-example)
4.6. [Fetch your index sequences](#fetch-your-index-sequences)
4.6.1. [Mixed-length primer](#mixed-length-primer)
4.7. [Compile the sample sheet](#compile-the-sample-sheet)
5. [Run `bcl-convert`](#run-bcl-convert)
6. [Examine output file size](#examine-output-file-size)
7. [Upload data onto AWS](#upload-data-onto-aws)
8. [Terminate AWS instance](#terminate-aws-instance)

## Get a Linux machine

Demultiplexing Novogene data requires modest computation. An AWS instance with 32GB memory and 1TB storage is recommended. The following steps assume your operating system is Ubuntu/Debian.

>Alternatively a Windows WSL2 will suffice so that you can do everything locally, and I think with some modification you can also use MacBook, though I haven't tested it.

## Install `bcl-convert` on the Linux machine

Download from Illumina [BCL Convert Software Downloads](https://support.illumina.com/sequencing/sequencing_software/bcl-convert/downloads.html) by selecting the latest CentOS version (here we use version 4.3.6 as an example). You will need to log in and enter a machine serial number (just use `LH00328`) to download it. Put the installer onto your Linux machine (for example with the name `bcl-convert-4.3.6.rpm`)

Convert its format, and install by running (make sure `alien` is installed).

```bash
sudo alien bcl-convert-4.3.6.rpm
sudo dpkg -i bcl-convert_4.3.6-3_amd64.deb
```

Test installment by running:

```bash
bcl-convert --version
```

## Get data from Novogene

Do this as soon as you receive email. Data expires after 7 days.

On your Linux machine, run the `wget` command given in the email. You will get the following directory structure:

```
├── usftp21.novogene.com
│   ├── <run_name>.BCL.tar  # the sequencing data
│   ├── MD5.txt
│   ├── Report_<a_random_string>.zip  # Novogene QC report
│   └── checkSize.xls
└── wget-log
```

Untar the sequencing data archive with

```
tar xf <run_name>.BCL.tar
```

and you will get a Illumina sequencer output data directory named by your run.

These steps can take several hours.

## Prepare sample sheet

An example sample sheet that can be used as a template is like this:

```
[Header],,,,
FileFormatVersion,2,,,
RunName,some_name,,,  # Name of your run
InstrumentPlatform,Novogene,,,
InstrumentType,NovaSeq X Plus,,,
,,,,
[Reads],,,,
Read1Cycles,151,,,
Read2Cycles,151,,,
Index1Cycles,10,,,
Index2Cycles,24,,,
,,,,
[BCLConvert_Settings],,,,
SoftwareVersion,4.3.6,,,
OverrideCycles,Y151;I10;N14I10;Y151,,,  # Cycle configuration of your run
,,,,
[BCLConvert_Data],,,,  # everything underneath
Lane,Sample_ID,index,index2,
4,sample_1,TAAGGCGAAT,ACCTAAGCCT,
4,sample_2,CGTACTAGAT,ACCTAAGCCT,
4,sample_3,AGGCAGAAAT,ACCTAAGCCT,
```

I am marking lines you most likely need to customize with a # symbol. Note that the # is not needed in the actual sample sheet. Now I will explain the content in the sample sheet.

### Name your experiment

Name properly to make sample sheet sharing or archiving easier. 

### Know your Instrument

In house, we use MiSeq the most. For Novogene, our sequencer model is NovaSeq X Plus.

### Figure out which lane is used

In our current setting, you are allocated only one lane on the NovaSeq flow cell. As Illumina index is independent across lanes, sample sheet requires lane-index combination for demultiplexing.

Inspect content of `usftp21.novogene.com/<run_name>/Data/Intensities/BaseCalls/`. You should see only one sub-folder, and from its name you can tell the lane number (e.g. `L001`, meaning your sample is on Lane 1).

### Decide cycle configuration

Illumina sequencing is composed of multiple rounds of cycles. Each cycle generates readout for one base, and thus cycle length is equivalent to sequence length. Depending on sequencing protocols, there can be these 4 rounds:

1. Read 1 round
2. Index 1 (aka i7) round
3. Index 2 (aka i5) round
4. Read 2 round

and each round can have difference cycle number. The cycle configuration of your run can be found in `usftp21.novogene.com/<run_name>/RunInfo.xml`, which says:

```
...
                <Reads>
                        <Read Number="1" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
                        <Read Number="2" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="N"/>
                        <Read Number="3" NumCycles="24" IsIndexedRead="Y" IsReverseComplement="Y"/>
                        <Read Number="4" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
                </Reads>
...
```

This means your are using **paired-end 151bp sequencing with dual indexing**. Your index 1 round has 10 cycles and index 2 round has 24 cycles.

However, our Illumina index sequences are either 8bp or 10bp. That's why we need to tell the software which cycles in the index rounds actually correspond to our index.

There is a specific syntax to follow to "tell the software". Here are two examples:

1. If your index sequences are 8bp, write in sample sheet
   ```
   OverrideCycles,Y151;I8N2;N16I8;Y151,,, 
   ```

2. If your index sequences are 10bp, write

   ```
   OverrideCycles,Y151;I10;N14I10;Y151,,,
   ```

#### Another example

Another example with a slightly different `usftp21.novogene.com/<run_name>/RunInfo.xml` (`NumCycles="10"` for `Read Number="3"`):

```
...
                <Reads>
                        <Read Number="1" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
                        <Read Number="2" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="N"/>
                        <Read Number="3" NumCycles="10" IsIndexedRead="Y" IsReverseComplement="Y"/>
                        <Read Number="4" NumCycles="151" IsIndexedRead="N" IsReverseComplement="N"/>
                </Reads>
...
```

1. If your index sequences are 8bp, write in sample sheet

   ```
   OverrideCycles,Y151;I8N2;N2I8;Y151,,, 
   ```

2. If your index sequences are 10bp, write

   ```
   OverrideCycles,Y151;I10;I10;Y151,,,
   ```

3. If 8bp and 10bp are mixed, also write

   ```
   OverrideCycles,Y151;I10;I10;Y151,,,
   ```

   And refer to [mixed -length primer](#mixed-length-primer) about index sequence modification.

### Fetch your index sequences

Assemble i7 (index) and i5 (index 2) sequences based on your plate layout during library construction.

> Note: with `bcl-convert`, we always write the index sequences in forward orientation regardless of sequencers. This behavior is different from `bcl2fastq`.

> It is tedious but important to distinguish these following 5 kinds of index sequences:
>
> - <u>bases in adapter</u>
> - <u>bases in forward orientation</u>
> - <u>bases in reverse orientation</u>
> - <u>bases in sample sheet</u>
> - <u>bases output by demultiplexing software</u> (e.g., `bcl2fastq` or `bcl-convert`).
>
> For i7, <u>bases in adapter</u> is always the rc (reverse complement) of <u>bases in sample sheet</u>.
>
> For i5, <u>bases in adapter</u> is equal to <u>bases in forward orientation</u> and the rc of <u>bases in reverse orientation</u>. The bases in sample sheet depends on sequencing instrument configuration ("index first" or "read first"), demultiplexing software, and sample sheet version.
>
> Therefore, using `bcl-convert` and sample sheet v2, <u>i7 in sample sheet</u> should be the rc of <u>i7 in adapter</u>, and <u>i5 in sample sheet</u> should be equal to <u>i5 in adapter</u>.
>
> For complete information on i5 orientation, refer to [Illumina](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/AdapterSeq/Overview.htm).

#### Mixed-length primer

If both 8bp and 10bp indices are used in the sample run, add `AT` to the end of i7 (`[i7]AT`) and `AC` to the start of i5 (`AC[i5]`). These are the last two bases of P7 (after taking reverse complement) and P5, respectively.

#### Compile the sample sheet

After you figure out the above items, you can compile your sample sheet by modifying from the template sample sheet. Given our common Novogene setup, you will just have to change the run name, cycle configuration, and index data.

## Run `bcl-convert`

Increase operating system's limit on opened files by running

```bash
ulimit -n 1000000
```

Then run:

```bash
bcl-convert \
    --bcl-input-directory <parent_dir>/usftp21.novogene.com/<run_name> \
    --output-directory <absolute_path_to_fastq_output_dir> \
    --sample-sheet <parent_dir>/usftp21.novogene.com/<run_name>/SampleSheet.csv \
    --force \
    --bcl-num-parallel-tiles 1
```

Replace the path arguments according to your machine. The execution can take an hour.

## Examine output file size

Check fastq file size by running

```bash
ls <fastq_output_dir> -lahS
```

File size is a reasonable approximation for sequencing depth, so by looking at file sizes you can confirm a couple of things:

1. Do most of you samples have nonzero sequencing depth?

   Why are there samples with zero sequencing depth? Are those the problematic samples during library construction? Are the index sequences in sample sheet actually incorrect? Or it's just that some samples got unlucky during sequencing.

2. Does it look like that your samples have a relatively uniform sequencing depth? By random chance, some samples will get more sequences than others, but if the difference is too drastic or there are too many outliers, try to be more cautious with your library normalization (qPCR or Qubit) next time.

3. Are there many reads where demultiplexing failed? These reads will go into the undetermined files. It is fine to have big the undetermined files (very likely they are the biggest, way bigger than the next biggest files). You are good as long as it does not account for a unreasonably huge proportion of all reads (like 60 or 70%).

## Upload data onto AWS

The fastq files will be mostly what we need in the future, but for the sake of data archiving, upload all data except for `.cbcl` files (since there are too many of these small files, very unfriendly to AWS) onto your AWS bucket.

```bash
aws s3 sync . s3://<your_bucket>/<directory_name_you_like>  --exclude "*.cbcl"
```

(assuming both the Novogene data and fastq output are in the current directory)

## Terminate AWS instance

If you opened an AWS instance for the demultiplexing, now you can terminate the instance.

---

The end
