# PoreCamp

## File formats

Oxford Nanopore are very bad at releasing official definitions of file formats, therefore much guess work is involved.

Most of the early ONT data was released from the SQK-MAP-005 kits - this includes 
[the MARC data](http://f1000research.com/articles/4-1075/v1), [Mick's B fragilis dataset](http://gigadb.org/dataset/100177) and [Nick's first E coli dataset](http://gigadb.org/dataset/100102).  These data were encoded in what can be best described as FAST5 v.1.0 (ONT don't actually assign version numbers!)

Then SQK-MAP-006 came along, which was a major chemistry change that increased throughput.  A major change is that metrichor has switched from a 5mer model to a 6mer model.   [Nick has also released E coli SQK-MAP-006 data](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/), and because he was very quick to do this, the files are still in FAST5 v.1.0

However, in November 2015, ONT released a new file format, which we can call FAST5 v1.1.  The major difference is that the template and complement FASTQ and events data have been moved to a new group within the FAST5 file, separate to the 2D data.  This actually makes logical sense, but can make data analysis difficult.

We will work here with SQK-MAP-006 data from Nick (so V1.0)

We will also hopefully work with data generated on the course (so V1.1)

#### Major file format differences

The major difference is where the template and complement data are.  In version 1.0 they are all in a group called Basecall_2D_000; however, in v1.1 they have been moved to Basecall_1D_000

FAST5 v1.0
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**2D**_000/BaseCalled_template/
* /Analyses/Basecall_**2D**_000/BaseCalled_complement/

FAST5 v1.1
* /Analyses/Basecall_**2D**_000/BaseCalled_2D/
* /Analyses/Basecall_**1D**_000/BaseCalled_template/
* /Analyses/Basecall_**1D**_000/BaseCalled_complement/



## File handling

MAP006-1 contains some SQK-MAP-006 data from the Loman lab.  There are some ~69000 files in the directory therefore ls can be uninformative.  However, we can use find to reveal directory structure

```sh
find MAP006-1/ -type d
```
```
MAP006-1/
MAP006-1/logs
MAP006-1/MAP006-1_downloads
MAP006-1/MAP006-1_downloads/fail
MAP006-1/MAP006-1_downloads/pass
```

The top-level directory contains fast5 files that have yet to be basecalled
The MAP006-1_downloads directory contains fast5 processed by metrichor
The pass and fail directories contain files that metrichor has determined have passed/failed

## Using HDF5 command line tools

h5ls and h5dump can be quite useful.

h5ls reveals the structure of fast5 files.  Adding the -r flag makes this recursive:

```sh
h5ls MAP006-1/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch480_file17_strand.fast5
```
```
Analyses                 Group
Sequences                Group
UniqueGlobalKey          Group
```

 Adding the -r flag makes this recursive
 ```sh
 h5ls -r MAP006-1/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch480_file17_strand.fast5
 ```
 ```
/                        Group
/Analyses                Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/abasic_detection Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/hairpin_detection Group
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_16 Group
/Analyses/EventDetection_000/Reads/Read_16/Events Dataset {1614/Inf}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can compare this to a base-called file

```
h5ls -r MAP006-1/MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand.fast5
```
```
/                        Group
/Analyses                Group
/Analyses/Basecall_2D_000 Group
/Analyses/Basecall_2D_000/BaseCalled_2D Group
/Analyses/Basecall_2D_000/BaseCalled_2D/Alignment Dataset {8394}
/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_complement Group
/Analyses/Basecall_2D_000/BaseCalled_complement/Events Dataset {6562}
/Analyses/Basecall_2D_000/BaseCalled_complement/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_complement/Model Dataset {4096}
/Analyses/Basecall_2D_000/BaseCalled_template Group
/Analyses/Basecall_2D_000/BaseCalled_template/Events Dataset {7052}
/Analyses/Basecall_2D_000/BaseCalled_template/Fastq Dataset {SCALAR}
/Analyses/Basecall_2D_000/BaseCalled_template/Model Dataset {4096}
/Analyses/Basecall_2D_000/Configuration Group
/Analyses/Basecall_2D_000/Configuration/aggregator Group
/Analyses/Basecall_2D_000/Configuration/basecall_1d Group
/Analyses/Basecall_2D_000/Configuration/basecall_2d Group
/Analyses/Basecall_2D_000/Configuration/calibration_strand Group
/Analyses/Basecall_2D_000/Configuration/general Group
/Analyses/Basecall_2D_000/Configuration/hairpin_align Group
/Analyses/Basecall_2D_000/Configuration/post_processing Group
/Analyses/Basecall_2D_000/Configuration/post_processing.3000Hz Group
/Analyses/Basecall_2D_000/Configuration/post_processing.5000Hz Group
/Analyses/Basecall_2D_000/Configuration/recipes Group
/Analyses/Basecall_2D_000/Configuration/split_hairpin Group
/Analyses/Basecall_2D_000/HairpinAlign Group
/Analyses/Basecall_2D_000/HairpinAlign/Alignment Dataset {6072}
/Analyses/Basecall_2D_000/InputEvents Soft Link {Analyses/EventDetection_000/Reads/Read_23/Events}
/Analyses/Basecall_2D_000/Log Dataset {SCALAR}
/Analyses/Basecall_2D_000/Summary Group
/Analyses/Basecall_2D_000/Summary/basecall_1d_complement Group
/Analyses/Basecall_2D_000/Summary/basecall_1d_template Group
/Analyses/Basecall_2D_000/Summary/basecall_2d Group
/Analyses/Basecall_2D_000/Summary/hairpin_align Group
/Analyses/Basecall_2D_000/Summary/post_process_complement Group
/Analyses/Basecall_2D_000/Summary/post_process_template Group
/Analyses/Basecall_2D_000/Summary/split_hairpin Group
/Analyses/Calibration_Strand_000 Group
/Analyses/Calibration_Strand_000/Configuration Group
/Analyses/Calibration_Strand_000/Configuration/aggregator Group
/Analyses/Calibration_Strand_000/Configuration/basecall_1d Group
/Analyses/Calibration_Strand_000/Configuration/basecall_2d Group
/Analyses/Calibration_Strand_000/Configuration/calibration_strand Group
/Analyses/Calibration_Strand_000/Configuration/general Group
/Analyses/Calibration_Strand_000/Configuration/hairpin_align Group
/Analyses/Calibration_Strand_000/Configuration/post_processing Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.3000Hz Group
/Analyses/Calibration_Strand_000/Configuration/post_processing.5000Hz Group
/Analyses/Calibration_Strand_000/Configuration/recipes Group
/Analyses/Calibration_Strand_000/Configuration/split_hairpin Group
/Analyses/Calibration_Strand_000/InputEvents Soft Link {Analyses/EventDetection_000/Reads/Read_23/Events}
/Analyses/Calibration_Strand_000/Log Dataset {SCALAR}
/Analyses/Calibration_Strand_000/Summary Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_2d Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_complement Group
/Analyses/Calibration_Strand_000/Summary/calibration_strand_template Group
/Analyses/EventDetection_000 Group
/Analyses/EventDetection_000/Configuration Group
/Analyses/EventDetection_000/Configuration/abasic_detection Group
/Analyses/EventDetection_000/Configuration/event_detection Group
/Analyses/EventDetection_000/Configuration/hairpin_detection Group
/Analyses/EventDetection_000/Reads Group
/Analyses/EventDetection_000/Reads/Read_23 Group
/Analyses/EventDetection_000/Reads/Read_23/Events Dataset {13677}
/Sequences               Group
/Sequences/Meta          Group
/UniqueGlobalKey         Group
/UniqueGlobalKey/channel_id Group
/UniqueGlobalKey/context_tags Group
/UniqueGlobalKey/tracking_id Group
```

We can also use h5ls -d to extract specific datasets from the files:
```sh
 h5ls -d MAP006-1/MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand.fast5/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq
 ```
 ```
 Fastq                    Dataset {SCALAR}
    Data:
        (0) "@72e97803-06e3-44db-8000-287bb520a673_Basecall_2D_000_2d LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand\nGTCTCTGTCTGTATCTGATATTGCTGTAACT
   ...
   <snip>
   ...
   &&'''(&&%&%'''&&'&&''&'&&(&(%('&%&%&'&&%%%\n"
```

Unsurprisingly h5dump dumps the entire file to STDOUT

```sh
h5dump MAP006-1/MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand.fast5
```

## Browsing HDF5 files

Any HDF5 file can be opened using hdfview and browsed/edited in a GUI

```sh
hdfview MAP006-1/MAP006-1_downloads/pass/LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch150_file24_strand.fast5 &
```

## Basic manipulation in poRe

poRe is a library for R available from [SourceForge](https://sourceforge.net/projects/rpore/) and published in [bioinformatics](http://bioinformatics.oxfordjournals.org/content/31/1/114).  poRe is incredibly simple to install and relies simply on R 3.0 or above and a few additional libraries.

The poRe library is set up to read v1.1 data by default, and offers users parameters to enable reading of v1.0 data.  Let's start it up.

```sh
R
```
```R
library(poRe)
```

#### FASTQ

We'll see how to extract FASTQ from entire directories below, but here are some exemples of single file analysis

```R
f5 <- "MAP006-1/MAP006-1_downloads/pass//LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch300_file44_strand.fast5"
get_fastq(f5)
```

This returns a list with a single fastq datasets, "2D".  However, where are template and complement?  Unfortunately the defaults are set to v1.1 and this file is v1.0.  We can use the path.t and path.c parameters to tell poRe where the template and complement data are

```R
get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", path.c="/Analyses/Basecall_2D_000/")
fqs <- get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", path.c="/Analyses/Basecall_2D_000/")
names(fqs)
```

If we don't want to extract all 3, we can choose which to extract using the "which" argument

```R
get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", which="template")
```

Working with lists is easy in R, and if you want the FASTQ as a string:

```R
fq <- get_fastq(f5, path.t="/Analyses/Basecall_2D_000/", which="all")
fq$template
```

From here you can see that it's incredibly simple to write a FASTQ extraction script:

```R
# path to fast5 files
f5dir <- "test"
# get vector of all fast5 files
f5files <- dir(f5dir, pattern="\\.fast5$", full.names = TRUE)
# iterate over files
for (f5 in f5files) {
    # extract 2D fastq
    fq <- get_fastq(f5, which="2D")
    # check fq is a list and contains 2D
    if (typeof(fq) == "list" && exists("2D", where=fq)) {
        # cat to "" (STDOUT but could be the name of a file
        cat(fq[["2D"]], file = "", sep = "\n", fill = FALSE)
    }
}
```

However, we have scripts that do this already, so there is no need to write your own!

#### FASTA data

Simply use get_fasta instead of get_fastq :-)

#### Events data

We extract events data with the get.events function.  As events data are co-located with FASTQ in the FAST5 file, then for Nick's SQK-MAP-006 data we need to give it the paths again.

get.events again returns a list, wth the template and complement events as data frames

```R
f5 <- "MAP006-1/MAP006-1_downloads/pass//LomanLabz_PC_Ecoli_K12_MG1655_20150924_MAP006_1_5005_1_ch300_file44_strand.fast5"
ev <- get.events(f5, path.t = "/Analyses/Basecall_2D_000/", path.c = "/Analyses/Basecall_2D_000/")
names(ev)
head(ev$template)
head(ev$complement)
```

* start: time in seconds
* mean: the mean signal of the event
* stdev: the standard deviation when sampling
* length: length of the event
* model_state: the predicted kmer
* model_level: the value from the model
* move: how many moves the model has made in that step

#### Extracting the model

Extracting the model itself in poRe is also easy, but relies on a little legacy code that needs to be updated... next version ;-)

```R
mods <- get.models(f5, tmodel = "/Analyses/Basecall_2D_000/BaseCalled_template/Model", 
                        cmodel = "/Analyses/Basecall_2D_000/BaseCalled_complement/Model")
head(mods$template)
head(mods$complement)
```

Here we can see the models that ONT use to turn events into kmers, but this isn't simply a case of looking up the closest match:

```R
head(ev$template)
# first event kmer is TGTTTC, event mean is 53.08718, model mean is 53.64642
#       mean    start      stdv      length model_state model_level move
#   53.08718 21803.40 0.7518740 0.007636122      TTTTTC    53.64642    0

mods$template[mods$template$kmer=="TTTTTC",]

# however the model level for this kmer is quite different:
#       kmer level_mean level_stdv  sd_mean  sd_stdv   weight
#     TTTTTC   42.37862   0.518937 0.700281 0.223579 1805.695

```

What gives?  Well, ONT apply three parameters to the model to "scale" the model to fit each read.  These parameters are drift, scale and shift and they are read dependent!

For this read, we have:
* drift	0.004207995
* scale	1.034605
* shift	9.801294

The drift parameter is applied per event and is scaled by the number of seconds since the start of the read.  So "drift * (time - min(time))" can be subtracted from the event mean, or added to the model mean.

Scale and shift are then applied in a classic linear model: (model.mean * scale) + shift.  In our case:

* (42.37862 * 1.034605) + 9.801294 = 53.64643

Which is the model value that shows up in the events table above. 

#### Extracting the model parameters

For v1.0:

```R
get.model.params(f5, tsum="/Analyses/Basecall_2D_000/Summary/basecall_1d_template", 
                     csum = "/Analyses/Basecall_2D_000/Summary/basecall_1d_complement")
```

For v1.1:
```R
get.model.params(f5)
```

#### Plotting Squiggles

Plotting squiggles works directly from the events data.  The parameters minseconds and maxseconds control which events to plot:

```R
plot.squiggle(ev$template)
plot.squiggle(ev$template, minseconds=5)
plot.squiggle(ev$complement)
```

## Run QC in poRe

If you ran pore_rt() during your nanopore run then you will have access to a metadata text file that can be used for run QC etc.  Otherwise we will have to create one (more below).  In the meantime, we have meta data for the pass and fail folders for Nick's SQK-MAP-006 run

```R
# load in pass data
pass <- read.table("pass.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(pass) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")

# load in the fail data
fail <- read.table("fail.meta.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

# set standard/expected column names
colnames(fail) <- c("filename","channel_num","read_num","read_start_time",
                    "status","tlen","clen","len2d",
                    "run_id","read_id","barcode","exp_start")
head(pass)
head(fail)
```

### yield

Yield over time can be plotted with plot.cumulative.yield

```R
yield.p <- plot.cumulative.yield(pass)
yield.f <- plot.cumulative.yield(fail)
```

The calculated cumulative yields are returned as data.frames

```R
head(yield.p)
head(yield.f)
```

To plot for a single barcode (Nick's SQK-MAP-006 data is not barcoded):

```R
barcode <- "no_barcode"
plot.cumulative.yield(pass[pass$barcode==barcode,])
```

### length histogram

We can see a histogram of read lengths

```R
plot.length.histogram(pass)
plot.length.histogram(fail)
```

There's quite a long faled template read!

```R
max(fail$tlen)
# [1] 379692
fail [fail$tlen==379692,]
```

Plotting in ggplot2 if you really want to

```R
m <- ggplot(pass, aes(x=len2d))
m + geom_histogram(binwidth=500)
```

### occupancy

The MinION flowcell is arranged into 512 channels in 4 blocks, and we can see the layout using poRe:

```R
show.layout()
```

We first need to calculate some statistics summarised by channel:

```R
# pass data
pass.s <- summarise.by.channel(pass)
head(pass.s)

fail.s <- summarise.by.channel(fail)
head(fail.s)
```

The rows of the result are the channel numbers, and the columns tell us how many channels appear in our summary data, and either the number (n) or sumulative length (l) of template, complement and 2d reads from each channel.

We can plot these:

```R
# pass

# the number of times the channel appears in any context
plot.channel.summary(pass.s)

# cumulative 2D length
plot.channel.summary(pass.s, report.col="l2d")

# fail

# the number of times the channel appears in any context
plot.channel.summary(fail.s)

# cumulative 2D length
plot.channel.summary(fail.s, report.col="l2d")

```

### Extracting meta-data from fast5
If you haven't run pore_rt(), then you can extract meta-data directly from the fast5 files.  This takes a long time as we have top open each fle and extract the attributes.  For demo purposes, we have selected 10% randomly of Nick's SQK-MAP-006 data in folder MAP006-1.10 

```R
pass <- read.meta.info("MAP006-1.10/MAP006-1_downloads/pass/", 
                        path.t="/Analyses/Basecall_2D_000/", 
                        path.c="/Analyses/Basecall_2D_000/", 
                        pass="P")
                        
fail <-  read.meta.info("MAP006-1.10/MAP006-1_downloads/fail/", 
                        path.t="/Analyses/Basecall_2D_000/", 
                        path.c="/Analyses/Basecall_2D_000/", 
                        pass="F")
                        
plot.cumulative.yield(pass)
plot.cumulative.yield(fail)
```

## Extracting FASTQ from the command-line


Command-line scripts for extracting FASTQ can be pulled from [github](https://github.com/mw55309/poRe_scripts).  For Nick's SQK-MAP-006 data we use the old_format scripts:

```sh
# 2D
./poRe_scripts/old_format/extract2D MAP006-1/MAP006-1_downloads/pass/

# template
./poRe_scripts/old_format/extractTemplate MAP006-1/MAP006-1_downloads/pass/

# complement
./poRe_scripts/old_format/extractComplement MAP006-1/MAP006-1_downloads/pass/
```

For any newer data, we can use the scripts in new_format


