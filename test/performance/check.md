# Performance testing

In this document we want to make sure we don't regress on speed.

On a `AMD Ryzen 7 3700X 8-Core Processor`:

2 cores:

```
        Command being timed: "bin/smoothxg -t 2 -g test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.seqwish.gfa -j 5k -e 5k -l 700,900,1100 -m test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.maf -C consensus,10,100:test/data/gi_568815592_32578768-32589835.txt:y,1000:test/data/gi_568815592_32578768-32589835.txt:n,10000 -o test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.gfa -r 12"
        User time (seconds): 5.75
        System time (seconds): 1.06
        Percent of CPU this job got: 26%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:25.45
```

Note the debug version is about 35s.

8 cores does not make much difference:

```
Command being timed: "bin/smoothxg -t 8 -g test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.seqwish.gfa -j 5k -e 5k -l 700,900,1100 -m test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.maf -C consensus,10,100:test/data/gi_568815592_32578768-32589835.txt:y,1000:test/data/gi_568815592_32578768-32589835.txt:n,10000 -o test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.gfa -r 12"
        User time (seconds): 12.42
        System time (seconds): 4.85
        Percent of CPU this job got: 73%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:23.42
```

Compiling with LTO creates a slightly faster runtime on 2 threads

```
        Command being timed: "bin/smoothxg -t 2 -g test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.seqwish.gfa -j 5k -e 5k -l 700,900,1100 -m test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.maf -C consensus,10,100:test/data/gi_568815592_32578768-32589835.txt:y,1000:test/data/gi_568815592_32578768-32589835.txt:n,10000 -o test/data/DRB1-3123.fa.gz.pggb-s3000-p70-n10-a70-K16-k8-w10000-j5000-e5000.smooth.gfa -r 12"
        User time (seconds): 5.43
        System time (seconds): 1.19
        Percent of CPU this job got: 26%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:25.01

```

Honoring -Ofast gives some speedup

```
        User time (seconds): 5.55
        System time (seconds): 1.02
        Percent of CPU this job got: 26%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:24.40
```

The static build with GNU Guix is same

```
        User time (seconds): 5.35
        System time (seconds): 1.12
        Percent of CPU this job got: 26%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:24.46
```
