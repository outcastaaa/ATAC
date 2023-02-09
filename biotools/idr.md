[详细参数](https://github.com/nboley/idr)  
[原理解读](https://mbd.baidu.com/ug_share/mbox/4a83aa9e65/share?product=smartapp&tk=d8a8687f3aa19279e39c5415ccedf749&share_url=https%3A%2F%2Fkfs479.smartapps.cn%2Fpages%2Fblogdetail%2Fblogdetail%3Fid%3D5398770%26_swebfr%3D1%26_swebFromHost%3Dbaiduboxapp&domain=mbd.baidu.com)
# 详细参数
```bash
$ idr -h
usage: idr [-h] --samples SAMPLES SAMPLES [--peak-list PEAK_LIST] [--input-file-type {narrowPeak,broadPeak,bed,gff}]
           [--rank RANK] [--output-file OUTPUT_FILE] [--output-file-type {narrowPeak,broadPeak,bed}]
           [--log-output-file LOG_OUTPUT_FILE] [--idr-threshold IDR_THRESHOLD]
           [--soft-idr-threshold SOFT_IDR_THRESHOLD] [--use-old-output-format] [--plot] [--use-nonoverlapping-peaks]
           [--peak-merge-method {sum,avg,min,max}] [--initial-mu INITIAL_MU] [--initial-sigma INITIAL_SIGMA]
           [--initial-rho INITIAL_RHO] [--initial-mix-param INITIAL_MIX_PARAM] [--fix-mu] [--fix-sigma]
           [--dont-filter-peaks-below-noise-mean] [--use-best-multisummit-IDR] [--allow-negative-scores]
           [--random-seed RANDOM_SEED] [--max-iter MAX_ITER] [--convergence-eps CONVERGENCE_EPS] [--only-merge-peaks]
           [--verbose] [--quiet] [--version]

Program: IDR (Irreproducible Discovery Rate)
Version: 2.0.3
Contact: Nathan Boley <npboley@gmail.com>

optional arguments:
  -h, --help            show this help message and exit
  --samples SAMPLES SAMPLES, -s SAMPLES SAMPLES
                        Files containing peaks and scores.
  --peak-list PEAK_LIST, -p PEAK_LIST
                        If provided, all peaks will be taken from this file.
  --input-file-type {narrowPeak,broadPeak,bed,gff}
                        File type of --samples and --peak-list.
  --rank RANK           Which column to use to rank peaks.
                        Options: signal.value p.value q.value columnIndex
                        Defaults:
                                narrowPeak/broadPeak: signal.value
                                bed: score
  --output-file OUTPUT_FILE, -o OUTPUT_FILE
                        File to write output to.
                        Default: idrValues.txt
  --output-file-type {narrowPeak,broadPeak,bed}
                        Output file type. Defaults to input file type when available, otherwise bed.
  --log-output-file LOG_OUTPUT_FILE, -l LOG_OUTPUT_FILE
                        File to write output to. Default: stderr
  --idr-threshold IDR_THRESHOLD, -i IDR_THRESHOLD
                        Only return peaks with a global idr threshold below this value.
                        Default: report all peaks
  --soft-idr-threshold SOFT_IDR_THRESHOLD
                        Report statistics for peaks with a global idr below this value but return all peaks with an idr below --idr.
                        Default: 0.05
  --use-old-output-format
                        Use old output format.
  --plot                Plot the results to [OFNAME].png
  --use-nonoverlapping-peaks
                        Use peaks without an overlapping match and set the value to 0.
  --peak-merge-method {sum,avg,min,max}
                        Which method to use for merging peaks.
                                Default: 'sum' for signal/score/column indexes, 'min' for p/q-value.
  --initial-mu INITIAL_MU
                        Initial value of mu. Default: 0.10
  --initial-sigma INITIAL_SIGMA
                        Initial value of sigma. Default: 1.00
  --initial-rho INITIAL_RHO
                        Initial value of rho. Default: 0.20
  --initial-mix-param INITIAL_MIX_PARAM
                        Initial value of the mixture params. Default: 0.50
  --fix-mu              Fix mu to the starting point and do not let it vary.
  --fix-sigma           Fix sigma to the starting point and do not let it vary.
  --dont-filter-peaks-below-noise-mean
                        Allow signal points that are below the noise mean (should only be used if you know what you are doing).
  --use-best-multisummit-IDR
                        Set the IDR value for a group of multi summit peaks (same chr/start/stop but different summit) to the best value across all peaks. This is a work around for peak callers that don't do a good job splitting scores across multi summit peaks.
  --allow-negative-scores
                        Allow negative values for scores. (should only be used if you know what you are doing)
  --random-seed RANDOM_SEED
                        The random seed value (sor braking ties). Default: 0
  --max-iter MAX_ITER   The maximum number of optimization iterations. Default: 3000
  --convergence-eps CONVERGENCE_EPS
                        The maximum change in parameter value changes for convergence. Default: 1.00e-06
  --only-merge-peaks    Only return the merged peak list.
  --verbose             Print out additional debug information
  --quiet               Don't print any status messages
  --version             show program's version number and exit
/home/linuxbrew/.linuxbrew/opt/python@3.9/bin/idr -h
```
