
# Sub-Script Description


## Review

Rhea_Chip inclusion of some auxiliary scripts, However, these scripts are available used independently.

## Description

###. bam_to_chip_gender.py : Use BAM file to determine gender
###. depthQC.py : For QC, Exon plot base on bam file
###. vcf_phasing.py : VCF phasing

## Usage
```
# bam_to_chip_gender.py
python bam_to_chip_gender.py <bamfile> <outfile>

#vcf_phasing.py
python vcf_phasing.py <bamfile> <VcfIn> <VcfOut>

# depthQC.py 
python depthQC.py [-h] [--version] --bam [bamfile] --config [config.json] -o [outdir]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir=OUTDIR
                        The output dir [ default:
                        /Users/huang/Desktop/Rhea_Chip/bin ]
  -b BAMFILE, --bam=BAMFILE
                        The bamfile (required)
  -c CONFIG, --config=CONFIG
                        The config file (required)
  -p, --plot            Plot exon depth graph

```

## License

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

